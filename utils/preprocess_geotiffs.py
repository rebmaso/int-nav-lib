import os
import subprocess
import glob

def preprocess_osgearth(input_dir, output_name, file_ext="tif", data_type="ortho"):
    """
    Preprocesses a directory of imagery/DEMs for osgEarth.
    
    Args:
        input_dir (str): Path to folder containing source files.
        output_name (str): Base name for the output files.
        file_ext (str): Extension of source files (tif, tiff, dt1, etc).
        data_type (str): "ortho" for imagery or "dem" for elevation.
    """
    
    # Setup paths
    vrt_file = os.path.join(input_dir, f"{output_name}.vrt")
    warped_vrt = os.path.join(input_dir, f"{output_name}_warped.vrt")
    processed_tif = os.path.join(input_dir, f"{output_name}_optimized.tif")
    mbtiles_file = os.path.join(input_dir, f"{output_name}.mbtiles")
    
    # 1. Build VRT (Virtual Dataset)
    # Combines all files in the directory into one virtual source
    print(f"--- Creating VRT from {input_dir} ---")
    search_path = os.path.join(input_dir, f"*.{file_ext}")
    files = glob.glob(search_path)
    if not files:
        print(f"No files found matching {search_path}")
        return

    subprocess.run(["gdalbuildvrt", vrt_file] + files)

    # 2. Warp to Geodetic (EPSG:4326)
    # osgEarth is fastest when data is already in the globe's native SRS
    print("--- Reprojecting to WGS84 (EPSG:4326) ---")
    subprocess.run(["gdalwarp", "-t_srs", "epsg:4326", "-of", "VRT", vrt_file, warped_vrt])

    # 3. Translate to Tiled GeoTIFF
    # -co TILED=YES makes random access much faster
    print("--- Creating Internal Tiles ---")
    translate_cmd = ["gdal_translate", "-of", "GTiff", "-co", "TILED=YES", "-co", "BIGTIFF=YES", warped_vrt, processed_tif]
    
    if data_type.lower() == "ortho":
        # Imagery gets JPEG compression to save space
        translate_cmd.extend(["-co", "COMPRESS=JPEG"])
    else:
        # DEMs (Elevation) get no lossy compression to preserve height accuracy
        print("--- DEM Mode: Compression Disabled ---")
        
    subprocess.run(translate_cmd)

    # 4. Build Overviews (Pyramids)
    print("--- Generating Overviews ---")
    resampling = "average" if data_type.lower() == "ortho" else "bilinear"
    subprocess.run(["gdaladdo", "-r", resampling, processed_tif])

    # 5. Convert to MBTiles using osgearth_conv
    print("--- Converting to MBTiles (osgearth_conv) ---")
    driver_in = "gdalimage" if data_type.lower() == "ortho" else "gdalelevation"
    driver_out = "mbtilesimage" if data_type.lower() == "ortho" else "mbtileselevation"
    out_format = "jpg" if data_type.lower() == "ortho" else "tiff"
    conv_cmd = [
        "osgearth_conv",
        "--in", "driver", driver_in, "--in", "url", processed_tif,
        "--out", "driver", driver_out, "--out", "filename", mbtiles_file, 
        "--out", "format", out_format,
        # "--threads", "2",
        # "--max-level", "17"
    ]

    subprocess.run(conv_cmd)

    print(f"\nDone! Final output: {mbtiles_file}")

# --- USER PARAMETERS ---
if __name__ == "__main__":
    params = {
        "input_dir": "/home/tommaso/Datasets/swisstopo/ch.swisstopo.swissimage-dop10-wADEOh8l",   # Path to your files
        "output_name": "dataset",     # Final filename
        "file_ext": "tif",               # Extension to look for
        "data_type": "ortho"             # "ortho" or "dem"
    }

    preprocess_osgearth(**params)
