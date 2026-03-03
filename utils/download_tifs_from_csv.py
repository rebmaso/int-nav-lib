import os
import pandas as pd
import requests
from urllib.parse import urlparse

def download_tifs_from_csv(csv_path, column_idx, download_folder):
    # 1. Create the directory if it doesn't exist
    if not os.path.exists(download_folder):
        os.makedirs(download_folder)
        print(f"Created directory: {download_folder}")

    # 2. Load the CSV
    try:
        df = pd.read_csv(csv_path)
    except Exception as e:
        print(f"Error reading CSV: {e}")
        return

    # 3. Iterate through the URLs
    for index, row in df.iterrows():
        url = row[column_idx]
        
        # Basic check to ensure it's a string and looks like a TIF
        if isinstance(url, str) and url.lower().endswith(('.tif', '.tiff')):
            # Extract filename from URL
            filename = os.path.basename(urlparse(url).path)
            save_path = os.path.join(download_folder, filename)
            temp_save_path = save_path + ".part"

            # Skip if file already exists
            if os.path.exists(save_path):
                print(f"Skipping already downloaded file: {filename}")
                continue

            try:
                print(f"Downloading: {filename}...")
                response = requests.get(url, stream=True, timeout=60) # Increased timeout
                response.raise_for_status() # Check for HTTP errors

                total_size = int(response.headers.get('content-length', 0))

                with open(temp_save_path, 'wb') as f:
                    downloaded_size = 0
                    for chunk in response.iter_content(chunk_size=8192):
                        f.write(chunk)
                        downloaded_size += len(chunk)
                
                if total_size != 0 and downloaded_size < total_size:
                    raise IOError(f"Incomplete download. Expected {total_size}, got {downloaded_size}")

                os.rename(temp_save_path, save_path)
                print(f"Successfully saved to {save_path}")

            except Exception as e:
                print(f"Failed to download {url}. Error: {e}")
                # Clean up partial file if it exists
                if os.path.exists(temp_save_path):
                    os.remove(temp_save_path)
        else:
            print(f"Skipping non-TIF URL at row {index}: {url}")

# --- Configuration ---
CSV_FILE = '/home/tommaso/Downloads/ch.swisstopo.swissimage-dop10-wADEOh8l.csv'          # Path to your CSV
COLUMN_IDX = 0      # The header name where URLs are stored
DESTINATION = '/home/tommaso/Downloads/ch.swisstopo.swissimage-dop10-wADEOh8l' # Local folder name

if __name__ == "__main__":
    download_tifs_from_csv(CSV_FILE, COLUMN_IDX, DESTINATION)