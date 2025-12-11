#!/bin/sh

# This script builds the IntNavLib library and all associated applications.
#
# Usage:
# ./build.sh [--clean | --coverage] [BuildType]
#
# Options:
#   --clean      Remove all build folders before building
#   --coverage   Build with code coverage enabled
#   
#   BuildType    "Debug" or "Release" (default: Release)
#
# The script will:
# 1. Build IntNavLib with the specified build type.
# 2. Install it to $HOME/local_install.
# 3. Build all applications in the 'apps' directory, linking against the
#    just-installed IntNavLib.

set -e

# --- Configuration ---
BUILD_TYPE="Release"
CLEAN_BUILD=false
BUILD_COVERAGE=false

# Parse args
for arg in "$@"; do
  case "$arg" in
    --clean)
      CLEAN_BUILD=true
      ;;
    --coverage)
      BUILD_COVERAGE=true
      ;;
    [Dd]ebug)
      BUILD_TYPE="Debug"
      ;;
    [Rr]elease)
      BUILD_TYPE="Release"
      ;;
  esac
done

INSTALL_PREFIX="$HOME/local_install"
REPO_ROOT="$(git rev-parse --show-toplevel)"

echo "================================================="
echo "IntNavLib Full Build Script"
echo "================================================="
echo "Build Type:           $BUILD_TYPE"
echo "Installation Prefix:  $INSTALL_PREFIX"
echo "Repository Root:      $REPO_ROOT"
echo "Clean Build:          $CLEAN_BUILD"
echo "Check Coverage        $BUILD_COVERAGE"
echo "================================================="

# --- Clean step ---
if [ "$CLEAN_BUILD" = true ]; then
  echo
  echo "--> Cleaning all build directories..."
  rm -rf "$REPO_ROOT/build"
  find "$REPO_ROOT/apps" -type d -name build -exec rm -rf {} +
  find "$REPO_ROOT/apps" -type d -name install -exec rm -rf {} +
  echo "--> Clean complete."
  echo
fi

# --- 1. Build and install IntNavLib library ---
echo "--> Building and installing IntNavLib library..."
echo

LIB_BUILD_DIR="$REPO_ROOT/build"
mkdir -p "$LIB_BUILD_DIR"
cd "$LIB_BUILD_DIR"

cmake .. -DCMAKE_BUILD_TYPE="$BUILD_TYPE" -DCMAKE_INSTALL_PREFIX="$INSTALL_PREFIX" -DBUILD_COVERAGE="$BUILD_COVERAGE" -DCMAKE_CXX_COMPILER_LAUNCHER=ccache

make -j$(nproc)
make install

echo
echo "--> IntNavLib installed successfully to $INSTALL_PREFIX"
echo

# --- 2. Build all applications ---
echo "--> Building all applications in 'apps' directory..."
echo

APPS_DIR="$REPO_ROOT/apps"

for app_path in "$APPS_DIR"/*; do
  if [ -d "$app_path" ]; then
    APP_NAME=$(basename "$app_path")
    echo "-------------------------------------------------"
    echo "--> Building application: $APP_NAME"
    echo "-------------------------------------------------"

    if [ -f "$app_path/package.xml" ]; then
      echo "--> Detected ROS 2 package. Using colcon."
      cd "$app_path"
      colcon build --install-base "install" \
        --cmake-args -DCMAKE_BUILD_TYPE="$BUILD_TYPE" -DCMAKE_PREFIX_PATH="$INSTALL_PREFIX" -DCMAKE_CXX_COMPILER_LAUNCHER=ccache
    else
      echo "--> Detected standard CMake project."
      APP_BUILD_DIR="$app_path/build"
      mkdir -p "$APP_BUILD_DIR"
      cd "$APP_BUILD_DIR"
      cmake .. -DCMAKE_BUILD_TYPE="$BUILD_TYPE" -DCMAKE_PREFIX_PATH="$INSTALL_PREFIX" -DCMAKE_CXX_COMPILER_LAUNCHER=ccache
      make -j$(nproc)
    fi
    echo "--> Finished building $APP_NAME."
    echo
  fi
done

echo "================================================="
echo "All builds completed successfully!"
echo "================================================="
