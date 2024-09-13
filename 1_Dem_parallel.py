# Afshin.Ghahramani@des.qld.gov.au
"""
This script automates the workflow of unzipping LAS files, generating individual DEMs using PDAL,
and merging them into a single GeoTIFF. Intermediate files are efficiently cleaned up after processing,
and parallel computing is employed for faster unzipping and DEM creation.
"""

# Afshin.Ghahramani@des.qld.gov.au
#PDAL is used in the generate_dem_with_pdal(las_file) function.
#GDAL is utilized indirectly through PDAL's "writers.gdal" stage and the rasterio library.

# Import necessary libraries
import os
import glob
import zipfile
import json
import subprocess
import rasterio
from rasterio.merge import merge
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm

# Define the directories for input data and intermediate processing
#input_directory = "C:\\Temp\\ell_Las\\"
input_directory = "C:\\Temp\\Rock_wood\\ell_Las\\"


#input_directory = "G:\\SCIENCE\\GullyScanning\\AERIAL_LiDAR data\\AerometrexLiDAR2023\\OP-005598_GullyStreambank\\Rookwood\\GDA2020_MGA_Z56_Rock_Wood_2308_Submission\\GDA2020_MGA_Z56_Rock_Wood_2308_Submission\\Las_1.4\\ell_las\\"

unzipped_directory = os.path.join(input_directory, "unzipped_files")

# Check if the directory for the unzipped LAS files exists. If not, create it.
if not os.path.exists(unzipped_directory):
    os.makedirs(unzipped_directory)


def unzip_las_files(zip_file):
    """
    Function to unzip LAS files.

    Args:
        zip_file (str): Path to the .zip file containing the LAS data.
    """
    # Extract all contents of the given zip file to the designated directory
    with zipfile.ZipFile(zip_file, 'r') as zip_ref:
        zip_ref.extractall(unzipped_directory)


def generate_dem_with_pdal(las_file):
    """
    Function to generate a Digital Elevation Model (DEM) from a LAS file using PDAL.

    Args:
        las_file (str): Path to the LAS file.

    Returns:
        str: Path to the generated DEM file in GeoTIFF format.
    """
    # Define the output DEM file path by replacing .las extension with _DEM.tif
    output_tif = las_file.replace('.las', '_DEM.tif')

    # Define the PDAL processing pipeline to convert LAS to DEM
    pipeline = {
        "pipeline": [
            {
                "type": "readers.las",  # Reader to load the LAS file
                "filename": las_file
            },
            {
                # Filter to retain only ground points (Classification code: 2)
                "type": "filters.range",
                "limits": "Classification[2:2]"
            },
            {
                # Use GDAL to generate a DEM with mean interpolation method
                "type": "writers.gdal",
                "gdaldriver": "GTiff",  # Specify output format as GeoTIFF
                "output_type": "mean",  # Use mean elevation for grid cells
                "resolution": 1.0,  # Define spatial resolution of the DEM
                "filename": output_tif
            }
        ]
    }

    # Convert the PDAL pipeline to a JSON string
    pipeline_str = json.dumps(pipeline)

    # Execute the PDAL pipeline using command-line
    subprocess.run(["pdal", "pipeline", "--stdin"], input=pipeline_str, text=True, check=True)

    # Return the path to the generated DEM
    return output_tif


def main():
    """
    Main function to handle the entire LAS processing workflow.
    """
    # Use parallel processing to unzip all LAS files for efficiency
    zip_files = glob.glob(os.path.join(input_directory, "*.zip"))
    with ProcessPoolExecutor() as executor:
        list(tqdm(executor.map(unzip_las_files, zip_files), total=len(zip_files), desc="Unzipping"))

    # Gather all unzipped LAS files for further processing
    las_files = glob.glob(os.path.join(unzipped_directory, "*.las"))

    # Use parallel processing to generate DEMs from each LAS file
    dem_files = []
    with ProcessPoolExecutor() as executor:
        futures = {executor.submit(generate_dem_with_pdal, las_file): las_file for las_file in las_files}
        for future in tqdm(as_completed(futures), total=len(las_files), desc="Processing LAS files with PDAL"):
            dem_files.append(future.result())

    # Merge the generated DEMs into a single combined DEM
    print("Merging DEMs...")

    # Open all DEMs and store them in a list
    src_files_to_mosaic = [rasterio.open(dem_file) for dem_file in dem_files]

    # Merge the DEMs
    merged_dem, out_trans = merge(src_files_to_mosaic)

    # Close individual DEMs after merging
    for src in src_files_to_mosaic:
        src.close()

    # Copy metadata from the first DEM file
    out_meta = src_files_to_mosaic[0].meta.copy()

    # Update the metadata with new dimensions, transform, and compression settings
    out_meta.update({
        "driver": "GTiff",
        "height": merged_dem.shape[1],
        "width": merged_dem.shape[2],
        "transform": out_trans,
        "compress": "lzw"
    })

    output_path = os.path.join(input_directory, "merged_DEM.tif")

    # Write the merged DEM to a GeoTIFF file
    with rasterio.open(output_path, 'w', **out_meta) as dest:
        dest.write(merged_dem)

    # Clean up: Remove individual LAS files and generated DEMs after processing
    for file_path in tqdm(las_files, desc="Cleaning up LAS files"):
        os.remove(file_path)
    for file_path in tqdm(dem_files, desc="Cleaning up DEM files"):
        os.remove(file_path)

    # Delete the directory used for the unzipped LAS files
    os.rmdir(unzipped_directory)

    print("Processing complete!")


# Check if this script is the main module being executed, and if so, run the main function
if __name__ == '__main__':
    main()
