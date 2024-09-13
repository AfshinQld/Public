import laspy  # Used for reading .las files, a common format for LiDAR data.
import numpy as np  # Provides support for efficient numerical computations.
import os  # Allows for operating system dependent functionality like reading file paths.

# This script processes LiDAR data in .las format, calculating the density of points per square meter across each file in a
# specified directory. It iterates through each .las file in the input directory, creates a 2D histogram representing the distribution
# of points within 1 square meter cells across the spatial extent of the data, and outputs the average point density per square meter
# for each file. This is useful for understanding the distribution and density of LiDAR points across a surveyed area.



# Specify the input directory where your .las files are located
input_directory = r"C:\Temp\Raglan\Single output\Southern site\\"

def count_points_per_square_meter(las_file):
    """
    Counts the number of points per square meter in a given .las file.

    Parameters:
    las_file (str): The path to the .las file.

    Returns:
    np.ndarray: A 2D histogram array of points per square meter.
    """
    las = laspy.read(las_file)
    x_points = las.x
    y_points = las.y
    x_min, x_max = np.min(x_points), np.max(x_points)
    y_min, y_max = np.min(y_points), np.max(y_points)
    cell_size = 1  # Defines the size of the square meter cell

    # Calculate the number of cells needed to cover the area
    x_cells = int(np.ceil((x_max - x_min) / cell_size))
    y_cells = int(np.ceil((y_max - y_min) / cell_size))

    # Create a 2D histogram of points within each cell
    histogram, _, _ = np.histogram2d(x_points, y_points, bins=[x_cells, y_cells], range=[[x_min, x_max], [y_min, y_max]])
    return histogram

def process_directory(input_directory):
    """
    Processes all .las files in the specified directory.

    Parameters:
    input_directory (str): The directory containing .las files.
    """
    # Iterate over all files in the input directory
    for filename in os.listdir(input_directory):
        if filename.endswith(".las"):
            las_file_path = os.path.join(input_directory, filename)
            print(f"Processing {filename}...")
            points_per_square_meter = count_points_per_square_meter(las_file_path)
            # Print results for each file
            print(f"Histogram of points per square meter for {filename}:")
            print(points_per_square_meter)
            print(f"Average points per square meter for {filename}: {np.mean(points_per_square_meter)}\n")

# Process all .las files in the specified directory
process_directory(input_directory)
