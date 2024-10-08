# Python code comparing:
#   • the blood vessel manual extraction using QuPath
#   • the blood vessel extraction using Matlab code 

# Files 
file_qupath = '/Users/coli171/Library/CloudStorage/OneDrive-TheUniversityofAuckland/Imaging/Test Matlab/Case7-DCISM.svs-points.tsv'
file_matlab = '/Users/coli171/Library/CloudStorage/OneDrive-TheUniversityofAuckland/Imaging/Test Matlab/Case7_centroidscoordinates.tsv'

# Libraries 
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

from matplotlib import cm
from matplotlib.ticker import LinearLocator 
from matplotlib.colors import LinearSegmentedColormap

from scipy.spatial import ConvexHull

from ImageAnalysisRunner import *

# Parameters
area_sample = 67000 # to verify
Normalise = 1 # normalisation of the coordinates 

# Plots

# Comparison between the Matlab detection and the manual detection  
# x_coordinates_qupath = XNodesCoordinates(file_qupath, Normalise)
# x_coordinates_matlab = XNodesCoordinates(file_matlab, Normalise)
# y_coordinates_qupath = YNodesCoordinates(file_qupath, Normalise)
# y_coordinates_matlab = YNodesCoordinates(file_matlab, Normalise)

# plt.scatter(x_coordinates_qupath, y_coordinates_qupath, c = 'b')
# plt.scatter(x_coordinates_matlab, y_coordinates_matlab, c = 'g', linestyle = '--')


array_harmonisedcoordinates_qupath, x_harmonisedcoordinates_qupath, y_harmonisedcoordinates_qupath = HarmonisingCoordinates(file_qupath)
array_harmonisedcoordinates_matlab, x_harmonisedcoordinates_matlab, y_harmonisedcoordinates_matlab = HarmonisingCoordinates(file_matlab)

# x_harmonisedcoordinates_matlab, y_harmonisedcoordinates_matlab, x_harmonisedcoordinates_qupath, y_harmonisedcoordinates_qupath = HarmonisingCoordinatesBetweenFiles(file_matlab, file_qupath)

plt.scatter(x_harmonisedcoordinates_qupath, y_harmonisedcoordinates_qupath, c = 'b')
plt.scatter(x_harmonisedcoordinates_matlab, y_harmonisedcoordinates_matlab, c = 'g')

plt.xlabel('x')
plt.ylabel('y')
plt.show()