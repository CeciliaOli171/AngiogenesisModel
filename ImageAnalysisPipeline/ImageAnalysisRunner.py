# python code for image analysis
# • reading QuPath file -> extract nodes coordinates
# • using coordinates to obtain surface area -> convex hull method 

# Libraries 
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import geojson 

from matplotlib import cm
from matplotlib.ticker import LinearLocator 
from matplotlib.colors import LinearSegmentedColormap

from scipy.spatial import ConvexHull

# Extract the data 
def NodesCoordinates(file_nodes):
    # initialisation of the list
    list_nodes = []

    # we open the file
    f = open(file_nodes, 'r')
    f.readline()

    for line in f:
        for x in line.split()[0:]:
            list_nodes.append(float(x))

    NumberNodes = int((len(list_nodes))/2)
    array_nodescoordinates = np.reshape(list_nodes, (NumberNodes,2))

    # we close the file 
    f.close()

    return array_nodescoordinates

def XNodesCoordinates(file_nodes, Normalise):
    x_nodescoordinates = []
    array_nodescoordinates = NodesCoordinates(file_nodes)
    TotalNbCells = array_nodescoordinates.shape[0]

    for k in range(TotalNbCells):
        if(Normalise == 1):
            norm_node = np.sqrt(array_nodescoordinates[k][0]**2 + array_nodescoordinates[k][1]**2)
        else:
            norm_node = 1  
        x_nodescoordinates.append(array_nodescoordinates[k][0]/norm_node)

    return x_nodescoordinates

def YNodesCoordinates(file_nodes, Normalise):
    y_nodescoordinates = []
    array_nodescoordinates = NodesCoordinates(file_nodes)
    TotalNbCells = array_nodescoordinates.shape[0]

    for k in range(TotalNbCells):
        if(Normalise == 1):
            norm_node = np.sqrt(array_nodescoordinates[k][0]**2 + array_nodescoordinates[k][1]**2)
        else:
            norm_node = 1 
        y_nodescoordinates.append(array_nodescoordinates[k][1]/norm_node)

    return y_nodescoordinates

def OriginNode(file_nodes):
    array_nodescoordinates = NodesCoordinates(file_nodes)
    TotalNbCells = array_nodescoordinates.shape[0]

    # initialisation of the origin point
    origin_node = array_nodescoordinates[0]

    for k in range(1, TotalNbCells):
        if(array_nodescoordinates[k][0] < origin_node[0]):
            origin_node[0] = array_nodescoordinates[k][0]
            origin_node[1] = array_nodescoordinates[k][1]
        if(array_nodescoordinates[k][0] == origin_node[0]):
            if(array_nodescoordinates[k][1] < origin_node[1]):
                origin_node[0] = array_nodescoordinates[k][0]
                origin_node[1] = array_nodescoordinates[k][1]

    return origin_node

def HarmonisingCoordinates(file_nodes):
    array_nodescoordinates = NodesCoordinates(file_nodes)
    TotalNbCells = array_nodescoordinates.shape[0]

    # initialisation
    array_harmonisednodescoordinates = []
    x_harmonisednodescoordinates = []
    y_harmonisednodescoordinates = []
    origin_node = OriginNode(file_nodes)

    for k in range(TotalNbCells):
        x_harmonised = array_nodescoordinates[k][0] - origin_node[0]
        y_harmonised = array_nodescoordinates[k][1] - origin_node[1]

        array_harmonisednodescoordinates.append([x_harmonised, y_harmonised])
        x_harmonisednodescoordinates.append(x_harmonised)
        y_harmonisednodescoordinates.append(y_harmonised)

    return array_harmonisednodescoordinates, x_harmonisednodescoordinates, y_harmonisednodescoordinates

def HarmonisingCoordinatesBetweenFiles(file_nodes_1, file_nodes_2):
    array_nodescoordinates_1 = NodesCoordinates(file_nodes_1)
    array_nodescoordinates_2 = NodesCoordinates(file_nodes_2)
    TotalNbCells_1 = array_nodescoordinates_1.shape[0]
    TotalNbCells_2 = array_nodescoordinates_2.shape[0]

    # initialisation
    x_harmonisednodescoordinates_1 = []
    y_harmonisednodescoordinates_1 = []
    origin_node_1 = OriginNode(file_nodes_1)
    x_harmonisednodescoordinates_2 = []
    y_harmonisednodescoordinates_2 = []
    origin_node_2 = OriginNode(file_nodes_2)

    for k in range(TotalNbCells_1):
        x_harmonised_1 = array_nodescoordinates_1[k][0] + np.abs(origin_node_1[0]-origin_node_2[0])
        y_harmonised_1 = array_nodescoordinates_1[k][1] + np.abs(origin_node_1[1]-origin_node_2[1])

        x_harmonisednodescoordinates_1.append(x_harmonised_1)
        y_harmonisednodescoordinates_1.append(y_harmonised_1)
        x_harmonisednodescoordinates_2.append(array_nodescoordinates_2[k][0])
        y_harmonisednodescoordinates_2.append(array_nodescoordinates_2[k][1])

    return x_harmonisednodescoordinates_1, y_harmonisednodescoordinates_1, x_harmonisednodescoordinates_2, y_harmonisednodescoordinates_2

# Quantitative Analysis 
# Number of Cells 
def NbCells(file_nodes):
    array_nodescoordinates = NodesCoordinates(file_nodes)
    TotalNbCells = array_nodescoordinates.shape[0]

    return TotalNbCells

# Study of the surface area : only for samples : not representative of the whole lesion 
def Area(file_nodes):
    # Method with the convex hull 
    array_nodescoordinates = NodesCoordinates(file_nodes)
    hull = ConvexHull(array_nodescoordinates)
    area = hull.area
    return area 

# Study of the density 
def DensityNodes(file_nodes, area_sample):
    area = Area(file_nodes)
    density = area/area_sample
    return density

def DensityCoordinates(coordinates):
    # area of the convex hull 
    hull = ConvexHull(coordinates)
    area = hull.area
    density = area
    return density

# Study of the center of gravity 
def CenterOfGravity(coordinates):
    # initialisation 
    centerofgravity = np.zeros(2)
    totalweight = 0

    NbPoints = len(coordinates)

    # total weight of annotation
    for k in range(NbPoints):
        xk_norm = np.sqrt(coordinates[k][0]**2 + coordinates[k][1]**2)
        totalweight += xk_norm

    # center of gravity 
    for k in range(NbPoints):
        xk_norm = np.sqrt(coordinates[k][0]**2 + coordinates[k][1]**2)
        centerofgravity[0] += (xk_norm/totalweight)*coordinates[k][0]
        centerofgravity[1] += (xk_norm/totalweight)*coordinates[k][1]

    return centerofgravity