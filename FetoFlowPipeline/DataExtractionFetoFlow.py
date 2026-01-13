# This is a class that analysis quantitatively the angiogenesis model
# you can use this class to:
# - read a file 
# - extract the information 

import numpy as np
import statistics as stats
import placentagen as pg
import pandas as pd

from scipy.spatial import ConvexHull

class DataExtractionFetoFlow:
    def read_csv_to_numpy(filename):
        # IA generated function
        # Read the CSV file into a pandas DataFrame
        df = pd.read_csv(filename)
        # Convert the DataFrame to a NumPy array
        numpy_array = df.to_numpy()
        return numpy_array

    def NodesDictionary(file_nodes, dim):
        # we open the file
        f = open(file_nodes, 'r')

        for line in f:
            pass
        last_line = line

        list_nodes = [float(x) for x in last_line.split()[1:]]

        NumberNodes = int((len(list_nodes))/dim)
        array_nodes = np.reshape(list_nodes, (NumberNodes,dim))

        # creation of the dictionary
        dict_nodes = {}
        for k in range(NumberNodes):
            if(dim == 2):
                array_nodes_dim = np.append(array_nodes[k], 0)
                dict_nodes[k] = array_nodes_dim
            else:
                dict_nodes[k] = array_nodes[k]

        # we close the file 
        f.close()

        return dict_nodes
    
    def NodesCoordinates(file_nodescoordinates, dim):
        # initialisation of the list
        list_nodescoordinates = []

        # we open the file
        f = open(file_nodescoordinates, 'r')

        for line in f:
            pass
        last_line = line

        list_nodescoordinates = [float(x) for x in last_line.split()[1:]]

        NumberNodes = int((len(list_nodescoordinates))/dim)
        array_nodescoordinates = np.reshape(list_nodescoordinates, (NumberNodes,dim))

        # we close the file 
        f.close()

        return array_nodescoordinates

    def ConnectivityArray(file_connectivity, dim):
        # we open the file
        f = open(file_connectivity, 'r')

        for line in f:
            pass
        last_line = line

        list_connectivity = [int(x) for x in last_line.split()[1:]]

        NumberElements = int((len(list_connectivity))/2)
        arraylist_connectivity = np.reshape(list_connectivity, (NumberElements,2))

        # creation of the list of tuples
        array_connectivity = []
        for k in range(1, NumberElements):
            tuple_connectivity = (int(arraylist_connectivity[k][1]), int(arraylist_connectivity[k][0]))
            array_connectivity.append(tuple_connectivity)

        # we close the file 
        f.close()

        return array_connectivity
    
    def read_viznodes_to_csv(file_nodes, main_file_path, dim):
        data_nodes = DataExtractionFetoFlow.NodesCoordinates(file_nodes, 2)
        output_path_nodes = main_file_path + "NodesCoordinatesAngiogenesisModel.csv"
        df = pd.DataFrame(data_nodes)
        df.to_csv(output_path_nodes, index=False)     
        
    def read_vizconnectivity_to_csv(file_connectivity, main_file_path, dim):
        data_connectivity = DataExtractionFetoFlow.ConnectivityArray(file_connectivity, 2)
        output_path_cellconnectivity = main_file_path + "CellConnectivityAngiogenesisModel.csv"
        df = pd.DataFrame(data_connectivity)
        df.to_csv(output_path_cellconnectivity, index=False)

