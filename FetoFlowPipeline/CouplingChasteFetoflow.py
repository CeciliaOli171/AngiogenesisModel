# This is a class that couples FetoFlow python library to Chaste c++ library 

import numpy as np
import statistics as stats
import SimpleITK as sitk
import re, os

from fetoflow import *

from collections import defaultdict
from typing import Iterable, Hashable, Tuple, List, Dict, Set, Optional
from scipy.spatial import ConvexHull


class CouplingChasteFetoflow:

    # this one should be inside Chaste 
    def ChasteDataForFetoflow(data_chaste):
        # converts c_vector to list
        return data_fetoflow

    # this should also be inside Chaste 
    def FetoflowDataForChaste(data_fetoflow):
        # converts list to c_vector 
        return data_chaste

    def NodesDictionary(nodes_chaste, dim):
        NumberNodes = int((len(nodes_chaste))/dim)
        array_nodes = np.reshape(nodes_chaste, (NumberNodes,dim))

        # creation of the dictionary
        dict_nodes = {}
        for k in range(NumberNodes):
            if(dim == 2):
                array_nodes_dim = np.append(array_nodes[k], 0)
                dict_nodes[k] = array_nodes_dim
            else:
                dict_nodes[k] = array_indexed_nodes[k]

        return dict_nodes
    
    def ConnectivityArray(connectivity_chaste, dim):
        NumberElements = int((len(connectivity_chaste))/2)
        arrayconnectivity = np.reshape(connectivity_chaste, (NumberElements,2))

        # creation of the list of tuples
        array_connectivity = []
        for k in range(1, NumberElements):
            tuple_connectivity = (int(arraylist_connectivity[k][1]), int(arraylist_connectivity[k][0]))
            array_connectivity.append(tuple_connectivity)

        return array_connectivity

    def CalculatesNodesEdgesForFetoflow(nodes_chaste, connectivity_chaste, dim):
        nodes = CouplingChasteFetoflow.NodesDictionary(nodes_chaste, dim)
        elements = CouplingChasteFetoflow.ConnectivityArray(connectivity_chaste, dim)
        return nodes, elements

    def CalculatesFlowsPressures(nodes, elements, dim, inlet_pressure, outlet_pressure, artery_radius, decay_factor, viscosity_type):
        bcs = generate_boundary_conditions(inlet_pressure = inlet_pressure, outlet_pressure = outlet_pressure, inlet_flow=None)
        G = create_geometry(nodes, elements, artery_radius, decay_factor, True, 5e-6, 1.46) 
        G = calculate_resistance(G, viscosity_model=viscosity_type)

        if SmallSystem:
            A,b,bc_export = create_small_matrices(G,bcs,branching_angles=False)
            p,q = solve_small_system(A,b,G,bc_export)
        else:
            A,b = create_matrices(G, bcs=bcs)
            p,q = solve_system(A, b, num_nodes=G.number_of_nodes(), num_edges=G.number_of_edges())

        pressures = [pressure[1] for pressure in p.items()]
        flows = [flow[1]*3600e9 for flow in q.items()] # we convert it into Chaste dimensions (mm^3.h)

        return pressures, flows

    def CalculatesWSS(flows, viscosity, artery_radius):
        wss = []
        for k in range(len(flows)):
            wss.append((4*viscosity*flows[k])/(np.pi*artery_radius**3))

        return wss

    def CouplingFunction(nodes_chaste, elements_chaste, dim, inlet_pressure, outlet_pressure, artery_radius, decay_factor, viscosity_type, viscosity):
        # we start by converting the nodes and elements to the correct form for Fetoflow
        nodes, edges = CouplingChasteFetoflow.CalculatesNodesEdgesForFetoflow(nodes_chaste, connectivity_chaste, dim)

        # we calculate the pressure and flow 
        pressures, flows = CouplingChasteFetoflow.CalculatesFlowsPressures(nodes, elements, dim, inlet_pressure, outlet_pressure, artery_radius, decay_factor, viscosity_type)

        # we calculate the wall shear stress 
        wss = CouplingChasteFetoflow.CalculatesWSS(flows, viscosity, artery_radius)

        # convert the results for Chaste 
        pressures_chaste = pressures
        flows_chaste = flows
        wss_chaste = wss

        return pressures_chaste, flows_chaste, wss_chaste