# This is a class that analysis quantitatively the angiogenesis model
# you can use this class to:
# - read a file 
# - extract the information 

import numpy as np
import statistics as stats
import placentagen as pg
import pandas as pd
import networkx as nx
import SimpleITK as sitk
import re, os

from fetoflow import *
import placentagen as pg

from collections import defaultdict
from typing import Iterable, Hashable, Tuple, List, Dict, Set, Optional
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

    def nodes_elements_calculation(main_pathway, InitialisationFiles, hpc, local, SeedNb, SourceNb, dim):
        if hpc:
            main_pathway_seedsource = main_pathway + str(SeedNb) + "Source" + str(SourceNb) + "/results_from_time_0/"
        elif local:
            main_pathway_seedsource = main_pathway + str(SeedNb) + "Source" + str(SourceNb) + "/"
        file_nodes = main_pathway_seedsource + "results.viznodes"
        file_connectivity = main_pathway_seedsource + "results.vizconnectivity"
        file_branches_number = main_pathway_seedsource + "results.vizbranchnumber"

        # FetoFlow 
        nodes = DataExtractionFetoFlow.NodesDictionary(file_nodes, dim)
        elements = DataExtractionFetoFlow.ConnectivityArray(file_connectivity, dim)

        # removing elements that appear twice
        Graphnx = nx.DiGraph()
        Graphnx.add_edges_from(elements)
        to_remove = set()
        for u, v in Graphnx.edges():
            if Graphnx.has_edge(v, u):
                to_remove.add((u, v))
                to_remove.add((v, u))
        Graphnx.remove_edges_from(to_remove)
        elements_reordered = Graphnx.edges()

        # re-order elements no necessary anymore because edges oriented directly inside Chaste
        # elems_reordered = DataExtractionFetoFlow.get_branches(nodes, elems, file_branches_number)
        # elements_reordered = DataExtractionFetoFlow.EdgesOrientation(elems_reordered)

        return nodes, elements_reordered

    def flow_pressure_calculation(main_pathway, nodes, elements, inlet_pressure, outlet_pressure, umbilical_artery_radius, decay_factor, viscosity_type, SmallSystem, hpc, local, SeedNb, SourceNb, dim):
        # function calculating the flow and pressure in a vascular network using FetoFlow
        if hpc:
            main_pathway_seedsource = main_pathway + str(SeedNb) + "Source" + str(SourceNb) + "/results_from_time_0/"
        elif local:
            main_pathway_seedsource = main_pathway + str(SeedNb) + "Source" + str(SourceNb) + "/"
        file_nodes = main_pathway_seedsource + "results.viznodes"
        
        # initialise geometry
        bcs = generate_boundary_conditions(inlet_pressure = inlet_pressure, outlet_pressure = outlet_pressure, inlet_flow=None)
        G = create_geometry(nodes, elements, umbilical_artery_radius, decay_factor, True, 5e-6, 1.46) 
        G = calculate_resistance(G, viscosity_model=viscosity_type)

        if SmallSystem:
            A,b,bc_export = create_small_matrices(G,bcs,branching_angles=False)
            p,q = solve_small_system(A,b,G,bc_export)
        else:
            A,b = create_matrices(G, bcs=bcs)
            p,q = solve_system(A, b, num_nodes=G.number_of_nodes(), num_edges=G.number_of_edges())

        pressures = pd.DataFrame([{"Node" : node, "Pressure" : pressure} for node,pressure in p.items()])
        pressures.to_csv(main_pathway_seedsource + "simulation_pressures.csv")
        flows = pd.DataFrame([{"Element" : element, "Flow" : flow*3600e9} for element,flow in q.items()])
        flows.to_csv(main_pathway_seedsource + "simulation_flows.csv")

        # extracts nodes coordinates, edges, pressure and flow in arrays 
        nodes_ps = DataExtractionFetoFlow.NodesCoordinates(file_nodes, dim)
        edges_ps = np.array([elem for elem in elements])
        pressure_ps = np.array([pressure[1] for pressure in p.items()])
        flow_ps = np.array([flow[1]*3600e9 for flow in q.items()]) # we convert it into Chaste dimensions (mm^3.h)

        return nodes_ps, edges_ps, pressure_ps, flow_ps

    def find_vesseltips(nodes, elements):
         # function selecting all the indices of the outlet edges for a graph defined by the nodes and elements lists 
        indices_vesseltips = []
        for k in range(1,len(nodes)):
            potential_vesseltip_indice = k # it is the indice of the node not the edge
            n = 0
            edge_indice = []
            m = 0
            for elem in elements:
                if(elem[0] == potential_vesseltip_indice):
                    n += 1
                    edge_indice.append(m)
                if(elem[1] == potential_vesseltip_indice):
                    n += 1
                    edge_indice.append(m)
                m +=1
            if(n == 1):
                indices_vesseltips.append(edge_indice[0])

        return indices_vesseltips
    
    def find_vesseltips_insidelesion(nodes, elements, ref_point):
        # function selecting all the indices of the outlet edges inside the lesion (located at x < ref_point) for a graph defined by the nodes and elements lists 
        indices_vesseltips = []
        for k in range(1,len(nodes)):
            potential_vesseltip_indice = k # it is the indice of the node not the edge
            n = 0
            edge_indice = []
            m = 0
            for elem in elements:
                if(elem[0] == potential_vesseltip_indice):
                    n += 1
                    edge_indice.append(m)
                if(elem[1] == potential_vesseltip_indice):
                    n += 1
                    edge_indice.append(m)
                m +=1
            if(n == 1):
                indices_vesseltips.append(k)
        
        indices_vesseltips_insidelesion = []
        for i in indices_vesseltips:
            x = nodes[i]
            if x[0] < ref_point:
                indices_vesseltips_insidelesion.append(i)

        return indices_vesseltips_insidelesion
    
    def find_vesseltips_outsidelesion(nodes, elements, ref_point):
        # function selecting all the indices of the outlet edges inside the lesion (located at x < ref_point) for a graph defined by the nodes and elements lists 
        indices_vesseltips = []
        for k in range(1,len(nodes)):
            potential_vesseltip_indice = k # it is the indice of the node not the edge
            n = 0
            edge_indice = []
            m = 0
            for elem in elements:
                if(elem[0] == potential_vesseltip_indice):
                    n += 1
                    edge_indice.append(m)
                if(elem[1] == potential_vesseltip_indice):
                    n += 1
                    edge_indice.append(m)
                m +=1
            if(n == 1):
                indices_vesseltips.append(k)
        
        indices_vesseltips_insidelesion = []
        for i in indices_vesseltips:
            x = nodes[i]
            if x[0] > ref_point:
                indices_vesseltips_insidelesion.append(i)

        return indices_vesseltips_insidelesion

    def BranchNumber(file_branchesnumber):
        # initialisation of the list
        list_branchesnumber = []

        # we open the file
        f = open(file_branchesnumber, 'r')

        for line in f:
            pass
        last_line = line

        list_branchesnumber = np.array([int(x) for x in last_line.split()[1:]])

        # we close the file 
        f.close()

        return list_branchesnumber

    def BranchesList(file_branchesnumber):
        branches_list = []
        branchnumbers_list = DataExtractionFetoFlow.BranchNumber(file_branchesnumber)

        maxbranchnumber = int(max(branchnumbers_list))
        for k in range(maxbranchnumber):
            branch = np.where(branchnumbers_list == k)
            # we add it to the single branch list 
            # we add the branch to the whole branches list
            branches_list.append(branch)

        return branches_list

    def get_branches(nodes, elems, file_branches_number):
        branches_list_connectivity = []
        branches_list = DataExtractionFetoFlow.BranchesList(file_branches_number)

        # for all nodes from the list of branches, we extract the connectivity list 
        for branch in branches_list:
            branch_connectivity = []
            for elem in elems:
                if(np.isin(elem[0],branch) and np.isin(elem[1],branch)):
                    branch_connectivity.append(elem)
            branches_list_connectivity.append(branch_connectivity)

        # we extract all the vessel tips of the network
        vesseltips_list = []
        for n in range(len(elems)):
            if(sum(n in s for s in elems) == 1):
                vesseltips_list.append(n)

        # we extract all the branching vessels of the network
        branchingvessels_list = []
        for n in range(len(elems)):
            if(sum(n in s for s in elems) > 2):
                branchingvessels_list.append(n)

        # we loop over the branches
        branchesconnectivity_reordered_list = []
        BranchTotalNb = len(branches_list)
        for k in range(BranchTotalNb):
            branchconnectivity_reordered = []
            branches_list_final, branches_list_connectivity_final, branchconnectivity_reordered_final, vesseltips_list = DataExtractionFetoFlow.BranchesConnectivity(branches_list[k], branches_list_connectivity[k], branchconnectivity_reordered, vesseltips_list, branchingvessels_list)
            branchesconnectivity_reordered_list.append(branchconnectivity_reordered_final)
        
        return branchesconnectivity_reordered_list

    def BranchesConnectivity(branch_list, branchconnectivity_list, branchconnectivity_reordered, vesseltips_list, branchingvessels_list):
        branch_list_final = branch_list
        branchconnectivity_list_final = branchconnectivity_list

        if not branch_list or not branchconnectivity_list:
            return branch_list_final, branchconnectivity_list_final, branchconnectivity_reordered, vesseltips_list, branchingvessels_list
        else:
            # case 0: branch with the INLET tip cell 
            if(np.isin(0, branch_list)):
                # we orient the other from this tip cell until the next branching cell:
                # we search for the edge containing the cell, we flip it if necessary 
                # if there is another branching cell, we repeat the process 
                cell_considered = 0
                branchconnectivity_reordered_part = []
                while(sum(cell_considered in s for s in branchconnectivity_list) == 2):
                    for elem in branchconnectivity_list:
                        if(elem[0] == cell_considered):
                            branchconnectivity_reordered_part.append((cell_considered, elem[0]))
                            branch_list.remove(cell_considered)
                            cell_considered = elem[0]
                            branchconnectivity_list.remove(elem)
                            break 
                        elif(elem[1] == cell_considered):
                            branchconnectivity_reordered_part.append((cell_considered, elem[1]))
                            branch_list.remove(cell_considered)
                            cell_considered = elem[1]
                            branchconnectivity_list.remove(elem)
                            break
                branchconnectivity_reordered.append(branchconnectivity_reordered_part)
                print(branchconnectivity_reordered)
                print(branch_list)
                print(branchconnectivity_list)
                
            # case 1: branch with an OUTLET tip cell 
            elif(any(x in vesseltips_list for x in branch_list)):
                # we orient the other towards the tip cell 
                cell_considered = next((x for x in vesseltips_list if x in branch_list), None)
                branchconnectivity_reordered_part = []
                while(sum(cell_considered in s for s in branchconnectivity_list) == 2):
                    for elem in branchconnectivity_list:
                        if(elem[0] == cell_considered):
                            branchconnectivity_reordered_part.append((elem[0], cell_considered))
                            branch_list.remove(cell_considered)
                            cell_considered = elem[0]
                            branchconnectivity_list.remove(elem)
                            break 
                        elif(elem[1] == cell_considered):
                            branchconnectivity_reordered_part.append((elem[1], cell_considered))
                            branch_list.remove(cell_considered)
                            cell_considered = elem[1]
                            branchconnectivity_list.remove(elem)
                            break
                branchconnectivity_reordered.append(branchconnectivity_reordered_part)

            # case 2: branch between two branching vessels 
            else:
                # we orient the other towards the next branching segment 
                cell_considered = next((x for x in branchingvessels_list if x in branch_list), None)
                branchconnectivity_reordered_part = []
                while(sum(cell_considered in s for s in branchconnectivity_list) == 2):
                    for elem in branchconnectivity_list:
                        if(elem[0] == cell_considered):
                            branchconnectivity_reordered_part.append((cell_considered, elem[0]))
                            branch_list.remove(cell_considered)
                            cell_considered = elem[0]
                            branchconnectivity_list.remove(elem)
                            break 
                        elif(elem[1] == cell_considered):
                            branchconnectivity_reordered_part.append((cell_considered, elem[1]))
                            branch_list.remove(cell_considered)
                            cell_considered = elem[1]
                            branchconnectivity_list.remove(elem)
                            break
                branchconnectivity_reordered.append(branchconnectivity_reordered_part)
        
        branch_list_final, branchconnectivity_list_final, branchconnectivity_reordered, vesseltips_list, branchingvessels_list = DataExtractionFetoFlow.BranchesConnectivity(branch_list, branchconnectivity_list, branchconnectivity_reordered, vesseltips_list, branchingvessels_list)

        return branch_list_final, branchconnectivity_list_final, branchconnectivity_reordered, vesseltips_list, branchingvessels_list


    def EdgesOrientation(elems):
        # if node connected to only one other node = edge 
        # if node connected to 2 other nodes = inside branch, must be the same direction as the others 
        # if node connected to 3 other nodes = branching node, probably different directions 
        elems_reordered = elems

        nodesinsidebranch = []
        for n in range(len(elems)):
            if(sum(n in s for s in elems) == 2):
                nodesinsidebranch.append(n) # we initialise the list of branching nodes indices with the inlet and outlet nodes

        for node in nodesinsidebranch:
            edgesinsidebranch0 = []
            edgesinsidebranch1 = []
            for e in elems_reordered:
                if(e[0] == node):
                    edgesinsidebranch0.append(e)
                if(e[1] == node):
                    edgesinsidebranch1.append(e)
            if(len(edgesinsidebranch0) == 0):
                neighbours = [n[0] for n in edgesinsidebranch1]
                for i in neighbours:
                    if(sum(i in s for s in elems) == 2):
                        elems_reordered.remove((i, node))
                        elems_reordered.append((node,i))
                        break
            if(len(edgesinsidebranch1) == 0):
                neighbours = [n[1] for n in edgesinsidebranch0]
                for i in neighbours:
                    if(sum(i in s for s in elems) == 2):
                        elems_reordered.remove((node, i))
                        elems_reordered.append((i,node))
                        break

        return elems_reordered

    # from copilto to represent thebox plots
    def style_boxplot(bp, edge_color):
        # Make interiors transparent and edges colored
        for box in bp["boxes"]:
            box.set(facecolor="none", edgecolor=edge_color, linewidth=1.8)
        for whisker in bp["whiskers"]:
            whisker.set(color=edge_color, linewidth=1.4)
        for cap in bp["caps"]:
            cap.set(color=edge_color, linewidth=1.4)
        for median in bp["medians"]:
            median.set(color=edge_color, linewidth=2.0)
        # Fliers (outliers)
        for flier in bp["fliers"]:
            flier.set(marker="o", markerfacecolor="none", markeredgecolor=edge_color, alpha=0.9)