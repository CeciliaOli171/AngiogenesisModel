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

    # from copilot
    def orient_downstream_by_bfs_no_reindex(
        edges_dir: Iterable[Tuple[Hashable, Hashable]],
        inlet: Hashable,
        *,
        weight: Optional[str] = None  # None: unweighted; else edge attribute name for Dijkstra on undirected view
    ) -> List[Tuple[Hashable, Hashable]]:
        """
        Force a downstream orientation from `inlet` without renumbering nodes.
        - Computes levels from the inlet on the UNDIRECTED view (so we don't get stuck by current directions).
        - Creates a strict total order: (distance, discovery_index, node_id).
        - Flips any edge that violates that order.

        Parameters
        ----------
        edges_dir : iterable[(u, v)]
            Your CURRENT directed connectivity (tuples).
        inlet : hashable
            The inlet node id.
        weight : str or None
            If you have edge weights (e.g., 'length'), compute distances by Dijkstra on the UNDIRECTED graph.

        Returns
        -------
        new_edges : list[(u, v)]
            Same node IDs, same connected pairs; only directions flipped to follow downstream.
        """
        # Build directed and undirected views
        DG = nx.DiGraph()
        DG.add_edges_from(edges_dir)
        if inlet not in DG:
            raise ValueError("Inlet node is not present in the graph.")

        UG = nx.Graph()
        UG.add_edges_from(DG.edges())

        # Compute distances from inlet on UNDIRECTED view (so we get a full outward notion)
        if weight is None:
            # Unweighted shortest path distances
            dist = nx.single_source_shortest_path_length(UG, inlet)
        else:
            # Weighted distances: construct a weighted undirected graph
            # Expect edges_dir to be (u, v) and weights stored elsewhere, or
            # pre-attach weights in a separate step before calling this function.
            # If you already have a Graph with weights, you can modify this section accordingly.
            # Fallback to treating edges as weight=1 if not present.
            # Here we assume weight attributes are on DG; copy to UG:
            for u, v in DG.edges():
                w = DG[u][v].get(weight, 1.0)
                UG.add_edge(u, v, **{weight: w})
            dist = nx.single_source_dijkstra_path_length(UG, inlet, weight=weight)

        # BFS discovery order on UNDIRECTED view to break ties
        discovery_index: Dict[Hashable, int] = {}
        idx = 0
        for n in nx.bfs_tree(UG, source=inlet):
            discovery_index[n] = idx
            idx += 1

        # For any nodes not reached by BFS (disconnected), we’ll place them after reachable ones
        for n in UG.nodes():
            if n not in discovery_index:
                discovery_index[n] = idx
                idx += 1

        # Build a strict rank key for every node
        # Use very large value for dist if unreachable from inlet
        def node_key(n):
            d = dist.get(n, float('inf'))
            di = discovery_index.get(n, float('inf'))
            return (d, di, n)

        # Flip edges that go "upstream" in this order
        new_edges: List[Tuple[Hashable, Hashable]] = []
        seen = set()
        for (u, v) in DG.edges():
            ku, kv = node_key(u), node_key(v)
            if ku <= kv:
                e = (u, v)
            else:
                e = (v, u)
            if e not in seen:
                # Optionally: if the reverse is already present, you can choose to drop it
                # seen.add((e[1], e[0]))   # uncomment to forbid both directions between same pair
                seen.add(e)
                new_edges.append(e)

        return new_edges

    def nodes_elements_calculation(main_pathway, InitialisationFiles, hpc, local, SeedNb, SourceNb, dim):
        if hpc:
            main_pathway_seedsource = main_pathway + str(SeedNb) + "Source" + str(SourceNb) + "/results_from_time_0/"
        elif local:
            main_pathway_seedsource = main_pathway + str(SeedNb) + "Source" + str(SourceNb) + "/"
        file_nodes = main_pathway_seedsource + "results.viznodes"
        file_connectivity = main_pathway_seedsource + "results.vizconnectivity"

        if InitialisationFiles:
            DataExtractionFetoFlow.read_viznodes_to_csv(file_nodes, main_pathway_seedsource, dim)
            DataExtractionFetoFlow.read_vizconnectivity_to_csv(file_connectivity, main_pathway_seedsource, dim)

        # Placentagen visualisation: check that the data is well organised for fetoflow
        if hpc: 
            nodes = DataExtractionFetoFlow.read_csv_to_numpy("/hpc/coli171/Results/PaperAngiogenesisModel2025/PaperModel2025Analysis2D/Figures/NodesCoordinatesAngiogenesisModel.csv")
        elif local:
            nodes = DataExtractionFetoFlow.read_csv_to_numpy(main_pathway_seedsource + "NodesCoordinatesAngiogenesisModel.csv")
        zeros_column = np.zeros((nodes.shape[0], 1))
        nodes_3d = np.hstack((nodes, zeros_column))
        indices = np.arange(len(nodes_3d)).reshape(-1, 1)
        indices = indices.astype(int)
        nodes_indexed = np.concatenate((indices, nodes_3d), axis=1)

        if hpc:
            elems = DataExtractionFetoFlow.read_csv_to_numpy("/hpc/coli171/Results/PaperAngiogenesisModel2025/PaperModel2025Analysis2D/Figures/CellConnectivityAngiogenesisModel.csv") # to modify for more than just one test
        elif local:
            elems = DataExtractionFetoFlow.read_csv_to_numpy(main_pathway_seedsource + "CellConnectivityAngiogenesisModel.csv")
        elems = elems.astype(int)
        elems = elems[1:]
        indices = np.arange(len(elems)).reshape(-1, 1)
        indices = indices.astype(int)
        elems = np.concatenate((indices, elems), axis=1)

        elems,branch_id,branch_start,branch_end,cycle_bool,seen_elements = pg.fix_elem_direction(nodes_indexed[0,1:4],elems,nodes_indexed)

        pg.export_ex_coords(nodes_indexed,'nodes', main_pathway_seedsource + 'nodes','exnode')
        pg.export_exelem_1d(elems, 'elems', main_pathway_seedsource + 'elems')
        pg.export_ip_coords(nodes_indexed,'nodes', main_pathway_seedsource + 'N_ipelem')
        pg.export_ipelem_1d(elems,'elems', main_pathway_seedsource + 'E_ipelem')

        df = pd.DataFrame(elems, columns=['Index', 'N1', 'N2'])
        df.to_csv(main_pathway_seedsource + "elems.csv", index=False, header=True)
        df = pd.DataFrame(nodes_indexed, columns=['Index', 'X', 'Y', 'Z'])
        df.to_csv(main_pathway_seedsource + "nodes.csv", index=False, header=True)

        # FetoFlow 
        nodes = DataExtractionFetoFlow.NodesDictionary(file_nodes, dim)

        # also not optimised at all but works for the moment
        elements = DataExtractionFetoFlow.ConnectivityArray(file_connectivity, dim)
        Graphnx = nx.DiGraph()
        Graphnx.add_edges_from(elements)

        # removing elements that appear twice
        to_remove = set()
        for u, v in Graphnx.edges():
            if Graphnx.has_edge(v, u):
                to_remove.add((u, v))
                to_remove.add((v, u))
        Graphnx.remove_edges_from(to_remove)
        elements = Graphnx.edges()

        edges = elements
        inlet = 0
        edges_downstream = DataExtractionFetoFlow.orient_downstream_by_bfs_no_reindex(edges, inlet, weight=None)
        elements_reordered = np.array(edges_downstream, dtype=object)

        return nodes, elements_reordered

    def flow_pressure_calculation(main_pathway, nodes, elements, inlet_pressure, outlet_pressure, umbilical_artery_radius, decay_factor, viscosity_type, SmallSystem, hpc, local, SeedNb, SourceNb, dim):
        if hpc:
            main_pathway_seedsource = main_pathway + str(SeedNb) + "Source" + str(SourceNb) + "/results_from_time_0/"
        elif local:
            main_pathway_seedsource = main_pathway + str(SeedNb) + "Source" + str(SourceNb) + "/"
        file_nodes = main_pathway_seedsource + "results.viznodes"
        
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
        flows = pd.DataFrame([{"Element" : element, "Flow" : flow} for element,flow in q.items()])
        flows.to_csv(main_pathway_seedsource + "simulation_flows.csv")

        nodes_ps = DataExtractionFetoFlow.NodesCoordinates(file_nodes,dim)
        edges_ps = np.array([element for element in elements])
        pressure_ps = np.array([pressure[1] for pressure in p.items()])
        flow_ps = np.array([flow[1]*3600e9 for flow in q.items()]) # we convert it into Chaste dimensions (mm^3.h)

        return nodes_ps, edges_ps, pressure_ps, flow_ps
    
    def find_vesseltips(nodes, elements):
        indices_vesseltips = []
        for k in range(len(nodes)):
            potential_vesseltip_indice = k
            n = 0
            for elem in elements:
                if(elem[0] == potential_vesseltip_indice):
                    n += 1
                if(elem[1] == potential_vesseltip_indice):
                    n += 1
            if(n == 1):
                indices_vesseltips.append(potential_vesseltip_indice)

        return indices_vesseltips
