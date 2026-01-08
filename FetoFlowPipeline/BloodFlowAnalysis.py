from fetoflow import *
import placentagen as pg
import pandas as pd 
import networkx as nx

import SimpleITK as sitk
import re, os
from collections import defaultdict

import polyscope as ps

import numpy as np
import matplotlib.pyplot as plt

from DataExtractionFetoFlow import DataExtractionFetoFlow
runner = DataExtractionFetoFlow

main_pathway = "/hpc/coli171/Results/PaperAngiogenesisModel2025/PaperModel2025Analysis2D/" + "CoupledModel2DAnalyticalApproxPde/CoupledModel2DAnalyticalApproxPdeSeed" + str(20) + "Source" + str(5) + "/results_from_time_0/"
#main_pathway = "/Users/coli171/Desktop/test/"
file_nodes = main_pathway + "results.viznodes"
file_connectivity = main_pathway + "results.vizconnectivity"
dim = 2


# Placentagen visualisation: check that the data is well organised for fetoflow
nodes = runner.read_csv_to_numpy("/hpc/coli171/Results/PaperAngiogenesisModel2025/PaperModel2025Analysis2D/Figures/NodesCoordinatesAngiogenesisModel.csv")
#nodes = runner.read_csv_to_numpy(main_pathway + "NodesCoordinatesAngiogenesisModel.csv")
zeros_column = np.zeros((nodes.shape[0], 1))
nodes_3d = np.hstack((nodes, zeros_column))
indices = np.arange(len(nodes_3d)).reshape(-1, 1)
indices = indices.astype(int)
nodes_indexed = np.concatenate((indices, nodes_3d), axis=1)

elems = runner.read_csv_to_numpy("/hpc/coli171/Results/PaperAngiogenesisModel2025/PaperModel2025Analysis2D/Figures/CellConnectivityAngiogenesisModel.csv")
#elems = runner.read_csv_to_numpy(main_pathway + "elems.csv")
elems = elems.astype(int)
elems = elems[1:]
indices = np.arange(len(elems)).reshape(-1, 1)
indices = indices.astype(int)
elems = np.concatenate((indices, elems), axis=1)

E_conn = pg.element_connectivity_1D(nodes_3d, elems)
terminal = pg.calc_terminal_branch(nodes_3d,elems)
elems,branch_id,branch_start,branch_end,cycle_bool,seen_elements = pg.fix_elem_direction(nodes_indexed[0,1:4],elems,nodes_indexed)
terminal = pg.calc_terminal_branch(nodes_3d,elems)

pg.export_ex_coords(nodes_indexed,'nodes', main_pathway + 'nodes','exnode')
pg.export_exelem_1d(elems, 'elems', main_pathway + 'elems')
pg.export_ip_coords(nodes_indexed,'nodes', main_pathway + 'N_ipelem')
pg.export_ipelem_1d(elems,'elems', main_pathway + 'E_ipelem')

df = pd.DataFrame(elems, columns=['Index', 'N1', 'N2'])
df.to_csv(main_pathway + "elems.csv", index=False, header=True)
df = pd.DataFrame(nodes_indexed, columns=['Index', 'X', 'Y', 'Z'])
df.to_csv(main_pathway + "nodes.csv", index=False, header=True)


# FetoFlow 
nodes = runner.NodesDictionary(file_nodes, dim)
elements = runner.ConnectivityArray(file_connectivity, dim)

inlet_pressure, outlet_pressure = 6650, 2660
bcs = generate_boundary_conditions(inlet_pressure = inlet_pressure, outlet_pressure = outlet_pressure, inlet_flow=None)

umbilical_artery_radius, decay_factor = 1.8 / 1000, 1.38 
umbilical_vein_radius, decay_factor_vein = 4.0 / 1000, 1.46
arteries_only = False 
viscosity_type = 'constant'

G = create_geometry(nodes, elements, umbilical_artery_radius, decay_factor, arteries_only, umbilical_vein_radius, decay_factor_vein)

G = calculate_resistance(G, viscosity_model=viscosity_type)

A,b,bc_export = create_small_matrices(G,bcs,branching_angles=False)
p,q = solve_small_system(A,b,G,bc_export)

pressures = pd.DataFrame([{"Node" : node, "Pressure" : pressure} for node,pressure in p.items()])
pressures.to_csv(main_pathway + "simulation_pressures.csv")
flows = pd.DataFrame([{"Element" : element, "Flow" : flow} for element,flow in q.items()])
flows.to_csv(main_pathway + "simulation_flows.csv")

nodes_ps = np.array([node for node in p.items()])
print(len(nodes_ps))
edges_ps = np.array([element for element in elements])
print(len(edges_ps))
pressure_ps = np.array([pressure[1] for pressure in p.items()])
print(len(pressure_ps))
flow_ps = np.array([flow[1] for flow in q.items()])
print(len(flow_ps))

# Polyscope 
ps.init("openGL3_egl")

# represent vascular tree 
tree = ps.register_curve_network("vascular tree", nodes_ps, edges_ps, color=[155/255,155/255,155/255])

# represent flow and pressure
# tree.add_scalar_quantity("flow", flow_ps, cmap='reds', enabled=True)
# tree.add_scalar_quantity("pressure (at node)", pressure_ps, cmap='jet', enabled=True)

# presentation settings
ps.set_ground_plane_mode("none")
ps.set_navigation_style("free")
ps.set_up_dir("y_up")
ps.set_front_dir("z_front")
ps.set_background_color([0,0,0])
ps.show()