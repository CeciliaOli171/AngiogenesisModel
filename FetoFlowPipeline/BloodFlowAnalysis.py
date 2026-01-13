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

# Settings
hpc = False
local = True

InitialisationFiles = True

SeedNb = 20
SourceNb = 10

if hpc:
    # example with analytical approximation of PDE of VEGF
    main_pathway = "/hpc/coli171/Results/PaperAngiogenesisModel2025/PaperModel2025Analysis2D/" + "CoupledModel2DAnalyticalApproxPde/CoupledModel2DAnalyticalApproxPdeSeed" + str(SeedNb) + "Source" + str(SourceNb) + "/results_from_time_0/"
    # example with constant VEGF
    # main_pathway = "/hpc/coli171/Results/PaperAngiogenesisModel2025/PaperModel2025Analysis2D/" + "CoupledModel2DConstant/CoupledModel2DConstantSeed" + str(20) + "Source" + str(5) + "/results_from_time_0/"
elif local:
    # example with analytical approximation of PDE of VEGF
    #main_pathway = "/Users/coli171/Desktop/test/AnalyticalApproxPde/Seed" + str(SeedNb) + "Source" + str(SourceNb) + "/"
    # example with constant VEGF
    main_pathway = "/Users/coli171/Desktop/test/Constant/"
file_nodes = main_pathway + "results.viznodes"
file_connectivity = main_pathway + "results.vizconnectivity"
dim = 2

if InitialisationFiles:
    runner.read_viznodes_to_csv(file_nodes, main_pathway, dim)
    runner.read_vizconnectivity_to_csv(file_connectivity, main_pathway, dim)


# Placentagen visualisation: check that the data is well organised for fetoflow
if hpc: 
    nodes = runner.read_csv_to_numpy("/hpc/coli171/Results/PaperAngiogenesisModel2025/PaperModel2025Analysis2D/Figures/NodesCoordinatesAngiogenesisModel.csv")
elif local:
    nodes = runner.read_csv_to_numpy(main_pathway + "NodesCoordinatesAngiogenesisModel.csv")
zeros_column = np.zeros((nodes.shape[0], 1))
nodes_3d = np.hstack((nodes, zeros_column))
indices = np.arange(len(nodes_3d)).reshape(-1, 1)
indices = indices.astype(int)
nodes_indexed = np.concatenate((indices, nodes_3d), axis=1)

if hpc:
    elems = runner.read_csv_to_numpy("/hpc/coli171/Results/PaperAngiogenesisModel2025/PaperModel2025Analysis2D/Figures/CellConnectivityAngiogenesisModel.csv") # to modify for more than just one test
elif local:
    elems = runner.read_csv_to_numpy(main_pathway + "CellConnectivityAngiogenesisModel.csv")
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

# at the inlet, we have an arterioles (pressure between 30 and 60 mmHg) and at the outlet, we have capillaries (average of 20 mmHg) -> Al-Nuaimi et al.
# pressure in Pa, 50 mmHg = 6666.1 Pa and 20 mmHg = 2666.45 Pa
# pressure difference between inlet and outlet of 30 mmHg -> Bazmara et al. but pressure difference of 9mmHg -> Chaplain et al.
inlet_pressure, outlet_pressure = 6666.1, 1333.22 
bcs = generate_boundary_conditions(inlet_pressure = inlet_pressure, outlet_pressure = outlet_pressure, inlet_flow=None)

umbilical_artery_radius, decay_factor = 5e-6, 1.38 # capillaries have a radius of approx. 5 micrometers
umbilical_vein_radius, decay_factor_vein = 5e-6, 1.46
arteries_only = True # very important otherwise creation of veins (multiply nodes by two)
viscosity_type = 'constant'

G = create_geometry(nodes, elements, umbilical_artery_radius, decay_factor, arteries_only, umbilical_vein_radius, decay_factor_vein) 

G = calculate_resistance(G, viscosity_model=viscosity_type)

# small system commands
# A,b,bc_export = create_small_matrices(G,bcs,branching_angles=False)
# p,q = solve_small_system(A,b,G,bc_export)

# full matrix commands 
A, b = create_matrices(G, bcs=bcs)
p,q = solve_system(A, b, num_nodes=G.number_of_nodes(), num_edges=G.number_of_edges())

pressures = pd.DataFrame([{"Node" : node, "Pressure" : pressure} for node,pressure in p.items()])
pressures.to_csv(main_pathway + "simulation_pressures.csv")
flows = pd.DataFrame([{"Element" : element, "Flow" : flow} for element,flow in q.items()])
flows.to_csv(main_pathway + "simulation_flows.csv")

#nodes_ps = np.array([node for node in nodes])
nodes_ps = runner.NodesCoordinates(file_nodes, dim)
edges_ps = np.array([element for element in elements])
pressure_ps = np.array([pressure[1] for pressure in p.items()])
flow_ps = np.array([flow[1] for flow in q.items()])

# Polyscope 
if hpc:
    ps.init("openGL3_egl")
elif local:
    ps.init()

# represent vascular tree 
tree = ps.register_curve_network("vascular tree", nodes_ps, edges_ps, color=[155/255,155/255,155/255])

# represent flow and pressure
tree.add_scalar_quantity("flow", flow_ps, defined_on='edges', cmap='reds', enabled=True)
tree.add_scalar_quantity("pressure", pressure_ps, cmap='jet', enabled=True)

# presentation settings
ps.set_ground_plane_mode("none")
ps.set_navigation_style("free")
ps.set_up_dir("y_up")
ps.set_front_dir("z_front")
ps.set_background_color([0,0,0])
ps.show()