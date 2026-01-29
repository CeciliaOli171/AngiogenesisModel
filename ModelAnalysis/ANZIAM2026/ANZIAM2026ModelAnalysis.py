## Libraries 
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import statistics as stats
import seaborn as sns
  
from matplotlib import cm
from matplotlib.ticker import LinearLocator 
from matplotlib.colors import LinearSegmentedColormap
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.patches as patches
from matplotlib.colors import LogNorm

from scipy.spatial import ConvexHull
from scipy.stats import linregress
import scipy.integrate as integrate

import pandas as pd

import os, sys
sys.path.insert(1, "/".join(os.path.realpath(__file__).split("/")[0:-2]))
from ParametersSensitivityRunner import ParametersSensitivityRunner

runner = ParametersSensitivityRunner

## Parameters 
Tcycle = 13.0
time = 2000
Tend = 2000
timestep = 120
TotalTestSourceNb = 10
TotalTestNb = 17
dim = 2
# According to Stratton et al., 2002: "87% were less than 1 cm wide and 47% were less than 5 mm wide." + "97% less than 1 cm deep, and 60% were less than 5 mm deep" [in 2D] i.e. since 1 CD = 10 micrometers = 1e-2 mm therefore 5 mm = 500 CD
# According to literature, for cancer, the tumour and the main blood vessel are separated by 40-140 micrometer, above that, necrosis and death of the tumour i.e. we can use that for the endometriotic model
ref_point = 30
ref_point_centre_lesion = 0 # do not really correspond to the centre of the lesion per se  
AreaPlane = ref_point*250

# VEGF concentration parameters
l_max = 150
c_0 = 0.1
Kc = 1.5e-2
lambdaSprout = 0.08

# Psprout parameters 
cmax = 0.8
cmin = 0.3
Pmax = 0.98
Pmin = 0.5

n = (1/np.log(cmax/cmin))*np.log((Pmax/Pmin)*(1-Pmin)/(1-Pmax))
K = cmax*((1-Pmax)/Pmax)**(1/n)

# Paths 
# HPC
main_file_path = "/hpc/coli171/Results/PaperAngiogenesisModel2025/PaperModel2025Analysis2D/" # 2D 
#main_file_path = "/hpc/coli171/Results/PaperAngiogenesisModel2025/PaperModel2025Analysis3D/" # 3D
# Local
#main_file_path = "/Users/coli171/Chaste/Output/PaperAngiogenesisModel2025/PaperModel2025Analysis2D/" # 2D
#main_file_path = "/Users/coli171/Chaste/Output/PaperAngiogenesisModel2025/PaperModel2025Analysis3D/" # 3D

# file_path_analyticalapproxpde = main_file_path + "CoupledModel2DAnalyticalApproxPde/CoupledModel2DAnalyticalApproxPdeSeed" + str(20) + "Source" + str(5)
# file_nodescoordinates = file_path_analyticalapproxpde + "/results_from_time_0/results.viznodes"
# file_branchesnumber = file_path_analyticalapproxpde + "/results_from_time_0/results.vizbranchnumber"
# file_anastomosis = file_path_analyticalapproxpde + "/results_from_time_0/results.vizanastomosis"
# file_cellmutation = file_path_analyticalapproxpde + "/results_from_time_0/results.vizmutationstates"
# file_cellconnectivity = file_path_analyticalapproxpde + "/results_from_time_0/results.vizconnectivity"
# # nodes coordinates
# data = runner.NodesCoordinates(file_nodescoordinates, 2)
# output_path = main_file_path + "Figures/NodesCoordinatesAngiogenesisModel.csv"
# df = pd.DataFrame(data)
# df.to_csv(output_path, index=False)
# # branches number
# data_branches = runner.BranchNumber(file_branchesnumber)
# output_path_branches = main_file_path + "Figures/BranchesNumberAngiogenesisModel.csv"
# df = pd.DataFrame(data_branches)
# df.to_csv(output_path_branches, index=False)
# # anastomosis 
# data_anastomosis = runner.AnastomosisTypes(file_anastomosis)
# output_path_anastomosis = main_file_path + "Figures/AnastomosisAngiogenesisModel.csv"
# df = pd.DataFrame(data_anastomosis)
# df.to_csv(output_path_anastomosis, index=False)
# # cell mutation
# data_cellmutation = runner.MutationStates(file_cellmutation)
# output_path_cellmutation = main_file_path + "Figures/CellMutationAngiogenesisModel.csv"
# df = pd.DataFrame(data_cellmutation)
# df.to_csv(output_path_cellmutation, index=False)
# # cell connectivity
# data_cellconnectivity = runner.ConnectivityArray(file_cellconnectivity, 2)
# output_path_cellconnectivity = main_file_path + "Figures/CellConnectivityAngiogenesisModel.csv"
# df = pd.DataFrame(data_cellconnectivity)
# df.to_csv(output_path_cellconnectivity, index=False)

## Settings
TestBaseline = False
GraphVEGFConcentration = False
GraphVEGFConcentrationKc = False
GraphVEGFGradient = False 
GraphVEGFGradientKc = False
GraphVEGFGradientXAxis = False

GraphPsprout = False
GraphPsproutSourceTerm = False
GraphPsproutPosition = False
GraphPsproutKc = False

GraphAnastomosis = False
GraphAnastomosisBranches = False
GraphAnastomosisVesselTips = False
GraphAnastomosisRatio = False
GraphAnastomosisRatioBranches = False
GraphTotalAnastomosisSource = False

GraphNbBranches = False
GraphNbVesselTipsTime = False
GraphNbBranchesNbVesselTips = False
GraphBarNbBranches = False

GraphExpectedNbBranches = False
GraphExpectedNbBranchesVEGF = False
GraphNbBranchesANDExpectedNbBranchesGlobalHypothesis = False
GraphNbBranchesANDExpectedNbBranchesLocalHypothesis = False
GraphNbBranchesANDExpectedNbBranchesSameGraph = False

GraphFurthestCell = False
GraphTimeReachingFurthestCell = False
GraphTimeNormFirstVesselTip = False

GraphNbCellsPlane = True
GraphCellDensityInsideLesion = False
GraphDensityCellInsideLesionComparedToTotalCells = False
GraphBarCellsInPlane = False
GraphErrorNbCellsPlanePositionLesion = False

GraphConvergenceTime = False
GraphFirstTimeCellReachingCentreLesion = False
GraphFirstTimeCellReachingLesion = False
GraphErrorFirstTimeCellReachingLesionPositionLesion = False
GraphVascularisation = False

## Functions 

def Psprout(x, c_max, Kc, GlobalOrLocal):
    if GlobalOrLocal == 1:
        c_0 = 0.1
        vegf = (c_max-c_0)*np.exp(-Kc*x)+c_0
        Psprout = lambdaSprout*vegf**n/(K**n + vegf**n)
    else:
        Psprout = lambdaSprout*c_max**n/(K**n + c_max**n)

    return Psprout

def Panastomosis(GlobalOrLocal, Pcomputation):
    if GlobalOrLocal == 1:
        Panastomosis = Pcomputation
    else:
        Panastomosis = Pcomputation

    return Panastomosis

def ExpectationSprouting(x, c_max, Kc, GlobalOrLocal):
    # case of Psprout depending on VEGF
    Sk = Psprout(l_max-0.3*x, c_max, Kc, GlobalOrLocal) # vessel segment division hyp
    #Sk = Psprout(l_max-0.3*0.5*x, c_max, Kc, GlobalOrLocal) # tip cell division hyp
    return Sk 

def ExpectationAnastomosis(GlobalOrLocal, Pcomputation):
    Ak = Panastomosis(GlobalOrLocal, Pcomputation) # vessel segment division hyp
    return Ak 

def ExpectationBranches(t, c_max, Kc, GlobalOrLocal, Pcomputation):
    v = 0.68
    Bn = 2*(1+ExpectationSprouting(v*t, c_max, Kc, GlobalOrLocal)-ExpectationAnastomosis(GlobalOrLocal, Pcomputation))**t - 1
    return Bn

## Plots

if GraphVEGFConcentration:
    fig = plt.subplots(figsize = (12,8), dpi = 300)
    x = np.linspace(0,150,1000)
    y1 = [(0.1 - c_0)*np.exp(-Kc*item)+c_0 for item in x]
    y2 = [(0.2 - c_0)*np.exp(-Kc*item)+c_0 for item in x]
    y5 = [(0.5 - c_0)*np.exp(-Kc*item)+c_0 for item in x]
    y10 = [(1.0 - c_0)*np.exp(-Kc*item)+c_0 for item in x]
    plt.plot(x, y1, label= r'$c_{max} = 0.1$', color = 'xkcd:pale blue', linewidth=5.0)
    plt.plot(x, y2, label= r'$c_{max} = 0.2$', color = 'xkcd:lightblue', linewidth=5.0)
    plt.plot(x, y5, label= r'$c_{max} = 0.5$', color = 'xkcd:blue', linewidth=5.0)
    plt.plot(x, y10, label= r'$c_{max} = 1.0$', color = 'xkcd:royal blue', linewidth=5.0)
    plt.xlabel(r'$x$', fontsize=26) 
    plt.ylabel(r'$c(x)$', fontsize=26)
    plt.xticks(fontsize = 24)
    plt.yticks(fontsize = 24)
    plt.legend(fontsize = 24)
    plt.show()
    plt.savefig(main_file_path + "Figures/SteadyStateCaseVEGFConcentration.png")

if GraphVEGFGradient:
    fig = plt.subplots(figsize = (12,8), dpi = 300)
    x = np.linspace(0,150,1000)
    y1 = [-Kc*(0.1 - c_0)*np.exp(-Kc*item) for item in x]
    y2 = [-Kc*(0.2 - c_0)*np.exp(-Kc*item) for item in x]
    y5 = [-Kc*(0.5 - c_0)*np.exp(-Kc*item) for item in x]
    y10 = [-Kc*(1.0 - c_0)*np.exp(-Kc*item) for item in x]
    plt.plot(x, y1, label= r'$c_{max} = 0.1$', color = 'xkcd:lavender', linewidth=5.0)
    plt.plot(x, y2, label= r'$c_{max} = 0.2$', color = 'xkcd:wisteria', linewidth=5.0)
    plt.plot(x, y5, label= r'$c_{max} = 0.5$', color='xkcd:purple', linewidth=5.0)
    plt.plot(x, y10, label=r'$c_{max} = 1.0$', color='xkcd:royal purple', linewidth=5.0)
    plt.xlabel(r'$x$', fontsize=26) 
    plt.ylabel(r'$\nabla c(x)[0]$', fontsize=26)
    plt.xticks(fontsize = 24)
    plt.yticks(fontsize = 24)
    plt.legend(fontsize = 24)

    plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
    ax = plt.gca()
    offset_text = ax.yaxis.get_offset_text()
    offset_text.set_position((0.0, 1.02))  # Move it to the top
    offset_text.set_horizontalalignment('left')
    offset_text.set_fontsize(24)

    plt.show()
    plt.savefig(main_file_path + "Figures/SteadyStateCaseVEGFGradient.png")

if GraphVEGFGradientXAxis:
    fig = plt.subplots(figsize = (12,8), dpi = 300)
    y = np.linspace(0,300,1000)
    x1 = [0 for item in y]
    x2 = [Kc*(0.2 - c_0)*np.exp(-Kc*item)+c_0 for item in y]
    x5 = [Kc*(0.5 - c_0)*np.exp(-Kc*item)+c_0 for item in y]
    x10 = [Kc*(1.0 - c_0)*np.exp(-Kc*item)+c_0 for item in y]
    #plt.plot(x1, y, label= r'$c_{max} = 0.1$', color = 'xkcd:pale blue')
    plt.plot(x2, y, label= r'$c_{max} = 0.2$', color = 'xkcd:lightblue')
    plt.plot(x5, y, label=r'$c_{max} = 0.5$', color='xkcd:blue')
    plt.plot(x10, y, label=r'$c_{max} = 1.0$', color='xkcd:royal blue')
    plt.ylabel(r'$x$', fontsize=18) 
    plt.xlabel(r'$|\nabla c(x)|$', fontsize=18)
    plt.xticks(fontsize = 16)
    plt.yticks(fontsize = 16)
    plt.legend(fontsize = 12)
    plt.show()

# Linear function and hill function (n = 1, 2 and 3)
if GraphPsproutSourceTerm:
    fig = plt.subplots(figsize = (12,8), dpi = 300)
    sourceterm = np.linspace(0,1.0,100)
    sourcetermticks = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
    Psprout = [lambdaSprout*(item)**n/(K**n + (item)**n) for item in sourceterm]
    plt.plot(sourceterm, Psprout, color='xkcd:orange', linewidth=5.0)
    plt.xlabel(r'$c$', fontsize=26) 
    plt.ylabel(r'$P_{sprout}$', fontsize=26)
    plt.xticks(sourcetermticks, fontsize = 24)
    plt.yticks([0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08], fontsize = 24)
    plt.show()
    plt.savefig(main_file_path + "Figures/SteadyStateCasePsproutSourceTerm.png")

if GraphPsproutPosition:
    fig = plt.subplots(figsize = (12,8), dpi = 300)
    x = np.linspace(0, 150, 1000)
    Psprout1 = [lambdaSprout*c_0**n/(K**n + c_0**n) for item in x]
    Psprout2 = [lambdaSprout*((0.2 - c_0)*np.exp(-Kc*item)+c_0)**n/(K**n + ((0.1 - c_0)*np.exp(-Kc*item)+c_0)**n) for item in x]
    Psprout5 = [lambdaSprout*((0.5 - c_0)*np.exp(-Kc*item)+c_0)**n/(K**n + ((0.5 - c_0)*np.exp(-Kc*item)+c_0)**n) for item in x]
    Psprout10 = [lambdaSprout*((1.0 - c_0)*np.exp(-Kc*item)+c_0)**n/(K**n + ((1.0 - c_0)*np.exp(-Kc*item)+c_0)**n) for item in x]
    plt.plot(x, Psprout1, label=r'$c_{max} = 0.1$', color='xkcd:beige', linewidth=5.0)
    plt.plot(x, Psprout2, label=r'$c_{max} = 0.2$', color='xkcd:tan', linewidth=5.0)
    plt.plot(x, Psprout5, label=r'$c_{max} = 0.5$', color='xkcd:light brown', linewidth=5.0)
    plt.plot(x, Psprout10, label=r'$c_{max} = 1.0$', color='xkcd:brown', linewidth=5.0)
    plt.xlabel(r'$x$', fontsize=26) 
    plt.ylabel(r'$P_{sprout}$', fontsize=26)
    plt.xticks(fontsize = 24)
    plt.yticks([0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08], fontsize = 24)
    plt.legend(fontsize = 24)
    plt.show()
    plt.savefig(main_file_path + "Figures/SteadyStateCasePsproutPosition.png")


# Number of branches depending on source term 
if GraphNbBranches:
    # plot 
    fig, ax = plt.subplots(figsize = (12,8), dpi = 300, layout='constrained')

    # variable 
    sourceterm = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]

    # loop over all the files 
    totalnumberbranches_analyticalapproxpde_average = np.zeros(10)
    totalnumberbranches_constant_average = np.zeros(10)

    for k in range(1, 11):
        totalnumberbranches_analyticalapproxpde = []
        totalnumberbranches_constant = []
        for j in [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]:
            file_path_analyticalapproxpde = main_file_path + "CoupledModel2DAnalyticalApproxPde/CoupledModel2DAnalyticalApproxPdeSeed" + str(j) + "Source" + str(k)
            file_path_constant = main_file_path + "CoupledModel2DConstant/CoupledModel2DConstantSeed" + str(j) + "Source" + str(k)

            file_branchesnumber_analyticalapproxpde = file_path_analyticalapproxpde + "/results_from_time_0/results.vizbranchnumber"
            file_branchesnumber_constant = file_path_constant + "/results_from_time_0/results.vizbranchnumber"
            file_anastomosis_analyticalapproxpde = file_path_analyticalapproxpde + "/results_from_time_0/results.vizanastomosis"
            file_anastomosis_constant = file_path_constant + "/results_from_time_0/results.vizanastomosis"

            totalnumberbranches_analyticalapproxpde.append(runner.TotalNumberBranches(file_branchesnumber_analyticalapproxpde, file_anastomosis_analyticalapproxpde))
            totalnumberbranches_constant.append(runner.TotalNumberBranches(file_branchesnumber_constant, file_anastomosis_constant))

        # scatter plot 
        ax.scatter([sourceterm[k-1] for i in range(10)], totalnumberbranches_analyticalapproxpde, alpha=0.75, color='#757575', marker='.', s = 180.0)
        ax.scatter([sourceterm[k-1] for i in range(10)], totalnumberbranches_constant, alpha=0.75, color='#c90000', marker='D', s = 35)

        # average of the results for one source term 
        totalnumberbranches_analyticalapproxpde_average[k-1] = stats.mean(totalnumberbranches_analyticalapproxpde)
        totalnumberbranches_constant_average[k-1] = stats.mean(totalnumberbranches_constant)
 
    ax.plot(sourceterm, totalnumberbranches_constant_average, marker = 'D', markersize = 6, label = 'ECM Hypothesis', color='#c90000')
    ax.plot(sourceterm, totalnumberbranches_analyticalapproxpde_average, label = 'Lesion Hypothesis', marker = '.', markersize = 15.0, color='#757575')

    ax.legend(loc='upper left', fontsize = 24)
    ax.set_xlabel(r'$c_{max}$', fontsize = 30)
    ax.set_ylabel('Number of Branches', fontsize = 26)
    ax.set_xticks(sourceterm)
    ax.tick_params(axis = 'both', labelsize = 24)
    plt.show()
    plt.savefig(main_file_path + "Figures/NbBranches.png")

if GraphAnastomosisVesselTips:
    # plot 
    fig, ax = plt.subplots(1, 2, figsize = (12,8), dpi = 300, layout='constrained')

    # colours and colourmaps
    cmap_analyticalapproxpde = LinearSegmentedColormap.from_list('red_transparent', [(1.0, 0.0, 0.0, 0.0), (0.75, 0.0, 0.0, 1.0)])
    cmap_constant = LinearSegmentedColormap.from_list('gold_transparent', [(1.0, 1.0, 1.0, 1.0), (0.0, 0.5, 0.5, 1.0)])

    # we represent the percentage of tip cells going through anastomosis at each time steps 
    # we have two different plots: one for the constant case and one for the steady-state of PDE case 
    # for each cmax, we represent the percentage of the ratio anastomosis/tip cells previous time for each time step (except t=0)
    # one cmax = one colour
    sourceterm = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    t = np.linspace(0, Tend, int(Tend/13)+1)

    # Storage for all data points across all (j, k)
    all_tip_cells_analyticalapproxpde = []
    all_anastomosis_analyticalapproxpde = []
    all_tip_cells_constant = []
    all_anastomosis_constant = []

    for k in range(1,11):
        anastomosis_analyticalapproxpde = []
        anastomosis_constant = []
        for j in [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]:
            file_path_analyticalapproxpde = main_file_path + "CoupledModel2DAnalyticalApproxPde/CoupledModel2DAnalyticalApproxPdeSeed" + str(j) + "Source" + str(k)
            file_path_constant = main_file_path + "CoupledModel2DConstant/CoupledModel2DConstantSeed" + str(j) + "Source" + str(k)

            file_anastomosis_analyticalapproxpde = file_path_analyticalapproxpde + "/results_from_time_0/results.vizanastomosis"
            file_anastomosis_constant = file_path_constant + "/results_from_time_0/results.vizanastomosis"

            file_cellmutation_analyticalapproxpde = file_path_analyticalapproxpde + "/results_from_time_0/results.vizmutationstates"
            file_cellmutation_constant = file_path_constant + "/results_from_time_0/results.vizmutationstates"

            tipcells_analyticalapproxpde = runner.VesselTipsPerTimeStep(file_cellmutation_analyticalapproxpde)
            anastomosis_analyticalapproxpde = runner.AnastomosisPerTimeStep(file_anastomosis_analyticalapproxpde)
            tipcells_constant = runner.VesselTipsPerTimeStep(file_cellmutation_constant)
            anastomosis_constant = runner.AnastomosisPerTimeStep(file_anastomosis_constant)

            all_tip_cells_analyticalapproxpde.extend(tipcells_analyticalapproxpde)
            all_anastomosis_analyticalapproxpde.extend(anastomosis_analyticalapproxpde)
            all_tip_cells_constant.extend(tipcells_constant)
            all_anastomosis_constant.extend(anastomosis_constant)

            # scatter plot 
            # if(j == 10):
            #     ax[1].scatter(tipcells_analyticalapproxpde, anastomosis_analyticalapproxpde, color = colors[k-1], alpha = 0.5, label = fr"$c_{{max}}={sourceterm[k-1]}$")
            #     ax[0].scatter(tipcells_constant, anastomosis_constant, color = colors[k-1], alpha = 0.5, label = fr"$c_{{max}}={sourceterm[k-1]}$")
            # else:
            #     ax[1].scatter(tipcells_analyticalapproxpde, anastomosis_analyticalapproxpde, color = colors[k-1], alpha = 0.5)
            #     ax[0].scatter(tipcells_constant, anastomosis_constant, color = colors[k-1], alpha = 0.5)

    all_tip_cells_constant = np.round(np.array(all_tip_cells_constant, dtype=np.float64)).astype(np.int64)
    all_anastomosis_constant = np.round(np.array(all_anastomosis_constant, dtype=np.float64)).astype(np.int64)
    all_tip_cells_analyticalapproxpde = np.round(np.array(all_tip_cells_analyticalapproxpde, dtype=np.float64)).astype(np.int64)
    all_anastomosis_analyticalapproxpde = np.round(np.array(all_anastomosis_analyticalapproxpde, dtype=np.float64)).astype(np.int64)

    idx_constant = np.argsort(all_tip_cells_constant)
    idx_analyticalapproxpde = np.argsort(all_tip_cells_analyticalapproxpde)

    all_tip_cells_sorted_constant = all_tip_cells_constant[idx_constant]
    all_anastomosis_sorted_constant = all_anastomosis_constant[idx_constant]
    all_tip_cells_sorted_analyticalapproxpde = all_tip_cells_analyticalapproxpde[idx_analyticalapproxpde]
    all_anastomosis_sorted_analyticalapproxpde = all_anastomosis_analyticalapproxpde[idx_analyticalapproxpde]

    ux_constant, inv_constant = np.unique(all_tip_cells_sorted_constant, return_inverse=True)
    means_constant = np.array([all_anastomosis_sorted_constant[inv_constant == i].mean() for i in range(len(ux_constant))])
    ux_analyticalapproxpde, inv_analyticalapproxpde = np.unique(all_tip_cells_sorted_analyticalapproxpde, return_inverse=True)
    means_analyticalapproxpde = np.array([all_anastomosis_sorted_analyticalapproxpde[inv_analyticalapproxpde == i].mean() for i in range(len(ux_analyticalapproxpde))])

    # ---- CONSTANT DATA: square-binned colormap ----
    # Define ranges like your hexbin 'extent'
    xmin0, xmax0 = 0, max(all_tip_cells_constant)
    ymin0, ymax0 = 0, max(all_anastomosis_constant)

    # Compute 2D histogram with square bins (gridsize -> bins)
    H0, xedges0, yedges0 = np.histogram2d(
        all_tip_cells_constant,
        all_anastomosis_constant,
        bins=(150, 14),                             # match gridsize=(30,9)
        range=[[xmin0, xmax0], [ymin0, ymax0]]
    )

    # Apply mincnt=1 by masking zeros
    H0_masked = np.ma.masked_where(H0 < 1, H0)

    # Draw squares via pcolormesh, log-normalized like bins='log'
    mesh0 = ax[0].pcolormesh(
        xedges0, yedges0, H0_masked.T,           # transpose to align orientation
        cmap=cmap_constant,
        norm=LogNorm(),                          # equivalent to bins='log'
        shading='auto'                           # fills full squares
    )

    # ---- ANALYTICAL APPROX PDE DATA: square-binned colormap ----
    xmin1, xmax1 = 0, max(all_tip_cells_analyticalapproxpde)
    ymin1, ymax1 = 0, max(all_anastomosis_analyticalapproxpde)

    H1, xedges1, yedges1 = np.histogram2d(
        all_tip_cells_analyticalapproxpde,
        all_anastomosis_analyticalapproxpde,
        bins=(50, 5),                            # match gridsize=(10,3)
        range=[[xmin1, xmax1], [ymin1, ymax1]]
    )

    H1_masked = np.ma.masked_where(H1 < 1, H1)

    mesh1 = ax[1].pcolormesh(
        xedges1, yedges1, H1_masked.T,
        cmap=cmap_analyticalapproxpde,
        norm=LogNorm(),
        shading='auto'
    )

    weights_constant = []
    for k in range(len(all_anastomosis_constant)):
        if (all_tip_cells_constant[k] <= 0):
            weights_constant.append(1)
        elif (all_anastomosis_constant[k] == 0):
            weights_constant.append(1)
        else:
            weights_constant.append(2)
    sqrt_weights_constant = np.sqrt(weights_constant)

    all_tip_cells_constant_weighted = all_tip_cells_constant * sqrt_weights_constant
    all_anastomosis_constant_weighted = all_anastomosis_constant * sqrt_weights_constant

    m_constant, p_constant, r_value, p_value, std_err = linregress(
        all_tip_cells_constant_weighted,
        all_anastomosis_constant_weighted
    )
    ax[0].plot(
        np.linspace(0, 150, 150),
        [m_constant * item + p_constant for item in np.linspace(0, 150, 150)],
        '--k', linewidth=5.0, label=fr'$P_{{a, ECM}} = {m_constant:.3g}$'
    )

    legend_box0 = patches.FancyBboxPatch(
        (0.05, 0.83), 0.3, 0.13, transform=ax[0].transAxes,
        boxstyle="round,pad=0.02", facecolor='white',
        edgecolor='black', linewidth=1, alpha=0.25
    )
    ax[0].add_patch(legend_box0)
    ax[0].text(0.06, 0.83, fr'$P_{{a, ECM}} = {m_constant:.3g}$',
            transform=ax[0].transAxes, fontsize=14, verticalalignment='center')
    ax[0].text(0.35, 0.92, 'Density',
            transform=ax[0].transAxes, fontsize=14,
            verticalalignment='center', horizontalalignment='right')

    cax0 = inset_axes(ax[0], width="100%", height="45%", loc='upper left',
                    bbox_to_anchor=(0.075, 0.85, 0.1, 0.1),
                    bbox_transform=ax[0].transAxes, borderpad=0)
    cb0 = fig.colorbar(mesh0, cax=cax0, orientation='horizontal')  # use mesh0
    cb0.ax.tick_params(labelsize=10)
    cb0.set_label('')

    legend_box1 = patches.FancyBboxPatch(
        (0.05, 0.87), 0.3, 0.09, transform=ax[1].transAxes,
        boxstyle="round,pad=0.02", facecolor='white',
        edgecolor='black', linewidth=1, alpha=0.25
    )
    ax[1].add_patch(legend_box1)
    ax[1].text(0.35, 0.92, 'Density',
            transform=ax[1].transAxes, fontsize=14,
            verticalalignment='center', horizontalalignment='right')

    cax1 = inset_axes(ax[1], width="100%", height="45%", loc='upper left',
                    bbox_to_anchor=(0.075, 0.85, 0.1, 0.1),
                    bbox_transform=ax[1].transAxes, borderpad=0)
    cb1 = fig.colorbar(mesh1, cax=cax1, orientation='horizontal')  # use mesh1
    cb1.ax.tick_params(labelsize=10)


    # weights_constant = []
    # for k in range(len(all_anastomosis_constant)):
    #     if(all_tip_cells_constant[k] <= 0):
    #         weights_constant.append(1)
    #     elif(all_anastomosis_constant[k] == 0):
    #         weights_constant.append(1)
    #     else:
    #         weights_constant.append(2)
    # sqrt_weights_constant = np.sqrt(weights_constant)

    # all_tip_cells_constant_weighted = all_tip_cells_constant * sqrt_weights_constant
    # all_anastomosis_constant_weighted = all_anastomosis_constant * sqrt_weights_constant

    # m_constant, p_constant, r_value, p_value, std_err = linregress(all_tip_cells_constant_weighted, all_anastomosis_constant_weighted)
    # ax[0].plot(np.linspace(0,150,150), [m_constant*item+p_constant for item in np.linspace(0,150,150)], '--k', linewidth = 5.0, label=fr'$P_{{a, ECM}} = {m_constant:.3g}$')

    # # m_analyticalapproxpde, p_analyticalapproxpde, r_value, p_value, std_err = linregress(all_tip_cells_analyticalapproxpde, all_anastomosis_analyticalapproxpde)
    # # ax[0].plot(np.linspace(0,50,50), [m_analyticalapproxpde*item+p_analyticalapproxpde for item in np.linspace(0,50,50)], '--k', linewidth = 5.0, label=fr'$a_{{lesion}} = {m_analyticalapproxpde}$') 
    
    # legend_box0 = patches.FancyBboxPatch((0.05, 0.83), 0.3, 0.13, transform=ax[0].transAxes, boxstyle="round,pad=0.02", facecolor='white', edgecolor='black', linewidth=1, alpha=0.25)
    # ax[0].add_patch(legend_box0)
    # ax[0].text(0.06, 0.83, fr'$P_{{a, ECM}} = {m_constant:.3g}$', transform=ax[0].transAxes, fontsize=14, verticalalignment='center')
    # ax[0].text(0.35, 0.92, 'Density', transform=ax[0].transAxes, fontsize=14, verticalalignment='center', horizontalalignment='right')
    # cax0 = inset_axes(ax[0], width="100%", height="45%", loc='upper left', bbox_to_anchor=(0.075, 0.85, 0.1, 0.1), bbox_transform=ax[0].transAxes, borderpad=0)
    # cb0 = fig.colorbar(hb0, cax=cax0, orientation='horizontal')
    # cb0.ax.tick_params(labelsize=10)
    # cb0.set_label('')

    # legend_box1 = patches.FancyBboxPatch((0.05, 0.87), 0.3, 0.09, transform=ax[1].transAxes, boxstyle="round,pad=0.02", facecolor='white', edgecolor='black', linewidth=1, alpha=0.25)
    # ax[1].add_patch(legend_box1)
    # ax[1].text(0.35, 0.92, 'Density', transform=ax[1].transAxes, fontsize=14, verticalalignment='center', horizontalalignment='right')
    # cax1 = inset_axes(ax[1], width="100%", height="45%", loc='upper left', bbox_to_anchor=(0.075, 0.85, 0.1, 0.1), bbox_transform=ax[1].transAxes, borderpad=0)
    # cb1 = fig.colorbar(hb1, cax=cax1, orientation='horizontal')
    # cb1.ax.tick_params(labelsize=10)
    # cb1.set_label('')

    ax[0].scatter(ux_constant.tolist(), means_constant, label = 'Mean anastomosis values', color = 'xkcd:black')
    #ax[1].plot(ux_analyticalapproxpde.tolist(), means_analyticalapproxpde, label = 'Mean anastomosis values', color = 'xkcd:black')

    ax[1].set_title('Lesion Hypothesis', fontsize = 24)
    ax[0].set_title('ECM Hypothesis', fontsize = 24)
    ax[1].set_xlabel('Number of Tip Cells', fontsize = 26)
    ax[0].set_xlabel('Number of Tip Cells', fontsize = 26)
    ax[0].set_ylabel('Anastomosis Events', fontsize = 26)
    ax[1].tick_params(axis = 'both', labelsize = 24)
    ax[0].tick_params(axis = 'both', labelsize = 24)
    plt.savefig(main_file_path + "Figures/AnastomosisVesselTips.png")


# Expected number of branches
if GraphExpectedNbBranches:
    fig = plt.subplots(figsize = (12,8), dpi = 300)

    t = np.linspace(0.0001, 3000, 50)
    sourceterm = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]

    y2 = [ExpectationBranches(item/Tcycle, 0.2, Kc) for item in t]
    y5 = [ExpectationBranches(item/Tcycle, 0.5, Kc) for item in t]
    y10 = [ExpectationBranches(item/Tcycle, 1.0, Kc) for item in t]

    plt.plot(t, y2, label = 'cmax = 0.2', color = 'xkcd:yellow', marker = 'x')
    plt.plot(t, y5, label = 'cmax = 0.5', color = 'xkcd:orange', marker = 'x')
    plt.plot(t, y10, label = 'cmax = 1.0', color = 'xkcd:red', marker = 'x')
    plt.xlabel(r'$t$', fontsize=18) 
    plt.ylabel('Expected Number of Branches', fontsize=18)
    plt.xticks(fontsize = 16)
    plt.yticks(fontsize = 16)
    plt.legend(fontsize = 12)
    plt.show()
    plt.savefig(main_file_path + "Figures/ExpectedNbBranches.png")

if GraphNbBranchesANDExpectedNbBranchesSameGraph:
    # plot 
    fig, ax = plt.subplots(1, 2, figsize = (12,8), dpi = 300, layout='constrained')

    # variable 
    sourceterm = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]

    # loop over all the files 
    totalnumberbranches_analyticalapproxpde_average = np.zeros(10)
    totalnumberbranches_constant_average = np.zeros(10)

    for k in range(1, 11):
        totalnumberbranches_analyticalapproxpde = []
        totalnumberbranches_constant = []
        for j in [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]:
            file_path_analyticalapproxpde = main_file_path + "CoupledModel2DAnalyticalApproxPde/CoupledModel2DAnalyticalApproxPdeSeed" + str(j) + "Source" + str(k)
            file_path_constant = main_file_path + "CoupledModel2DConstant/CoupledModel2DConstantSeed" + str(j) + "Source" + str(k)

            file_cellmutation_analyticalapproxpde = file_path_analyticalapproxpde + "/results_from_time_0/results.vizmutationstates"
            file_cellmutation_constant = file_path_constant + "/results_from_time_0/results.vizmutationstates"
            file_anastomosis_analyticalapproxpde = file_path_analyticalapproxpde + "/results_from_time_0/results.vizanastomosis"
            file_anastomosis_constant = file_path_constant + "/results_from_time_0/results.vizanastomosis"

            totalnumberbranches_analyticalapproxpde.append(runner.TotalNumberBranches(file_cellmutation_analyticalapproxpde, file_anastomosis_analyticalapproxpde))
            totalnumberbranches_constant.append(runner.TotalNumberBranches(file_cellmutation_constant, file_anastomosis_constant))

        # scatter plot 
        ax[0].scatter([sourceterm[k-1] for i in range(10)], totalnumberbranches_analyticalapproxpde, alpha=0.5, color='#C00000')
        ax[1].scatter([sourceterm[k-1] for i in range(10)], totalnumberbranches_constant, alpha=0.5, color='#008080')

        # average of the results for one source term 
        totalnumberbranches_analyticalapproxpde_average[k-1] = stats.mean(totalnumberbranches_analyticalapproxpde)
        totalnumberbranches_constant_average[k-1] = stats.mean(totalnumberbranches_constant)

    yGlobal = [ExpectationBranches(Tend/Tcycle, item, Kc, 0) for item in sourceterm]
    yLocal = [ExpectationBranches(Tend/Tcycle, item, Kc, 1) for item in sourceterm]
 
    ax[1].plot(sourceterm, yGlobal, label = 'Expected Number Branches (Global Hypothesis)', color = 'xkcd:black')
    ax[0].plot(sourceterm, yLocal, label = 'Expected Number Branches (Local Hypothesis)', color = 'xkcd:black')

    ax[1].plot(sourceterm, totalnumberbranches_constant_average, marker = 'D', markersize = 6, label = 'Global Hypothesis', color='#008080')
    ax[0].plot(sourceterm, totalnumberbranches_analyticalapproxpde_average, label = 'Local Hypothesis', marker = '.', markersize = 15.0, color='#C00000')

    ax[0].legend(loc='upper left', fontsize = 24)
    ax[1].legend(loc='upper left', fontsize = 24)
    ax[0].set_xlabel(r'$c_{max}$', fontsize = 26)
    ax[1].set_xlabel(r'$c_{max}$', fontsize = 26)
    ax[0].set_ylabel('Number of Branches', fontsize = 26)
    ax[0].set_xticks(sourceterm)
    ax[1].set_xticks(sourceterm)
    ax[0].tick_params(axis = 'both', labelsize = 24)
    ax[1].tick_params(axis = 'both', labelsize = 24)
    plt.savefig(main_file_path + "Figures/ExpectedNbBranchesBothHypotheses.png")

if GraphExpectedNbBranchesVEGF:
    fig = plt.subplots(figsize = (12,8), dpi = 300)
    sourceterm = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]

    # data from Zhang et al. for comparison 
    datapoints = [[0.7, 107], [1.77, 105], [2.05, 90], [1.3, 92], [2.25, 52], [1.8, 70], [1.13, 15], [2.3, 75], [2.25, 67], [1.95, 47], [1.7, 68], [1.25, 45], [1.3, 42], [0.48, 48], [0.25, 52], [0.49, 20], [0.4, 29], [1.0, 23], [1.1, 28], [0.9, 10], [0.6, 12], [0.85, 98], [0.7, 83], [0.8, 82], [0.68, 21], [0.27, 28], [0.98, 28], [0.47, 6], [0.4, 14], [0.32, 16], [0.9, 78], [0.8, 72], [0.92, 68], [0.77, 61], [0.94, 57], [0.7, 42], [0.58, 42], [0.57, 41]]
    microvesseldensity = [datapoints[k][0]/2.5 for k in range(len(datapoints))]
    VEGFgeneexpressionlinearregression = [20.34*4*item*2.5+30.90*4 for item in microvesseldensity]
    VEGFgeneexpression = [datapoints[k][1]*4 for k in range(len(datapoints))]

    Kc1 = 1.5e-2
    Kc3 = 2.2e-2
    Kc5 = 2.6e-2
    Kc10 = 3e-2

    TendExperimental = 3000
    nbbranches = [ExpectationBranches(TendExperimental/Tcycle, item, Kc1, 0) for item in sourceterm]
    nbbranches1 = [ExpectationBranches(TendExperimental/Tcycle, item, Kc1, 1) for item in sourceterm]
    nbbranches3 = [ExpectationBranches(TendExperimental/Tcycle, item, Kc3, 1) for item in sourceterm]
    nbbranches5 = [ExpectationBranches(TendExperimental/Tcycle, item, Kc5, 1) for item in sourceterm]
    nbbranches10 = [ExpectationBranches(TendExperimental/Tcycle, item, Kc10, 1) for item in sourceterm]

    plt.plot(sourceterm, nbbranches, color='xkcd:red', label = 'Global Hypothesis', linewidth = 5.0)
    # plt.plot(sourceterm, nbbranches1, color='xkcd:coral', label = 'Local Hypothesis (m = 1)', linewidth = 5.0)
    # plt.plot(sourceterm, nbbranches3, color='xkcd:pale red', label = 'Local Hypothesis (m = 3)', linewidth = 5.0)
    # plt.plot(sourceterm, nbbranches5, color='xkcd:red', label = 'Local Hypothesis (m = 5)', linewidth = 5.0)
    # plt.plot(sourceterm, nbbranches10, color='xkcd:dark red', label = 'Local Hypothesis (m = 10)', linewidth = 5.0)
    plt.scatter(microvesseldensity, VEGFgeneexpression, label = 'Experimental Results \n (Zhang et al., 2006)', color = 'g')
    plt.plot(microvesseldensity, VEGFgeneexpressionlinearregression, color = 'g', linewidth = 5.0)
    plt.xlabel(r'$c$', fontsize=26) 
    plt.ylabel('Number of Branches', fontsize=26)
    plt.xticks(fontsize = 24)
    plt.yticks(fontsize = 24)
    plt.legend(loc = 'lower right', fontsize = 24)
    plt.savefig(main_file_path + "Figures/ExpectedNbBranchesExperiment.png")


# Number of cells after a certain plane depending on source term for different Psprout 
if GraphNbCellsPlane:
    # plot 
    fig, ax = plt.subplots(figsize = (12,8), dpi = 300, layout='constrained')

    sourceterm = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]

    # loop over all the files 
    cellsafterplane_analyticalapproxpde_average = np.zeros(10)
    cellsafterplane_constant_average = np.zeros(10)
    for k in range(1,11):
        cellsafterplane_analyticalapproxpde = []
        cellsafterplane_constant = []
        for j in [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]:
            file_path_analyticalapproxpde = main_file_path + "CoupledModel2DAnalyticalApproxPde/CoupledModel2DAnalyticalApproxPdeSeed" + str(j) + "Source" + str(k)
            file_path_constant = main_file_path + "CoupledModel2DConstant/CoupledModel2DConstantSeed" + str(j) + "Source" + str(k)

            file_cellmutation_analyticalapproxpde = file_path_analyticalapproxpde + "/results_from_time_0/results.vizmutationstates"
            file_nodescoordinates_analyticalapproxpde = file_path_analyticalapproxpde + "/results_from_time_0/results.viznodes"
            file_cellmutation_constant = file_path_constant + "/results_from_time_0/results.vizmutationstates"
            file_nodescoordinates_constant = file_path_constant + "/results_from_time_0/results.viznodes"

            cellsafterplane_analyticalapproxpde.append(runner.NbCellsAfterPlane(file_nodescoordinates_analyticalapproxpde, ref_point, dim))
            cellsafterplane_constant.append(runner.NbCellsAfterPlane(file_nodescoordinates_constant, ref_point, dim))

        # scatter plot 
        ax.scatter([sourceterm[k-1] for i in range(10)], cellsafterplane_analyticalapproxpde, alpha=0.65, color='#757575')
        ax.scatter([sourceterm[k-1] for i in range(10)], cellsafterplane_constant, alpha = 0.65, color='#c90000')

        # average of the results for one source term 
        cellsafterplane_analyticalapproxpde_average[k-1] = stats.mean(cellsafterplane_analyticalapproxpde)
        cellsafterplane_constant_average[k-1] = stats.mean(cellsafterplane_constant)

    ax.plot(sourceterm, cellsafterplane_constant_average, marker = 'D', markersize = 6, label = 'ECM Hypothesis', color='#c90000')
    ax.plot(sourceterm, cellsafterplane_analyticalapproxpde_average, label = 'Lesion Hypothesis', marker = '.', markersize = 15.0, color='#757575')

    ax.set_xlabel(r'$c_{max}$', fontsize = 26)
    ax.set_ylabel('Vessel Segments Inside Lesion', fontsize = 26)
    ax.set_xticks(sourceterm)
    ax.legend(loc='upper left', fontsize = 24)
    ax.tick_params(axis = 'both', labelsize = 24)
    plt.show()
    plt.savefig(main_file_path + "Figures/CellsAfterPlane" + str(ref_point) + ".png")

if GraphDensityCellInsideLesionComparedToTotalCells:
    # plot 
    fig, ax = plt.subplots(figsize = (12,8), dpi = 300, layout='constrained')

    sourceterm = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]

    # loop over all the files 
    celldensitytotal_analyticalapproxpde_average = np.zeros(10)
    celldensitytotal_constant_average = np.zeros(10)
    for k in range(1,11):
        celldensitytotal_analyticalapproxpde = []
        celldensitytotal_constant = []
        for j in [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]:
            file_path_analyticalapproxpde = main_file_path + "CoupledModel2DAnalyticalApproxPde/CoupledModel2DAnalyticalApproxPdeSeed" + str(j) + "Source" + str(k)
            file_path_constant = main_file_path + "CoupledModel2DConstant/CoupledModel2DConstantSeed" + str(j) + "Source" + str(k)

            file_cellmutation_analyticalapproxpde = file_path_analyticalapproxpde + "/results_from_time_0/results.vizmutationstates"
            file_nodescoordinates_analyticalapproxpde = file_path_analyticalapproxpde + "/results_from_time_0/results.viznodes"
            file_cellmutation_constant = file_path_constant + "/results_from_time_0/results.vizmutationstates"
            file_nodescoordinates_constant = file_path_constant + "/results_from_time_0/results.viznodes"

            celldensitytotal_analyticalapproxpde.append(100*runner.NbCellsAfterPlane(file_nodescoordinates_analyticalapproxpde, ref_point, dim)/runner.TotalNumberCells(file_cellmutation_analyticalapproxpde))
            celldensitytotal_constant.append(100*runner.NbCellsAfterPlane(file_nodescoordinates_constant, ref_point, dim)/runner.TotalNumberCells(file_cellmutation_constant))

        # scatter plot 
        ax.scatter([sourceterm[k-1] for i in range(10)], celldensitytotal_analyticalapproxpde, alpha=0.65, color='#757575', marker='.', s = 180.0)
        ax.scatter([sourceterm[k-1] for i in range(10)], celldensitytotal_constant, alpha = 0.65, color='#c90000', marker='D', s = 35)

        # average of the results for one source term 
        celldensitytotal_analyticalapproxpde_average[k-1] = stats.mean(celldensitytotal_analyticalapproxpde)
        celldensitytotal_constant_average[k-1] = stats.mean(celldensitytotal_constant)


    ax.plot(sourceterm, celldensitytotal_constant_average, marker = 'D', markersize = 6, label = 'ECM Hypothesis', color='#c90000')
    ax.plot(sourceterm, celldensitytotal_analyticalapproxpde_average, label = 'Lesion Hypothesis', marker = '.', markersize = 15.0, color='#757575')

    ax.set_xlabel(r'$c_{max}$', fontsize = 30)
    ax.set_ylabel("Percentage of Networks' Cells \n Inside Lesion", fontsize = 24)
    ax.set_xticks(sourceterm)
    ax.legend(loc='upper left', fontsize = 24)
    ax.tick_params(axis = 'both', labelsize = 24)
    plt.savefig(main_file_path + "Figures/CellDensityInsideLesionComparedToTotalCells.png")

if GraphErrorNbCellsPlanePositionLesion:
    # plot 
    fig, ax = plt.subplots(figsize = (12,8), dpi = 300, layout='constrained')

    DistanceLesionVessel = [60, 40, 20]
    cellsinsidelesion_analyticalapproxpde = []
    cellsinsidelesion_constant = []

    # loop over all the files 
    for i in DistanceLesionVessel:
        celldensity_analyticalapproxpde = []
        celldensity_constant = []

        for k in [10]:
            for j in [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]:
                file_path_analyticalapproxpde = main_file_path + "CoupledModel2DAnalyticalApproxPde/CoupledModel2DAnalyticalApproxPdeSeed" + str(j) + "Source" + str(k)
                file_path_constant = main_file_path + "CoupledModel2DConstant/CoupledModel2DConstantSeed" + str(j) + "Source" + str(k)

                file_cellmutation_analyticalapproxpde = file_path_analyticalapproxpde + "/results_from_time_0/results.vizmutationstates"
                file_nodescoordinates_analyticalapproxpde = file_path_analyticalapproxpde + "/results_from_time_0/results.viznodes"
                file_cellmutation_constant = file_path_constant + "/results_from_time_0/results.vizmutationstates"
                file_nodescoordinates_constant = file_path_constant + "/results_from_time_0/results.viznodes"

                celldensity_analyticalapproxpde.append(runner.NbCellsAfterPlane(file_nodescoordinates_analyticalapproxpde, i, dim))
                celldensity_constant.append(runner.NbCellsAfterPlane(file_nodescoordinates_constant, i, dim))

        cellsinsidelesion_analyticalapproxpde.append(celldensity_analyticalapproxpde)
        cellsinsidelesion_constant.append(celldensity_constant)

    DistanceLesionVesselTicks = [90, 110, 130]
    width = 5
    position_constant = [item - width/2 for item in DistanceLesionVesselTicks]
    position_analyticalapproxpde = [item + width/2 for item in DistanceLesionVesselTicks]

    box_colors = ['#C00000', '#008080']
    bp_analyticalapproxpde = ax.violinplot(cellsinsidelesion_analyticalapproxpde, positions = position_analyticalapproxpde, widths=width, showmeans=True, showmedians=False, showextrema=True) 
    bp_constant = ax.violinplot(cellsinsidelesion_constant, positions = position_constant, widths=width, showmeans=True, showmedians=False, showextrema=True) 

    for i in range(len(position_analyticalapproxpde)):
        d = position_analyticalapproxpde[i]
        min_value = min(cellsinsidelesion_analyticalapproxpde[i])
        max_value = max(cellsinsidelesion_analyticalapproxpde[i])
        ax.vlines(d, min_value, max_value, linewidth = 2.0, color = box_colors[0])
    for i in range(len(position_constant)):
        d = position_constant[i]
        min_value = min(cellsinsidelesion_constant[i])
        max_value = max(cellsinsidelesion_constant[i])
        ax.vlines(d, min_value, max_value, linewidth = 2.0, color = box_colors[1])

    for pc in bp_analyticalapproxpde['bodies']:
        pc.set_facecolor(box_colors[0])
        pc.set_edgecolor(box_colors[0])
        pc.set_alpha(0.5)
        pc.set_linewidth(0)
    bp_analyticalapproxpde['cmeans'].set_color(box_colors[0])
    bp_analyticalapproxpde['cmeans'].set_linewidth(2)
    bp_analyticalapproxpde['cmins'].set_color(box_colors[0])
    bp_analyticalapproxpde['cmins'].set_linewidth(2)
    bp_analyticalapproxpde['cmaxes'].set_color(box_colors[0])
    bp_analyticalapproxpde['cmaxes'].set_linewidth(2)

    for pc in bp_constant['bodies']:
        pc.set_facecolor(box_colors[1])
        pc.set_edgecolor(box_colors[1])
        pc.set_alpha(0.5)
        pc.set_linewidth(0)
    bp_constant['cmeans'].set_color(box_colors[1])
    bp_constant['cmeans'].set_linewidth(2)
    bp_constant['cmins'].set_color(box_colors[1])
    bp_constant['cmins'].set_linewidth(2)
    bp_constant['cmaxes'].set_color(box_colors[1])
    bp_constant['cmaxes'].set_linewidth(2)

    ax.set_xlabel('Distance Lesion / Parent Vessel', fontsize = 26)
    ax.set_ylabel('Cell Density Inside Lesion', fontsize = 26)
    custom_lines = [plt.Line2D([0], [0], color=box_colors[1], lw=4), plt.Line2D([0], [0], color=box_colors[0], lw=4)]
    ax.legend(custom_lines, ['ECM Hypothesis', 'Lesion Hypothesis'], loc='upper right', fontsize = 24)
    ax.set_xticks(DistanceLesionVesselTicks)
    ax.tick_params(axis = 'both', labelsize = 24)
    plt.savefig(main_file_path + "Figures/ErrorCellDensityPositionLesion.png")


# Times
if GraphFirstTimeCellReachingLesion:
    # plot 
    fig, ax = plt.subplots(figsize = (12,8), dpi = 300, layout='constrained')

    sourceterm = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]

    # loop over all the files 
    timefirstreachinglesion_analyticalapproxpde_average = np.zeros(10)
    timefirstreachinglesion_constant_average = np.zeros(10)
    for k in range(1,11):
        timefirstreachinglesion_analyticalapproxpde = []
        timefirstreachinglesion_constant = []
        for j in [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]:
            file_path_analyticalapproxpde = main_file_path + "CoupledModel2DAnalyticalApproxPde/CoupledModel2DAnalyticalApproxPdeSeed" + str(j) + "Source" + str(k)
            file_path_constant = main_file_path + "CoupledModel2DConstant/CoupledModel2DConstantSeed" + str(j) + "Source" + str(k)

            file_cellmutation_analyticalapproxpde = file_path_analyticalapproxpde + "/results_from_time_0/results.vizmutationstates"
            file_nodescoordinates_analyticalapproxpde = file_path_analyticalapproxpde + "/results_from_time_0/results.viznodes"
            file_cellmutation_constant = file_path_constant + "/results_from_time_0/results.vizmutationstates"
            file_nodescoordinates_constant = file_path_constant + "/results_from_time_0/results.viznodes"

            timefirstreachinglesion_analyticalapproxpde.append(runner.TimeFirstReachingPlane(file_nodescoordinates_analyticalapproxpde, ref_point, dim))
            timefirstreachinglesion_constant.append(runner.TimeFirstReachingPlane(file_nodescoordinates_constant, ref_point, dim))

        # scatter plot 
        ax.scatter([sourceterm[k-1] for i in range(10)], timefirstreachinglesion_analyticalapproxpde, alpha=0.65, color='#C00000', marker='.', s = 180.0)
        ax.scatter([sourceterm[k-1] for i in range(10)], timefirstreachinglesion_constant, alpha=0.65, color='#008080', marker='D', s = 35)

        # average of the results for one source term 
        timefirstreachinglesion_analyticalapproxpde_average[k-1] = stats.mean(timefirstreachinglesion_analyticalapproxpde)
        timefirstreachinglesion_constant_average[k-1] = stats.mean(timefirstreachinglesion_constant)

    ax.plot(sourceterm, timefirstreachinglesion_constant_average, marker = 'D', markersize = 6, label = 'ECM Hypothesis', color='#008080')
    ax.plot(sourceterm, timefirstreachinglesion_analyticalapproxpde_average, label = 'Lesion Hypothesis', marker = '.', markersize = 15.0, color='#C00000')
 
    ax.legend(loc='upper left', fontsize = 24)
    ax.set_xlabel(r'$c_{max}$', fontsize = 30)
    ax.set_ylabel('Time For Network to Reach Lesion', fontsize = 26)
    ax.set_xticks(sourceterm)
    ax.set_yticks([1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000])
    ax.set_ylim(1160, 2015)
    ax.tick_params(axis = 'both', labelsize = 24)
    plt.show()
    plt.savefig(main_file_path + "Figures/TimeFirstCellReachingLesion.png")

if GraphVascularisation:
    # plot 
    fig, ax = plt.subplots(figsize = (12,8), dpi = 300, layout='constrained')

    # loop over all the files 
    vascularisation_analyticalapproxpde = []
    vascularisation_constant = []

    for k in range(1,11):
        time_constant = []
        time_analyticalapproxpde = []
        vascularisationvalue_constant = 0
        vascularisationvalue_analyticalapproxpde = 0
        for j in [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]:
            file_path_analyticalapproxpde = main_file_path + "CoupledModel2DAnalyticalApproxPde/CoupledModel2DAnalyticalApproxPdeSeed" + str(j) + "Source" + str(k)
            file_path_constant = main_file_path + "CoupledModel2DConstant/CoupledModel2DConstantSeed" + str(j) + "Source" + str(k)

            file_nodescoordinates_analyticalapproxpde = file_path_analyticalapproxpde + "/results_from_time_0/results.viznodes"
            file_nodescoordinates_constant = file_path_constant + "/results_from_time_0/results.viznodes"

            timeseed_analyticalapproxpde = runner.TimeFirstReachingPlane(file_nodescoordinates_analyticalapproxpde, ref_point, dim)
            timeseed_constant = runner.TimeFirstReachingPlane(file_nodescoordinates_constant, ref_point, dim)

            # remove the not vascularised network
            if timeseed_analyticalapproxpde < 1989.0:
                time_analyticalapproxpde.append(timeseed_analyticalapproxpde)
                vascularisationvalue_analyticalapproxpde += 1
            if timeseed_constant < 1989.0:
                time_constant.append(timeseed_constant)
                vascularisationvalue_constant += 1

        vascularisation_analyticalapproxpde.append(vascularisationvalue_analyticalapproxpde)
        vascularisation_constant.append(vascularisationvalue_constant)

    SourceTermTicks = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    barWidth = 0.025
    position_constant = [item - barWidth/2 for item in SourceTermTicks]
    position_analyticalapproxpde = [item + barWidth/2 for item in SourceTermTicks]

    bp_analyticalapproxpde = ax.bar(position_analyticalapproxpde, vascularisation_analyticalapproxpde, width=barWidth, color = '#757575', label = 'Lesion Hypothesis') 
    bp_constant = ax.bar(position_constant, vascularisation_constant, width=barWidth, color = '#c90000', label = 'ECM Hypothesis') 
   
    ax.set_xlabel(r'$c_{max}$', fontsize = 26)
    ax.set_ylabel('Number of Realisations with \n Lesion Vascularised', fontsize = 26)
    ax.legend(loc='upper left', fontsize = 24)
    ax.set_xticks(SourceTermTicks)
    ax.set_yticks([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
    ax.tick_params(axis = 'both', labelsize = 24)
    plt.savefig(main_file_path + "Figures/VascularisationBarPlot" + str(ref_point) + ".png")