## Libraries 
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import statistics as stats

from matplotlib import cm
from matplotlib.ticker import LinearLocator 
from matplotlib.colors import LinearSegmentedColormap

from scipy.spatial import ConvexHull
import scipy.integrate as integrate

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
ref_point = 40
ref_point_centre_lesion = 0 # do not really correspond to the centre of the lesion per se  
AreaPlane = ref_point*150

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

GraphNbBranches = False
GraphAnastomosis = False
GraphAnastomosisBranches = False
GraphAnastomosisTipCells = False
GraphAnastomosisRatio = False
GraphAnastomosisRatioBranches = False
GraphAnastomosisSource = False
GraphExpectedNbBranches = False
GraphExpectedNbBranchesVEGF = False
GraphNbBranchesANDExpectedNbBranchesGlobalHypothesis = False
GraphNbBranchesANDExpectedNbBranchesLocalHypothesis = False
GraphNbBranchesANDExpectedNbBranchesSameGraph = False
GraphBarNbBranches = False

GraphFurthestCell = False
GraphTimeReachingFurthestCell = False
GraphTimeNormFirstTipCell = False

GraphNbCellsPlane = False
GraphCellDensityInsideLesion = False
GraphDensityCellInsideLesionComparedToTotalCells = False
GraphBarCellsInPlane = False
GraphErrorNbCellsPlanePositionLesion = True

GraphConvergenceTime = False
GraphFirstTimeCellReachingCentreLesion = False
GraphFirstTimeCellReachingLesion = False
GraphErrorFirstTimeCellReachingLesionPositionLesion = True
GraphErrorVascularisationPercentagePositionLesion = True

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

if GraphVEGFConcentrationKc:
    fig = plt.subplots(figsize = (12,8), dpi = 300)
    x = np.linspace(0,150,1000)
    C = 1.0
    Kc1 = 1.5e-2
    Kc3 = 2.2e-2
    Kc5 = 2.6e-2
    Kc10 = 3e-2
    y1 = [(C - c_0)*np.exp(-Kc1*item)+c_0 for item in x]
    y3 = [(C - c_0)*np.exp(-Kc3*item)+c_0 for item in x]
    y5 = [(C - c_0)*np.exp(-Kc5*item)+c_0 for item in x]
    y10 = [(C - c_0)*np.exp(-Kc10*item)+c_0 for item in x]
    plt.plot(x, y1, label= r'$m = 1$', color = 'xkcd:pale blue', linewidth = 5.0)
    plt.plot(x, y3, label= r'$m = 3$', color = 'xkcd:lightblue', linewidth = 5.0)
    plt.plot(x, y5, label= r'$m = 5$', color = 'xkcd:blue', linewidth = 5.0)
    plt.plot(x, y10, label= r'$m = 10$', color = 'xkcd:royal blue', linewidth = 5.0)
    plt.xlabel(r'$x$', fontsize=26) 
    plt.ylabel(r'$c(x)$', fontsize=26)
    plt.xticks(fontsize = 24)
    plt.yticks([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0], fontsize = 24)
    plt.legend(fontsize = 24)
    plt.show()
    plt.savefig(main_file_path + "Figures/VEGFConcentrationKc.png")

if GraphVEGFGradientKc:
    fig = plt.subplots(figsize = (12,8), dpi = 300)
    x = np.linspace(0,150,1000)
    C = 1.0
    Kc1 = 1.5e-2
    Kc3 = 2.2e-2
    Kc5 = 2.6e-2
    Kc10 = 3e-2
    y1 = [-Kc1*(C - c_0)*np.exp(-Kc1*item) for item in x]
    y3 = [-Kc3*(C - c_0)*np.exp(-Kc3*item) for item in x]
    y5 = [-Kc5*(C - c_0)*np.exp(-Kc5*item) for item in x]
    y10 = [-Kc10*(C - c_0)*np.exp(-Kc10*item) for item in x]
    plt.plot(x, y1, label= r'$m = 1$', color = 'xkcd:lavender', linewidth = 5.0)
    plt.plot(x, y3, label= r'$m = 3$', color = 'xkcd:wisteria', linewidth = 5.0)
    plt.plot(x, y5, label= r'$m = 5$', color = 'xkcd:purple', linewidth = 5.0)
    plt.plot(x, y10, label= r'$m = 10$', color = 'xkcd:royal purple', linewidth = 5.0)
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
    plt.savefig(main_file_path + "Figures/VEGFGradientKc.png")


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
if GraphPsprout:
    x = np.linspace(0,1,100)
    y = [0.2*(item**n/(item**n + K**n)) for item in x]
    #y0 = [0.2*item for item in x] # linear function
    plt.plot(x, y, label='hill function (n=%.2f, K=%.2f)' %(n,K), color='xkcd:red')
    # plt.plot(x, y0, label='linear function', color='xkcd:brown green')
    plt.xlabel(r'$c$', fontsize=24) 
    plt.ylabel(r'$P_{sprout}$', fontsize=24)
    plt.xticks(fontsize = 24)
    plt.yticks(fontsize = 24)
    plt.legend(fontsize = 24)
    plt.show()

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

if GraphPsproutKc:
    fig = plt.subplots(figsize = (12,8), dpi = 300)
    x = np.linspace(0, 150, 1000)

    C = 1.0
    Kc1 = 1.5e-2
    Kc3 = 2.2e-2
    Kc5 = 2.6e-2
    Kc10 = 3e-2

    Psprout1 = [lambdaSprout*((C - c_0)*np.exp(-Kc1*item)+c_0)**n/(K**n + ((C - c_0)*np.exp(-Kc1*item)+c_0)**n) for item in x]
    Psprout3 = [lambdaSprout*((C - c_0)*np.exp(-Kc3*item)+c_0)**n/(K**n + ((C - c_0)*np.exp(-Kc3*item)+c_0)**n) for item in x]
    Psprout5 = [lambdaSprout*((C - c_0)*np.exp(-Kc5*item)+c_0)**n/(K**n + ((C - c_0)*np.exp(-Kc5*item)+c_0)**n) for item in x]
    Psprout10 = [lambdaSprout*((C - c_0)*np.exp(-Kc10*item)+c_0)**n/(K**n + ((C - c_0)*np.exp(-Kc10*item)+c_0)**n) for item in x]

    plt.plot(x, Psprout1, label=r'$m = 1$', color='xkcd:beige', linewidth=5.0)
    plt.plot(x, Psprout3, label=r'$m = 3$', color='xkcd:tan', linewidth=5.0)
    plt.plot(x, Psprout5, label=r'$m = 5$', color='xkcd:light brown', linewidth=5.0)
    plt.plot(x, Psprout10, label=r'$m = 10$', color='xkcd:brown', linewidth=5.0)
    plt.xlabel(r'$x$', fontsize=26) 
    plt.ylabel(r'$P_{sprout}$', fontsize=26)
    plt.xticks(fontsize = 24)
    plt.yticks([0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08], fontsize = 24)
    plt.legend(fontsize = 24)
    plt.show()
    plt.savefig(main_file_path + "Figures/PsproutKc.png")

# Time convergence depending on Psprout and method (PDE or analytical)
if GraphConvergenceTime:
    PsproutFunctions = ["Linear \n Function", "Hill Function \n (n=1)", "Hill Function \n (n=2)", "Hill Function \n (n=3)"]
    ConvergenceTimePDE = [711, 848, 1007, 821]
    ConvergenceTimeAnalytical = [13, 14, 15, 16]
    barWidth = 0.5
    fig = plt.subplots(figsize = (12,8), dpi = 300)
    PsproutFunctionsPDE = np.arange(len(PsproutFunctions))
    PsproutFunctionsAnalytical = [elem + barWidth for elem in PsproutFunctionsPDE]
    plt.bar(PsproutFunctionsPDE, ConvergenceTimePDE, color = 'xkcd:green', width = barWidth, label = 'PDE Solution')
    plt.bar(PsproutFunctionsAnalytical, ConvergenceTimeAnalytical, color = 'xkcd:olive', width = barWidth, label = 'Analytical Approximation')
    plt.xticks(PsproutFunctionsPDE+barWidth/2, PsproutFunctions)
    plt.ylabel('Convergence Time')
    plt.legend()
    plt.show()


# Number of branches depending on source term for the different Psprout 
if GraphNbBranches:
    # plot 
    fig, ax = plt.subplots(figsize = (12,8), dpi = 300, layout='constrained')

    # variable 
    sourceterm = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]

    # loop over all the files 
    totalnumberbranches_analyticalapproxpde_average = np.zeros(10)
    totalnumberbranches_constant_average = np.zeros(10)

    error_analyticalapproxpde_lower = []
    error_analyticalapproxpde_upper = []
    error_constant_lower = []
    error_constant_upper = []

    for k in range(1, 11):
        totalnumberbranches_analyticalapproxpde = []
        totalnumberbranches_constant = []
        for j in [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]:
            file_path_analyticalapproxpde = main_file_path + "CoupledModel2DAnalyticalApproxPde/CoupledModel2DAnalyticalApproxPdeSeed" + str(j) + "Source" + str(k)
            file_path_constant = main_file_path + "CoupledModel2DConstant/CoupledModel2DConstantSeed" + str(j) + "Source" + str(k)

            file_cellmutation_analyticalapproxpde = file_path_analyticalapproxpde + "/results_from_time_0/results.vizmutationstates"
            file_cellmutation_constant = file_path_constant + "/results_from_time_0/results.vizmutationstates"

            totalnumberbranches_analyticalapproxpde.append(runner.TotalNumberBranches(file_cellmutation_analyticalapproxpde))
            totalnumberbranches_constant.append(runner.TotalNumberBranches(file_cellmutation_constant))

        # scatter plot 
        ax.scatter([sourceterm[k-1] for i in range(10)], totalnumberbranches_analyticalapproxpde, alpha=0.5, color='#8c0000')
        ax.scatter([sourceterm[k-1] for i in range(10)], totalnumberbranches_constant, alpha=0.5, color='#ff0004')

        # average of the results for one source term 
        totalnumberbranches_analyticalapproxpde_average[k-1] = stats.mean(totalnumberbranches_analyticalapproxpde)
        totalnumberbranches_constant_average[k-1] = stats.mean(totalnumberbranches_constant)

        low_error_analyticalapproxpde = stats.mean(totalnumberbranches_analyticalapproxpde)-max(min(totalnumberbranches_analyticalapproxpde), 1e-2)
        high_error_analyticalapproxpde = max(totalnumberbranches_analyticalapproxpde)-stats.mean(totalnumberbranches_analyticalapproxpde)
        low_error_constant = stats.mean(totalnumberbranches_analyticalapproxpde)-max(min(totalnumberbranches_constant),1e-2)
        high_error_constant = max(totalnumberbranches_constant)-stats.mean(totalnumberbranches_constant)

        error_analyticalapproxpde_lower.append(low_error_analyticalapproxpde)
        error_analyticalapproxpde_upper.append(high_error_analyticalapproxpde)
        error_constant_lower.append(low_error_constant)
        error_constant_upper.append(high_error_constant)

    error_analyticalapproxpde = np.array([error_analyticalapproxpde_lower, error_analyticalapproxpde_upper])
    error_constant = np.array([error_constant_lower, error_constant_upper])
 
    ax.plot(sourceterm, totalnumberbranches_constant_average, marker = 'D', markersize = 6, label = 'Global Hypothesis', color='#ff0004')
    ax.plot(sourceterm, totalnumberbranches_analyticalapproxpde_average, label = 'Local Hypothesis', marker = '.', markersize = 15.0, color='#8c0000')
    #ax.errorbar(sourceterm, totalnumberbranches_constant_average, yerr = error_constant, marker = 'D', markersize = 6, label = 'Global Hypothesis', color='xkcd:brown')
    #ax.errorbar(sourceterm, totalnumberbranches_analyticalapproxpde_average, yerr = error_analyticalapproxpde, label = 'Local Hypothesis', marker = '.', markersize = 15.0, color='xkcd:tan')
    ax.legend(loc='upper left', fontsize = 24)

    ax.set_xlabel(r'$c_{max}$', fontsize = 26)
    ax.set_ylabel('Number of Branches', fontsize = 26)
    ax.set_xticks(sourceterm)
    ax.tick_params(axis = 'both', labelsize = 24)
    plt.show()
    plt.savefig(main_file_path + "Figures/NbBranches.png")

def Psprout(x, c_max, Kc, GlobalOrLocal):
    if GlobalOrLocal == 1:
        vegf = (c_max-c_0)*np.exp(-Kc*x)+c_0
        Psprout = lambdaSprout*vegf**n/(K**n + vegf**n)
    else:
        Psprout = lambdaSprout*c_max**n/(K**n + c_max**n)

    return Psprout

def Panastomosis(x, c_max, Kc, GlobalOrLocal, a):
    n_anastomosis = n
    K_anastomosis = K

    if GlobalOrLocal == 1:
        lambdaAnastomosis = a*lambdaSprout
        vegf = (c_max-c_0)*np.exp(-Kc*x)+c_0
        Panastomosis = lambdaAnastomosis*vegf**n_anastomosis/(K_anastomosis**n_anastomosis + vegf**n_anastomosis)
        #Panastomosis = a
    else:
        lambdaAnastomosis = a*lambdaSprout
        Panastomosis = lambdaAnastomosis*c_max**n_anastomosis/(K_anastomosis**n_anastomosis + c_max**n_anastomosis)
        #Panastomosis = a

    return Panastomosis

def ExpectationSprouting(t, c_max, Kc, GlobalOrLocal):
    # case of Psprout depending on VEGF
    Sk = Psprout(l_max-0.3*t, c_max, Kc, GlobalOrLocal)
    #Sk = Psprout(l_max-0.3*0.5*t, c_max, Kc, GlobalOrLocal)
    return Sk 

def ExpectationAnastomosis(t, c_max, Kc, GlobalOrLocal, a):
    Ak = Panastomosis(l_max-0.3*t, c_max, Kc, GlobalOrLocal, a)
    #Ak = Panastomosis(l_max-0.3*0.5*t, c_max, Kc, GlobalOrLocal, a)
    return Ak 

def ExpectationBranches(t, c_max, Kc, GlobalOrLocal, a):
    #Bn = 2*(1+ExpectationSprouting(t, c_max, Kc, GlobalOrLocal)-ExpectationAnastomosis(t, c_max, Kc, GlobalOrLocal, a))**t - 1
    Bn = 2*(1+(1-a)*ExpectationSprouting(t, c_max, Kc, GlobalOrLocal))**t - 1
    return Bn

if GraphAnastomosis:
    # plot 
    fig, ax = plt.subplots(1, 2, figsize = (12,8), dpi = 300, layout='constrained')

    # we represent the percentage of tip cells going through anastomosis at each time steps 
    # we have two different plots: one for the constant case and one for the steady-state of PDE case 
    # for each cmax, we represent the percentage of the ratio anastomosis/tip cells previous time for each time step (except t=0)
    # one cmax = one colour
    sourceterm = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    colors = ["xkcd:tomato red", "xkcd:coral", "xkcd:orange", "xkcd:sunshine yellow", "xkcd:lime", "xkcd:bright turquoise", "xkcd:bright sky blue", "xkcd:royal blue", "xkcd:bright lavender", "xkcd:magenta"]
    t = np.linspace(0, Tend, int(Tend/13)+1)

    for k in range(1,11):
        anastomosis_analyticalapproxpde = []
        anastomosis_constant = []
        for j in [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]:
            file_path_analyticalapproxpde = main_file_path + "CoupledModel2DAnalyticalApproxPde/CoupledModel2DAnalyticalApproxPdeSeed" + str(j) + "Source" + str(k)
            file_path_constant = main_file_path + "CoupledModel2DConstant/CoupledModel2DConstantSeed" + str(j) + "Source" + str(k)

            file_anastomosis_analyticalapproxpde = file_path_analyticalapproxpde + "/results_from_time_0/results.vizanastomosis"
            file_anastomosis_constant = file_path_constant + "/results_from_time_0/results.vizanastomosis"

            anastomosis_analyticalapproxpde = runner.AnastomosisPerTimeStep(file_anastomosis_analyticalapproxpde)
            anastomosis_constant = runner.AnastomosisPerTimeStep(file_anastomosis_constant)

            if(j == 10):
                ax[0].scatter(t, anastomosis_analyticalapproxpde, color = colors[k-1], alpha = 0.5, label = f"cmax={sourceterm[k-1]}")
                ax[1].scatter(t, anastomosis_constant, color = colors[k-1], alpha = 0.5, label = f"cmax={sourceterm[k-1]}")
            else:
                ax[0].scatter(t, anastomosis_analyticalapproxpde, color = colors[k-1], alpha = 0.5)
                ax[1].scatter(t, anastomosis_constant, color = colors[k-1], alpha = 0.5)
    
    ax[0].set_title('Local Hypothesis', fontsize = 24)
    ax[1].set_title('Global Hypothesis', fontsize = 24)
    ax[0].set_xlabel(r'$t$', fontsize = 26)
    ax[1].set_xlabel(r'$t$', fontsize = 26)
    ax[0].set_ylabel('Anastomosis Events', fontsize = 26)
    ax[0].legend(loc='upper left', fontsize = 12)
    # ax[1].legend(loc='upper left', fontsize = 24)
    ax[0].tick_params(axis = 'both', labelsize = 24)
    ax[1].tick_params(axis = 'both', labelsize = 24)
    plt.savefig(main_file_path + "Figures/Anastomosis.png")

if GraphAnastomosisBranches:
    # plot 
    fig, ax = plt.subplots(1, 2, figsize = (12,8), dpi = 300, layout='constrained')

    # we represent the percentage of tip cells going through anastomosis at each time steps 
    # we have two different plots: one for the constant case and one for the steady-state of PDE case 
    # for each cmax, we represent the percentage of the ratio anastomosis/tip cells previous time for each time step (except t=0)
    # one cmax = one colour
    sourceterm = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    colors = ["xkcd:tomato red", "xkcd:coral", "xkcd:orange", "xkcd:sunshine yellow", "xkcd:lime", "xkcd:bright turquoise", "xkcd:bright sky blue", "xkcd:royal blue", "xkcd:bright lavender", "xkcd:magenta"]
    t = np.linspace(0, Tend, int(Tend/13)+1)

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

            branches_analyticalapproxpde = runner.BranchesPerTimeStep(file_cellmutation_analyticalapproxpde)
            anastomosis_analyticalapproxpde = runner.AnastomosisPerTimeStep(file_anastomosis_analyticalapproxpde)
            branches_constant = runner.BranchesPerTimeStep(file_cellmutation_constant)
            anastomosis_constant = runner.AnastomosisPerTimeStep(file_anastomosis_constant)

            if(j == 10):
                ax[0].scatter(branches_analyticalapproxpde, anastomosis_analyticalapproxpde, color = colors[k-1], alpha = 0.5, label = fr"$c_{{\max}}={sourceterm[k-1]}$")
                ax[1].scatter(branches_constant, anastomosis_constant, color = colors[k-1], alpha = 0.5, label = fr"$c_{{\max}}={sourceterm[k-1]}$")
            else:
                ax[0].scatter(branches_analyticalapproxpde, anastomosis_analyticalapproxpde, color = colors[k-1], alpha = 0.5)
                ax[1].scatter(branches_constant, anastomosis_constant, color = colors[k-1], alpha = 0.5)

    coeff_analyticalapproxpde = 13/140
    coeff_constant = 120/750
    ax[0].plot(np.linspace(0,150,150), [coeff_analyticalapproxpde*item for item in np.linspace(0,150,150)], '--k', linewidth = 5.0, label= fr'$a_{{local}} = {coeff_analyticalapproxpde}$')
    ax[1].plot(np.linspace(0,750,750), [coeff_constant*item for item in np.linspace(0,750,750)], '--k', linewidth = 5.0, label= fr'$a_{{global}} = {coeff_constant}$')
    
    ax[0].set_title('Local Hypothesis', fontsize = 24)
    ax[1].set_title('Global Hypothesis', fontsize = 24)
    ax[0].set_xlabel('Nb of Branches', fontsize = 26)
    ax[1].set_xlabel('Number of Branches', fontsize = 26)
    ax[0].set_ylabel('Anastomosis Events', fontsize = 26)
    ax[0].legend(loc='upper left', fontsize = 12)
    ax[1].legend(loc='upper left', fontsize = 12)
    ax[0].tick_params(axis = 'both', labelsize = 24)
    ax[1].tick_params(axis = 'both', labelsize = 24)
    plt.savefig(main_file_path + "Figures/AnastomosisBranches.png")

if GraphAnastomosisTipCells:
    # plot 
    fig, ax = plt.subplots(1, 2, figsize = (12,8), dpi = 300, layout='constrained')

    # we represent the percentage of tip cells going through anastomosis at each time steps 
    # we have two different plots: one for the constant case and one for the steady-state of PDE case 
    # for each cmax, we represent the percentage of the ratio anastomosis/tip cells previous time for each time step (except t=0)
    # one cmax = one colour
    sourceterm = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    colors = ["xkcd:tomato red", "xkcd:coral", "xkcd:orange", "xkcd:sunshine yellow", "xkcd:lime", "xkcd:bright turquoise", "xkcd:bright sky blue", "xkcd:royal blue", "xkcd:bright lavender", "xkcd:magenta"]
    t = np.linspace(0, Tend, int(Tend/13)+1)

    list_m_analyticalapproxpde = []
    list_m_constant = []

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

            tipcells_analyticalapproxpde = runner.TipCellsPerTimeStep(file_cellmutation_analyticalapproxpde)
            anastomosis_analyticalapproxpde = runner.AnastomosisPerTimeStep(file_anastomosis_analyticalapproxpde)
            tipcells_constant = runner.TipCellsPerTimeStep(file_cellmutation_constant)
            anastomosis_constant = runner.AnastomosisPerTimeStep(file_anastomosis_constant)

            if(j == 10):
                ax[0].scatter(tipcells_analyticalapproxpde, anastomosis_analyticalapproxpde, color = colors[k-1], alpha = 0.5, label = f"cmax={sourceterm[k-1]}")
                ax[1].scatter(tipcells_constant, anastomosis_constant, color = colors[k-1], alpha = 0.5, label = f"cmax={sourceterm[k-1]}")
            else:
                ax[0].scatter(tipcells_analyticalapproxpde, anastomosis_analyticalapproxpde, color = colors[k-1], alpha = 0.5)
                ax[1].scatter(tipcells_constant, anastomosis_constant, color = colors[k-1], alpha = 0.5)

    ax[0].plot(np.linspace(0,50,50), [0.31*item for item in np.linspace(0,50,50)], '--k', linewidth = 5.0, label=r'$a_{local} = 0.31$')
    ax[1].plot(np.linspace(0,150,150), [0.52*item for item in np.linspace(0,150,150)], '--k', linewidth = 5.0, label=r'$a_{global} = 0.52$') 
    
    ax[0].set_title('Local Hypothesis', fontsize = 24)
    ax[1].set_title('Global Hypothesis', fontsize = 24)
    ax[0].set_xlabel('Nb of Tip Cells', fontsize = 26)
    ax[1].set_xlabel('Nb of Tip Cells', fontsize = 26)
    ax[0].set_ylabel('Anastomosis Events', fontsize = 26)
    ax[0].legend(loc='upper left', fontsize = 12)
    ax[1].legend(loc='upper left', fontsize = 12)
    ax[0].tick_params(axis = 'both', labelsize = 24)
    ax[1].tick_params(axis = 'both', labelsize = 24)
    plt.savefig(main_file_path + "Figures/AnastomosisTipCells.png")

if GraphAnastomosisRatio:
    # plot 
    fig, ax = plt.subplots(1, 2, figsize = (12,8), dpi = 300, layout='constrained')

    # we represent the percentage of tip cells going through anastomosis at each time steps 
    # we have two different plots: one for the constant case and one for the steady-state of PDE case 
    # for each cmax, we represent the percentage of the ratio anastomosis/tip cells previous time for each time step (except t=0)
    # one cmax = one colour
    sourceterm = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    colors = ["xkcd:tomato red", "xkcd:coral", "xkcd:orange", "xkcd:sunshine yellow", "xkcd:lime", "xkcd:bright turquoise", "xkcd:bright sky blue", "xkcd:royal blue", "xkcd:bright lavender", "xkcd:magenta"]
    t = np.linspace(0, Tend, int(Tend/13))

    for k in range(1,11):
        anastomosisratio_analyticalapproxpde = []
        anastomosisratio_constant = []
        for j in [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]:
            file_path_analyticalapproxpde = main_file_path + "CoupledModel2DAnalyticalApproxPde/CoupledModel2DAnalyticalApproxPdeSeed" + str(j) + "Source" + str(k)
            file_path_constant = main_file_path + "CoupledModel2DConstant/CoupledModel2DConstantSeed" + str(j) + "Source" + str(k)

            file_cellmutation_analyticalapproxpde = file_path_analyticalapproxpde + "/results_from_time_0/results.vizmutationstates"
            file_anastomosis_analyticalapproxpde = file_path_analyticalapproxpde + "/results_from_time_0/results.vizanastomosis"
            file_cellmutation_constant = file_path_constant + "/results_from_time_0/results.vizmutationstates"
            file_anastomosis_constant = file_path_constant + "/results_from_time_0/results.vizanastomosis"

            anastomosisratio_analyticalapproxpde = runner.RatioAnastomosisPerTimeStep(file_anastomosis_analyticalapproxpde, file_cellmutation_analyticalapproxpde)
            anastomosisratio_constant = runner.RatioAnastomosisPerTimeStep(file_anastomosis_constant, file_cellmutation_constant)

            if(j == 10):
                ax[0].scatter(t, anastomosisratio_analyticalapproxpde, color = colors[k-1], alpha = 0.5, label = f"cmax={sourceterm[k-1]}")
                ax[1].scatter(t, anastomosisratio_constant, color = colors[k-1], alpha = 0.5, label = f"cmax={sourceterm[k-1]}")
            else:
                ax[0].scatter(t, anastomosisratio_analyticalapproxpde, color = colors[k-1], alpha = 0.5)
                ax[1].scatter(t, anastomosisratio_constant, color = colors[k-1], alpha = 0.5)
    
    ax[0].set_title('Local Hypothesis', fontsize = 24)
    ax[1].set_title('Global Hypothesis', fontsize = 24)
    ax[0].set_xlabel(r'$t$', fontsize = 26)
    ax[1].set_xlabel(r'$t$', fontsize = 26)
    ax[0].set_ylabel('Anastomosis Ratio', fontsize = 26)
    ax[0].legend(loc='upper left', fontsize = 12)
    # ax[1].legend(loc='upper left', fontsize = 24)
    ax[0].tick_params(axis = 'both', labelsize = 24)
    ax[1].tick_params(axis = 'both', labelsize = 24)
    plt.savefig(main_file_path + "Figures/AnastomosisRatio.png")

if GraphAnastomosisRatioBranches:
    # plot 
    fig, ax = plt.subplots(1, 2, figsize = (12,8), dpi = 300, layout='constrained')

    # we represent the percentage of tip cells going through anastomosis at each time steps 
    # we have two different plots: one for the constant case and one for the steady-state of PDE case 
    # for each cmax, we represent the percentage of the ratio anastomosis/tip cells previous time for each time step (except t=0)
    # one cmax = one colour
    sourceterm = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    colors = ["xkcd:tomato red", "xkcd:coral", "xkcd:orange", "xkcd:sunshine yellow", "xkcd:lime", "xkcd:bright turquoise", "xkcd:bright sky blue", "xkcd:royal blue", "xkcd:bright lavender", "xkcd:magenta"]
    t = np.linspace(0, Tend, 20)

    for k in range(1,11):
        anastomosisratio_analyticalapproxpde = []
        anastomosisratio_constant = []
        for j in [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]:
            file_path_analyticalapproxpde = main_file_path + "CoupledModel2DAnalyticalApproxPde/CoupledModel2DAnalyticalApproxPdeSeed" + str(j) + "Source" + str(k)
            file_path_constant = main_file_path + "CoupledModel2DConstant/CoupledModel2DConstantSeed" + str(j) + "Source" + str(k)

            file_cellmutation_analyticalapproxpde = file_path_analyticalapproxpde + "/results_from_time_0/results.vizmutationstates"
            file_anastomosis_analyticalapproxpde = file_path_analyticalapproxpde + "/results_from_time_0/results.vizanastomosis"
            file_cellmutation_constant = file_path_constant + "/results_from_time_0/results.vizmutationstates"
            file_anastomosis_constant = file_path_constant + "/results_from_time_0/results.vizanastomosis"

            tipcells_analyticalapproxpde = runner.TipCellsPerTimeStep(file_cellmutation_analyticalapproxpde)
            tipcells_constant = runner.TipCellsPerTimeStep(file_cellmutation_constant)
            anastomosisratio_analyticalapproxpde = runner.RatioAnastomosisPerTimeStep(file_anastomosis_analyticalapproxpde, file_cellmutation_analyticalapproxpde)
            anastomosisratio_constant = runner.RatioAnastomosisPerTimeStep(file_anastomosis_constant, file_cellmutation_constant)

            if(j == 10):
                ax[0].scatter(tipcells_analyticalapproxpde, anastomosisratio_analyticalapproxpde, color = colors[k-1], alpha = 0.5, label = f"cmax={sourceterm[k-1]}")
                ax[1].scatter(tipcells_analyticalapproxpde, anastomosisratio_constant, color = colors[k-1], alpha = 0.5, label = f"cmax={sourceterm[k-1]}")
            else:
                ax[0].scatter(tipcells_analyticalapproxpde, anastomosisratio_analyticalapproxpde, color = colors[k-1], alpha = 0.5)
                ax[1].scatter(tipcells_analyticalapproxpde, anastomosisratio_constant, color = colors[k-1], alpha = 0.5)
                
    ax[0].set_title('Local Hypothesis', fontsize = 24)
    ax[1].set_title('Global Hypothesis', fontsize = 24)
    ax[0].set_xlabel('Number Branches', fontsize = 26)
    ax[1].set_xlabel('Number Branches', fontsize = 26)
    ax[0].set_ylabel('Anastomosis Ratio', fontsize = 26)
    ax[1].legend(loc='upper right', fontsize = 12)
    ax[0].tick_params(axis = 'both', labelsize = 24)
    ax[1].tick_params(axis = 'both', labelsize = 24)
    plt.savefig(main_file_path + "Figures/AnastomosisRatioNbBranches.png")

if GraphAnastomosisSource:
    # plot 
    fig, ax = plt.subplots(1, 2, figsize = (12,8), dpi = 300, layout='constrained')

    # we represent the percentage of tip cells going through anastomosis at each time steps 
    # we have two different plots: one for the constant case and one for the steady-state of PDE case 
    # for each cmax, we represent the percentage of the ratio anastomosis/tip cells previous time for each time step (except t=0)
    # one cmax = one colour
    sourceterm = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    colors = ["xkcd:tomato red", "xkcd:coral", "xkcd:orange", "xkcd:sunshine yellow", "xkcd:lime", "xkcd:bright turquoise", "xkcd:bright sky blue", "xkcd:royal blue", "xkcd:bright lavender", "xkcd:magenta"]
    t = np.linspace(0, Tend, 20)

    for k in range(1,11):
        anastomosisratio_analyticalapproxpde = []
        anastomosisratio_constant = []
        for j in [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]:
            file_path_analyticalapproxpde = main_file_path + "CoupledModel2DAnalyticalApproxPde/CoupledModel2DAnalyticalApproxPdeSeed" + str(j) + "Source" + str(k)
            file_path_constant = main_file_path + "CoupledModel2DConstant/CoupledModel2DConstantSeed" + str(j) + "Source" + str(k)

            file_cellmutation_analyticalapproxpde = file_path_analyticalapproxpde + "/results_from_time_0/results.vizmutationstates"
            file_anastomosis_analyticalapproxpde = file_path_analyticalapproxpde + "/results_from_time_0/results.vizanastomosis"
            file_cellmutation_constant = file_path_constant + "/results_from_time_0/results.vizmutationstates"
            file_anastomosis_constant = file_path_constant + "/results_from_time_0/results.vizanastomosis"

            tipcells_analyticalapproxpde = runner.TipCellsPerTimeStep(file_cellmutation_analyticalapproxpde)
            tipcells_constant = runner.TipCellsPerTimeStep(file_cellmutation_constant)
            anastomosisratio_analyticalapproxpde = runner.RatioAnastomosisPerTimeStep(file_anastomosis_analyticalapproxpde, file_cellmutation_analyticalapproxpde)
            anastomosisratio_constant = runner.RatioAnastomosisPerTimeStep(file_anastomosis_constant, file_cellmutation_constant)

            if(j == 10):
                ax[0].scatter([sourceterm[k-1] for j in range(len(anastomosisratio_analyticalapproxpde))], anastomosisratio_analyticalapproxpde, color = colors[k-1], alpha = 0.5, label = f"cmax={sourceterm[k-1]}")
                ax[1].scatter([sourceterm[k-1] for j in range(len(anastomosisratio_constant))], anastomosisratio_constant, color = colors[k-1], alpha = 0.5, label = f"cmax={sourceterm[k-1]}")
            else:
                ax[0].scatter([sourceterm[k-1] for j in range(len(anastomosisratio_analyticalapproxpde))], anastomosisratio_analyticalapproxpde, color = colors[k-1], alpha = 0.5)
                ax[1].scatter([sourceterm[k-1] for j in range(len(anastomosisratio_constant))], anastomosisratio_constant, color = colors[k-1], alpha = 0.5)
                
    ax[0].set_title('Local Hypothesis', fontsize = 24)
    ax[1].set_title('Global Hypothesis', fontsize = 24)
    ax[0].set_xlabel(r'$c_{max}$', fontsize = 26)
    ax[1].set_xlabel(r'$c_{max}$', fontsize = 26)
    ax[0].set_ylabel('Anastomosis Ratio', fontsize = 26)
    ax[0].legend(loc='upper left', fontsize = 12)
    ax[0].tick_params(axis = 'both', labelsize = 24)
    ax[1].tick_params(axis = 'both', labelsize = 24)
    plt.savefig(main_file_path + "Figures/AnastomosisRatioSource.png")

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

if GraphNbBranchesANDExpectedNbBranchesGlobalHypothesis:
    # plot 
    fig, ax = plt.subplots(figsize = (12,8), dpi = 300, layout='constrained')

    # variable 
    sourceterm = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]

    # loop over all the files 
    totalnumberbranches_constant_average = np.zeros(10)

    for k in range(1, 11):
        totalnumberbranches_constant = []
        for j in [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]:
            file_path_constant = main_file_path + "CoupledModel2DConstant/CoupledModel2DConstantSeed" + str(j) + "Source" + str(k)
            file_cellmutation_constant = file_path_constant + "/results_from_time_0/results.vizmutationstates"

            totalnumberbranches_constant.append(runner.TotalNumberBranches(file_cellmutation_constant))

        # scatter plot 
        ax.scatter([sourceterm[k-1] for i in range(10)], totalnumberbranches_constant, alpha=0.5, color='xkcd:grey') #color='#ff0004'

        # average of the results for one source term 
        totalnumberbranches_constant_average[k-1] = stats.mean(totalnumberbranches_constant)

    # yGlobal4 = [ExpectationBranches(Tend/Tcycle, item, Kc, 0, 0.04) for item in sourceterm]
    # yGlobal41 = [ExpectationBranches(Tend/Tcycle, item, Kc, 0, 0.041) for item in sourceterm]
    # yGlobal42 = [ExpectationBranches(Tend/Tcycle, item, Kc, 0, 0.042) for item in sourceterm]
    # yGlobal44 = [ExpectationBranches(Tend/Tcycle, item, Kc, 0, 0.044) for item in sourceterm]
    # yGlobal46 = [ExpectationBranches(Tend/Tcycle, item, Kc, 0, 0.046) for item in sourceterm]
    # yGlobal48 = [ExpectationBranches(Tend/Tcycle, item, Kc, 0, 0.048) for item in sourceterm]
    # yGlobal5 = [ExpectationBranches(Tend/Tcycle, item, Kc, 0, 0.05) for item in sourceterm]

    aglobal = ((120/750)+1)/2
    yGlobal = [ExpectationBranches(Tend/Tcycle, item, Kc, 0, aglobal) for item in sourceterm]
 
    ax.plot(sourceterm, totalnumberbranches_constant_average, marker = 'D', markersize = 6, label = 'Global Hypothesis', color= 'xkcd:black') #color='#ff0004'
    ax.plot(sourceterm, yGlobal, label = fr'$a_{{global}} = {aglobal}$', color = 'xkcd:pink')

    ax.legend(loc='upper left', fontsize = 24)
    ax.set_xlabel(r'$c_{max}$', fontsize = 26)
    ax.set_ylabel('Number of Branches', fontsize = 26)
    ax.set_xticks(sourceterm)
    ax.tick_params(axis = 'both', labelsize = 24)
    plt.savefig(main_file_path + "Figures/ExpectedNbBranchesGlobalHypothesis.png")

if GraphNbBranchesANDExpectedNbBranchesLocalHypothesis:
    # plot 
    fig, ax = plt.subplots(figsize = (12,8), dpi = 300, layout='constrained')

    # variable 
    sourceterm = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]

    # loop over all the files 
    totalnumberbranches_analyticalapproxpde_average = np.zeros(10)

    for k in range(1, 11):
        totalnumberbranches_analyticalapproxpde = []
        for j in [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]:
            file_path_analyticalapproxpde = main_file_path + "CoupledModel2DAnalyticalApproxPde/CoupledModel2DAnalyticalApproxPdeSeed" + str(j) + "Source" + str(k)
            file_cellmutation_analyticalapproxpde = file_path_analyticalapproxpde + "/results_from_time_0/results.vizmutationstates"

            totalnumberbranches_analyticalapproxpde.append(runner.TotalNumberBranches(file_cellmutation_analyticalapproxpde))        

        # scatter plot 
        ax.scatter([sourceterm[k-1] for i in range(10)], totalnumberbranches_analyticalapproxpde, alpha=0.5, color='#8c0000')

        # average of the results for one source term 
        totalnumberbranches_analyticalapproxpde_average[k-1] = stats.mean(totalnumberbranches_analyticalapproxpde)

    ax.plot(sourceterm, totalnumberbranches_analyticalapproxpde_average, label = 'Local Hypothesis', marker = '.', markersize = 15.0, color='#8c0000')

    alocal = ((13/140)+1)/2
    yLocal = [ExpectationBranches(Tend/Tcycle, item, Kc, 1, alocal) for item in sourceterm]
    ax.plot(sourceterm, yLocal, label = fr'$a_{{local}} = {alocal}$', color = 'xkcd:pink')

    ax.legend(loc='upper left', fontsize = 24)
    ax.set_xlabel(r'$c_{max}$', fontsize = 26)
    ax.set_ylabel('Number of Branches', fontsize = 26)
    ax.set_xticks(sourceterm)
    ax.tick_params(axis = 'both', labelsize = 24)
    plt.savefig(main_file_path + "Figures/ExpectedNbBranchesLocalHypothesis.png")

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

            totalnumberbranches_analyticalapproxpde.append(runner.TotalNumberBranches(file_cellmutation_analyticalapproxpde))
            totalnumberbranches_constant.append(runner.TotalNumberBranches(file_cellmutation_constant))

        # scatter plot 
        ax[0].scatter([sourceterm[k-1] for i in range(10)], totalnumberbranches_analyticalapproxpde, alpha=0.5, color='#8c0000')
        ax[1].scatter([sourceterm[k-1] for i in range(10)], totalnumberbranches_constant, alpha=0.5, color='#ff0004')

        # average of the results for one source term 
        totalnumberbranches_analyticalapproxpde_average[k-1] = stats.mean(totalnumberbranches_analyticalapproxpde)
        totalnumberbranches_constant_average[k-1] = stats.mean(totalnumberbranches_constant)

    yGlobal = [ExpectationBranches(Tend/Tcycle, item, Kc, 0) for item in sourceterm]
    yLocal = [ExpectationBranches(Tend/Tcycle, item, Kc, 1) for item in sourceterm]
 
    ax[1].plot(sourceterm, yGlobal, label = 'Expected Number Branches (Global Hypothesis)', color = 'xkcd:black')
    ax[0].plot(sourceterm, yLocal, label = 'Expected Number Branches (Local Hypothesis)', color = 'xkcd:black')

    ax[1].plot(sourceterm, totalnumberbranches_constant_average, marker = 'D', markersize = 6, label = 'Global Hypothesis', color='#ff0004')
    ax[0].plot(sourceterm, totalnumberbranches_analyticalapproxpde_average, label = 'Local Hypothesis', marker = '.', markersize = 15.0, color='#8c0000')

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

    error_analyticalapproxpde_lower = []
    error_analyticalapproxpde_upper = []
    error_constant_lower = []
    error_constant_upper = []

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
        ax.scatter([sourceterm[k-1] for i in range(10)], cellsafterplane_analyticalapproxpde, alpha=0.5, color='#8c0000')
        ax.scatter([sourceterm[k-1] for i in range(10)], cellsafterplane_constant, alpha = 0.5, color='#ff0004')

        # average of the results for one source term 
        cellsafterplane_analyticalapproxpde_average[k-1] = stats.mean(cellsafterplane_analyticalapproxpde)
        cellsafterplane_constant_average[k-1] = stats.mean(cellsafterplane_constant)

        low_error_analyticalapproxpde = stats.mean(cellsafterplane_analyticalapproxpde)-max(min(cellsafterplane_analyticalapproxpde), 1e-2)
        high_error_analyticalapproxpde = max(cellsafterplane_analyticalapproxpde)-stats.mean(cellsafterplane_analyticalapproxpde)
        low_error_constant = stats.mean(cellsafterplane_analyticalapproxpde)-max(min(cellsafterplane_constant),1e-2)
        high_error_constant = max(cellsafterplane_constant)-stats.mean(cellsafterplane_constant)

        error_analyticalapproxpde_lower.append(low_error_analyticalapproxpde)
        error_analyticalapproxpde_upper.append(high_error_analyticalapproxpde)
        error_constant_lower.append(low_error_constant)
        error_constant_upper.append(high_error_constant)

    error_analyticalapproxpde = np.array([error_analyticalapproxpde_lower, error_analyticalapproxpde_upper])
    error_constant = np.array([error_constant_lower, error_constant_upper])

    ax.plot(sourceterm, cellsafterplane_constant_average, marker = 'D', markersize = 6, label = 'Global Hypothesis', color='#ff0004')
    ax.plot(sourceterm, cellsafterplane_analyticalapproxpde_average, label = 'Local Hypothesis', marker = '.', markersize = 15.0, color='#8c0000')
    #ax.errorbar(sourceterm, cellsafterplane_constant_average, yerr = error_constant, marker = 'D', markersize = 6, label = 'Global Hypothesis', color='xkcd:brown')
    #ax.errorbar(sourceterm, cellsafterplane_analyticalapproxpde_average, yerr = error_analyticalapproxpde, label = 'Local Hypothesis', marker = '.', markersize = 15.0, color='xkcd:tan')
    ax.set_xlabel(r'$c_{max}$', fontsize = 26)
    ax.set_ylabel('Cells In Lesion', fontsize = 26)
    ax.set_xticks(sourceterm)
    ax.legend(loc='upper left', fontsize = 24)
    ax.tick_params(axis = 'both', labelsize = 24)
    plt.show()
    plt.savefig(main_file_path + "Figures/CellsAfterPlane.png")

if GraphCellDensityInsideLesion:
    # plot 
    fig, ax = plt.subplots(figsize = (12,8), dpi = 300, layout='constrained')

    sourceterm = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]

    # loop over all the files 
    celldensity_analyticalapproxpde_average = np.zeros(10)
    celldensity_constant_average = np.zeros(10)
    for k in range(1,11):
        celldensity_analyticalapproxpde = []
        celldensity_constant = []
        for j in [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]:
            file_path_analyticalapproxpde = main_file_path + "CoupledModel2DAnalyticalApproxPde/CoupledModel2DAnalyticalApproxPdeSeed" + str(j) + "Source" + str(k)
            file_path_constant = main_file_path + "CoupledModel2DConstant/CoupledModel2DConstantSeed" + str(j) + "Source" + str(k)

            file_cellmutation_analyticalapproxpde = file_path_analyticalapproxpde + "/results_from_time_0/results.vizmutationstates"
            file_nodescoordinates_analyticalapproxpde = file_path_analyticalapproxpde + "/results_from_time_0/results.viznodes"
            file_cellmutation_constant = file_path_constant + "/results_from_time_0/results.vizmutationstates"
            file_nodescoordinates_constant = file_path_constant + "/results_from_time_0/results.viznodes"

            celldensity_analyticalapproxpde.append(runner.NbCellsAfterPlane(file_nodescoordinates_analyticalapproxpde, ref_point, dim)/AreaPlane)
            celldensity_constant.append(runner.NbCellsAfterPlane(file_nodescoordinates_constant, ref_point, dim)/AreaPlane)

        # scatter plot 
        ax.scatter([sourceterm[k-1] for i in range(10)], celldensity_analyticalapproxpde, alpha=0.5, color='#8c0000')
        ax.scatter([sourceterm[k-1] for i in range(10)], celldensity_constant, alpha = 0.5, color='#ff0004')

        # average of the results for one source term 
        celldensity_analyticalapproxpde_average[k-1] = stats.mean(celldensity_analyticalapproxpde)
        celldensity_constant_average[k-1] = stats.mean(celldensity_constant)

    ax.plot(sourceterm, celldensity_constant_average, marker = 'D', markersize = 6, label = 'Global Hypothesis', color='#ff0004')
    ax.plot(sourceterm, celldensity_analyticalapproxpde_average, label = 'Local Hypothesis', marker = '.', markersize = 15.0, color='#8c0000')

    ax.set_xlabel(r'$c_{max}$', fontsize = 26)
    ax.set_ylabel('Cell Density Inside Lesion', fontsize = 26)
    ax.set_xticks(sourceterm)
    ax.legend(loc='upper left', fontsize = 24)
    ax.tick_params(axis = 'both', labelsize = 24)
    plt.show()
    plt.savefig(main_file_path + "Figures/CellDensity.png")

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
        ax.scatter([sourceterm[k-1] for i in range(10)], celldensitytotal_analyticalapproxpde, alpha=0.5, color='#8c0000')
        ax.scatter([sourceterm[k-1] for i in range(10)], celldensitytotal_constant, alpha = 0.5, color='#ff0004')

        # average of the results for one source term 
        celldensitytotal_analyticalapproxpde_average[k-1] = stats.mean(celldensitytotal_analyticalapproxpde)
        celldensitytotal_constant_average[k-1] = stats.mean(celldensitytotal_constant)


    ax.plot(sourceterm, celldensitytotal_constant_average, marker = 'D', markersize = 6, label = 'Global Hypothesis', color='#ff0004')
    ax.plot(sourceterm, celldensitytotal_analyticalapproxpde_average, label = 'Local Hypothesis', marker = '.', markersize = 15.0, color='#8c0000')

    ax.set_xlabel(r'$c_{max}$', fontsize = 26)
    ax.set_ylabel('Cell Density Inside Lesion Compared to Total Network', fontsize = 26)
    ax.set_xticks(sourceterm)
    ax.legend(loc='upper left', fontsize = 24)
    ax.tick_params(axis = 'both', labelsize = 24)
    plt.savefig(main_file_path + "Figures/CellDensityInsideLesionComparedToTotalCells.png")

if GraphErrorNbCellsPlanePositionLesion:
    # plot 
    fig, ax = plt.subplots(figsize = (12,8), dpi = 300, layout='constrained')

    sourceterm = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]

    error_analyticalapproxpde_lower = []
    error_analyticalapproxpde_upper = []
    error_constant_lower = []
    error_constant_upper = []

    # loop over all the files 
    celldensity_analyticalapproxpde_average = np.zeros(10)
    celldensity_constant_average = np.zeros(10)
    for k in range(1,11):
        celldensity_analyticalapproxpde = []
        celldensity_constant = []
        for j in [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]:
            file_path_analyticalapproxpde = main_file_path + "CoupledModel2DAnalyticalApproxPde/CoupledModel2DAnalyticalApproxPdeSeed" + str(j) + "Source" + str(k)
            file_path_constant = main_file_path + "CoupledModel2DConstant/CoupledModel2DConstantSeed" + str(j) + "Source" + str(k)

            file_cellmutation_analyticalapproxpde = file_path_analyticalapproxpde + "/results_from_time_0/results.vizmutationstates"
            file_nodescoordinates_analyticalapproxpde = file_path_analyticalapproxpde + "/results_from_time_0/results.viznodes"
            file_cellmutation_constant = file_path_constant + "/results_from_time_0/results.vizmutationstates"
            file_nodescoordinates_constant = file_path_constant + "/results_from_time_0/results.viznodes"

            celldensity_analyticalapproxpde.append(runner.NbCellsAfterPlane(file_nodescoordinates_analyticalapproxpde, ref_point, dim)/AreaPlane)
            celldensity_constant.append(runner.NbCellsAfterPlane(file_nodescoordinates_constant, ref_point, dim)/AreaPlane)

        # scatter plot 
        ax.scatter([sourceterm[k-1] for i in range(10)], celldensity_analyticalapproxpde, alpha=0.5, color='#8c0000')
        ax.scatter([sourceterm[k-1] for i in range(10)], celldensity_constant, alpha = 0.5, color='#ff0004')

        # average of the results for one source term 
        celldensity_analyticalapproxpde_average[k-1] = stats.mean(celldensity_analyticalapproxpde)
        celldensity_constant_average[k-1] = stats.mean(celldensity_constant)

        low_error_analyticalapproxpde = stats.mean(celldensity_analyticalapproxpde)-max(min(celldensity_analyticalapproxpde), 1e-2)
        high_error_analyticalapproxpde = max(celldensity_analyticalapproxpde)-stats.mean(celldensity_analyticalapproxpde)
        low_error_constant = stats.mean(celldensity_analyticalapproxpde)-max(min(celldensity_constant),1e-2)
        high_error_constant = max(celldensity_constant)-stats.mean(celldensity_constant)

        error_analyticalapproxpde_lower.append(low_error_analyticalapproxpde)
        error_analyticalapproxpde_upper.append(high_error_analyticalapproxpde)
        error_constant_lower.append(low_error_constant)
        error_constant_upper.append(high_error_constant)

    error_analyticalapproxpde = np.array([error_analyticalapproxpde_lower, error_analyticalapproxpde_upper])
    error_constant = np.array([error_constant_lower, error_constant_upper])

    ax.plot(sourceterm, celldensity_constant_average, marker = 'D', markersize = 6, label = 'Global Hypothesis', color='#ff0004')
    ax.plot(sourceterm, celldensity_analyticalapproxpde_average, label = 'Local Hypothesis', marker = '.', markersize = 15.0, color='#8c0000')
    #ax.errorbar(sourceterm, celldensity_constant_average, yerr = error_constant, marker = 'D', markersize = 6, label = 'Global Hypothesis', color='xkcd:brown')
    #ax.errorbar(sourceterm, celldensity_analyticalapproxpde_average, yerr = error_analyticalapproxpde, label = 'Local Hypothesis', marker = '.', markersize = 15.0, color='xkcd:tan')
    ax.set_xlabel(r'$c_{max}$', fontsize = 26)
    ax.set_ylabel('Cell Density Inside Lesion', fontsize = 26)
    ax.set_xticks(sourceterm)
    ax.legend(loc='upper left', fontsize = 24)
    ax.tick_params(axis = 'both', labelsize = 24)
    plt.show()
    plt.savefig(main_file_path + "Figures/CellDensity.png")

# Furthest cell depending on source term for chemotactic force  
if GraphFurthestCell:
    # plot 
    fig, ax = plt.subplots(figsize = (12,8), dpi = 300, layout='constrained')

    sourceterm = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]

    error_analyticalapproxpde_lower = []
    error_analyticalapproxpde_upper = []
    error_constant_lower = []
    error_constant_upper = []

    # loop over all the files 
    furthestcell_analyticalapproxpde_average = np.zeros(10)
    furthestcell_constant_average = np.zeros(10)
    for k in range(1,11):
        furthestcell_analyticalapproxpde = []
        furthestcell_constant = []
        for j in [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]:
            file_path_analyticalapproxpde = main_file_path + "CoupledModel2DAnalyticalApproxPde/CoupledModel2DAnalyticalApproxPdeSeed" + str(j) + "Source" + str(k)
            file_path_constant = main_file_path + "CoupledModel2DConstant/CoupledModel2DConstantSeed" + str(j) + "Source" + str(k)

            file_cellmutation_analyticalapproxpde = file_path_analyticalapproxpde + "/results_from_time_0/results.vizmutationstates"
            file_nodescoordinates_analyticalapproxpde = file_path_analyticalapproxpde + "/results_from_time_0/results.viznodes"
            file_cellmutation_constant = file_path_constant + "/results_from_time_0/results.vizmutationstates"
            file_nodescoordinates_constant = file_path_constant + "/results_from_time_0/results.viznodes"

            furthestcell_analyticalapproxpde.append(runner.NormFurthestCell(file_nodescoordinates_analyticalapproxpde, file_cellmutation_analyticalapproxpde, dim))
            furthestcell_constant.append(runner.NormFurthestCell(file_nodescoordinates_constant, file_cellmutation_constant, dim))

        # scatter plot 
        ax.scatter([sourceterm[k-1] for i in range(10)], furthestcell_analyticalapproxpde, alpha=0.5, color='#8c0000')
        ax.scatter([sourceterm[k-1] for i in range(10)], furthestcell_constant, alpha=0.5, color='#ff0004')

        # average of the results for one source term 
        furthestcell_analyticalapproxpde_average[k-1] = stats.mean(furthestcell_analyticalapproxpde)
        furthestcell_constant_average[k-1] = stats.mean(furthestcell_constant)

        low_error_analyticalapproxpde = stats.mean(furthestcell_analyticalapproxpde)-min(furthestcell_analyticalapproxpde)
        high_error_analyticalapproxpde = max(furthestcell_analyticalapproxpde)-stats.mean(furthestcell_analyticalapproxpde)
        low_error_constant = stats.mean(furthestcell_analyticalapproxpde)-min(furthestcell_constant)
        high_error_constant = max(furthestcell_constant)-stats.mean(furthestcell_constant)

        error_analyticalapproxpde_lower.append(low_error_analyticalapproxpde)
        error_analyticalapproxpde_upper.append(high_error_analyticalapproxpde)
        error_constant_lower.append(low_error_constant)
        error_constant_upper.append(high_error_constant)

    error_analyticalapproxpde = np.array([error_analyticalapproxpde_lower, error_analyticalapproxpde_upper])
    error_constant = np.array([error_constant_lower, error_constant_upper])

    ax.plot(sourceterm, furthestcell_constant_average, marker = 'D', markersize = 6, label = 'Global Hypothesis', color='#ff0004')
    ax.plot(sourceterm, furthestcell_analyticalapproxpde_average, label = 'Local Hypothesis', marker = '.', markersize = 15.0, color='#8c0000')
    #ax.errorbar(sourceterm, furthestcell_constant_average, yerr = error_constant, marker = 'D', markersize = 6, label = 'Global Hypothesis', color='xkcd:brown')
    #ax.errorbar(sourceterm, furthestcell_analyticalapproxpde_average, yerr = error_analyticalapproxpde, label = 'Local Hypothesis', marker = '.', markersize = 15.0, color='xkcd:tan')
    ax.legend(loc='upper right', fontsize = 24)
    ax.set_xlabel(r'$c_{max}$', fontsize = 26)
    ax.set_ylabel('Network Length', fontsize = 26)
    ax.set_xticks(sourceterm)
    ax.tick_params(axis = 'both', labelsize = 24)
    plt.show()
    plt.savefig(main_file_path + "Figures/FurthestCell.png")

if GraphFirstTimeCellReachingCentreLesion:
    # plot 
    fig, ax = plt.subplots(figsize = (12,8), dpi = 300, layout='constrained')

    sourceterm = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]

    error_analyticalapproxpde_lower = []
    error_analyticalapproxpde_upper = []
    error_constant_lower = []
    error_constant_upper = []

    # loop over all the files 
    timefirstreachingcentrelesion_analyticalapproxpde_average = np.zeros(10)
    timefirstreachingcentrelesion_constant_average = np.zeros(10)
    for k in range(1,11):
        timefirstreachingcentrelesion_analyticalapproxpde = []
        timefirstreachingcentrelesion_constant = []
        for j in [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]:
            file_path_analyticalapproxpde = main_file_path + "CoupledModel2DAnalyticalApproxPde/CoupledModel2DAnalyticalApproxPdeSeed" + str(j) + "Source" + str(k)
            file_path_constant = main_file_path + "CoupledModel2DConstant/CoupledModel2DConstantSeed" + str(j) + "Source" + str(k)

            file_cellmutation_analyticalapproxpde = file_path_analyticalapproxpde + "/results_from_time_0/results.vizmutationstates"
            file_nodescoordinates_analyticalapproxpde = file_path_analyticalapproxpde + "/results_from_time_0/results.viznodes"
            file_cellmutation_constant = file_path_constant + "/results_from_time_0/results.vizmutationstates"
            file_nodescoordinates_constant = file_path_constant + "/results_from_time_0/results.viznodes"

            timefirstreachingcentrelesion_analyticalapproxpde.append(runner.TimeFirstReachingPlane(file_nodescoordinates_analyticalapproxpde, ref_point_centre_lesion, dim))
            timefirstreachingcentrelesion_constant.append(runner.TimeFirstReachingPlane(file_nodescoordinates_constant, ref_point_centre_lesion, dim))

        # scatter plot 
        ax.scatter([sourceterm[k-1] for i in range(10)], timefirstreachingcentrelesion_analyticalapproxpde, alpha=0.5, color='#8c0000')
        ax.scatter([sourceterm[k-1] for i in range(10)], timefirstreachingcentrelesion_constant, alpha=0.5, color='#ff0004')

        # average of the results for one source term 
        timefirstreachingcentrelesion_analyticalapproxpde_average[k-1] = stats.mean(timefirstreachingcentrelesion_analyticalapproxpde)
        timefirstreachingcentrelesion_constant_average[k-1] = stats.mean(timefirstreachingcentrelesion_constant)

        low_error_analyticalapproxpde = stats.mean(timefirstreachingcentrelesion_analyticalapproxpde)-max(min(timefirstreachingcentrelesion_analyticalapproxpde), 1e-9)
        high_error_analyticalapproxpde = max(timefirstreachingcentrelesion_analyticalapproxpde)-stats.mean(timefirstreachingcentrelesion_analyticalapproxpde)
        low_error_constant = stats.mean(timefirstreachingcentrelesion_analyticalapproxpde)-max(min(timefirstreachingcentrelesion_constant),1e-9)
        high_error_constant = max(timefirstreachingcentrelesion_constant)-stats.mean(timefirstreachingcentrelesion_constant)

        error_analyticalapproxpde_lower.append(low_error_analyticalapproxpde)
        error_analyticalapproxpde_upper.append(high_error_analyticalapproxpde)
        error_constant_lower.append(low_error_constant)
        error_constant_upper.append(high_error_constant)

    error_analyticalapproxpde = np.array([error_analyticalapproxpde_lower, error_analyticalapproxpde_upper])
    error_constant = np.array([error_constant_lower, error_constant_upper])

    ax.plot(sourceterm, timefirstreachingcentrelesion_constant_average, marker = 'D', markersize = 6, label = 'Global Hypothesis', color='#ff0004')
    ax.plot(sourceterm, timefirstreachingcentrelesion_analyticalapproxpde_average, label = 'Local Hypothesis', marker = '.', markersize = 15.0, color='#8c0000')
    #ax.errorbar(sourceterm, timefirstreachingcentrelesion_constant_average, yerr = error_constant, marker = 'D', markersize = 6, label = 'Global Hypothesis', color='xkcd:brown')
    #ax.errorbar(sourceterm, timefirstreachingcentrelesion_analyticalapproxpde_average, yerr = error_analyticalapproxpde, label = 'Local Hypothesis', marker = '.', markersize = 15.0, color='xkcd:tan')
    ax.legend(loc='upper right', fontsize = 24)
    ax.set_xlabel(r'$c_{max}$', fontsize = 26)
    ax.set_ylabel('Time to Reach Centre of Lesion', fontsize = 26)
    ax.set_xticks(sourceterm)
    ax.tick_params(axis = 'both', labelsize = 26)
    plt.show()
    plt.savefig(main_file_path + "Figures/TimeFirstCellReachingCentreLesion.png")

if GraphFirstTimeCellReachingLesion:
    # plot 
    fig, ax = plt.subplots(figsize = (12,8), dpi = 300, layout='constrained')

    sourceterm = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]

    error_analyticalapproxpde_lower = []
    error_analyticalapproxpde_upper = []
    error_constant_lower = []
    error_constant_upper = []

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

            # check that all simulations reach the lesion 
            if(runner.TimeFirstReachingPlane(file_nodescoordinates_analyticalapproxpde, ref_point, dim) > 1900):
                print("Approx PDE")
                print(k)
                print(j)
            if(runner.TimeFirstReachingPlane(file_nodescoordinates_constant, ref_point, dim) > 1900):
                print("Constant")
                print(k)
                print(j)

        # scatter plot 
        ax.scatter([sourceterm[k-1] for i in range(10)], timefirstreachinglesion_analyticalapproxpde, alpha=0.5, color='#8c0000')
        ax.scatter([sourceterm[k-1] for i in range(10)], timefirstreachinglesion_constant, alpha=0.5, color='#ff0004')

        # average of the results for one source term 
        timefirstreachinglesion_analyticalapproxpde_average[k-1] = stats.mean(timefirstreachinglesion_analyticalapproxpde)
        timefirstreachinglesion_constant_average[k-1] = stats.mean(timefirstreachinglesion_constant)

    ax.plot(sourceterm, timefirstreachinglesion_constant_average, marker = 'D', markersize = 6, label = 'Global Hypothesis', color='#ff0004')
    ax.plot(sourceterm, timefirstreachinglesion_analyticalapproxpde_average, label = 'Local Hypothesis', marker = '.', markersize = 15.0, color='#8c0000')
 
    ax.legend(loc='upper left', fontsize = 24)
    ax.set_xlabel(r'$c_{max}$', fontsize = 26)
    ax.set_ylabel('Time to Reach Lesion', fontsize = 26)
    ax.set_xticks(sourceterm)
    ax.tick_params(axis = 'both', labelsize = 24)
    plt.show()
    plt.savefig(main_file_path + "Figures/TimeFirstCellReachingLesion.png")

if GraphTimeReachingFurthestCell:
    # new function finding the furthest cell and keeping the time it reached the furthest point 
    # scatter plot + attributing colour depending on the source (put this in legend)

    # plot 
    fig, ax = plt.subplots(figsize = (12,8), dpi = 300, layout='constrained')

    sourceterm = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    t_analyticalapproxpde = []
    normfurthestcell_analyticalapproxpde = []
    t_constant = []
    normfurthestcell_constant = []

    for k in range(1,11):
        for j in [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]:
            file_path_analyticalapproxpde = main_file_path + "CoupledModel2DAnalyticalApproxPde/CoupledModel2DAnalyticalApproxPdeSeed" + str(j) + "Source" + str(k)
            file_path_constant = main_file_path + "CoupledModel2DConstant/CoupledModel2DConstantSeed" + str(j) + "Source" + str(k)

            file_cellmutation_analyticalapproxpde = file_path_analyticalapproxpde + "/results_from_time_0/results.vizmutationstates"
            file_nodescoordinates_analyticalapproxpde = file_path_analyticalapproxpde + "/results_from_time_0/results.viznodes"
            file_cellmutation_constant = file_path_constant + "/results_from_time_0/results.vizmutationstates"
            file_nodescoordinates_constant = file_path_constant + "/results_from_time_0/results.viznodes"

            [valuenormfurthestcell_analyticalapproxpde, valuet_analyticalapproxpde] = runner.NormAndTimeFurthestCell(file_nodescoordinates_analyticalapproxpde, file_cellmutation_analyticalapproxpde, dim)
            [valuenormfurthestcell_constant, valuet_constant] = runner.NormAndTimeFurthestCell(file_nodescoordinates_constant, file_cellmutation_constant, dim)

            normfurthestcell_analyticalapproxpde.append(valuenormfurthestcell_analyticalapproxpde)
            t_analyticalapproxpde.append(valuet_analyticalapproxpde)
            normfurthestcell_constant.append(valuenormfurthestcell_constant)
            t_constant.append(valuet_constant)

    # scatter plot 
    ax.scatter(t_analyticalapproxpde, normfurthestcell_analyticalapproxpde, alpha=0.5, color='#8c0000')
    ax.scatter(t_constant, normfurthestcell_constant, alpha=0.5, color='#ff0004')

    #ax.legend(loc='upper right', fontsize = 24)
    ax.set_xlabel(r'$t$', fontsize = 26)
    ax.set_ylabel('Norm Furthest Cell', fontsize = 26)
    ax.tick_params(axis = 'both', labelsize = 24)
    plt.show()
    plt.savefig(main_file_path + "Figures/TimeReachingFurthestCell.png")

if GraphTimeNormFirstTipCell:
    # plot 
    fig, ax = plt.subplots(figsize = (12,8), dpi = 300, layout='constrained')

    sourceterm = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    t = np.linspace(0, Tend, 21)

    normfirsttipcell_analyticalapproxpde = []
    normfirsttipcell_constant = []

    for k in range(1,11):
        for j in [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]:
            file_path_analyticalapproxpde = main_file_path + "CoupledModel2DAnalyticalApproxPde/CoupledModel2DAnalyticalApproxPdeSeed" + str(j) + "Source" + str(k)
            file_path_constant = main_file_path + "CoupledModel2DConstant/CoupledModel2DConstantSeed" + str(j) + "Source" + str(k)

            file_cellmutation_analyticalapproxpde = file_path_analyticalapproxpde + "/results_from_time_0/results.vizmutationstates"
            file_nodescoordinates_analyticalapproxpde = file_path_analyticalapproxpde + "/results_from_time_0/results.viznodes"
            file_cellmutation_constant = file_path_constant + "/results_from_time_0/results.vizmutationstates"
            file_nodescoordinates_constant = file_path_constant + "/results_from_time_0/results.viznodes"

            valuenormfirsttipcell_analyticalapproxpde = runner.NormFirstTipCell(file_nodescoordinates_analyticalapproxpde, file_cellmutation_analyticalapproxpde, dim)
            valuenormfirsttipcell_constant = runner.NormFirstTipCell(file_nodescoordinates_constant, file_cellmutation_constant, dim)

            # scatter plot 
            ax.plot(t, valuenormfirsttipcell_analyticalapproxpde, color='#8c0000')
            ax.plot(t, valuenormfirsttipcell_constant, color='#ff0004')

    #ax.legend(loc='upper right', fontsize = 24)
    ax.set_xlabel(r'$t$', fontsize = 26)
    ax.set_ylabel('Norm First Tip Cell', fontsize = 26)
    ax.tick_params(axis = 'both', labelsize = 24)
    plt.show()
    plt.savefig(main_file_path + "Figures/TimeReachingFirstTipCell.png")