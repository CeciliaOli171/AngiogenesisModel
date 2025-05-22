## Libraries 
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import statistics as stats

from matplotlib import cm
from matplotlib.ticker import LinearLocator 
from matplotlib.colors import LinearSegmentedColormap

from scipy.spatial import ConvexHull

import os, sys
sys.path.insert(1, "/".join(os.path.realpath(__file__).split("/")[0:-2]))
from ParametersSensitivityRunner import ParametersSensitivityRunner

runner = ParametersSensitivityRunner

## Parameters 
Tcycle = 13.0
time = 4000
Tend = 4000
TotalTestSourceNb = 10
TotalTestNb = 17
dim = 2
ref_point = 15

# VEGF concentration parameters
l_max = 220
c_0 = 0.1
Kc = 2e-2
lambdaSprout = 0.024

# Psprout parameters 
cmax = 0.8
cmin = 0.3
Pmax = 0.98
Pmin = 0.5

n = (1/np.log(cmax/cmin))*np.log((Pmax/Pmin)*(1-Pmin)/(1-Pmax))
K = cmax*((1-Pmax)/Pmax)**(1/n)

main_file_path = "/Users/coli171/Chaste/Output/PaperAngiogenesisModel2025/PaperModel2025Analysis2D/" # 2D
#main_file_path = "/Users/coli171/Chaste/Output/PaperAngiogenesisModel2025/PaperModel2025Analysis3D/" # 3D

## Settings
TestBaseline = False
GraphVEGFConcentration = False
GraphVEGFConcentrationKc = False
GraphVEGFGradient = False 
GraphVEGFGradientXAxis = False
GraphPsprout = False
GraphPsproutSourceTerm = False
GraphPsproutPosition = False
GraphNbBranches = False
GraphExpectedNbBranches = False
GraphExpectedNbBranchesVEGF = True
GraphBarNbBranches = False
GraphFurthestCell = False
GraphConvergenceTime = False
GraphNbCellsPlane = False
GraphBarCellsInPlane = False
GraphFirstTimeCellReachingLesion = False

## Plots

if GraphVEGFConcentration:
    fig = plt.subplots(figsize = (12,8), dpi = 300)
    x = np.linspace(-200,200,1000)
    y1 = [(0.1 - c_0)*np.exp(-Kc*np.abs(item))+c_0 for item in x]
    y2 = [(0.2 - c_0)*np.exp(-Kc*np.abs(item))+c_0 for item in x]
    y5 = [(0.5 - c_0)*np.exp(-Kc*np.abs(item))+c_0 for item in x]
    y10 = [(1.0 - c_0)*np.exp(-Kc*np.abs(item))+c_0 for item in x]
    plt.plot(x, y1, label= r'$c_{max} = 0.1$', color = 'xkcd:pale blue')
    plt.plot(x, y2, label= r'$c_{max} = 0.2$', color = 'xkcd:lightblue')
    plt.plot(x, y5, label= r'$c_{max} = 0.5$', color = 'xkcd:blue')
    plt.plot(x, y10, label= r'$c_{max} = 1.0$', color = 'xkcd:royal blue')
    plt.xlabel(r'$x$', fontsize=18) 
    plt.ylabel(r'$c(x)$', fontsize=18)
    plt.xticks(fontsize = 16)
    plt.yticks(fontsize = 16)
    plt.legend(fontsize = 12)
    plt.show()

if GraphVEGFConcentrationKc:
    fig = plt.subplots(figsize = (12,8), dpi = 300)
    x = np.linspace(0,200,1000)
    C = 1.0
    Kc1 = 1e-2
    Kc2 = 2e-2
    Kc5 = 5e-2
    Kc10 = 1e-1
    y1 = [(C - c_0)*np.exp(-Kc1*item)+c_0 for item in x]
    y2 = [(C - c_0)*np.exp(-Kc2*item)+c_0 for item in x]
    y5 = [(C - c_0)*np.exp(-Kc5*item)+c_0 for item in x]
    y10 = [(C - c_0)*np.exp(-Kc10*item)+c_0 for item in x]
    plt.plot(x, y1, label= r'$K_c = 1e-2$', color = 'xkcd:pale blue')
    plt.plot(x, y2, label= r'$K_c = 2e-2$', color = 'xkcd:lightblue')
    plt.plot(x, y5, label= r'$K_c = 5e-2$', color = 'xkcd:blue')
    plt.plot(x, y10, label= r'$K_c = 1e-1$', color = 'xkcd:royal blue')
    plt.xlabel(r'$x$', fontsize=18) 
    plt.ylabel(r'$c(x)$', fontsize=18)
    plt.xticks(fontsize = 16)
    plt.yticks([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0], fontsize = 16)
    plt.legend(fontsize = 12)
    plt.show()

if GraphVEGFGradient:
    fig = plt.subplots(figsize = (12,8), dpi = 300)
    x = np.linspace(0,300,1000)
    y1 = [0.1*np.sqrt(10-8.3)/np.sqrt(200)*(np.exp(-np.sqrt(10-8.3)/np.sqrt(200)*item)+np.exp(np.sqrt(10-8.3)/np.sqrt(200)*(item-2*100)))/(1-np.exp(-2*100*np.sqrt(10-8.3)/np.sqrt(200))) for item in x]
    y5 = [0.5*np.sqrt(10-8.3)/np.sqrt(200)*(np.exp(-np.sqrt(10-8.3)/np.sqrt(200)*item)+np.exp(np.sqrt(10-8.3)/np.sqrt(200)*(item-2*100)))/(1-np.exp(-2*100*np.sqrt(10-8.3)/np.sqrt(200))) for item in x]
    y10 = [1.0*np.sqrt(10-8.3)/np.sqrt(200)*(np.exp(-np.sqrt(10-8.3)/np.sqrt(200)*item)+np.exp(np.sqrt(10-8.3)/np.sqrt(200)*(item-2*100)))/(1-np.exp(-2*100*np.sqrt(10-8.3)/np.sqrt(200))) for item in x]
    plt.plot(x, y1, label= r'$c_{max} = 0.1$', color = 'xkcd:lightblue')
    plt.plot(x, y5, label=r'$c_{max} = 0.5$', color='xkcd:blue')
    plt.plot(x, y10, label=r'$c_{max} = 1.0$', color='xkcd:royal blue')
    plt.ylabel(r'$x$', fontsize=18) 
    plt.xlabel(r'$|\nabla c(x)|$', fontsize=18)
    plt.xticks(fontsize = 16)
    plt.yticks(fontsize = 16)
    plt.legend(fontsize = 12)
    plt.show()

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
    plt.plot(sourceterm, Psprout, label='hill function', color='xkcd:brown')
    plt.xlabel(r'$c$', fontsize=18) 
    plt.ylabel(r'$P_{sprout}$', fontsize=18)
    plt.xticks(sourcetermticks, fontsize = 16)
    plt.yticks([0.0, 0.01, 0.02, 0.03, 0.04], fontsize = 16)
    plt.legend(fontsize = 12)
    plt.show()

if GraphPsproutPosition:
    fig = plt.subplots(figsize = (12,8), dpi = 300)
    x = np.linspace(0, 300, 1000)
    Psprout1 = [0.02*0.1**n/(K**n + 0.1**n) for item in x]
    Psprout2 = [0.02*((0.2 - c_0)*np.exp(-Kc*item)+c_0)**n/(K**n + ((0.1 - c_0)*np.exp(-Kc*item)+c_0)**n) for item in x]
    Psprout5 = [0.02*((0.5 - c_0)*np.exp(-Kc*item)+c_0)**n/(K**n + ((0.5 - c_0)*np.exp(-Kc*item)+c_0)**n) for item in x]
    Psprout10 = [0.02*((1.0 - c_0)*np.exp(-Kc*item)+c_0)**n/(K**n + ((1.0 - c_0)*np.exp(-Kc*item)+c_0)**n) for item in x]
    #plt.plot(x, Psprout1, label=r'$c_{max} = 0.1$', color='xkcd:beige')
    plt.plot(x, Psprout2, label=r'$c_{max} = 0.2$', color='xkcd:tan')
    plt.plot(x, Psprout5, label=r'$c_{max} = 0.5$', color='xkcd:light brown')
    plt.plot(x, Psprout10, label=r'$c_{max} = 1.0$', color='xkcd:brown')
    plt.xlabel(r'$x$', fontsize=18) 
    plt.ylabel(r'$P_{sprout}$', fontsize=18)
    plt.xticks(fontsize = 16)
    plt.yticks([0.0, 0.01, 0.02, 0.03, 0.04], fontsize = 16)
    plt.legend(fontsize = 12)
    plt.show()

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
        for j in range(4, 21):
            file_path_analyticalapproxpde = main_file_path + "CoupledModel2DAnalyticalApproxPde/CoupledModel2DAnalyticalApproxPdeSeed" + str(j) + "Source" + str(k)
            file_path_constant = main_file_path + "CoupledModel2DConstant/CoupledModel2DConstantSeed" + str(j) + "Source" + str(k)

            file_cellmutation_analyticalapproxpde = file_path_analyticalapproxpde + "/results_from_time_0/results.vizmutationstates"
            file_cellmutation_constant = file_path_constant + "/results_from_time_0/results.vizmutationstates"

            totalnumberbranches_analyticalapproxpde.append(runner.TotalNumberBranches(file_cellmutation_analyticalapproxpde))
            totalnumberbranches_constant.append(runner.TotalNumberBranches(file_cellmutation_constant))

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

    print(error_analyticalapproxpde)
    print(error_constant)

    # plot 
    fig, ax = plt.subplots(figsize = (12,8), dpi = 300, layout='constrained')
 
    # ax.plot(sourceterm, totalnumberbranches_constant_average, marker = 'D', markersize = 6, label = 'Constant VEGF', color='xkcd:brown')
    # ax.plot(sourceterm, totalnumberbranches_analyticalapproxpde_average, label = 'Lesion Up-regulation', marker = '.', markersize = 15.0, color='xkcd:tan')
    ax.errorbar(sourceterm, totalnumberbranches_constant_average, yerr = error_constant, marker = 'D', markersize = 6, label = 'Constant VEGF', color='xkcd:brown')
    ax.errorbar(sourceterm, totalnumberbranches_analyticalapproxpde_average, yerr = error_analyticalapproxpde, label = 'Lesion Up-regulation', marker = '.', markersize = 15.0, color='xkcd:tan')
    ax.legend(loc='upper left', fontsize = 8)

    ax.set_xlabel(r'$c_{max}$', fontsize = 15)
    ax.set_ylabel('Number of Branches', fontsize = 15)
    ax.set_xticks(sourceterm)
    ax.tick_params(axis = 'both', labelsize = 15)
    plt.show()

def ExpectationTC(t, c_max, Kc):
    # case of Psprout depending on VEGF
    # Sk = (1/(0.85*t))*(0.02*Tcycle/(t*Kc))*(np.log(K**n + ((c_max-c_0)*np.exp(-t*Kc*(220-t/Tcycle)/Tcycle))**n)-np.log(K**n + ((c_max-c_0)*np.exp(-t*Kc*(220-0.15*t/Tcycle)/Tcycle))**n)) # case of integral: NOT WORKING
    Sk = lambdaSprout*((c_max - c_0)*np.exp(-Kc*(l_max-t/Tcycle))+c_0)**n/(K**n + ((c_max - c_0)*np.exp(-Kc*(l_max-t/Tcycle))+c_0)**n) # Hill function

    # constant case of Psprout = 0.02
    #Sk = 0.02
    return Sk 

def ExpectationBranches(t, c_max):
    Bn = 2*(1+ExpectationTC(t, c_max))**n - 1
    return Bn

if GraphExpectedNbBranches:
    fig = plt.subplots(figsize = (12,8), dpi = 300)

    x = np.linspace(0, l_max, 1000)
    t = np.linspace(0, 3000, 50)
    sourceterm = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]

    y = [2*(1+0.02)**(item/Tcycle)-1 for item in t]
    y2 = [2*(1+ExpectationTC(item, 0.2))**(item/Tcycle)-1 for item in t]
    y5 = [2*(1+ExpectationTC(item, 0.5))**(item/Tcycle)-1 for item in t]
    y10 = [2*(1+ExpectationTC(item, 1.0))**(item/Tcycle)-1 for item in t]

    # plot just Sk
    # y = [0.02 for item in t]
    # y2 = [ExpectationTC(item, 0.2) for item in t]
    # y5 = [ExpectationTC(item, 0.5) for item in t]
    # y10 = [ExpectationTC(item, 1.0) for item in t]

    plt.plot(t, y, label = 'constant background', color = 'xkcd:black')
    plt.plot(t, y2, label = 'cmax = 0.2', color = 'xkcd:yellow', marker = 'x')
    plt.plot(t, y5, label = 'cmax = 0.5', color = 'xkcd:orange', marker = 'x')
    plt.plot(t, y10, label = 'cmax = 1.0', color = 'xkcd:red', marker = 'x')
    plt.xlabel(r'$t$', fontsize=18) 
    plt.ylabel('Expected Number of Branches', fontsize=18)
    plt.xticks(fontsize = 16)
    plt.yticks(fontsize = 16)
    plt.legend(fontsize = 12)
    plt.show()

if GraphExpectedNbBranchesVEGF:
    fig = plt.subplots(figsize = (12,8), dpi = 300)

    x = np.linspace(0, l_max, 1000)
    t = np.linspace(0, 3000, 50)
    sourceterm = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]

    # data from Zhang et al. for comparaison 
    datapoints = [[0.7, 107], [1.77, 105], [2.05, 90], [1.3, 92], [2.25, 52], [1.8, 70], [1.13, 15], [2.3, 75], [2.25, 67], [1.95, 47], [1.7, 68], [1.25, 45], [1.3, 42], [0.48, 48], [0.25, 52], [0.49, 20], [0.4, 29], [1.0, 23], [1.1, 28], [0.9, 10], [0.6, 12], [0.85, 98], [0.7, 83], [0.8, 82], [0.68, 21], [0.27, 28], [0.98, 28], [0.47, 6], [0.4, 14], [0.32, 16], [0.9, 78], [0.8, 72], [0.92, 68], [0.77, 61], [0.94, 57], [0.7, 42], [0.58, 42], [0.57, 41]]
    microvesseldensity = [datapoints[k][0]/2.5 for k in range(len(datapoints))]
    VEGFgeneexpressionlinearregression = [20.34*4*item*2.5+30.90*4 for item in microvesseldensity]
    VEGFgeneexpression = [datapoints[k][1]*4 for k in range(len(datapoints))]

    nbbranches12 = [2*(1+ExpectationTC(3000, v, 1e-2))**(3000/Tcycle)-1 for v in sourceterm]
    nbbranches22 = [2*(1+ExpectationTC(3000, v, 2e-2))**(3000/Tcycle)-1 for v in sourceterm]
    nbbranches42 = [2*(1+ExpectationTC(3000, v, 4e-2))**(3000/Tcycle)-1 for v in sourceterm]
    nbbranches11 = [2*(1+ExpectationTC(3000, v, 1e-1))**(3000/Tcycle)-1 for v in sourceterm]

    plt.plot(sourceterm, nbbranches12, color='xkcd:pale blue', label = 'Kc = 1e-2')
    plt.plot(sourceterm, nbbranches22, color='xkcd:lightblue', label = 'Kc = 2e-2')
    plt.plot(sourceterm, nbbranches42, color='xkcd:blue', label = 'Kc = 4e-2')
    plt.plot(sourceterm, nbbranches11, color='xkcd:royal blue', label = 'Kc = 1e-1')
    plt.scatter(microvesseldensity, VEGFgeneexpression, color = 'g')
    plt.plot(microvesseldensity, VEGFgeneexpressionlinearregression, color = 'g')
    plt.xlabel(r'$t$', fontsize=18) 
    #plt.ylabel('Expected Number of Branches', fontsize=18)
    plt.xticks(fontsize = 16)
    plt.yticks(fontsize = 16)
    plt.legend(fontsize = 12)
    plt.show()

# Number of cells after a certain plane depending on source term for different Psprout 
if GraphNbCellsPlane:
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
        for j in range(4,21):
            file_path_analyticalapproxpde = main_file_path + "CoupledModel2DAnalyticalApproxPde/CoupledModel2DAnalyticalApproxPdeSeed" + str(j) + "Source" + str(k)
            file_path_constant = main_file_path + "CoupledModel2DConstant/CoupledModel2DConstantSeed" + str(j) + "Source" + str(k)

            file_cellmutation_analyticalapproxpde = file_path_analyticalapproxpde + "/results_from_time_0/results.vizmutationstates"
            file_nodescoordinates_analyticalapproxpde = file_path_analyticalapproxpde + "/results_from_time_0/results.viznodes"
            file_cellmutation_constant = file_path_constant + "/results_from_time_0/results.vizmutationstates"
            file_nodescoordinates_constant = file_path_constant + "/results_from_time_0/results.viznodes"

            cellsafterplane_analyticalapproxpde.append(runner.NbCellsAfterPlane(file_nodescoordinates_analyticalapproxpde, ref_point, dim))
            cellsafterplane_constant.append(runner.NbCellsAfterPlane(file_nodescoordinates_constant, ref_point, dim))

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

    print(error_analyticalapproxpde)
    print(error_constant)

    # plot 
    fig, ax = plt.subplots(figsize = (12,8), dpi = 300, layout='constrained')

    # ax.plot(sourceterm, cellsafterplane_constant_average, marker = 'D', markersize = 6, label = 'Constant VEGF', color='xkcd:brown')
    # ax.plot(sourceterm, cellsafterplane_analyticalapproxpde_average, label = 'Lesion Up-regulation', marker = '.', markersize = 15.0, color='xkcd:tan')
    ax.errorbar(sourceterm, cellsafterplane_constant_average, yerr = error_constant, marker = 'D', markersize = 6, label = 'Constant VEGF', color='xkcd:brown')
    ax.errorbar(sourceterm, cellsafterplane_analyticalapproxpde_average, yerr = error_analyticalapproxpde, label = 'Lesion Up-regulation', marker = '.', markersize = 15.0, color='xkcd:tan')
    ax.set_xlabel(r'$c_{max}$', fontsize = 15)
    ax.set_ylabel('Cells In Lesion', fontsize = 15)
    ax.set_xticks(sourceterm)
    ax.legend(loc='upper left', fontsize = 8)
    ax.tick_params(axis = 'both', labelsize = 15)
    plt.show()

# Furthest cell depending on source term for chemotactic force  
if GraphFurthestCell:
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
        for j in range(4,21):
            file_path_analyticalapproxpde = main_file_path + "CoupledModel2DAnalyticalApproxPde/CoupledModel2DAnalyticalApproxPdeSeed" + str(j) + "Source" + str(k)
            file_path_constant = main_file_path + "CoupledModel2DConstant/CoupledModel2DConstantSeed" + str(j) + "Source" + str(k)

            file_cellmutation_analyticalapproxpde = file_path_analyticalapproxpde + "/results_from_time_0/results.vizmutationstates"
            file_nodescoordinates_analyticalapproxpde = file_path_analyticalapproxpde + "/results_from_time_0/results.viznodes"
            file_cellmutation_constant = file_path_constant + "/results_from_time_0/results.vizmutationstates"
            file_nodescoordinates_constant = file_path_constant + "/results_from_time_0/results.viznodes"

            furthestcell_analyticalapproxpde.append(runner.NormFurthestCell(file_nodescoordinates_analyticalapproxpde, file_cellmutation_analyticalapproxpde, dim))
            furthestcell_constant.append(runner.NormFurthestCell(file_nodescoordinates_constant, file_cellmutation_constant, dim))

        # average of the results for one source term 
        print(furthestcell_constant)
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

    print(error_analyticalapproxpde)
    print(error_constant)

    # plot 
    fig, ax = plt.subplots(figsize = (12,8), dpi = 300, layout='constrained')

    # ax.plot(sourceterm, furthestcell_constant_average, marker = 'D', markersize = 6, label = 'Constant VEGF', color='xkcd:brown')
    # ax.plot(sourceterm, furthestcell_analyticalapproxpde_average, label = 'Lesion Up-regulation', marker = '.', markersize = 15.0, color='xkcd:tan')
    ax.errorbar(sourceterm, furthestcell_constant_average, yerr = error_constant, marker = 'D', markersize = 6, label = 'Constant VEGF', color='xkcd:brown')
    ax.errorbar(sourceterm, furthestcell_analyticalapproxpde_average, yerr = error_analyticalapproxpde, label = 'Lesion Up-regulation', marker = '.', markersize = 15.0, color='xkcd:tan')
    ax.legend(loc='upper right', fontsize = 8)
    ax.set_xlabel(r'$c_{max}$', fontsize = 15)
    ax.set_ylabel('Network Length', fontsize = 15)
    ax.set_xticks(sourceterm)
    ax.tick_params(axis = 'both', labelsize = 15)
    plt.show()

if GraphFirstTimeCellReachingLesion:
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
        for j in range(4,21):
            file_path_analyticalapproxpde = main_file_path + "CoupledModel2DAnalyticalApproxPde/CoupledModel2DAnalyticalApproxPdeSeed" + str(j) + "Source" + str(k)
            file_path_constant = main_file_path + "CoupledModel2DConstant/CoupledModel2DConstantSeed" + str(j) + "Source" + str(k)

            file_cellmutation_analyticalapproxpde = file_path_analyticalapproxpde + "/results_from_time_0/results.vizmutationstates"
            file_nodescoordinates_analyticalapproxpde = file_path_analyticalapproxpde + "/results_from_time_0/results.viznodes"
            file_cellmutation_constant = file_path_constant + "/results_from_time_0/results.vizmutationstates"
            file_nodescoordinates_constant = file_path_constant + "/results_from_time_0/results.viznodes"

            timefirstreachinglesion_analyticalapproxpde.append(runner.TimeFirstReachingPlane(file_nodescoordinates_analyticalapproxpde, ref_point, dim))
            timefirstreachinglesion_constant.append(runner.TimeFirstReachingPlane(file_nodescoordinates_constant, ref_point, dim))

        # average of the results for one source term 
        timefirstreachinglesion_analyticalapproxpde_average[k-1] = stats.mean(timefirstreachinglesion_analyticalapproxpde)
        timefirstreachinglesion_constant_average[k-1] = stats.mean(timefirstreachinglesion_constant)

        low_error_analyticalapproxpde = stats.mean(timefirstreachinglesion_analyticalapproxpde)-max(min(timefirstreachinglesion_analyticalapproxpde), 1e-9)
        high_error_analyticalapproxpde = max(timefirstreachinglesion_analyticalapproxpde)-stats.mean(timefirstreachinglesion_analyticalapproxpde)
        low_error_constant = stats.mean(timefirstreachinglesion_analyticalapproxpde)-max(min(timefirstreachinglesion_constant),1e-9)
        high_error_constant = max(timefirstreachinglesion_constant)-stats.mean(timefirstreachinglesion_constant)

        error_analyticalapproxpde_lower.append(low_error_analyticalapproxpde)
        error_analyticalapproxpde_upper.append(high_error_analyticalapproxpde)
        error_constant_lower.append(low_error_constant)
        error_constant_upper.append(high_error_constant)

    error_analyticalapproxpde = np.array([error_analyticalapproxpde_lower, error_analyticalapproxpde_upper])
    error_constant = np.array([error_constant_lower, error_constant_upper])

    # plot 
    fig, ax = plt.subplots(figsize = (12,8), dpi = 300, layout='constrained')

    # ax.plot(sourceterm, timefirstreachinglesion_constant_average, marker = 'D', markersize = 6, label = 'Constant VEGF', color='xkcd:brown')
    # ax.plot(sourceterm, timefirstreachinglesion_analyticalapproxpde_average, label = 'Lesion Up-regulation', marker = '.', markersize = 15.0, color='xkcd:tan')
    ax.errorbar(sourceterm, timefirstreachinglesion_constant_average, yerr = error_constant, marker = 'D', markersize = 6, label = 'Constant VEGF', color='xkcd:brown')
    ax.errorbar(sourceterm, timefirstreachinglesion_analyticalapproxpde_average, yerr = error_analyticalapproxpde, label = 'Lesion Up-regulation', marker = '.', markersize = 15.0, color='xkcd:tan')
    ax.legend(loc='upper right', fontsize = 8)
    ax.set_xlabel(r'$c_{max}$', fontsize = 15)
    ax.set_ylabel('Time to Reach Lesion', fontsize = 15)
    ax.set_xticks(sourceterm)
    ax.tick_params(axis = 'both', labelsize = 15)
    plt.show()