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
time = 300
Tend = 300
TotalTestSourceNb = 10
TotalTestNb = 20
dim = 2
ref_point = 15

# Psprout parameters 
cmax = 0.8
c0 = 0.5
Pmax = 0.9
Pmin = 0.5

n = (1/np.log(cmax/c0))*np.log((Pmax/Pmin)*(1-Pmin)/(1-Pmax))
K = cmax*((1-Pmax)/Pmax)**(1/n)

main_file_path = "/Users/coli171/Chaste/Output/PaperAngiogenesisModel2025/PaperModel2025Analysis2D/" # 2D
#main_file_path = "/Users/coli171/Chaste/Output/PaperAngiogenesisModel2025/PaperModel2025Analysis3D/" # 3D

## Settings
TestBaseline = False
GraphVEGFConcentration = False
GraphVEGFGradient = False 
GraphVEGFGradientXAxis = False
GraphPsprout = True
GraphPsproutSourceTerm = False
GraphNbBranches = False
GraphBarNbBranches = False
GraphFurthestCell = False
GraphConvergenceTime = False
GraphNbCellsPlane = False
GraphBarCellsInPlane = False
GraphFirstTimeCellReachingLesion = False

## Plots

if GraphVEGFConcentration:
    fig = plt.subplots(figsize = (12,8), dpi = 300)
    x = np.linspace(0,100,1000)
    y1 = [0.1*(np.exp(-np.sqrt(10-8.3)/np.sqrt(200)*item)-np.exp(np.sqrt(10-8.3)/np.sqrt(200)*(item-2*100)))/(1-np.exp(-2*100*np.sqrt(10-8.3)/np.sqrt(200))) for item in x]
    y5 = [0.5*(np.exp(-np.sqrt(10-8.3)/np.sqrt(200)*item)-np.exp(np.sqrt(10-8.3)/np.sqrt(200)*(item-2*100)))/(1-np.exp(-2*100*np.sqrt(10-8.3)/np.sqrt(200))) for item in x]
    y10 = [1.0*(np.exp(-np.sqrt(10-8.3)/np.sqrt(200)*item)-np.exp(np.sqrt(10-8.3)/np.sqrt(200)*(item-2*100)))/(1-np.exp(-2*100*np.sqrt(10-8.3)/np.sqrt(200))) for item in x]
    plt.plot(x, y1, label= r'$c_{max} = 0.1$', color = 'xkcd:lightblue')
    plt.plot(x, y5, label=r'$c_{max} = 0.5$', color='xkcd:blue')
    plt.plot(x, y10, label=r'$c_{max} = 1.0$', color='xkcd:royal blue')
    plt.xlabel(r'$x$', fontsize=18) 
    plt.ylabel(r'$c(x)$', fontsize=18)
    plt.xticks(fontsize = 16)
    plt.yticks(fontsize = 16)
    plt.legend(fontsize = 18)
    plt.show()

if GraphVEGFGradient:
    fig = plt.subplots(figsize = (12,8), dpi = 300)
    x = np.linspace(0,100,1000)
    ybaseline = [0.01 for item in x]
    y1 = [0.1*np.sqrt(10-8.3)/np.sqrt(200)*(np.exp(-np.sqrt(10-8.3)/np.sqrt(200)*item)+np.exp(np.sqrt(10-8.3)/np.sqrt(200)*(item-2*100)))/(1-np.exp(-2*100*np.sqrt(10-8.3)/np.sqrt(200))) for item in x]
    y5 = [0.5*np.sqrt(10-8.3)/np.sqrt(200)*(np.exp(-np.sqrt(10-8.3)/np.sqrt(200)*item)+np.exp(np.sqrt(10-8.3)/np.sqrt(200)*(item-2*100)))/(1-np.exp(-2*100*np.sqrt(10-8.3)/np.sqrt(200))) for item in x]
    y10 = [1.0*np.sqrt(10-8.3)/np.sqrt(200)*(np.exp(-np.sqrt(10-8.3)/np.sqrt(200)*item)+np.exp(np.sqrt(10-8.3)/np.sqrt(200)*(item-2*100)))/(1-np.exp(-2*100*np.sqrt(10-8.3)/np.sqrt(200))) for item in x]
    plt.plot(x, ybaseline, label='Constant VEGF', color='xkcd:black')
    plt.plot(x, y1, label= r'$c_{max} = 0.1$', color = 'xkcd:lightblue')
    plt.plot(x, y5, label=r'$c_{max} = 0.5$', color='xkcd:blue')
    # x5average = [stats.mean(x5) for item in y]
    # plt.plot(x5average, y)
    plt.plot(x, y10, label=r'$c_{max} = 1.0$', color='xkcd:royal blue')
    plt.ylabel(r'$x$', fontsize=18) 
    plt.xlabel(r'$|\nabla c(x)|$', fontsize=18)
    plt.xticks(fontsize = 16)
    plt.yticks(fontsize = 16)
    plt.legend(fontsize = 12)
    plt.show()

if GraphVEGFGradientXAxis:
    fig = plt.subplots(figsize = (12,8), dpi = 300)
    y = np.linspace(0,100,1000)
    xbaseline = [0.01 for item in y]
    x1 = [0.1*np.sqrt(10-8.3)/np.sqrt(200)*(np.exp(-np.sqrt(10-8.3)/np.sqrt(200)*item)+np.exp(np.sqrt(10-8.3)/np.sqrt(200)*(item-2*100)))/(1-np.exp(-2*100*np.sqrt(10-8.3)/np.sqrt(200))) for item in y]
    x5 = [0.5*np.sqrt(10-8.3)/np.sqrt(200)*(np.exp(-np.sqrt(10-8.3)/np.sqrt(200)*item)+np.exp(np.sqrt(10-8.3)/np.sqrt(200)*(item-2*100)))/(1-np.exp(-2*100*np.sqrt(10-8.3)/np.sqrt(200))) for item in y]
    x10 = [1.0*np.sqrt(10-8.3)/np.sqrt(200)*(np.exp(-np.sqrt(10-8.3)/np.sqrt(200)*item)+np.exp(np.sqrt(10-8.3)/np.sqrt(200)*(item-2*100)))/(1-np.exp(-2*100*np.sqrt(10-8.3)/np.sqrt(200))) for item in y]
    plt.plot(xbaseline, y, label='Constant VEGF', color='xkcd:black')
    plt.plot(x1, y, label= r'$c_{max} = 0.1$', color = 'xkcd:lightblue')
    plt.plot(x5, y, label=r'$c_{max} = 0.5$', color='xkcd:blue')
    # x5average = [stats.mean(x5) for item in y]
    # plt.plot(x5average, y)
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
    # y1 = [0.2*(item**1/(item**1 + 0.3**1)) for item in x]
    # y2 = [0.2*(item**2/(item**2 + 0.24**2)) for item in x]
    # y3 = [0.2*(item**3/(item**3 + 0.23**3)) for item in x]
    plt.plot(x, y, label='hill function (n=%.2f, K=%.2f)' %(n,K), color='xkcd:red')
    # plt.plot(x, y0, label='linear function', color='xkcd:brown green')
    # plt.plot(x, y1, label='hill function (n=1, K=0.3)', color='xkcd:tan')
    # plt.plot(x, y2, label='hill function (n=2, K=0.24)', color='xkcd:light brown')
    # plt.plot(x, y3, label='hill function (n=3, K=0.23)', color='xkcd:brown')
    plt.xlabel(r'$c$', fontsize=24) 
    plt.ylabel(r'$P_{sprout}$', fontsize=24)
    plt.xticks(fontsize = 24)
    plt.yticks(fontsize = 24)
    plt.legend(fontsize = 24)
    plt.show()

if GraphPsproutSourceTerm:
    fig = plt.subplots(figsize = (12,8), dpi = 300)
    sourceterm = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    #Psprout0 = [0.2*item*(np.exp(-np.sqrt(10-8.3)/np.sqrt(200)*15)-np.exp(np.sqrt(10-8.3)/np.sqrt(200)*(15-2*100)))/(1-np.exp(-2*100*np.sqrt(10-8.3)/np.sqrt(200))) for item in sourceterm]
    Psprout0 = [0.2*(item*(np.exp(-np.sqrt(10-8.3)/np.sqrt(200)*15)-np.exp(np.sqrt(10-8.3)/np.sqrt(200)*(15-2*100)))/(1-np.exp(-2*100*np.sqrt(10-8.3)/np.sqrt(200))))**0.66/(8.17**0.66 + (item*(np.exp(-np.sqrt(10-8.3)/np.sqrt(200)*15)-np.exp(np.sqrt(10-8.3)/np.sqrt(200)*(15-2*100)))/(1-np.exp(-2*100*np.sqrt(10-8.3)/np.sqrt(200))))**0.66) for item in sourceterm]
    Psprout1 = [0.2*(item*(np.exp(-np.sqrt(10-8.3)/np.sqrt(200)*15)-np.exp(np.sqrt(10-8.3)/np.sqrt(200)*(15-2*100)))/(1-np.exp(-2*100*np.sqrt(10-8.3)/np.sqrt(200))))**1/(0.3**1 + (item*(np.exp(-np.sqrt(10-8.3)/np.sqrt(200)*15)-np.exp(np.sqrt(10-8.3)/np.sqrt(200)*(15-2*100)))/(1-np.exp(-2*100*np.sqrt(10-8.3)/np.sqrt(200))))**1) for item in sourceterm]
    Psprout2 = [0.2*(item*(np.exp(-np.sqrt(10-8.3)/np.sqrt(200)*15)-np.exp(np.sqrt(10-8.3)/np.sqrt(200)*(15-2*100)))/(1-np.exp(-2*100*np.sqrt(10-8.3)/np.sqrt(200))))**2/(0.23**2 + (item*(np.exp(-np.sqrt(10-8.3)/np.sqrt(200)*15)-np.exp(np.sqrt(10-8.3)/np.sqrt(200)*(15-2*100)))/(1-np.exp(-2*100*np.sqrt(10-8.3)/np.sqrt(200))))**2) for item in sourceterm]
    Psprout3 = [0.2*(item*(np.exp(-np.sqrt(10-8.3)/np.sqrt(200)*15)-np.exp(np.sqrt(10-8.3)/np.sqrt(200)*(15-2*100)))/(1-np.exp(-2*100*np.sqrt(10-8.3)/np.sqrt(200))))**3/(0.24**3 + (item*(np.exp(-np.sqrt(10-8.3)/np.sqrt(200)*15)-np.exp(np.sqrt(10-8.3)/np.sqrt(200)*(15-2*100)))/(1-np.exp(-2*100*np.sqrt(10-8.3)/np.sqrt(200))))**3) for item in sourceterm]
    plt.plot(sourceterm, Psprout0, label='Hill Function (n=0.66, K=8.17)', color = 'xkcd:green brown')
    plt.plot(sourceterm, Psprout1, label='Hill Function (n=1, K=0.3)', color='xkcd:tan')
    plt.plot(sourceterm, Psprout2, label='Hill Function (n=2, K=0.23)', color='xkcd:light brown')
    plt.plot(sourceterm, Psprout3, label='Hill Function (n=3, K=0.24)', color='xkcd:brown')
    plt.xlabel(r'$c_{max}$', fontsize=18) 
    plt.ylabel(r'$P_{sprout}$', fontsize=18)
    plt.xticks(fontsize = 16)
    plt.yticks(fontsize = 16)
    plt.legend(fontsize = 8)
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

    # plot 
    fig, ax = plt.subplots(figsize = (12,8), dpi = 300, layout='constrained')
 
    ax.plot(sourceterm, totalnumberbranches_constant_average, marker = 'D', markersize = 6, label = 'Constant VEGF', color='xkcd:brown')
    ax.plot(sourceterm, totalnumberbranches_analyticalapproxpde_average, label = r'$P_{sprout}$', marker = '.', markersize = 15.0, color='xkcd:tan')
    ax.legend(loc='upper left', fontsize = 8)

    ax.set_xlabel(r'$c_{max}$', fontsize = 15)
    ax.set_ylabel('Number of Branches', fontsize = 15)
    ax.set_xticks(sourceterm)
    ax.set_yticks([0, 5, 10, 15, 20])
    ax.tick_params(axis = 'both', labelsize = 15)
    plt.show()

# Number of cells after a certain plane depending on source term for different Psprout 
if GraphNbCellsPlane:
    sourceterm = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]

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


    # plot 
    fig, ax = plt.subplots(figsize = (12,8), dpi = 300, layout='constrained')


    ax.plot(sourceterm, cellsafterplane_constant_average, marker = 'D', markersize = 6, label = 'Constant VEGF', color='xkcd:brown')
    ax.plot(sourceterm, cellsafterplane_analyticalapproxpde_average, label = r'$P_{sprout}$', marker = '.', markersize = 15.0, color='xkcd:tan')

    ax.set_xlabel(r'$c_{max}$', fontsize = 15)
    ax.set_ylabel('Cells In Lesion', fontsize = 15)
    ax.set_xticks(sourceterm)
    ax.set_yticks([10, 20, 30, 40, 50, 60])
    ax.legend(loc='lower right', fontsize = 8)
    ax.tick_params(axis = 'both', labelsize = 15)
    plt.show()

# Furthest cell depending on source term for chemotactic force  
if GraphFurthestCell:
    sourceterm = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]

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

            furthestcell_analyticalapproxpde.append(runner.NormFurthestTipCell(file_nodescoordinates_analyticalapproxpde, file_cellmutation_analyticalapproxpde, dim))
            furthestcell_constant.append(runner.NormFurthestTipCell(file_nodescoordinates_constant, file_cellmutation_constant, dim))

        # average of the results for one source term 
        furthestcell_analyticalapproxpde_average[k-1] = stats.mean(furthestcell_analyticalapproxpde)
        furthestcell_constant_average[k-1] = stats.mean(furthestcell_constant)

    # plot 
    fig, ax = plt.subplots(figsize = (12,8), dpi = 300, layout='constrained')


    ax.plot(sourceterm, furthestcell_constant_average, marker = 'D', markersize = 6, label = 'Constant VEGF', color='xkcd:brown')
    ax.plot(sourceterm, furthestcell_analyticalapproxpde_average, label = r'$P_{sprout}$', marker = '.', markersize = 15.0, color='xkcd:tan')
    ax.legend(loc='upper right', fontsize = 8)

    ax.set_xlabel(r'$c_{max}$', fontsize = 15)
    ax.set_ylabel('Network Length', fontsize = 15)
    ax.set_xticks(sourceterm)
    ax.tick_params(axis = 'both', labelsize = 15)
    plt.show()

if GraphFirstTimeCellReachingLesion:
    sourceterm = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]

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

    # plot 
    fig, ax = plt.subplots(figsize = (12,8), dpi = 300, layout='constrained')

    ax.plot(sourceterm, timefirstreachinglesion_constant_average, marker = 'D', markersize = 6, label = 'Constant VEGF', color='xkcd:brown')
    ax.plot(sourceterm, timefirstreachinglesion_analyticalapproxpde_average, label = r'$P_{sprout}$', marker = '.', markersize = 15.0, color='xkcd:tan')
    ax.legend(loc='upper right', fontsize = 8)

    ax.set_xlabel(r'$c_{max}$', fontsize = 15)
    ax.set_ylabel('Time to Reach Lesion', fontsize = 15)
    ax.set_xticks(sourceterm)
    ax.tick_params(axis = 'both', labelsize = 15)
    plt.show()