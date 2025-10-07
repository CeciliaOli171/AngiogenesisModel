## Libraries 
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import statistics as stats

from matplotlib import cm
from matplotlib.ticker import LinearLocator 
from matplotlib.colors import LinearSegmentedColormap

from scipy.spatial import ConvexHull
from ParametersSensitivityRunner import ParametersSensitivityRunner

runner = ParametersSensitivityRunner

## Parameters 
Tcycle = 13.0
time = 300
Tend = 500
TotalTestSourceNb = 10
TotalTestNb = 20
dim = 2
ref_point = 15

## Settings
TestBaseline = False
GraphVEGFConcentration = False
GraphVEGFGradient = True
GraphVEGFGradientXAxis = False
GraphPsprout = False
GraphPsproutSourceTerm = False
GraphNbBranches = False
GraphBarNbBranches = False
GraphFurthestCell = False
GraphConvergenceTime = False
GraphNbCellsPlane = False
GraphBarCellsInPlane = False
GraphFirstTimeCellReachingLesion = False 

## Plots

if TestBaseline:
    baseline_furthestcell_average = 0
    baseline_furthestcell = []
    for k in range(1,21):
        file_path = "/Users/coli171/Chaste/Output/BaselineANZIAM/AngioModelSeed" + str(k) 
        file_cellmutation = file_path + "/results_from_time_0/results.vizmutationstates"
        file_nodescoordinates = file_path + "/results_from_time_0/results.viznodes"

        baseline_furthestcell.append(runner.NormFurthestVesselTip(file_nodescoordinates, file_cellmutation, dim))

        # average of the results for one source term 
        baseline_furthestcell_average = stats.mean(baseline_furthestcell)

    baseline_nbcellsafterplane_average = 0
    baseline_nbcellsafterplane = []
    for k in range(1,21):
        file_path = "/Users/coli171/Chaste/Output/BaselineANZIAM/AngioModelSeed" + str(k) 
        file_cellmutation = file_path + "/results_from_time_0/results.vizmutationstates"
        file_nodescoordinates = file_path + "/results_from_time_0/results.viznodes"

        baseline_nbcellsafterplane.append(runner.NbCellsAfterPlane(file_nodescoordinates, ref_point, dim))

        # average of the results for one source term 
        baseline_nbcellsafterplane_average = stats.mean(baseline_nbcellsafterplane)

    baseline_timefirstreachinglesion_average = 0
    baseline_timefirstreachinglesion = []
    for k in range(1,21):
        file_path = "/Users/coli171/Chaste/Output/BaselineANZIAM/AngioModelSeed" + str(k) 
        file_cellmutation = file_path + "/results_from_time_0/results.vizmutationstates"
        file_nodescoordinates = file_path + "/results_from_time_0/results.viznodes"

        baseline_timefirstreachinglesion.append(runner.TimeFirstReachingPlane(file_nodescoordinates, ref_point, dim))

        # average of the results for one source term 
        baseline_timefirstreachinglesion_average = stats.mean(baseline_timefirstreachinglesion)

    print(baseline_furthestcell_average)
    print(baseline_nbcellsafterplane_average)
    print(baseline_timefirstreachinglesion_average)

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
    y = [0.2*item for item in x]
    y1 = [0.2*(item**1/(item**1 + 0.3**1)) for item in x]
    y2 = [0.2*(item**2/(item**2 + 0.24**2)) for item in x]
    y3 = [0.2*(item**3/(item**3 + 0.23**3)) for item in x]
    plt.plot(x, y, label= 'linear function', color='xkcd:brown green')
    plt.plot(x, y1, label='hill function (n=1, K=0.3)', color='xkcd:tan')
    plt.plot(x, y2, label='hill function (n=2, K=0.24)', color='xkcd:light brown')
    plt.plot(x, y3, label='hill function (n=3, K=0.23)', color='xkcd:brown')
    plt.xlabel(r'$c$', fontsize=24) 
    plt.ylabel(r'$P_{sprout}$', fontsize=24)
    plt.xticks(fontsize = 24)
    plt.yticks(fontsize = 24)
    plt.legend(fontsize = 24)
    plt.show()

if GraphPsproutSourceTerm:
    fig = plt.subplots(figsize = (12,8), dpi = 300)
    sourceterm = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    Psprout0 = [0.2*item*(np.exp(-np.sqrt(10-8.3)/np.sqrt(200)*15)-np.exp(np.sqrt(10-8.3)/np.sqrt(200)*(15-2*100)))/(1-np.exp(-2*100*np.sqrt(10-8.3)/np.sqrt(200))) for item in sourceterm]
    Psprout1 = [0.2*(item*(np.exp(-np.sqrt(10-8.3)/np.sqrt(200)*15)-np.exp(np.sqrt(10-8.3)/np.sqrt(200)*(15-2*100)))/(1-np.exp(-2*100*np.sqrt(10-8.3)/np.sqrt(200))))**1/(0.3**1 + (item*(np.exp(-np.sqrt(10-8.3)/np.sqrt(200)*15)-np.exp(np.sqrt(10-8.3)/np.sqrt(200)*(15-2*100)))/(1-np.exp(-2*100*np.sqrt(10-8.3)/np.sqrt(200))))**1) for item in sourceterm]
    Psprout2 = [0.2*(item*(np.exp(-np.sqrt(10-8.3)/np.sqrt(200)*15)-np.exp(np.sqrt(10-8.3)/np.sqrt(200)*(15-2*100)))/(1-np.exp(-2*100*np.sqrt(10-8.3)/np.sqrt(200))))**2/(0.23**2 + (item*(np.exp(-np.sqrt(10-8.3)/np.sqrt(200)*15)-np.exp(np.sqrt(10-8.3)/np.sqrt(200)*(15-2*100)))/(1-np.exp(-2*100*np.sqrt(10-8.3)/np.sqrt(200))))**2) for item in sourceterm]
    Psprout3 = [0.2*(item*(np.exp(-np.sqrt(10-8.3)/np.sqrt(200)*15)-np.exp(np.sqrt(10-8.3)/np.sqrt(200)*(15-2*100)))/(1-np.exp(-2*100*np.sqrt(10-8.3)/np.sqrt(200))))**3/(0.24**3 + (item*(np.exp(-np.sqrt(10-8.3)/np.sqrt(200)*15)-np.exp(np.sqrt(10-8.3)/np.sqrt(200)*(15-2*100)))/(1-np.exp(-2*100*np.sqrt(10-8.3)/np.sqrt(200))))**3) for item in sourceterm]
    plt.plot(sourceterm, Psprout0, label='Linear Function', color = 'xkcd:green brown')
    plt.plot(sourceterm, Psprout1, label='Hill Function (n=1)', color='xkcd:tan')
    plt.plot(sourceterm, Psprout2, label='Hill Function (n=2)', color='xkcd:light brown')
    plt.plot(sourceterm, Psprout3, label='Hill Function (n=3)', color='xkcd:brown')
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
if GraphBarNbBranches:
    # Linear function 
    # loop over all the files 
    linearfunction_totalnumberbranches = []
    for k in range(1, 11):
        for j in range(1, 51):
            file_path = "/Users/coli171/Chaste/Output/ModelWithPdes2DTestAnalyticalSeedSource/Test0/ModelWithPdes2DTest0AnalyticalSeed" + str(j) + "Source" + str(k)
            file_cellmutation = file_path + "/results_from_time_0/results.vizmutationstates"

            linearfunction_totalnumberbranches.append(runner.TotalNumberBranches(file_cellmutation))

    # Hill function (n=1)
    # loop over all the files 
    hillfunction1_totalnumberbranches = []
    for k in range(1, 11):
        for j in range(1, 51):
            file_path = "/Users/coli171/Chaste/Output/ModelWithPdes2DTestAnalyticalSeedSource/Test1/ModelWithPdes2DTest1AnalyticalSeed" + str(j) + "Source" + str(k)
            file_cellmutation = file_path + "/results_from_time_0/results.vizmutationstates"

            hillfunction1_totalnumberbranches.append(runner.TotalNumberBranches(file_cellmutation))


    # Hill function (n=2)
    # loop over all the files 
    hillfunction2_totalnumberbranches = []
    for k in range(1, 11):
        for j in range(1, 51):
            file_path = "/Users/coli171/Chaste/Output/ModelWithPdes2DTestAnalyticalSeedSource/Test2/ModelWithPdes2DTest2AnalyticalSeed" + str(j) + "Source" + str(k)
            file_cellmutation = file_path + "/results_from_time_0/results.vizmutationstates"

            hillfunction2_totalnumberbranches.append(runner.TotalNumberBranches(file_cellmutation))


    # Hill function (n=3)
    # loop over all the files 
    hillfunction3_totalnumberbranches = []
    for k in range(1, 11):
        for j in range(1, 51):
            file_path = "/Users/coli171/Chaste/Output/ModelWithPdes2DTestAnalyticalSeedSource/Test3/ModelWithPdes2DTest3AnalyticalSeed" + str(j) + "Source" + str(k)
            file_cellmutation = file_path + "/results_from_time_0/results.vizmutationstates"

            hillfunction3_totalnumberbranches.append(runner.TotalNumberBranches(file_cellmutation))


    # plot 
    fig = plt.subplots(figsize = (12,8), dpi = 300)
    PsproutFunctions = ["Linear \n Function", "Hill Function \n (n=1)", "Hill Function \n (n=2)", "Hill Function \n (n=3)"]
    BranchNumber = [linearfunction_totalnumberbranches, hillfunction1_totalnumberbranches, hillfunction2_totalnumberbranches, hillfunction3_totalnumberbranches]
    PsproutFunctionsFormula = np.arange(len(PsproutFunctions))
    box_colors = ['xkcd:brown green', 'xkcd:tan', 'xkcd:light brown', 'xkcd:brown']
    bp = plt.boxplot(BranchNumber, labels=PsproutFunctions, patch_artist=True, manage_ticks = True, medianprops=dict(color='black', linewidth=2), whiskerprops=dict(color='black', linewidth=1), flierprops=dict(marker='o', color='black', markersize=5)) 
    for i, patch in enumerate(bp['boxes']):  # `plt.gca().artists` gives the box patches
        patch.set_facecolor(box_colors[i])  # Set the facecolor for each box
    plt.ylabel('Number of Branches')
    plt.show()

if GraphNbBranches:
    sourceterm = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]

    # Linear function 
    # loop over all the files 
    linearfunction_totalnumberbranches_average = np.zeros(10)
    for k in range(1, 11):
        linearfunction_totalnumberbranches = []
        for j in range(1, 51):
            file_path = "/Users/coli171/Chaste/Output/ModelWithPdes2DTestAnalyticalSeedSource/Test0/ModelWithPdes2DTest0AnalyticalSeed" + str(j) + "Source" + str(k)
            file_cellmutation = file_path + "/results_from_time_0/results.vizmutationstates"

            linearfunction_totalnumberbranches.append(runner.TotalNumberBranches(file_cellmutation))

        # average of the results for one source term 
        linearfunction_totalnumberbranches_average[k-1] = stats.mean(linearfunction_totalnumberbranches)
    

    # Hill function (n=1)
    # loop over all the files 
    hillfunction1_totalnumberbranches_average = np.zeros(10)
    for k in range(1, 11):
        hillfunction1_totalnumberbranches = []
        for j in range(1, 11):
            file_path = "/Users/coli171/Chaste/Output/ModelWithPdes2DTestAnalyticalSeedSource/Test1/ModelWithPdes2DTest1AnalyticalSeed" + str(j) + "Source" + str(k)
            file_cellmutation = file_path + "/results_from_time_0/results.vizmutationstates"

            hillfunction1_totalnumberbranches.append(runner.TotalNumberBranches(file_cellmutation))

        # average of the results for one source term 
        hillfunction1_totalnumberbranches_average[k-1] = stats.mean(hillfunction1_totalnumberbranches)

    # Hill function (n=2)
    # loop over all the files 
    hillfunction2_totalnumberbranches_average = np.zeros(10)
    for k in range(1, 11):
        hillfunction2_totalnumberbranches = []
        for j in range(1, 11):
            file_path = "/Users/coli171/Chaste/Output/ModelWithPdes2DTestAnalyticalSeedSource/Test2/ModelWithPdes2DTest2AnalyticalSeed" + str(j) + "Source" + str(k)
            file_cellmutation = file_path + "/results_from_time_0/results.vizmutationstates"

            hillfunction2_totalnumberbranches.append(runner.TotalNumberBranches(file_cellmutation))

        # average of the results for one source term 
        hillfunction2_totalnumberbranches_average[k-1] = stats.mean(hillfunction2_totalnumberbranches)

    # Hill function (n=3)
    # loop over all the files 
    hillfunction3_totalnumberbranches_average = np.zeros(10)
    for k in range(1, 11):
        hillfunction3_totalnumberbranches = []
        for j in range(1, 11):
            file_path = "/Users/coli171/Chaste/Output/ModelWithPdes2DTestAnalyticalSeedSource/Test3/ModelWithPdes2DTest3AnalyticalSeed" + str(j) + "Source" + str(k)
            file_cellmutation = file_path + "/results_from_time_0/results.vizmutationstates"

            hillfunction3_totalnumberbranches.append(runner.TotalNumberBranches(file_cellmutation))

        # average of the results for one source term 
        hillfunction3_totalnumberbranches_average[k-1] += stats.mean(hillfunction3_totalnumberbranches)

    # plot 
    fig, ax = plt.subplots(figsize = (12,8), dpi = 300, layout='constrained')

    # ANZIAM Presentation
    # ax.plot(0.5, 17, marker = 'D', markerfacecolor='xkcd:black', markeredgecolor='xkcd:black', markersize = 10, label = 'Constant VEGF')
    # ax.plot(sourceterm, linearfunction_totalnumberbranches_average, label = 'Linear Function', marker = '.', markersize = 15.0, color='xkcd:brown green')
    # ax.plot(sourceterm, hillfunction1_totalnumberbranches_average, label = 'Hill Function (n=1)', marker = '.', markersize = 15.0, color='xkcd:tan')
    # ax.plot(sourceterm, hillfunction2_totalnumberbranches_average, label = 'Hill Function (n=2)', marker = '.', markersize = 15.0, color='xkcd:light brown')
    # ax.plot(sourceterm, hillfunction3_totalnumberbranches_average, label = 'Hill Function (n=3)', marker = '.', markersize = 15.0, color='xkcd:brown')
    # ax.legend(loc='upper left', fontsize = 8, ncol = 2)

    # ABI Poster 
    ax.plot(0.5, 17, marker = 'D', markerfacecolor='xkcd:red', markeredgecolor='xkcd:red', markersize = 6, label = 'Constant VEGF')
    ax.plot(sourceterm, hillfunction1_totalnumberbranches_average, label = 'Hill Function (n=1)', marker = '.', markersize = 15.0, color='xkcd:tan')
    ax.legend(loc='lower right', fontsize = 8)

    ax.set_xlabel(r'$c_{max}$', fontsize = 15)
    ax.set_ylabel('Number of Branches', fontsize = 15)
    ax.set_xticks(sourceterm)
    ax.set_yticks([0, 5, 10, 15, 20])
    ax.tick_params(axis = 'both', labelsize = 15)
    plt.show()

# Number of cells after a certain plane depending on source term for different Psprout 
if GraphNbCellsPlane:
    sourceterm = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]

    # Linear function 
    # loop over all the files 
    linearfunction_cellsafterplane_average = np.zeros(10)
    for k in range(1,11):
        linearfunction_cellsafterplane = []
        for j in range(1,51):
            file_path = "/Users/coli171/Chaste/Output/ModelWithPdes2DTestAnalyticalSeedSource/Test0/ModelWithPdes2DTest0AnalyticalSeed" + str(j) + "Source" + str(k)
            file_cellmutation = file_path + "/results_from_time_0/results.vizmutationstates"
            file_nodescoordinates = file_path + "/results_from_time_0/results.viznodes"

            linearfunction_cellsafterplane.append(runner.NbCellsAfterPlane(file_nodescoordinates, ref_point, dim))

        # average of the results for one source term 
        linearfunction_cellsafterplane_average[k-1] = stats.mean(linearfunction_cellsafterplane)
    

    # Hill function (n=1)
    # loop over all the files 
    hillfunction1_cellsafterplane_average = np.zeros(10)
    for k in range(1,11):
        hillfunction1_cellsafterplane = []
        for j in range(1,51):
            file_path = "/Users/coli171/Chaste/Output/ModelWithPdes2DTestAnalyticalSeedSource/Test1/ModelWithPdes2DTest1AnalyticalSeed" + str(j) + "Source" + str(k)
            file_cellmutation = file_path + "/results_from_time_0/results.vizmutationstates"
            file_nodescoordinates = file_path + "/results_from_time_0/results.viznodes"

            hillfunction1_cellsafterplane.append(runner.NbCellsAfterPlane(file_nodescoordinates, ref_point, dim))

        # average of the results for one source term 
        hillfunction1_cellsafterplane_average[k-1] = stats.mean(hillfunction1_cellsafterplane)

    # Hill function (n=2)
    # loop over all the files 
    hillfunction2_cellsafterplane_average = np.zeros(10)
    for k in range(1,11):
        hillfunction2_cellsafterplane = []
        for j in range(1,51):
            file_path = "/Users/coli171/Chaste/Output/ModelWithPdes2DTestAnalyticalSeedSource/Test2/ModelWithPdes2DTest2AnalyticalSeed" + str(j) + "Source" + str(k)
            file_cellmutation = file_path + "/results_from_time_0/results.vizmutationstates"
            file_nodescoordinates = file_path + "/results_from_time_0/results.viznodes"

            hillfunction2_cellsafterplane.append(runner.NbCellsAfterPlane(file_nodescoordinates, ref_point, dim))

        # average of the results for one source term 
        hillfunction2_cellsafterplane_average[k-1] = stats.mean(hillfunction2_cellsafterplane) 

    # Hill function (n=3)
    # loop over all the files 
    hillfunction3_cellsafterplane_average = np.zeros(10)
    for k in range(1,11):
        hillfunction3_cellsafterplane = []
        for j in range(1,51):
            file_path = "/Users/coli171/Chaste/Output/ModelWithPdes2DTestAnalyticalSeedSource/Test3/ModelWithPdes2DTest3AnalyticalSeed" + str(j) + "Source" + str(k)
            file_cellmutation = file_path + "/results_from_time_0/results.vizmutationstates"
            file_nodescoordinates = file_path + "/results_from_time_0/results.viznodes"

            hillfunction3_cellsafterplane.append(runner.NbCellsAfterPlane(file_nodescoordinates, ref_point, dim))

        # average of the results for one source term 
        hillfunction3_cellsafterplane_average[k-1] = stats.mean(hillfunction3_cellsafterplane)

    # plot 
    fig, ax = plt.subplots(figsize = (12,8), dpi = 300, layout='constrained')

    #ANZIAM
    # ax.plot(0.5, 34, marker = 'D', markerfacecolor='xkcd:black', markeredgecolor='xkcd:black', markersize = 6, label = 'Constant VEGF')
    # ax.plot(sourceterm, linearfunction_cellsafterplane_average, label = 'Linear Function', marker = '.', markersize = 15.0, color='xkcd:green brown')
    # ax.plot(sourceterm, hillfunction1_cellsafterplane_average, label = 'Hill Function (n=1)', marker = '.', markersize = 15.0, color='xkcd:tan')
    # ax.plot(sourceterm, hillfunction2_cellsafterplane_average, label = 'Hill Function (n=2)', marker = '.', markersize = 15.0, color='xkcd:light brown')
    # ax.plot(sourceterm, hillfunction3_cellsafterplane_average, label = 'Hill Function (n=3)', marker = '.', markersize = 15.0, color='xkcd:brown')

    #ABI Poster
    ax.plot(0.5, 34, marker = 'D', markerfacecolor='xkcd:red', markeredgecolor='xkcd:red', markersize = 6, label = 'Constant VEGF')
    ax.plot(sourceterm, hillfunction1_cellsafterplane_average, label = 'Hill Function (n=1)', marker = '.', markersize = 15.0, color='xkcd:tan')

    ax.set_xlabel(r'$c_{max}$', fontsize = 15)
    ax.set_ylabel('Cells In Lesion', fontsize = 15)
    ax.set_xticks(sourceterm)
    ax.set_yticks([15, 20, 25, 30, 35])
    ax.legend(loc='lower right', fontsize = 8)
    ax.tick_params(axis = 'both', labelsize = 15)
    plt.show()

# Randomness checking of cells in lesion 
if GraphBarCellsInPlane:
    # Linear function 
    # loop over all the files 
    linearfunction_cellsafterplane = []
    for k in range(1,11):
        for j in range(1,51):
            file_path = "/Users/coli171/Chaste/Output/ModelWithPdes2DTestAnalyticalSeedSource/Test0/ModelWithPdes2DTest0AnalyticalSeed" + str(j) + "Source" + str(k)
            file_cellmutation = file_path + "/results_from_time_0/results.vizmutationstates"
            file_nodescoordinates = file_path + "/results_from_time_0/results.viznodes"

            linearfunction_cellsafterplane.append(runner.NbCellsAfterPlane(file_nodescoordinates, ref_point, dim))
    

    # Hill function (n=1)
    # loop over all the files 
    hillfunction1_cellsafterplane = []
    for k in range(1,11):
        for j in range(1,51):
            file_path = "/Users/coli171/Chaste/Output/ModelWithPdes2DTestAnalyticalSeedSource/Test1/ModelWithPdes2DTest1AnalyticalSeed" + str(j) + "Source" + str(k)
            file_cellmutation = file_path + "/results_from_time_0/results.vizmutationstates"
            file_nodescoordinates = file_path + "/results_from_time_0/results.viznodes"

            hillfunction1_cellsafterplane.append(runner.NbCellsAfterPlane(file_nodescoordinates, ref_point, dim))


    # Hill function (n=2)
    # loop over all the files 
    hillfunction2_cellsafterplane = []
    for k in range(1,11):
        for j in range(1,51):
            file_path = "/Users/coli171/Chaste/Output/ModelWithPdes2DTestAnalyticalSeedSource/Test2/ModelWithPdes2DTest2AnalyticalSeed" + str(j) + "Source" + str(k)
            file_cellmutation = file_path + "/results_from_time_0/results.vizmutationstates"
            file_nodescoordinates = file_path + "/results_from_time_0/results.viznodes"

            hillfunction2_cellsafterplane.append(runner.NbCellsAfterPlane(file_nodescoordinates, ref_point, dim))


    # Hill function (n=3)
    # loop over all the files 
    hillfunction3_cellsafterplane = []
    for k in range(1,11):
        for j in range(1,51):
            file_path = "/Users/coli171/Chaste/Output/ModelWithPdes2DTestAnalyticalSeedSource/Test3/ModelWithPdes2DTest3AnalyticalSeed" + str(j) + "Source" + str(k)
            file_cellmutation = file_path + "/results_from_time_0/results.vizmutationstates"
            file_nodescoordinates = file_path + "/results_from_time_0/results.viznodes"

            hillfunction3_cellsafterplane.append(runner.NbCellsAfterPlane(file_nodescoordinates, ref_point, dim))


    # plot 
    fig = plt.subplots(figsize = (12,8), dpi = 300)
    PsproutFunctions = ["Linear \n Function", "Hill Function \n (n=1)", "Hill Function \n (n=2)", "Hill Function \n (n=3)"]
    CellsInPlane = [linearfunction_cellsafterplane, hillfunction1_cellsafterplane, hillfunction2_cellsafterplane, hillfunction3_cellsafterplane]
    PsproutFunctionsFormula = np.arange(len(PsproutFunctions))
    box_colors = ['xkcd:brown green', 'xkcd:tan', 'xkcd:light brown', 'xkcd:brown']
    bp = plt.boxplot(CellsInPlane, labels=PsproutFunctions, patch_artist=True, manage_ticks = True, medianprops=dict(color='black', linewidth=2), whiskerprops=dict(color='black', linewidth=1), flierprops=dict(marker='o', color='black', markersize=5)) 
    for i, patch in enumerate(bp['boxes']):  # `plt.gca().artists` gives the box patches
        patch.set_facecolor(box_colors[i])  # Set the facecolor for each box
    plt.ylabel('Cells In Plane')
    plt.show()

# Furthest cell depending on source term for chemotactic force  
if GraphFurthestCell:
    sourceterm = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]

    # Linear function 
    # loop over all the files 
    linearfunction_furthestcell_average = np.zeros(10)
    for k in range(1,11):
        linearfunction_furthestcell = []
        for j in range(1,51):
            file_path = "/Users/coli171/Chaste/Output/ModelWithPdes2DTestAnalyticalSeedSource/Test0/ModelWithPdes2DTest0AnalyticalSeed" + str(j) + "Source" + str(k)
            file_cellmutation = file_path + "/results_from_time_0/results.vizmutationstates"
            file_nodescoordinates = file_path + "/results_from_time_0/results.viznodes"

            linearfunction_furthestcell.append(runner.NormFurthestVesselTip(file_nodescoordinates, file_cellmutation, dim))

        # average of the results for one source term 
        linearfunction_furthestcell_average[k-1] = stats.mean(linearfunction_furthestcell)
    

    # Hill function (n=1)
    # loop over all the files 
    hillfunction1_furthestcell_average = np.zeros(10)
    for k in range(1,11):
        hillfunction1_furthestcell = []
        for j in range(1,51):
            file_path = "/Users/coli171/Chaste/Output/ModelWithPdes2DTestAnalyticalSeedSource/Test1/ModelWithPdes2DTest1AnalyticalSeed" + str(j) + "Source" + str(k)
            file_cellmutation = file_path + "/results_from_time_0/results.vizmutationstates"
            file_nodescoordinates = file_path + "/results_from_time_0/results.viznodes"

            hillfunction1_furthestcell.append(runner.NormFurthestVesselTip(file_nodescoordinates, file_cellmutation, dim))

        # average of the results for one source term 
        hillfunction1_furthestcell_average[k-1] = stats.mean(hillfunction1_furthestcell)

    # Hill function (n=2)
    # loop over all the files 
    hillfunction2_furthestcell_average = np.zeros(10)
    for k in range(1,11):
        hillfunction2_furthestcell = []
        for j in range(1,51):
            file_path = "/Users/coli171/Chaste/Output/ModelWithPdes2DTestAnalyticalSeedSource/Test2/ModelWithPdes2DTest2AnalyticalSeed" + str(j) + "Source" + str(k)
            file_cellmutation = file_path + "/results_from_time_0/results.vizmutationstates"
            file_nodescoordinates = file_path + "/results_from_time_0/results.viznodes"

            hillfunction2_furthestcell.append(runner.NormFurthestVesselTip(file_nodescoordinates, file_cellmutation, dim))

        # average of the results for one source term 
        hillfunction2_furthestcell_average[k-1] = stats.mean(hillfunction2_furthestcell)

    # Hill function (n=3)
    # loop over all the files 
    hillfunction3_furthestcell_average = np.zeros(10)
    for k in range(1,11):
        hillfunction3_furthestcell = []
        for j in range(1,51):
            file_path = "/Users/coli171/Chaste/Output/ModelWithPdes2DTestAnalyticalSeedSource/Test3/ModelWithPdes2DTest3AnalyticalSeed" + str(j) + "Source" + str(k)
            file_cellmutation = file_path + "/results_from_time_0/results.vizmutationstates"
            file_nodescoordinates = file_path + "/results_from_time_0/results.viznodes"

            hillfunction3_furthestcell.append(runner.NormFurthestVesselTip(file_nodescoordinates, file_cellmutation, dim))

        # average of the results for one source term 
        hillfunction3_furthestcell_average[k-1] = stats.mean(hillfunction3_furthestcell)

    # plot 
    fig, ax = plt.subplots(figsize = (12,8), dpi = 300, layout='constrained')

    # ANZIAM Presentation
    # ax.plot(0.5, 22, marker = 'D', markerfacecolor='xkcd:black', markeredgecolor='xkcd:black', markersize = 10, label = 'Constant VEGF')
    # ax.plot(sourceterm, linearfunction_furthestcell_average, label = 'Linear Function', marker = '.', markersize = 15.0, color='xkcd:green brown')
    # ax.plot(sourceterm, hillfunction1_furthestcell_average, label = 'Hill Function (n=1)', marker = '.', markersize = 15.0, color='xkcd:tan')
    # ax.plot(sourceterm, hillfunction2_furthestcell_average, label = 'Hill Function (n=2)', marker = '.', markersize = 15.0, color='xkcd:light brown')
    # ax.plot(sourceterm, hillfunction3_furthestcell_average, label = 'Hill Function (n=3)', marker = '.', markersize = 15.0, color='xkcd:brown')
    #ax.legend(loc='upper left', fontsize = 24)

    # ABI Poster
    ax.plot(0.5, 22, marker = 'D', markerfacecolor='xkcd:red', markeredgecolor='xkcd:red', markersize = 6, label = 'Constant VEGF')
    ax.plot(sourceterm, hillfunction1_furthestcell_average, label = 'Hill Function (n=1)', marker = '.', markersize = 15.0, color='xkcd:tan')
    ax.legend(loc='upper right', fontsize = 8)

    ax.set_xlabel(r'$c_{max}$', fontsize = 15)
    ax.set_ylabel('Network Length', fontsize = 15)
    ax.set_xticks(sourceterm)
    ax.tick_params(axis = 'both', labelsize = 15)
    plt.show()

if GraphFirstTimeCellReachingLesion:
    sourceterm = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]

    # Linear function 
    # loop over all the files 
    linearfunction_timefirstreachinglesion_average = np.zeros(10)
    for k in range(1,11):
        linearfunction_timefirstreachinglesion = []
        for j in range(1,51):
            file_path = "/Users/coli171/Chaste/Output/ModelWithPdes2DTestAnalyticalSeedSource/Test0/ModelWithPdes2DTest0AnalyticalSeed" + str(j) + "Source" + str(k)
            file_cellmutation = file_path + "/results_from_time_0/results.vizmutationstates"
            file_nodescoordinates = file_path + "/results_from_time_0/results.viznodes"

            linearfunction_timefirstreachinglesion.append(runner.TimeFirstReachingPlane(file_nodescoordinates, ref_point, dim))

        # average of the results for one source term 
        linearfunction_timefirstreachinglesion_average[k-1] = stats.mean(linearfunction_timefirstreachinglesion)
    

    # Hill function (n=1)
    # loop over all the files 
    hillfunction1_timefirstreachinglesion_average = np.zeros(10)
    for k in range(1,11):
        hillfunction1_timefirstreachinglesion = []
        for j in range(1,51):
            file_path = "/Users/coli171/Chaste/Output/ModelWithPdes2DTestAnalyticalSeedSource/Test1/ModelWithPdes2DTest1AnalyticalSeed" + str(j) + "Source" + str(k)
            file_cellmutation = file_path + "/results_from_time_0/results.vizmutationstates"
            file_nodescoordinates = file_path + "/results_from_time_0/results.viznodes"

            hillfunction1_timefirstreachinglesion.append(runner.TimeFirstReachingPlane(file_nodescoordinates, ref_point, dim))

        # average of the results for one source term 
        hillfunction1_timefirstreachinglesion_average[k-1] = stats.mean(hillfunction1_timefirstreachinglesion)

    # Hill function (n=2)
    # loop over all the files 
    hillfunction2_timefirstreachinglesion_average = np.zeros(10)
    for k in range(1,11):
        hillfunction2_timefirstreachinglesion = []
        for j in range(1,51):
            file_path = "/Users/coli171/Chaste/Output/ModelWithPdes2DTestAnalyticalSeedSource/Test2/ModelWithPdes2DTest2AnalyticalSeed" + str(j) + "Source" + str(k)
            file_cellmutation = file_path + "/results_from_time_0/results.vizmutationstates"
            file_nodescoordinates = file_path + "/results_from_time_0/results.viznodes"

            hillfunction2_timefirstreachinglesion.append(runner.TimeFirstReachingPlane(file_nodescoordinates, ref_point, dim))

        # average of the results for one source term 
        hillfunction2_timefirstreachinglesion_average[k-1] = stats.mean(hillfunction2_timefirstreachinglesion)

    # Hill function (n=3)
    # loop over all the files 
    hillfunction3_timefirstreachinglesion_average = np.zeros(10)
    for k in range(1,11):
        hillfunction3_timefirstreachinglesion = []
        for j in range(1,51):
            file_path = "/Users/coli171/Chaste/Output/ModelWithPdes2DTestAnalyticalSeedSource/Test3/ModelWithPdes2DTest3AnalyticalSeed" + str(j) + "Source" + str(k)
            file_cellmutation = file_path + "/results_from_time_0/results.vizmutationstates"
            file_nodescoordinates = file_path + "/results_from_time_0/results.viznodes"

            hillfunction3_timefirstreachinglesion.append(runner.TimeFirstReachingPlane(file_nodescoordinates, ref_point, dim))

        # average of the results for one source term 
        hillfunction3_timefirstreachinglesion_average[k-1] = stats.mean(hillfunction3_timefirstreachinglesion)

    # plot 
    fig, ax = plt.subplots(figsize = (12,8), dpi = 300, layout='constrained')

    # ANZIAM Presentation
    # ax.plot(0.5, 136, marker = 'D', markerfacecolor='xkcd:black', markeredgecolor='xkcd:black', markersize = 6,label = 'Constant VEGF')
    # ax.plot(sourceterm, linearfunction_timefirstreachinglesion_average, label = 'Linear Function', marker = '.', markersize = 15.0, color='xkcd:green brown')
    # ax.plot(sourceterm, hillfunction1_timefirstreachinglesion_average, label = 'Hill Function (n=1)', marker = '.', markersize = 15.0, color='xkcd:tan')
    # ax.plot(sourceterm, hillfunction2_timefirstreachinglesion_average, label = 'Hill Function (n=2)', marker = '.', markersize = 15.0, color='xkcd:light brown')
    # ax.plot(sourceterm, hillfunction3_timefirstreachinglesion_average, label = 'Hill Function (n=3)', marker = '.', markersize = 15.0, color='xkcd:brown')
    # ax.legend(loc='upper left', fontsize = 24)

    # ABI Poster
    ax.plot(0.5, 136, marker = 'D', markerfacecolor='xkcd:red', markeredgecolor='xkcd:red', markersize = 6, label = 'Constant VEGF')
    ax.plot(sourceterm, hillfunction1_timefirstreachinglesion_average, label = 'Hill Function (n=1)', marker = '.', markersize = 15.0, color='xkcd:tan')
    ax.legend(loc='upper right', fontsize = 8)

    ax.set_xlabel(r'$c_{max}$', fontsize = 15)
    ax.set_ylabel('Time to Reach Lesion', fontsize = 15)
    ax.set_xticks(sourceterm)
    ax.tick_params(axis = 'both', labelsize = 15)
    plt.show()