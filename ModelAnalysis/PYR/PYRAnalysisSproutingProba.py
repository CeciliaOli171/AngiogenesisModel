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
TotalTestNb = 10 # 9 for Psprout and 10 for the time step study
plane_normal_vector = [18, 0, 0, 0]
t = np.linspace(0, 300, 301)

## Settings 
StudyTimeStep = False

GraphNbBranches = False
GraphNbCells = False
GraphNbTipStalkCells = False
GraphPlaneNbCells = False 
GraphPlaneNbCellsArea = False
GraphPlaneNbCellsTimeStep = False
GraphPlaneNbCellsAreaTimeStep = False

BarPlotNbBranches = False
BarPlotNbCells = False

GraphDivisionTime = False
SmallValuesPsprout = False
BigValuesPsprout = False

ComputationalError = False

GraphAreaVolume = False
GraphPlaneArea = False 
GraphMinNbConsecutiveBranches = False
GraphMaxNbConsecutiveBranches = False
GraphTortuosity = False
GraphTortuosityNbBranches = True
GraphAverageClosestNeighbour = False 
GraphAverageClosestNeighbourTimeStep = False

## Initialisation
y_totalnumbercells_average = np.zeros(20)
y_numbercellsplane_average = np.zeros(20)
y_totalnumbertipcells_average = np.zeros(20)
y_totalnumberstalkcells_average = np.zeros(20)
y_totalnumbertipcellsmax = np.zeros(20)
y_totalnumberstalkcellsmax = np.zeros(20)
y_totalnumbertipcellsmin = np.zeros(20)
y_totalnumberstalkcellsmin = np.zeros(20)
y_totalnumberbranches_average = np.zeros(20)
y_areaplane_average = np.zeros(20)
y_area_average = np.zeros(20)
y_volume_average = np.zeros(20)
y_minconsecutivebranches_average = np.zeros(20)
y_maxconsecutivebranches_average = np.zeros(20)
y_computationalnumberofbranches_average = np.zeros(20)
y_expectednumberofbranches_average = np.zeros(20)
y_computationalnumberofbranchesmax = np.zeros(20)
y_computationalnumberofbranchesmin = np.zeros(20)
y_expectednumberofbranchesmax = np.zeros(20)
y_formulaexpectednumberofbranches_average = np.zeros(20)
y_tortuosity_average = np.zeros(20)
y_neighbouringdistance_average = np.zeros(20)

for j in range(1,10):
    x1 = []
    y_totalnumbercells1 = []
    y_numbercellsplane1 = []
    y_totalnumbertipcells1 = []
    y_totalnumberstalkcells1 = []
    y_totalnumberbranches1 = []
    y_areaplane1 = []
    y_area1 = []
    y_volume1 = []
    y_minconsecutivebranches1 = []
    y_maxconsecutivebranches1 = []
    y_computationalnumberofbranches1 = []
    y_expectednumberofbranches1 = []
    y_formulaexpectednumberofbranches1 = []
    y_tortuosity1 = []
    y_neighbouringdistance1 = []

    x2 = []
    y_totalnumbercells2 = []
    y_numbercellsplane2 = []
    y_totalnumbertipcells2 = []
    y_totalnumberstalkcells2 = []
    y_totalnumberbranches2 = []
    y_areaplane2 = []
    y_area2 = []
    y_volume2 = []
    y_minconsecutivebranches2 = []
    y_maxconsecutivebranches2 = []
    y_computationalnumberofbranches2 = []
    y_expectednumberofbranches2 = []
    y_formulaexpectednumberofbranches2 = []
    y_tortuosity2 = []
    y_neighbouringdistance2 = []

    for k in range(0,10):
        file_cellmutation = '/Users/coli171/Chaste/Output/TestParametersSensitivity/TESTPSPROUT_2NDVERSION/testPsprout0' + str(k) + "/test" + str(j) + "/results_from_time_0/results.vizmutationstates"
        file_nodescoordinates = '/Users/coli171/Chaste/Output/TestParametersSensitivity/TESTPSPROUT_2NDVERSION/testPsprout0' + str(k) + "/test" + str(j) + "/results_from_time_0/results.viznodes"
        file_branchesnumber = '/Users/coli171/Chaste/Output/TestParametersSensitivity/TESTPSPROUT_2NDVERSION/testPsprout0' + str(k) + "/test" + str(j) + "/results_from_time_0/results.vizbranchesnumber"
        file_consecutivebranches = '/Users/coli171/Chaste/Output/TestParametersSensitivity/TESTPSPROUT_2NDVERSION/testPsprout0' + str(k) + "/test" + str(j) + "/results_from_time_0/results.vizconsecutivebranches"
        file_birthtimecells = '/Users/coli171/Chaste/Output/TestParametersSensitivity/TESTPSPROUT_2NDVERSION/testPsprout0' + str(k) + "/test" + str(j) + "/results_from_time_0/results.vizbirthtimecell"
        file_tortuosity = '/Users/coli171/Chaste/Output/TestParametersSensitivity/TESTPSPROUT_2NDVERSION/testPsprout0' + str(k) + "/test" + str(j) + "/results_from_time_0/results.viztortuosity"
        
        Psprout = k*10**(-2)
        list_consecutivebranches = runner.ConsecutiveBranches(file_consecutivebranches)

        computationalnumberofbranches, expectednumberofbranches, formulaexpectednumberofbranches = runner.ComparisonConsecutiveBranches(file_cellmutation, Psprout, Tcycle)

        x1.append(Psprout)
        y_totalnumbercells1.append(runner.TotalNumberCells(file_cellmutation))
        y_numbercellsplane1.append(runner.NbCellsPlane(file_nodescoordinates, plane_normal_vector))
        y_totalnumbertipcells1.append(runner.TotalNumberVesselTips(file_cellmutation))
        y_totalnumberstalkcells1.append(runner.TotalNumberStalkCells(file_cellmutation))
        y_totalnumberbranches1.append(runner.TotalNumberBranches(file_cellmutation))

        y_areaplane1.append(runner.AreaCellsPlane(file_nodescoordinates, plane_normal_vector))
        y_area1.append(runner.AverageArea(file_nodescoordinates, file_cellmutation))
        y_volume1.append(runner.AverageVolume(file_nodescoordinates, file_cellmutation))

        y_minconsecutivebranches1.append(runner.MinimumConsecutiveBranches(list_consecutivebranches))
        y_maxconsecutivebranches1.append(runner.MaximumConsecutiveBranches(list_consecutivebranches))
        y_computationalnumberofbranches1.append(computationalnumberofbranches)
        y_expectednumberofbranches1.append(expectednumberofbranches)
        y_formulaexpectednumberofbranches1.append(formulaexpectednumberofbranches)

        y_tortuosity1.append(runner.Tortuosity(file_tortuosity))
        #y_neighbouringdistance1.append(runner.ClosestNeighbouringDistanceInPlane(file_nodescoordinates, plane_normal_vector))

    for k in range(0,10):
        file_cellmutation = "/Users/coli171/Chaste/Output/TestParametersSensitivity/TESTPSPROUT_2NDVERSION/testPsprout1" + str(k) + "/test" + str(j) + "/results_from_time_0/results.vizmutationstates"
        file_nodescoordinates = "/Users/coli171/Chaste/Output/TestParametersSensitivity/TESTPSPROUT_2NDVERSION/testPsprout1" + str(k) + "/test" + str(j) + "/results_from_time_0/results.viznodes"
        file_branchesnumber = "/Users/coli171/Chaste/Output/TestParametersSensitivity/TESTPSPROUT_2NDVERSION/testPsprout1" + str(k) + "/test" + str(j) + "/results_from_time_0/results.vizbranchesnumber"
        file_consecutivebranches = "/Users/coli171/Chaste/Output/TestParametersSensitivity/TESTPSPROUT_2NDVERSION/testPsprout1" + str(k) + "/test" + str(j) + "/results_from_time_0/results.vizconsecutivebranches"
        file_birthtimecells = "/Users/coli171/Chaste/Output/TestParametersSensitivity/TESTPSPROUT_2NDVERSION/testPsprout1" + str(k) + "/test" + str(j) + "/results_from_time_0/results.vizbirthtimecell"
        file_tortuosity = "/Users/coli171/Chaste/Output/TestParametersSensitivity/TESTPSPROUT_2NDVERSION/testPsprout1" + str(k) + "/test" + str(j) + "/results_from_time_0/results.viztortuosity"

        Psprout = k*10**(-2) + 0.1
        list_consecutivebranches = runner.ConsecutiveBranches(file_consecutivebranches)
        
        computationalnumberofbranches, expectednumberofbranches, formulaexpectednumberofbranches = runner.ComparisonConsecutiveBranches(file_cellmutation, Psprout, Tcycle)

        x2.append(Psprout)
        y_totalnumbercells2.append(runner.TotalNumberCells(file_cellmutation))
        y_numbercellsplane2.append(runner.NbCellsPlane(file_nodescoordinates, plane_normal_vector))
        y_totalnumbertipcells2.append(runner.TotalNumberVesselTips(file_cellmutation))
        y_totalnumberstalkcells2.append(runner.TotalNumberStalkCells(file_cellmutation))
        y_totalnumberbranches2.append(runner.TotalNumberBranches(file_cellmutation))

        y_areaplane2.append(runner.AreaCellsPlane(file_nodescoordinates, plane_normal_vector))
        y_area2.append(runner.AverageArea(file_nodescoordinates, file_cellmutation))
        y_volume2.append(runner.AverageVolume(file_nodescoordinates, file_cellmutation))

        y_minconsecutivebranches2.append(runner.MinimumConsecutiveBranches(list_consecutivebranches))
        y_maxconsecutivebranches2.append(runner.MaximumConsecutiveBranches(list_consecutivebranches))
        y_computationalnumberofbranches2.append(computationalnumberofbranches)
        y_expectednumberofbranches2.append(expectednumberofbranches)
        y_formulaexpectednumberofbranches2.append(formulaexpectednumberofbranches)

        y_tortuosity2.append(runner.Tortuosity(file_tortuosity))
        #y_neighbouringdistance2.append(runner.ClosestNeighbouringDistanceInPlane(file_nodescoordinates, plane_normal_vector))

    x = x1 + x2
    y_totalnumbercells = y_totalnumbercells1 + y_totalnumbercells2
    y_numbercellsplane = y_numbercellsplane1 + y_numbercellsplane2
    y_totalnumbertipcells = y_totalnumbertipcells1 + y_totalnumbertipcells2
    y_totalnumberstalkcells = y_totalnumberstalkcells1 + y_totalnumberstalkcells2
    y_totalnumberbranches = y_totalnumberbranches1 + y_totalnumberbranches2
    y_areaplane = y_areaplane1 + y_areaplane2
    y_area = y_area1 + y_area2
    y_volume = y_volume1 + y_volume2
    y_minconsecutivebranches = y_minconsecutivebranches1 + y_minconsecutivebranches2
    y_maxconsecutivebranches = y_maxconsecutivebranches1 + y_maxconsecutivebranches2
    y_computationalnumberofbranches = y_computationalnumberofbranches1 + y_computationalnumberofbranches2
    y_expectednumberofbranches = y_expectednumberofbranches1 + y_expectednumberofbranches2
    y_formulaexpectednumberofbranches = y_formulaexpectednumberofbranches1 + y_formulaexpectednumberofbranches2
    y_tortuosity = y_tortuosity1 + y_tortuosity2
    #y_neighbouringdistance = y_neighbouringdistance1 + y_neighbouringdistance2

    y_totalnumbercells_average= runner.AverageFunction(y_totalnumbercells_average, y_totalnumbercells, TotalTestNb)
    y_numbercellsplane_average = runner.AverageFunction(y_numbercellsplane_average, y_numbercellsplane, TotalTestNb)
    y_totalnumbertipcells_average = runner.AverageFunction(y_totalnumbertipcells_average, y_totalnumbertipcells, TotalTestNb)
    # y_totalnumberstalkcellsmax = runner.MaxFunction(y_totalnumberstalkcellsmax, y_totalnumberstalkcells)
    # y_totalnumbertipcellsmax = runner.MaxFunction(y_totalnumbertipcellsmax, y_totalnumbertipcells)
    # y_totalnumberstalkcellsmin = runner.MinFunction(y_totalnumberstalkcellsmin, y_totalnumberstalkcells)
    # y_totalnumbertipcellsmin = runner.MinFunction(y_totalnumbertipcellsmin, y_totalnumbertipcells)
    y_totalnumberstalkcells_average = runner.AverageFunction(y_totalnumberstalkcells_average, y_totalnumberstalkcells, TotalTestNb)
    y_totalnumberbranches_average = runner.AverageFunction(y_totalnumberbranches_average, y_totalnumberbranches, TotalTestNb)
    y_areaplane_average = runner.AverageFunction(y_areaplane_average, y_areaplane, TotalTestNb)
    y_area_average = runner.AverageFunction(y_area_average, y_area, TotalTestNb)
    y_volume_average = runner.AverageFunction(y_volume_average, y_volume, TotalTestNb)
    y_minconsecutivebranches_average = runner.AverageFunction(y_minconsecutivebranches_average, y_minconsecutivebranches, TotalTestNb)
    y_maxconsecutivebranches_average = runner.AverageFunction(y_maxconsecutivebranches_average, y_maxconsecutivebranches, TotalTestNb)
    y_computationalnumberofbranches_average = runner.AverageFunction(y_computationalnumberofbranches_average, y_computationalnumberofbranches, TotalTestNb)
    y_expectednumberofbranches_average = runner.AverageFunction(y_expectednumberofbranches_average, y_expectednumberofbranches, TotalTestNb)
    # y_computationalnumberofbranchesmax = runner.MaxFunction(y_computationalnumberofbranchesmax, y_computationalnumberofbranches)
    # y_computationalnumberofbranchesmin = runner.MinFunction(y_computationalnumberofbranchesmin, y_computationalnumberofbranches)
    # y_expectednumberofbranchesmax = runner.MaxFunction(y_expectednumberofbranchesmax, y_expectednumberofbranches)
    y_formulaexpectednumberofbranches_average = runner.AverageFunction(y_formulaexpectednumberofbranches_average, y_formulaexpectednumberofbranches, TotalTestNb)
    y_tortuosity_average = runner.AverageFunction(y_tortuosity_average, y_tortuosity, TotalTestNb)
    #y_neighbouringdistance_average = runner.AverageFunction(y_neighbouringdistance_average, y_neighbouringdistance, TotalTestNb)

    # Line plot of number of branches 
    if GraphNbBranches:
        plt.scatter(x, y_expectednumberofbranches, color = 'r', alpha = 0.5)
        plt.scatter(x, y_computationalnumberofbranches, color = 'b', alpha = 0.5)
        if(j == 9):
            plt.scatter(x, y_computationalnumberofbranches, color = 'xkcd:blue', alpha = 0.7, label = 'Number of Branches for one Realisation')
        else:
            plt.scatter(x, y_computationalnumberofbranches, color = 'xkcd:blue', alpha = 0.7)

    if GraphTortuosityNbBranches:
        if(j == 9):
            plt.scatter(y_computationalnumberofbranches, y_tortuosity, color = 'xkcd:turquoise', alpha = 0.7, label = 'Tortuosity for one Realisation')
        else:
            plt.scatter(y_computationalnumberofbranches, y_tortuosity, color = 'xkcd:turquoise', alpha = 0.7)

    # Line plot of number of tip and stalk cells 
    if GraphNbTipStalkCells:
        if(j == 9):
            plt.scatter(x, y_totalnumbertipcells, color = 'xkcd:red', alpha = 0.7, label = 'Number of Tip Cells for one Realisation')
            plt.scatter(x, y_totalnumberstalkcells, color = 'xkcd:orange', alpha = 0.7, label = 'Number of Stalk Cells for one Realisation')
        else:
            plt.scatter(x, y_totalnumbertipcells, color = 'xkcd:red', alpha = 0.7)
            plt.scatter(x, y_totalnumberstalkcells, color = 'xkcd:orange', alpha = 0.7)

    if GraphNbCells:
        if(j == 9):
            plt.scatter(x, y_totalnumbercells, color = 'xkcd:red', alpha = 0.7, label = 'Total Number of Cells for one Realisation')
        else:
            plt.scatter(x, y_totalnumbercells, color = 'xkcd:red', alpha = 0.7)

    if GraphPlaneNbCells:
        if(j == 9):
            plt.scatter(x, y_numbercellsplane, color = 'xkcd:green', alpha = 0.7, label = 'Cell Density in Plane for one Realisation')
        else:
            plt.scatter(x, y_numbercellsplane, color = 'xkcd:green', alpha = 0.7)

    # if GraphAverageClosestNeighbour:
    #     if(j == 9):
    #         plt.scatter(x, y_neighbouringdistance, color = 'xkcd:turquoise', alpha = 0.7, label = 'Minimal Neighbouring Distance in Plane for one Realisation')
    #     else:
    #         plt.scatter(x, y_neighbouringdistance, color = 'xkcd:turquoise', alpha = 0.7)

    if GraphPlaneArea:
        if(j == 9):
            plt.scatter(x, y_areaplane, color = 'xkcd:green', alpha = 0.7, label = 'Number of Tip Cells for one Realisation')
        else:
            plt.scatter(x, y_areaplane, color = 'xkcd:green', alpha = 0.7)

    # Line plot of tortuosity
    if GraphTortuosity:
        if(j == 9):
            plt.scatter(x, y_tortuosity, color = 'xkcd:turquoise', alpha = 0.7, label = 'Tortuosity for one Realisation')
        else: 
            plt.scatter(x, y_tortuosity, color = 'xkcd:turquoise', alpha = 0.7)

if ComputationalError:
    y_computationalnumberofbranches_error = np.zeros((2,20))
    y_totalnumberstalkcells_error = np.zeros((2,20))
    y_totalnumbertipcells_error = np.zeros((2,20))

    for k in range(20):
        y_computationalnumberofbranches_error[0][k] = y_computationalnumberofbranches_average[k] - y_computationalnumberofbranchesmin[k]
        y_computationalnumberofbranches_error[1][k] = y_computationalnumberofbranchesmax[k] - y_computationalnumberofbranches_average[k] 

        y_totalnumberstalkcells_error[0][k] = y_totalnumberstalkcells_average[k] - y_totalnumberstalkcellsmin[k]
        y_totalnumberstalkcells_error[1][k] = y_totalnumberstalkcellsmax[k] - y_totalnumberstalkcells_average[k]

        y_totalnumbertipcells_error[0][k] = y_totalnumbertipcells_average[k] - y_totalnumbertipcellsmin[k]
        y_totalnumbertipcells_error[1][k] = y_totalnumbertipcellsmax[k] - y_totalnumbertipcells_average[k]

        if(y_computationalnumberofbranches_error[1][k] < 1e-14):
            y_computationalnumberofbranches_error[1][k] = 0

        if(y_totalnumberstalkcells_error[1][k] < 1e-14):
            y_totalnumberstalkcells_error[1][k] = 0

        if(y_totalnumbertipcells_error[1][k] < 1e-14):
            y_totalnumbertipcells_error[1][k] = 0

if StudyTimeStep:
    y_planenbcellstimestep_average = np.zeros(301)
    y_averageclosestneighbourtimestep_average = np.zeros(301)
    y_planenbcellstimestep_firstaverage = []
    y_averageclosestneighbourtimestep_firstaverage = []
    y_planenbcellstimestep_secondaverage = []
    y_averageclosestneighbourtimestep_secondaverage = []

    list_color = ['xkcd:lightblue', 'xkcd:turquoise', 'xkcd:blue', 'xkcd:magenta', 'xkcd:pink', 'xkcd:gold', 'xkcd:orange', 'xkcd:red', 'xkcd:green', 'xkcd:lightgreen']

    for k in range(1,TotalTestNb+1):
        file_path = '/Users/coli171/Chaste/Output/TestParametersSensitivity/TestAreaVolumeTortuosity/Test' + str(k) + '330'
        file_nodescoordinates = file_path + "/results_from_time_0/results.viznodes"

        y_planenbcellstimestep = []
        y_averageclosestneighbourtimestep = []

        y_planenbcellstimestep = runner.NbCellsInPlaneAtEachTimeStep(file_nodescoordinates,plane_normal_vector)
        #y_averageclosestneighbourtimestep = runner.ClosestNeighbouringDistanceInPlaneAtEachTimeStep(file_nodescoordinates, plane_normal_vector)

        # Scatter plot of number of cells in a plane for each time step
        if GraphPlaneNbCellsTimeStep:
            if(k == 9):
                t_reduced = [t[k] for k in range(200, 301)]
                y_planenbcellstimestep_reduced = [y_planenbcellstimestep[k] for k in range(200, 301)]
                plt.scatter(t_reduced, y_planenbcellstimestep_reduced, color = list_color[k-1], alpha = 0.7, label = 'Number of Cells in a Plane for one Realisation')
            else: 
                t_reduced = [t[k] for k in range(200, 301)]
                y_planenbcellstimestep_reduced = [y_planenbcellstimestep[k] for k in range(200, 301)]
                plt.scatter(t_reduced, y_planenbcellstimestep_reduced, color = list_color[k-1], alpha = 0.7)

        # Scatter plot of neighbouring distance in a plane for each time step
        if GraphAverageClosestNeighbourTimeStep:
            if(j == 9):
                plt.scatter(t, y_averageclosestneighbourtimestep, color = 'xkcd:turquoise', alpha = 0.7, label = 'Average Neighbouring Distance for one Realisation')
            else: 
                plt.scatter(t, y_averageclosestneighbourtimestep, color = 'xkcd:turquoise', alpha = 0.7)


## PLOT THE GRAPHS ##

# BRANCHES 
# Bar Plot #
# expected and computational number of branches : only for the case of Psprout variation 
if BarPlotNbBranches:
    barWidth = 0.005
    fig = plt.subplots(figsize = (12,8), dpi = 300)
    x_expected = x
    x_computational = [elem + barWidth for elem in x]
    plt.bar(x_expected, y_expectednumberofbranches_average, color = 'r', width = barWidth, edgecolor = 'grey', label = 'Expected Number of Branches')
    plt.bar(x_computational, y_computationalnumberofbranches_average, color = 'b', width = barWidth, edgecolor = 'grey', label = 'Computational Number of Branches')
    plt.xticks([r for r in x])
    plt.ylabel('Number of Branches')

# Line Plot #
if GraphNbBranches:
    # formula plot
    plt.plot(x, y_formulaexpectednumberofbranches_average, label = r'$2 \cdot (1+P_{sprout})^{\lfloor \frac{T_{end}}{T_{cycle}} \rfloor} - 1$', color = 'xkcd:green')

    # max plots
    #plt.plot(x, y_expectednumberofbranchesmax, color = 'r', label = 'Max Expected Number of Branches')
    #plt.plot(x, y_computationalnumberofbranchesmax, color = 'b', label = 'Max Computational Number of Branches')

    # average plots 
    plt.plot(x, y_expectednumberofbranches_average, color = 'xkcd:orange', linestyle = '--', label = 'Average Expected Number of Branches')
    plt.plot(x, y_computationalnumberofbranches_average, color = 'xkcd:blue', linestyle = '--', label = 'Average Number of Branches for 10 Realisations', linewidth = 5.0)

    if ComputationalError:
        plt.errorbar(x, y_computationalnumberofbranches_average, y_computationalnumberofbranches_error, color = 'b', linestyle = 'None', marker = 'o')

    plt.xticks([0.00, 0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18], fontsize = 28)
    plt.yticks(fontsize=28)
    plt.ylabel('Number of Branches', fontsize = 28)
    plt.xlabel(r'$P_{sprout}$', fontsize = 32)

# Line plot depending on time of division #
if GraphDivisionTime:
    x1 = []
    x2 = []
    y_totalnumberbranches1 = []
    y_totalnumberbranches2 = []
    y_timefirstbranch1 = []
    y_timefirstbranch2 = []

    y_timefirstbranch1_total = []
    y_totalnumberbranches1_total = []
    y_timefirstbranch2_total = []
    y_totalnumberbranches2_total = []

    i = 0
    a = 0
    list_color_2 = ['xkcd:chartreuse', 'xkcd:blue', 'xkcd:aquamarine', 'xkcd:azure',  'xkcd:turquoise'] # 'xkcd:green', 'xkcd:lime', 'xkcd:lightgreen', 'xkcd:yellowgreen', 'xkcd:lightblue',
    list_color_1 = ['xkcd:blue', 'xkcd:magenta', 'xkcd:orange', 'xkcd:red', 'xkcd:green'] # 'xkcd:yellow', 'xkcd:gold', 'xkcd:wheat', 'xkcd:beige', 'xkcd:pink', 

    if SmallValuesPsprout:
        for k in range(1, 10, 2):
            Psprout1 = k*10**(-2) 
            x1.append(Psprout1)

            y_timefirstbranch1 = []
            y_totalnumberbranches1 = []

            for j in range(1,10):
                file_branchesnumber1 = "/Users/coli171/Chaste/Output/TestParametersSensitivity/TESTPSPROUT_2NDVERSION/testPsprout0" + str(k) + "/test" + str(j) + "/results_from_time_0/results.vizbranchesnumber"
                file_cellmutation1 = "/Users/coli171/Chaste/Output/TestParametersSensitivity/TESTPSPROUT_2NDVERSION/testPsprout0" + str(k) + "/test" + str(j) + "/results_from_time_0/results.vizmutationstates"

                y_timefirstbranch1.append(runner.TimeFirstBranch(file_branchesnumber1))
                y_totalnumberbranches1.append(runner.TotalNumberBranches(file_cellmutation1))

            plt.scatter(y_timefirstbranch1, y_totalnumberbranches1, color = list_color_1[i], alpha = 0.7, s = 200.0, label = r'$P_{sprout} =$ %1.2f'%Psprout1)

            i += 1

    if BigValuesPsprout:
        for k in range(0, 10, 2):
            Psprout2 = k*10**(-2) + 0.1
            x2.append(Psprout2)

            y_timefirstbranch2 = []
            y_totalnumberbranches2 = []

            for j in range(1,10):
                file_branchesnumber2 = "/Users/coli171/Chaste/Output/TestParametersSensitivity/TESTPSPROUT_2NDVERSION/testPsprout1" + str(k) + "/test" + str(j) + "/results_from_time_0/results.vizbranchesnumber"
                file_cellmutation2 = "/Users/coli171/Chaste/Output/TestParametersSensitivity/TESTPSPROUT_2NDVERSION/testPsprout1" + str(k) + "/test" + str(j) + "/results_from_time_0/results.vizmutationstates"

                y_timefirstbranch2.append(runner.TimeFirstBranch(file_branchesnumber2))
                y_totalnumberbranches2.append(runner.TotalNumberBranches(file_cellmutation2))
                #y_timefirstbranch2_total.append(TimeFirstBranch(file_branchesnumber2))
                #y_totalnumberbranches2_total.append(TotalNumberBranches(file_cellmutation2))

            plt.scatter(y_timefirstbranch2, y_totalnumberbranches2, color = list_color_1[a], alpha = 0.8, s = 200.0, label = r'$P_{sprout} =$ %1.2f'%Psprout2)

            b2, a2 = np.polyfit(y_timefirstbranch2, y_totalnumberbranches2, deg = 1)
            xseq = np.linspace(0, 300, num = 100)

            if(k == 8):
                plt.plot(xseq , a2+b2*xseq, color=list_color_1[a], linestyle = '--', linewidth = 5.0, label = 'Linear regression', alpha = 0.4)
            else:
                plt.plot(xseq , a2+b2*xseq, color=list_color_1[a], linestyle = '--', linewidth = 5.0, alpha = 0.4)

            a += 1

    plt.xlim(-10, 300)
    plt.ylim(-5, 175)
    plt.xticks(fontsize = 28)
    plt.yticks(fontsize = 28)
    plt.xlabel('Time of First Sprouting Division', fontsize = 28)
    plt.ylabel('Number of Branches', fontsize = 28)



# CELL TYPES
# Bar Plot #
if BarPlotNbCells:
    barWidth = 0.003
    fig = plt.subplots(figsize = (12,8))
    x_tip = x 
    x_stalk = [elem + barWidth for elem in x_tip]
    x_total = [elem + barWidth for elem in x_stalk]
    plt.bar(x_tip, y_totalnumbertipcells_average, color = 'b', width = barWidth/3, edgecolor = 'grey', label = 'Number of Tip Cells')
    plt.bar(x_stalk, y_totalnumberstalkcells_average, color = 'g', width = barWidth/3, edgecolor = 'grey', label = 'Number of Stalk Cells')
    plt.bar(x_total, y_totalnumbercells_average, color = 'r', width = barWidth/3, edgecolor = 'grey', label = 'Total Number of Cells')
    plt.xticks([0.00, 0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18], fontsize = 18)
    plt.yticks(fontsize=18)
    plt.ylabel('Number of Cells', fontsize = 18)
    plt.xlabel(r'$P_{sprout}$', fontsize = 18)


# Line Plot #
if GraphNbTipStalkCells:
    plt.plot(x, y_totalnumbertipcells_average, color = 'xkcd:red', linestyle = '--', label = 'Average Number of Tip Cells', linewidth = 5.0)
    plt.plot(x, y_totalnumberstalkcells_average, color = 'xkcd:orange', linestyle = '--', label = 'Average Number of Stalk Cells', linewidth = 5.0)

    # max plots 
    # plt.plot(x, y_totalnumbertipcellsmax, color = 'b', label = 'Maximum Number of Tip Cells')
    # plt.plot(x, y_totalnumberstalkcellsmax, color = 'g', label = 'Maximum Number of Stalk Cells')

    if ComputationalError:
        plt.errorbar(x, y_totalnumbertipcells_average, y_totalnumbertipcells_error, color = 'b', linestyle = 'None', marker = 'o')
        plt.errorbar(x, y_totalnumberstalkcells_average, y_totalnumberstalkcells_error, color = 'g', linestyle = 'None', marker = 'o')

    plt.xticks([0.00, 0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18], fontsize = 28)
    plt.yticks(fontsize=28)
    plt.ylabel('Number of Cells', fontsize = 28)
    plt.xlabel(r'$P_{sprout}$', fontsize = 32)

if GraphNbCells:
    plt.plot(x, y_totalnumbercells_average, color = 'xkcd:red', linestyle = '--', label = 'Average Total Number of Cells', linewidth = 5.0)
    plt.xticks([0.00, 0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18], fontsize = 28)
    plt.yticks(fontsize=28)
    plt.ylabel('Number of Cells', fontsize = 28)
    plt.xlabel(r'$P_{sprout}$', fontsize = 32)

if GraphPlaneNbCells:
    x_graphplane = [0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14]
    y_numbercellsplane_average_graphplane = [y_numbercellsplane_average[k] for k in range(7)]
    plt.plot(x_graphplane, y_numbercellsplane_average_graphplane, color = 'xkcd:darkgreen', linestyle = '-', marker = '.', label = 'Average Cell Density in Plane', linewidth = 3.0, markersize = 12.0)
    plt.xticks([0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14], fontsize = 28)
    plt.yticks(fontsize = 28)
    plt.ylabel('Cell Density in Plane', fontsize = 28)
    plt.xlabel(r'$P_{sprout}$', fontsize = 32)

if GraphPlaneNbCellsArea:
    x_graphplane = [0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14]
    y_numbercellsplane_average_graphplane = [y_numbercellsplane_average[k] for k in range(7)]
    y_areaplane_average_graphplane = [y_areaplane_average[k] for k in range(7)]

    fig, ax1 = plt.subplots()
    ax1.set_xlabel(r'$P_{sprout}$', fontsize = 32)
    ax1.tick_params(axis='x', labelsize = 18)
    color1 = 'tab:green'
    ax1.set_ylabel('Number of Cells in Plane', fontsize = 28, color = color1)
    ax1.plot(x_graphplane, y_numbercellsplane_average_graphplane, color = color1, linestyle = '-', marker = '.', linewidth = 3.0, markersize = 12.0)
    ax1.tick_params(axis='y', labelcolor=color1, labelsize = 18)
    ax2 = ax1.twinx()
    color2 = 'tab:red'
    ax2.set_ylabel('Area of Plane', fontsize = 28, color = color2)
    ax2.plot(x_graphplane, y_areaplane_average_graphplane, color = color2, linestyle = '-', marker = '.', linewidth = 3.0, markersize = 12.0)
    ax2.tick_params(axis='y', labelcolor=color2, labelsize = 18)
    fig.tight_layout()

# GEOMETRY 
## AREA AND VOLUME ## 
if GraphAreaVolume:
    plt.plot(x, y_area_average, linestyle = '-', marker = 'x', color = 'xkcd:olive', label = 'Area', linewidth = 5.0, markersize = 15.0, markeredgewidth = 3.0)
    plt.plot(x, y_volume_average, linestyle = '-', marker = 'x', color = 'xkcd:green', label = 'Volume', linewidth = 5.0, markersize = 15.0, markeredgewidth = 3.0)
    plt.xlabel(r'$P_{sprout}$', fontsize = 32)
    plt.xticks([0.00, 0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18], fontsize = 28)
    plt.yticks(fontsize = 28)

if GraphPlaneArea:
    x_graphplane = [0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14]
    y_areaplane_average_graphplane = [y_areaplane_average[k] for k in range(7)]
    plt.plot(x_graphplane, y_areaplane_average_graphplane, linestyle = '-', marker = 'x', color = 'xkcd:olive', label = 'Area', linewidth = 5.0, markersize = 15.0, markeredgewidth = 3.0)
    plt.xlabel(r'$P_{sprout}$', fontsize = 32)
    plt.xticks([0.08, 0.1, 0.12, 0.14], fontsize = 28)
    plt.yticks([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10], fontsize = 28)
    plt.yscale('log')

# minimum number of consecutive branches 
if GraphMinNbConsecutiveBranches:
    plt.ylabel('Minimum number of consecutive branches')
    plt.plot(x, y_minconsecutivebranches_average)
    plt.xlabel('Psprout')

# maximum number of consecutive branches 
if GraphMaxNbConsecutiveBranches:
    plt.ylabel('Maximum number of consecutive branches')
    plt.plot(x, y_maxconsecutivebranches_average)
    plt.xlabel('Psprout')

## TORTUOSITY ## 
if GraphTortuosity:
    plt.ylabel('Tortuosity', fontsize = 28)
    plt.plot(x, y_tortuosity_average, linestyle = '--', color = 'xkcd:teal', label = 'Average Tortuosity', linewidth = 5.0)
    plt.xlabel(r'$P_{sprout}$', fontsize = 32)
    plt.xticks([0.00, 0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18], fontsize = 28)
    plt.yticks(fontsize = 28)

if GraphTortuosityNbBranches:
    plt.ylabel('Tortuosity', fontsize = 28)
    #plt.plot(y_computationalnumberofbranches_average, y_tortuosity_average, color = 'xkcd:teal', label = 'Average Tortuosity', linestyle = '-', marker='.', linewidth = 3.0, markersize = 13.0)
    plt.xlabel('Total Number of Branches', fontsize = 32)
    plt.xticks(fontsize = 28)
    plt.yticks(fontsize = 28)

# Neighbouring Distance
if GraphAverageClosestNeighbour:
    x_graphplane = [0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14]
    y_neighbouringdistance_average_graphplane = [y_neighbouringdistance_average[k]for k in range(7)]
    plt.plot(x_graphplane, y_neighbouringdistance_average_graphplane, linestyle = '-', marker = '.', markersize = 15.0, color = 'xkcd:teal', linewidth = 3.0, label = 'Average Minimal Neighbouring Distance in Plane')
    plt.xlabel(r'$P_{sprout}$', fontsize = 32)
    plt.xticks([0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14], fontsize = 28)
    plt.ylabel('Distance to Neighbour', fontsize = 28)
    plt.yticks(fontsize = 28)


# Time dependent 
if GraphPlaneNbCellsTimeStep:
    plt.ylabel('Number of Cells in Plane', fontsize = 28)
    #plt.plot(t, y_planenbcellstimestep_average, linestyle = '--', label='Average Number of Cells', color = 'xkcd:black', linewidth = 5.0)
    plt.xlabel('t', fontsize = 32)
    #plt.ylim([0, 1.4])
    plt.xticks(fontsize = 28)
    plt.yticks(fontsize = 28)
    #plt.yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0], fontsize = 28)

if GraphAverageClosestNeighbourTimeStep:
    plt.ylabel('Distance to Neighbour', fontsize = 28)
    plt.plot(t, y_averageclosestneighbourtimestep_average, linestyle = '--', color = 'xkcd:teal', linewidth = 5.0)
    plt.xlabel('t', fontsize = 32)
    plt.xticks([0, 50, 100, 150, 200, 250, 300], fontsize = 28)
    plt.yticks(fontsize = 28)


## PLOT ## 
# plt.legend(loc = 'upper left', fontsize = 24, scatterpoints = 1)
plt.legend(loc = 'upper right', fontsize = 22, ncol = 1)
# plt.legend(loc = 'lower left', fontsize = 18)
plt.show()