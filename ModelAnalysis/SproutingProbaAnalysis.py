# Libraries 
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

from matplotlib import cm
from matplotlib.ticker import LinearLocator 
from matplotlib.colors import LinearSegmentedColormap

from scipy.spatial import ConvexHull

import ParametersSensitivityRunner

# Parameters 
Tcycle = 13.0
time = 300

y_totalnumbercells_average = np.zeros(20)
y_totalnumbertipcells_average = np.zeros(20)
y_totalnumberstalkcells_average = np.zeros(20)
y_totalnumbertipcellsmax = np.zeros(20)
y_totalnumberstalkcellsmax = np.zeros(20)
y_totalnumbertipcellsmin = np.zeros(20)
y_totalnumberstalkcellsmin = np.zeros(20)
y_totalnumberbranches_average = np.zeros(20)
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

for j in range(1,10):

    x1 = []
    y_totalnumbercells1 = []
    y_totalnumbertipcells1 = []
    y_totalnumberstalkcells1 = []
    y_totalnumberbranches1 = []
    y_area1 = []
    y_volume1 = []
    y_minconsecutivebranches1 = []
    y_maxconsecutivebranches1 = []
    y_computationalnumberofbranches1 = []
    y_expectednumberofbranches1 = []
    y_formulaexpectednumberofbranches1 = []
    y_tortuosity1 = []

    x2 = []
    y_totalnumbercells2 = []
    y_totalnumbertipcells2 = []
    y_totalnumberstalkcells2 = []
    y_totalnumberbranches2 = []
    y_area2 = []
    y_volume2 = []
    y_minconsecutivebranches2 = []
    y_maxconsecutivebranches2 = []
    y_computationalnumberofbranches2 = []
    y_expectednumberofbranches2 = []
    y_formulaexpectednumberofbranches2 = []
    y_tortuosity2 = []

    for k in range(10):
        file_cellmutation = '/Users/coli171/Chaste/Output/TestParametersSensitivity/TESTPSPROUT_2NDVERSION/testPsprout0' + str(k) + "/test" + str(j) + "/results_from_time_0/results.vizmutationstates"
        file_nodescoordinates = '/Users/coli171/Chaste/Output/TestParametersSensitivity/TESTPSPROUT_2NDVERSION/testPsprout0' + str(k) + "/test" + str(j) + "/results_from_time_0/results.viznodes"
        file_branchesnumber = '/Users/coli171/Chaste/Output/TestParametersSensitivity/TESTPSPROUT_2NDVERSION/testPsprout0' + str(k) + "/test" + str(j) + "/results_from_time_0/results.vizbranchesnumber"
        file_consecutivebranches = '/Users/coli171/Chaste/Output/TestParametersSensitivity/TESTPSPROUT_2NDVERSION/testPsprout0' + str(k) + "/test" + str(j) + "/results_from_time_0/results.vizconsecutivebranches"
        file_birthtimecells = '/Users/coli171/Chaste/Output/TestParametersSensitivity/TESTPSPROUT_2NDVERSION/testPsprout0' + str(k) + "/test" + str(j) + "/results_from_time_0/results.vizbirthtimecell"
        file_tortuosity = '/Users/coli171/Chaste/Output/TestParametersSensitivity/TESTPSPROUT_2NDVERSION/testPsprout0' + str(k) + "/test" + str(j) + "/results_from_time_0/results.viztortuosity"
        
        Psprout = k*10**(-2)
        list_consecutivebranches = ConsecutiveBranches(file_consecutivebranches)

        computationalnumberofbranches, expectednumberofbranches, formulaexpectednumberofbranches = ComparisonConsecutiveBranches(file_cellmutation, Psprout, Tcycle)

        x1.append(Psprout)
        y_totalnumbercells1.append(TotalNumberCells(file_cellmutation))
        y_totalnumbertipcells1.append(TotalNumberTipCells(file_cellmutation))
        y_totalnumberstalkcells1.append(TotalNumberStalkCells(file_cellmutation))
        y_totalnumberbranches1.append(TotalNumberBranches(file_cellmutation))

        y_area1.append(AverageArea(file_nodescoordinates, file_cellmutation))
        y_volume1.append(AverageVolume(file_nodescoordinates, file_cellmutation))

        y_minconsecutivebranches1.append(MinimumConsecutiveBranches(list_consecutivebranches))
        y_maxconsecutivebranches1.append(MaximumConsecutiveBranches(list_consecutivebranches))
        y_computationalnumberofbranches1.append(computationalnumberofbranches)
        y_expectednumberofbranches1.append(expectednumberofbranches)
        y_formulaexpectednumberofbranches1.append(formulaexpectednumberofbranches)

        y_tortuosity1.append(Tortuosity(file_tortuosity))

    for k in range(10):
        file_cellmutation = "/Users/coli171/Chaste/Output/TestParametersSensitivity/TESTPSPROUT_2NDVERSION/testPsprout1" + str(k) + "/test" + str(j) + "/results_from_time_0/results.vizmutationstates"
        file_nodescoordinates = "/Users/coli171/Chaste/Output/TestParametersSensitivity/TESTPSPROUT_2NDVERSION/testPsprout1" + str(k) + "/test" + str(j) + "/results_from_time_0/results.viznodes"
        file_branchesnumber = "/Users/coli171/Chaste/Output/TestParametersSensitivity/TESTPSPROUT_2NDVERSION/testPsprout1" + str(k) + "/test" + str(j) + "/results_from_time_0/results.vizbranchesnumber"
        file_consecutivebranches = "/Users/coli171/Chaste/Output/TestParametersSensitivity/TESTPSPROUT_2NDVERSION/testPsprout1" + str(k) + "/test" + str(j) + "/results_from_time_0/results.vizconsecutivebranches"
        file_birthtimecells = "/Users/coli171/Chaste/Output/TestParametersSensitivity/TESTPSPROUT_2NDVERSION/testPsprout1" + str(k) + "/test" + str(j) + "/results_from_time_0/results.vizbirthtimecell"
        file_tortuosity = "/Users/coli171/Chaste/Output/TestParametersSensitivity/TESTPSPROUT_2NDVERSION/testPsprout1" + str(k) + "/test" + str(j) + "/results_from_time_0/results.viztortuosity"

        Psprout = k*10**(-2) + 0.1
        list_consecutivebranches = ConsecutiveBranches(file_consecutivebranches)
        
        computationalnumberofbranches, expectednumberofbranches, formulaexpectednumberofbranches = ComparisonConsecutiveBranches(file_cellmutation, Psprout, Tcycle)

        x2.append(Psprout)
        y_totalnumbercells2.append(TotalNumberCells(file_cellmutation))
        y_totalnumbertipcells2.append(TotalNumberTipCells(file_cellmutation))
        y_totalnumberstalkcells2.append(TotalNumberStalkCells(file_cellmutation))
        y_totalnumberbranches2.append(TotalNumberBranches(file_cellmutation))

        y_area2.append(AverageArea(file_nodescoordinates, file_cellmutation))
        y_volume2.append(AverageVolume(file_nodescoordinates, file_cellmutation))

        y_minconsecutivebranches2.append(MinimumConsecutiveBranches(list_consecutivebranches))
        y_maxconsecutivebranches2.append(MaximumConsecutiveBranches(list_consecutivebranches))
        y_computationalnumberofbranches2.append(computationalnumberofbranches)
        y_expectednumberofbranches2.append(expectednumberofbranches)
        y_formulaexpectednumberofbranches2.append(formulaexpectednumberofbranches)

        y_tortuosity2.append(Tortuosity(file_tortuosity))

    x = x1 + x2
    y_totalnumbercells = y_totalnumbercells1 + y_totalnumbercells2
    y_totalnumbertipcells = y_totalnumbertipcells1 + y_totalnumbertipcells2
    y_totalnumberstalkcells = y_totalnumberstalkcells1 + y_totalnumberstalkcells2
    y_totalnumberbranches = y_totalnumberbranches1 + y_totalnumberbranches2
    y_area = y_area1 + y_area2
    y_volume = y_volume1 + y_volume2
    y_minconsecutivebranches = y_minconsecutivebranches1 + y_minconsecutivebranches2
    y_maxconsecutivebranches = y_maxconsecutivebranches1 + y_maxconsecutivebranches2
    y_computationalnumberofbranches = y_computationalnumberofbranches1 + y_computationalnumberofbranches2
    y_expectednumberofbranches = y_expectednumberofbranches1 + y_expectednumberofbranches2
    y_formulaexpectednumberofbranches = y_formulaexpectednumberofbranches1 + y_formulaexpectednumberofbranches2
    y_tortuosity = y_tortuosity1 + y_tortuosity2

    y_totalnumbercells_average= AverageFunction(y_totalnumbercells_average, y_totalnumbercells)
    y_totalnumbertipcells_average = AverageFunction(y_totalnumbertipcells_average, y_totalnumbertipcells)
    y_totalnumberstalkcellsmax = MaxFunction(y_totalnumberstalkcellsmax, y_totalnumberstalkcells)
    y_totalnumbertipcellsmax = MaxFunction(y_totalnumbertipcellsmax, y_totalnumbertipcells)
    y_totalnumberstalkcellsmin = MinFunction(y_totalnumberstalkcellsmin, y_totalnumberstalkcells)
    y_totalnumbertipcellsmin = MinFunction(y_totalnumbertipcellsmin, y_totalnumbertipcells)
    y_totalnumberstalkcells_average = AverageFunction(y_totalnumberstalkcells_average, y_totalnumberstalkcells)
    y_totalnumberbranches_average = AverageFunction(y_totalnumberbranches_average, y_totalnumberbranches)
    y_area_average = AverageFunction(y_area_average, y_area)
    y_volume_average = AverageFunction(y_volume_average, y_volume)
    y_minconsecutivebranches_average = AverageFunction(y_minconsecutivebranches_average, y_minconsecutivebranches)
    y_maxconsecutivebranches_average = AverageFunction(y_maxconsecutivebranches_average, y_maxconsecutivebranches)
    y_computationalnumberofbranches_average = AverageFunction(y_computationalnumberofbranches_average, y_computationalnumberofbranches)
    y_expectednumberofbranches_average = AverageFunction(y_expectednumberofbranches_average, y_expectednumberofbranches)
    y_computationalnumberofbranchesmax = MaxFunction(y_computationalnumberofbranchesmax, y_computationalnumberofbranches)
    y_computationalnumberofbranchesmin = MinFunction(y_computationalnumberofbranchesmin, y_computationalnumberofbranches)
    y_expectednumberofbranchesmax = MaxFunction(y_expectednumberofbranchesmax, y_expectednumberofbranches)
    y_formulaexpectednumberofbranches_average = AverageFunction(y_formulaexpectednumberofbranches_average, y_formulaexpectednumberofbranches)
    y_tortuosity_average = AverageFunction(y_tortuosity_average, y_tortuosity)


    # Line plot of number of branches 
    # plt.scatter(x, y_expectednumberofbranches, color = 'r', alpha = 0.5)
    # plt.scatter(x, y_computationalnumberofbranches, color = 'b', alpha = 0.5)
    # if(j == 9):
    #     plt.scatter(x, y_computationalnumberofbranches, color = 'xkcd:blue', alpha = 0.7, label = 'Number of Branches for one Realisation')
    # else:
    #     plt.scatter(x, y_computationalnumberofbranches, color = 'xkcd:blue', alpha = 0.7)

    # Line plot of number of cells 
    # if(j == 9):
    #     plt.scatter(x, y_totalnumbertipcells, color = 'xkcd:red', alpha = 0.7, label = 'Number of Tip Cells for one Realisation')
    #     plt.scatter(x, y_totalnumberstalkcells, color = 'xkcd:orange', alpha = 0.7, label = 'Number of Stalk Cells for one Realisation')
    # else:
    #     plt.scatter(x, y_totalnumbertipcells, color = 'xkcd:red', alpha = 0.7)
    #     plt.scatter(x, y_totalnumberstalkcells, color = 'xkcd:orange', alpha = 0.7)

    # Line plot of tortuosity
    # if(j == 9):
    #     plt.scatter(x, y_tortuosity, color = 'xkcd:turquoise', alpha = 0.7, label = 'Tortuosity for one Realisation')
    # else: 
    #     plt.scatter(x, y_tortuosity, color = 'xkcd:turquoise', alpha = 0.7)

# y_computationalnumberofbranches_error = np.zeros((2,20))
# y_totalnumberstalkcells_error = np.zeros((2,20))
# y_totalnumbertipcells_error = np.zeros((2,20))

# for k in range(20):
#     y_computationalnumberofbranches_error[0][k] = y_computationalnumberofbranches_average[k] - y_computationalnumberofbranchesmin[k]
#     y_computationalnumberofbranches_error[1][k] = y_computationalnumberofbranchesmax[k] - y_computationalnumberofbranches_average[k] 

#     y_totalnumberstalkcells_error[0][k] = y_totalnumberstalkcells_average[k] - y_totalnumberstalkcellsmin[k]
#     y_totalnumberstalkcells_error[1][k] = y_totalnumberstalkcellsmax[k] - y_totalnumberstalkcells_average[k]

#     y_totalnumbertipcells_error[0][k] = y_totalnumbertipcells_average[k] - y_totalnumbertipcellsmin[k]
#     y_totalnumbertipcells_error[1][k] = y_totalnumbertipcellsmax[k] - y_totalnumbertipcells_average[k]

#     if(y_computationalnumberofbranches_error[1][k] < 1e-14):
#         y_computationalnumberofbranches_error[1][k] = 0

#     if(y_totalnumberstalkcells_error[1][k] < 1e-14):
#         y_totalnumberstalkcells_error[1][k] = 0

#     if(y_totalnumbertipcells_error[1][k] < 1e-14):
#         y_totalnumbertipcells_error[1][k] = 0


## PLOT THE GRAPHS ##

# BRANCHES 
# Bar Plot #
# expected and computational number of branches : only for the case of Psprout variation 
# barWidth = 0.005
# fig = plt.subplots(figsize = (12,8), dpi = 300)
# x_expected = x
# x_computational = [elem + barWidth for elem in x]
# plt.bar(x_expected, y_expectednumberofbranches_average, color = 'r', width = barWidth, edgecolor = 'grey', label = 'Expected Number of Branches')
# plt.bar(x_computational, y_computationalnumberofbranches_average, color = 'b', width = barWidth, edgecolor = 'grey', label = 'Computational Number of Branches')
# plt.xticks([r for r in x])
# plt.ylabel('Number of Branches')

# Line Plot #
# plt.plot(x, y_formulaexpectednumberofbranches_average, label = r'$2 \cdot (1+P_{sprout})^{\lfloor \frac{T_{end}}{T_{cycle}} \rfloor} - 1$', color = 'xkcd:green')
# plt.plot(x, y_expectednumberofbranchesmax, color = 'r', label = 'Max Expected Number of Branches')
# plt.plot(x, y_computationalnumberofbranchesmax, color = 'b', label = 'Max Computational Number of Branches')
# plt.plot(x, y_expectednumberofbranches_average, color = 'xkcd:orange', linestyle = '--', label = 'Average Expected Number of Branches')
# plt.plot(x, y_computationalnumberofbranches_average, color = 'xkcd:blue', linestyle = '--', label = 'Average Number of Branches for 10 Realisations', linewidth = 5.0)
# plt.plot(x, y_formulaexpectednumberofbranches_average, label = 'Expected Number of Branches', color = 'xkcd:green', linewidth = 5.0)
# plt.errorbar(x, y_computationalnumberofbranches_average, y_computationalnumberofbranches_error, color = 'b', linestyle = 'None', marker = 'o')
# plt.xticks([0.00, 0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18], fontsize = 28)
# plt.yticks(fontsize=28)
# plt.ylabel('Number of Branches', fontsize = 28)
# plt.xlabel(r'$P_{sprout}$', fontsize = 32)

# Line plot depending on time of division #
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


# for k in range(1, 10, 2):
#     Psprout1 = k*10**(-2) 
#     x1.append(Psprout1)

#     y_timefirstbranch1 = []
#     y_totalnumberbranches1 = []

#     for j in range(1,10):
#         file_branchesnumber1 = "/Users/coli171/Chaste/Output/TestParametersSensitivity/TESTPSPROUT_2NDVERSION/testPsprout0" + str(k) + "/test" + str(j) + "/results_from_time_0/results.vizbranchesnumber"
#         file_cellmutation1 = "/Users/coli171/Chaste/Output/TestParametersSensitivity/TESTPSPROUT_2NDVERSION/testPsprout0" + str(k) + "/test" + str(j) + "/results_from_time_0/results.vizmutationstates"

#         y_timefirstbranch1.append(TimeFirstBranch(file_branchesnumber1))
#         y_totalnumberbranches1.append(TotalNumberBranches(file_cellmutation1))

#     plt.scatter(y_timefirstbranch1, y_totalnumberbranches1, color = list_color_1[i], alpha = 0.7, s = 200.0, label = r'$P_{sprout} =$ %1.2f'%Psprout1)

#     i += 1

for k in range(0, 10, 2):
    Psprout2 = k*10**(-2) + 0.1
    x2.append(Psprout2)

    y_timefirstbranch2 = []
    y_totalnumberbranches2 = []

    for j in range(1,10):
        file_branchesnumber2 = "/Users/coli171/Chaste/Output/TestParametersSensitivity/TESTPSPROUT_2NDVERSION/testPsprout1" + str(k) + "/test" + str(j) + "/results_from_time_0/results.vizbranchesnumber"
        file_cellmutation2 = "/Users/coli171/Chaste/Output/TestParametersSensitivity/TESTPSPROUT_2NDVERSION/testPsprout1" + str(k) + "/test" + str(j) + "/results_from_time_0/results.vizmutationstates"

        y_timefirstbranch2.append(TimeFirstBranch(file_branchesnumber2))
        y_totalnumberbranches2.append(TotalNumberBranches(file_cellmutation2))
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
# barWidth = 0.003
# fig = plt.subplots(figsize = (12,8))
# x_tip = x 
# x_stalk = [elem + barWidth for elem in x_tip]
# x_total = [elem + barWidth for elem in x_stalk]
# plt.bar(x_tip, y_totalnumbertipcells_average, color = 'b', width = barWidth/3, edgecolor = 'grey', label = 'Number of Tip Cells')
# plt.bar(x_stalk, y_totalnumberstalkcells_average, color = 'g', width = barWidth/3, edgecolor = 'grey', label = 'Number of Stalk Cells')
# plt.bar(x_total, y_totalnumbercells_average, color = 'r', width = barWidth/3, edgecolor = 'grey', label = 'Total Number of Cells')
# plt.xticks([0.00, 0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18], fontsize = 18)
# plt.yticks(fontsize=18)
# plt.ylabel('Number of Cells', fontsize = 18)
# plt.xlabel(r'$P_{sprout}$', fontsize = 18)


# Line Plot #
# plt.plot(x, y_totalnumbertipcellsmax, color = 'b', label = 'Maximum Number of Tip Cells')
# plt.plot(x, y_totalnumberstalkcellsmax, color = 'g', label = 'Maximum Number of Stalk Cells')
# plt.plot(x, y_totalnumbertipcells_average, color = 'xkcd:red', linestyle = '--', label = 'Average Number of Tip Cells', linewidth = 5.0)
# plt.plot(x, y_totalnumberstalkcells_average, color = 'xkcd:orange', linestyle = '--', label = 'Average Number of Stalk Cells', linewidth = 5.0)
# plt.errorbar(x, y_totalnumbertipcells_average, y_totalnumbertipcells_error, color = 'b', linestyle = 'None', marker = 'o')
# plt.errorbar(x, y_totalnumberstalkcells_average, y_totalnumberstalkcells_error, color = 'g', linestyle = 'None', marker = 'o')
# plt.xticks([0.00, 0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18], fontsize = 28)
# plt.yticks(fontsize=28)
# plt.ylabel('Number of Cells', fontsize = 28)
# plt.xlabel(r'$P_{sprout}$', fontsize = 32)


# GEOMETRY 
## AREA AND VOLUME ## 
# plt.plot(x, y_area_average, linestyle = '-', marker = 'x', color = 'xkcd:olive', label = 'Area', linewidth = 5.0, markersize = 15.0, markeredgewidth = 3.0)
# plt.plot(x, y_volume_average, linestyle = '-', marker = 'x', color = 'xkcd:green', label = 'Volume', linewidth = 5.0, markersize = 15.0, markeredgewidth = 3.0)
# plt.xlabel(r'$P_{sprout}$', fontsize = 32)
# plt.xticks([0.00, 0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18], fontsize = 28)
# plt.yticks(fontsize = 28)

# minimum number of consecutive branches 
# plt.ylabel('Minimum number of consecutive branches')
# plt.plot(x, y_minconsecutivebranches_average)
# plt.xlabel('Psprout')

# maximum number of consecutive branches 
# plt.ylabel('Maximum number of consecutive branches')
# plt.plot(x, y_maxconsecutivebranches_average)
# plt.xlabel('Psprout')

## TORTUOSITY ## 
# plt.ylabel('Tortuosity', fontsize = 28)
# plt.plot(x, y_tortuosity_average, linestyle = '--', color = 'xkcd:teal', label = 'Average Tortuosity', linewidth = 5.0)
# plt.xlabel(r'$P_{sprout}$', fontsize = 32)
# plt.xticks([0.00, 0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18], fontsize = 28)
# plt.yticks(fontsize = 28)

## PLOT ## 
# plt.legend(loc = 'upper left', fontsize = 24, scatterpoints = 1)
plt.legend(loc = 'upper right', fontsize = 24, ncol = 1)
# plt.legend(loc = 'lower left', fontsize = 18)
plt.show()