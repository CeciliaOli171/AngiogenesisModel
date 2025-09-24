## Libraries and Runners
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
file_path = '/Users/coli171/Chaste/Output/TestParametersSensitivity/TestAreaVolumeTortuosity/Test'

## Settings 
MultipleRandomSeedTest = True
RandomSeed = 5 # if study for one random seed 

# no constraints !
Graph3D = False
AreaGraph3D = False 
VolumeGraph3D = False 
TortuosityGraph3D = False 
AreaColormap = False
VolumeColormap = False
TortuosityColormap = False

# constraint one parameter and study the two others 
TestCounter2DSigma = 2
TestCounter2DOmegap = 2
TestCounter2DOmegaa = 4
Graph2DSigmaOmegap = False
Graph2DSigmaOmegaa = False
Graph2DOmegapOmegaa= False
TotalNbCellsGraph2DSigmaOmegap = False
TotalNbCellsGraph2DSigmaOmegaa = False
TotalNbCellsGraph2DOmegapOmegaa = False
AreaGraph2DSigmaOmegap = False
AreaGraph2DSigmaOmegaa = False
AreaGraph2DOmegapOmegaa= False
VolumeGraph2DSigmaOmegap = False
VolumeGraph2DSigmaOmegaa = False
VolumeGraph2DOmegapOmegaa= False
TortuosityGraph2DSigmaOmegap = False
TortuosityGraph2DSigmaOmegaa = False
TortuosityGraph2DOmegapOmegaa= False

# constraint two parameters and study the third 
TestCounter1DSigma = 3
TestCounter1DOmegap = 3
TestCounter1DOmegaa = 1
AreaVolumeGraph1DSigma = False
AreaVolumeGraph1DOmegap= False
AreaVolumeGraph1DOmegaa = False
NbCellsPlane1DSigma = False
NbCellsPlane1DOmegap = False
NbCellsPlane1DOmegaa = False
NeighbouringDistancePlane1DSigma = False
NeighbouringDistancePlane1DOmegap = False
NeighbouringDistancePlane1DOmegaa = False 
GraphNbCellsAreaPlane1DSigma = True
GraphNbCellsAreaPlane1DOmegap = False
GraphNbCellsAreaPlane1DOmegaa = False 

## Parameters  
Tcycle = 13.0
time = 300
TotalTestNb = 10
plane_normal_vector = [18, 0, 0, 0]

list_sigma = np.array([1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10])
list_omegap = np.array([1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10])
list_omegaa = np.array([1.0, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8])

## Graphs 

## 3D GRAPHS ##

# Initialisation of the lists
x_sigma = []
x_omegap = []
x_omegaa = []
y_totalnumbercells = []
y_totalnumbertipcells = []
y_totalnumberstalkcells = []
y_totalnumberbranches = []
y_area = []
y_volume = []
y_tortuosity = []

list_color_area = []
list_color_volume = []
list_color_tortuosity = []

if Graph3D: 
    if MultipleRandomSeedTest:
        for i in range(0,8,2):
            for j in range(2,11,2):
                for k in range(2,11,2):
                    x_omegaa.append(10**(-i))
                    x_omegap.append(10**(-j))
                    x_sigma.append(10**(-k))

                    for l in range(1,11,1):
                        # initialisation 
                        value_area_multiple_seeds = []
                        value_volume_multiple_seeds = []
                        value_tortuosity_multiple_seeds = []

                        # we read the files 
                        main_file_path = file_path + str(l) + str(k) + str(j) + str(i) + "/results_from_time_0/"

                        file_cellmutation = main_file_path + "results.vizmutationstates"
                        file_nodescoordinates = main_file_path + "results.viznodes"
                        file_branchesnumber = main_file_path + "results.vizbranchesnumber"
                        file_consecutivebranches = main_file_path + "results.vizconsecutivebranches"
                        file_birthtimecells = main_file_path + "results.vizbirthtimecell"
                        file_tortuosity = main_file_path + "results.viztortuosity"

                        # we collect the results
                        list_consecutivebranches = runner.ConsecutiveBranches(file_consecutivebranches)

                        y_totalnumbercells.append(runner.TotalNumberCells(file_cellmutation))
                        y_totalnumbertipcells.append(runner.TotalNumberVesselTips(file_cellmutation))
                        y_totalnumberstalkcells.append(runner.TotalNumberStalkCells(file_cellmutation))
                        y_totalnumberbranches.append(runner.TotalNumberBranches(file_cellmutation))

                        value_area_seed = runner.AverageArea(file_nodescoordinates, file_cellmutation)
                        value_volume_seed = runner.AverageVolume(file_nodescoordinates, file_cellmutation)
                        value_tortuosity_seed = runner.Tortuosity(file_tortuosity)

                        value_area_multiple_seeds.append(value_area_seed)
                        value_volume_multiple_seeds.append(value_volume_seed)
                        value_tortuosity_multiple_seeds.append(value_tortuosity_seed)

                    value_area = stats.mean(value_area_multiple_seeds)
                    value_volume = stats.mean(value_volume_multiple_seeds)
                    value_tortuosity = stats.mean(value_tortuosity_multiple_seeds)

                    y_area.append(value_area)
                    y_volume.append(value_volume)
                    y_tortuosity.append(value_tortuosity)

    else:
        for i in range(0,8,2):
            for j in range(2,11,2):
                for k in range(2,11,2):
                    x_omegaa.append(10**(-i))
                    x_omegap.append(10**(-j))
                    x_sigma.append(10**(-k))

                    # we read the files 
                    main_file_path = file_path + str(RandomSeed) + str(k) + str(j) + str(i) + "/results_from_time_0/"

                    file_cellmutation = main_file_path + "results.vizmutationstates"
                    file_nodescoordinates = main_file_path + "results.viznodes"
                    file_branchesnumber = main_file_path + "results.vizbranchesnumber"
                    file_consecutivebranches = main_file_path + "results.vizconsecutivebranches"
                    file_birthtimecells = main_file_path + "results.vizbirthtimecell"
                    file_tortuosity = main_file_path + "results.viztortuosity"

                    # we collect the results
                    list_consecutivebranches = runner.ConsecutiveBranches(file_consecutivebranches)

                    y_totalnumbercells.append(runner.TotalNumberCells(file_cellmutation))
                    y_totalnumbertipcells.append(runner.TotalNumberVesselTips(file_cellmutation))
                    y_totalnumberstalkcells.append(runner.TotalNumberStalkCells(file_cellmutation))
                    y_totalnumberbranches.append(runner.TotalNumberBranches(file_cellmutation))

                    value_area = runner.AverageArea(file_nodescoordinates, file_cellmutation)
                    value_volume = runner.AverageVolume(file_nodescoordinates, file_cellmutation)
                    value_tortuosity = runner.Tortuosity(file_tortuosity)

                    y_area.append(value_area)
                    y_volume.append(value_volume)
                    y_tortuosity.append(value_tortuosity)


    for value_area in y_area:
        if(value_area < 126):
            list_color_area.append('xkcd:lightgreen')
        elif(126 <= value_area < 160):
            list_color_area.append('xkcd:green')
        elif(160 <= value_area < 190):
            list_color_area.append('xkcd:lightblue')
        elif(190 <= value_area < 210):
            list_color_area.append('xkcd:blue')
        elif(210 <= value_area < 340):
            list_color_area.append('xkcd:orange')
        else:
            list_color_area.append('xkcd:red')

    for value_volume in y_volume:
        if(value_volume < 200):
            list_color_volume.append('xkcd:lightgreen')
        elif(200 <= value_volume < 210):
            list_color_volume.append('xkcd:green')
        elif(360 <= value_volume < 420):
            list_color_volume.append('xkcd:lightblue')
        elif(450 <= value_volume < 575):
            list_color_volume.append('xkcd:blue')
        elif(790 <= value_volume < 1000):
            list_color_volume.append('xkcd:orange')
        else:
            list_color_volume.append('xkcd:red')

    for value_tortuosity in y_tortuosity:
        if(value_tortuosity < 1.2):
            list_color_tortuosity.append('xkcd:lightgreen')
        elif(1.2 <= value_tortuosity < 1.4):
            list_color_tortuosity.append('xkcd:green')
        elif(1.4 <= value_tortuosity < 1.6):
            list_color_tortuosity.append('xkcd:lightblue')
        elif(1.6 <= value_tortuosity < 1.8):
            list_color_tortuosity.append('xkcd:blue')
        elif(1.8 <= value_tortuosity < 2.0):
            list_color_tortuosity.append('xkcd:orange')
        else:
            list_color_tortuosity.append('xkcd:red')

# Colormaps
if TortuosityColormap:
    list_color = ['xkcd:lightgreen', 'xkcd:white', 'xkcd:green', 'xkcd:white', 'xkcd:lightblue', 'xkcd:white', 'xkcd:blue', 'xkcd:white', 'xkcd:orange', 'xkcd:white', 'xkcd:red']
    fig, ax = plt.subplots(figsize=(1,6), layout ='constrained', dpi=300)
    bounds_tortuosity = [1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2]
    cmap_tortuosity = mpl.colors.ListedColormap(list_color)
    norm_tortuosity = mpl.colors.BoundaryNorm(bounds_tortuosity, cmap_tortuosity.N)
    fig.colorbar(mpl.cm.ScalarMappable(cmap=cmap_tortuosity, norm=norm_tortuosity), cax=ax, orientation = 'vertical', spacing = 'proportional')
    ax.tick_params(labelsize=10)

if AreaColormap:
    list_color = ['xkcd:lightgreen', 'xkcd:white', 'xkcd:green', 'xkcd:white', 'xkcd:lightblue', 'xkcd:white', 'xkcd:blue', 'xkcd:white', 'xkcd:orange', 'xkcd:white', 'xkcd:red']
    fig, ax = plt.subplots(figsize=(1,6), layout ='constrained', dpi=300)
    bounds_area = [100, 120, 160, 190, 210, 340, 410]
    cmap_area = mpl.colors.ListedColormap(list_color)
    norm_area = mpl.colors.BoundaryNorm(bounds_area, cmap_area.N)
    fig.colorbar(mpl.cm.ScalarMappable(cmap=cmap_area, norm=norm_area), cax=ax, orientation = 'vertical', spacing = 'proportional')
    ax.tick_params(labelsize=10)

if VolumeColormap:
    list_color = ['xkcd:lightgreen', 'xkcd:green', 'xkcd:white', 'xkcd:lightblue', 'xkcd:white', 'xkcd:blue', 'xkcd:white', 'xkcd:orange', 'xkcd:white', 'xkcd:red']
    fig, ax = plt.subplots(figsize=(1,6), layout ='constrained', dpi=300)
    bounds_volume = [100, 200, 360, 420, 450, 575, 790, 1000, 1150]
    cmap_volume = mpl.colors.ListedColormap(list_color)
    norm_volume = mpl.colors.BoundaryNorm(bounds_volume, cmap_volume.N)
    fig.colorbar(mpl.cm.ScalarMappable(cmap=cmap_volume, norm=norm_volume), cax=ax, orientation = 'vertical', spacing = 'proportional')
    ax.tick_params(labelsize=8)

# Tortuosity
if TortuosityGraph3D:
    fig = plt.figure(figsize=(5,4), dpi = 300)
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(np.log10(x_sigma), np.log10(x_omegap), np.log10(x_omegaa), s=[20*item for item in y_tortuosity], c=list_color_tortuosity)
    ax.set_xlabel(r'$\sigma$', fontsize = 14)
    ax.set_ylabel(r'$\omega_p$', fontsize = 14)
    ax.set_zlabel(r'$\omega_a$', fontsize = 14)
    ax.set_xticks([-10, -8, -6, -4, -2], labels = ['1e-10', '1e-8', '1e-6', '1e-4', '1e-2'])
    ax.set_yticks([-10, -8, -6, -4, -2], labels = ['1e-10', '1e-8', '1e-6', '1e-4', '1e-2'])
    ax.set_zticks([-6, -4, -2], labels = ['1e-6', '1e-4', '1e-2'])
    ax.set_title('Tortuosity', loc = 'left', fontsize = 14)

# Area
if AreaGraph3D:
    fig = plt.figure(figsize=(5,4), dpi = 300)
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(np.log10(x_sigma), np.log10(x_omegap), np.log10(x_omegaa), s=[0.5*item for item in y_area], c=list_color_area)
    ax.set_xlabel(r'$\sigma$', fontsize = 14)
    ax.set_ylabel(r'$\omega_p$', fontsize = 14)
    ax.set_zlabel(r'$\omega_a$', fontsize = 14)
    ax.set_xticks([-10, -8, -6, -4, -2], labels = ['1e-10', '1e-8', '1e-6', '1e-4', '1e-2'])
    ax.set_yticks([-10, -8, -6, -4, -2], labels = ['1e-10', '1e-8', '1e-6', '1e-4', '1e-2'])
    ax.set_zticks([-6, -4, -2], labels = ['1e-6', '1e-4', '1e-2'])
    ax.set_title('Area', loc = 'left', fontsize = 14)

# Volume
if VolumeGraph3D:
    fig = plt.figure(figsize=(5,4), dpi = 300)
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(np.log10(x_sigma), np.log10(x_omegap), np.log10(x_omegaa), s=[0.5*item for item in y_volume], c=list_color_volume)
    ax.set_xlabel(r'$\sigma$', fontsize = 14)
    ax.set_ylabel(r'$\omega_p$', fontsize = 14)
    ax.set_zlabel(r'$\omega_a$', fontsize = 14)
    ax.set_xticks([-10, -8, -6, -4, -2], labels = ['1e-10', '1e-8', '1e-6', '1e-4', '1e-2'])
    ax.set_yticks([-10, -8, -6, -4, -2], labels = ['1e-10', '1e-8', '1e-6', '1e-4', '1e-2'])
    ax.set_zticks([-6, -4, -2], labels = ['1e-6', '1e-4', '1e-2'])
    ax.set_title('Volume', loc = 'left', fontsize = 14)


## 1D GRAPHS ##

# Initilisation of the lists 
x_sigma = []
x_omegap = []
x_omegaa = []
y_area_sigma = []
y_volume_sigma = []
y_area_omegap = []
y_volume_omegap = []
y_area_omegaa = []
y_volume_omegaa = []

if AreaVolumeGraph1DSigma:
    for k in range(2,11):
        x_sigma.append(10**(-k))

        # we read the files 
        main_file_path = file_path + str(k) + str(TestCounter1DOmegap) + str(TestCounter1DOmegaa) + "/results_from_time_0/"

        file_cellmutation = main_file_path + "results.vizmutationstates"
        file_nodescoordinates = main_file_path + "results.viznodes"

        value_area = runner.AverageArea(file_nodescoordinates, file_cellmutation)
        value_volume = runner.AverageVolume(file_nodescoordinates, file_cellmutation)

        y_area_sigma.append(value_area)
        y_volume_sigma.append(value_volume)

    plt.plot(x_sigma, y_area_sigma, 'rx-', label = 'Area')
    plt.plot(x_sigma, y_volume_sigma, 'gx-', label = 'Volume')
    plt.xscale('log')
    plt.xlabel(r'$\sigma$', fontsize = 18)
    plt.xticks(fontsize = 18)
    plt.yticks(fontsize = 18)

# Area (red) and Volume (green) for omegap
if AreaVolumeGraph1DOmegap:
    for k in range(2,11):
        x_omegap.append(10**(-k))

        # we read the files 
        main_file_path = file_path + str(TestCounter1DSigma) + str(k) + str(TestCounter1DOmegaa) + "/results_from_time_0/"

        file_cellmutation = main_file_path + "results.vizmutationstates"
        file_nodescoordinates = main_file_path + "results.viznodes"

        value_area = runner.AverageArea(file_nodescoordinates, file_cellmutation)
        value_volume = runner.AverageVolume(file_nodescoordinates, file_cellmutation)

        y_area_omegap.append(value_area)
        y_volume_omegap.append(value_volume)

    plt.plot(x_omegap, y_area_omegap, 'rx-', label = 'Area')
    plt.plot(x_omegap, y_volume_omegap, 'gx-', label = 'Volume')
    plt.xscale('log')
    plt.xlabel(r'$\omega_p$', fontsize = 18)
    plt.xticks(fontsize = 18)
    plt.yticks(fontsize = 18)

# Area (red) and Volume (green) for omegaa 
if AreaVolumeGraph1DOmegaa:
    for k in range(0,9):
        x_omegaa.append(10**(-k))

        # we read the files 
        main_file_path = file_path + str(TestCounter1DSigma) + str(TestCounter1DOmegap) + str(k) + "/results_from_time_0/"

        file_cellmutation = main_file_path + "results.vizmutationstates"
        file_nodescoordinates = main_file_path + "results.viznodes"

        value_area = runner.AverageArea(file_nodescoordinates, file_cellmutation)
        value_volume = runner.AverageVolume(file_nodescoordinates, file_cellmutation)

        y_area_omegaa.append(value_area)
        y_volume_omegaa.append(value_volume)

    plt.plot(x_omegaa, y_area_omegaa, 'rx-', label = 'Area')
    plt.plot(x_omegaa, y_volume_omegaa, 'gx-', label = 'Volume')
    plt.xscale('log')
    plt.xlabel(r'$\omega_a$', fontsize = 18)
    plt.xticks(fontsize = 18)
    plt.yticks(fontsize = 18)

if NbCellsPlane1DSigma:
    y_numbercellsplane_average = np.zeros(9)
    for l in range(1,TotalTestNb+1,1):
        y_numbercellsplane = []
        for k in range(2,11,1):
            file_nodescoordinates = file_path + str(l) + str(k) + str(TestCounter1DOmegap) + str(TestCounter1DOmegaa) + "/results_from_time_0/" + "results.viznodes"
            y_numbercellsplane.append(runner.NbCellsPlane(file_nodescoordinates, plane_normal_vector))
        
        if(l == 9):
            plt.scatter(list_sigma, y_numbercellsplane, color = 'xkcd:green', alpha = 0.7, label = 'Cell Density in Plane for one Realisation')
        else:
            plt.scatter(list_sigma, y_numbercellsplane, color = 'xkcd:green', alpha = 0.7)

        y_numbercellsplane_average = runner.AverageFunction(y_numbercellsplane_average,y_numbercellsplane,TotalTestNb)

    plt.plot(list_sigma, y_numbercellsplane_average, color = 'xkcd:darkgreen', linestyle = '--', label = 'Average Cell Density in Plane', linewidth = 5.0)
    plt.xscale('log')
    plt.ylim([0.0, 1.4])
    plt.yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0], fontsize=28)
    plt.ylabel('Cell Density in Plane', fontsize = 28)
    plt.xlabel(r'$\sigma$', fontsize = 32)

if NbCellsPlane1DOmegap:
    y_numbercellsplane_average = np.zeros(9)
    for l in range(1,TotalTestNb+1,1):
        y_numbercellsplane = []
        for k in range(2,11,1):
            file_nodescoordinates = file_path + str(l) + str(TestCounter1DSigma) + str(k) + str(TestCounter1DOmegaa) + "/results_from_time_0/" + "results.viznodes"
            y_numbercellsplane.append(runner.NbCellsPlane(file_nodescoordinates, plane_normal_vector))
        
        if(l == 9):
            plt.scatter(list_omegap, y_numbercellsplane, color = 'xkcd:green', alpha = 0.7, label = 'Cell Density in Plane for one Realisation')
        else:
            plt.scatter(list_omegap, y_numbercellsplane, color = 'xkcd:green', alpha = 0.7)

        y_numbercellsplane_average = runner.AverageFunction(y_numbercellsplane_average,y_numbercellsplane,TotalTestNb)

    plt.plot(list_omegap, y_numbercellsplane_average, color = 'xkcd:darkgreen', linestyle = '--', label = 'Average Cell Density in Plane', linewidth = 5.0)
    plt.xscale('log')
    plt.ylim([0.0, 1.4])
    plt.yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0], fontsize=28)
    plt.ylabel('Cell Density in Plane', fontsize = 28)
    plt.xlabel(r'$\omega_p$', fontsize = 32)

if NbCellsPlane1DOmegaa:
    y_numbercellsplane_average = np.zeros(9)
    for l in range(1,TotalTestNb+1,1):
        y_numbercellsplane = []
        for k in range(0,9,1):
            file_nodescoordinates = file_path + str(l) + str(TestCounter1DSigma) + str(TestCounter1DOmegap) + str(k) + "/results_from_time_0/" + "results.viznodes"
            y_numbercellsplane.append(runner.NbCellsPlane(file_nodescoordinates, plane_normal_vector))
        
        if(l == 9):
            plt.scatter(list_omegaa, y_numbercellsplane, color = 'xkcd:green', alpha = 0.7, label = 'Cell Density in Plane for one Realisation')
        else:
            plt.scatter(list_omegaa, y_numbercellsplane, color = 'xkcd:green', alpha = 0.7)

        y_numbercellsplane_average = runner.AverageFunction(y_numbercellsplane_average,y_numbercellsplane,TotalTestNb)

    plt.plot(list_omegaa, y_numbercellsplane_average, color = 'xkcd:darkgreen', linestyle = '--', label = 'Average Cell Density in Plane', linewidth = 5.0)
    plt.xscale('log')
    plt.ylim([0.0, 1.4])
    plt.yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0], fontsize=28)
    plt.ylabel('Cell Density in Plane', fontsize = 28)
    plt.xlabel(r'$\omega_a$', fontsize = 32)


if NeighbouringDistancePlane1DSigma:
    y_neighbouringdistanceplane_average = np.zeros(9)
    for l in range(1,TotalTestNb+1,1):
        y_neighbouringdistanceplane = []
        for k in range(2,11,1):
            file_nodescoordinates = file_path + str(l) + str(k) + str(TestCounter1DOmegap) + str(TestCounter1DOmegaa) + "/results_from_time_0/" + "results.viznodes"
            y_neighbouringdistanceplane.append(runner.ClosestNeighbouringDistanceInPlane(file_nodescoordinates, plane_normal_vector))
        
        if(l == 9):
            plt.scatter(list_sigma, y_neighbouringdistanceplane, color = 'xkcd:green', alpha = 0.7, label = 'Minimal Neighbouring Distance in Plane for one Realisation')
        else:
            plt.scatter(list_sigma, y_neighbouringdistanceplane, color = 'xkcd:green', alpha = 0.7)

        y_neighbouringdistanceplane_average = runner.AverageFunction(y_neighbouringdistanceplane_average,y_neighbouringdistanceplane,TotalTestNb)

    plt.plot(list_sigma, y_neighbouringdistanceplane_average, color = 'xkcd:darkgreen', linestyle = '--', label = 'Average Minimal Neighbouring Distance in Plane', linewidth = 5.0)
    plt.xscale('log')
    plt.yticks(fontsize=28)
    plt.ylabel('Minimal Neighbouring Distance in Plane', fontsize = 28)
    plt.xlabel(r'$\sigma$', fontsize = 32)

if NeighbouringDistancePlane1DOmegap:
    y_neighbouringdistanceplane_average = np.zeros(9)
    for l in range(1,TotalTestNb+1,1):
        y_neighbouringdistanceplane = []
        for k in range(2,11,1):
            file_nodescoordinates = file_path + str(l) + str(TestCounter1DSigma) + str(k) + str(TestCounter1DOmegaa) + "/results_from_time_0/" + "results.viznodes"
            y_neighbouringdistanceplane.append(runner.NbCellsPlane(file_nodescoordinates, plane_normal_vector))
        
        if(l == 9):
            plt.scatter(list_omegap, y_neighbouringdistanceplane, color = 'xkcd:green', alpha = 0.7, label = 'Minimal Neighbouring Distance in Plane for one Realisation')
        else:
            plt.scatter(list_omegap, y_neighbouringdistanceplane, color = 'xkcd:green', alpha = 0.7)

        y_neighbouringdistanceplane_average = runner.AverageFunction(y_neighbouringdistanceplane_average,y_neighbouringdistanceplane,TotalTestNb)

    plt.plot(list_omegap, y_neighbouringdistanceplane_average, color = 'xkcd:darkgreen', linestyle = '--', label = 'Average Minimal Neighbouring Distance in Plane', linewidth = 5.0)
    plt.xscale('log')
    plt.yticks(fontsize=28)
    plt.ylabel('Minimal Neighbouring Distance in Plane', fontsize = 28)
    plt.xlabel(r'$\omega_p$', fontsize = 32)

if NeighbouringDistancePlane1DOmegaa:
    y_neighbouringdistanceplane_average = np.zeros(9)
    for l in range(1,TotalTestNb+1,1):
        y_neighbouringdistanceplane = []
        for k in range(0,9,1):
            file_nodescoordinates = file_path + str(l) + str(TestCounter1DSigma) + str(TestCounter1DOmegap) + str(k) + "/results_from_time_0/" + "results.viznodes"
            y_neighbouringdistanceplane.append(runner.NbCellsPlane(file_nodescoordinates, plane_normal_vector))
        
        if(l == 9):
            plt.scatter(list_omegaa, y_neighbouringdistanceplane, color = 'xkcd:green', alpha = 0.7, label = 'Minimal Neighbouring Distance in Plane for one Realisation')
        else:
            plt.scatter(list_omegaa, y_neighbouringdistanceplane, color = 'xkcd:green', alpha = 0.7)

        y_neighbouringdistanceplane_average = runner.AverageFunction(y_neighbouringdistanceplane_average,y_neighbouringdistanceplane,TotalTestNb)

    plt.plot(list_omegaa, y_neighbouringdistanceplane_average, color = 'xkcd:darkgreen', linestyle = '--', label = 'Average Minimal Neighbouring Distance in Plane', linewidth = 5.0)
    plt.xscale('log')
    plt.yticks(fontsize=28)
    plt.ylabel('Minimal Neighbouring Distance in Plane', fontsize = 28)
    plt.xlabel(r'$\omega_a$', fontsize = 32)

if GraphNbCellsAreaPlane1DSigma:
    y_numbercellsplane_average = np.zeros(9)
    y_areaplane_average = np.zeros(9)
    for l in range(1,TotalTestNb+1,1):
        y_numbercellsplane = []
        y_areaplane = []
        for k in range(2,11,1):
            file_nodescoordinates = file_path + str(l) + str(k) + str(TestCounter1DOmegap) + str(TestCounter1DOmegaa) + "/results_from_time_0/" + "results.viznodes"
            y_numbercellsplane.append(runner.NbCellsPlane(file_nodescoordinates, plane_normal_vector))
            y_areaplane.append(runner.AreaCellsPlane(file_nodescoordinates, plane_normal_vector))
        
        # if(l == 9):
        #     plt.scatter(list_sigma, y_numbercellsplane, color = 'xkcd:green', alpha = 0.7, label = 'Cell Density in Plane for one Realisation')
        # else:
        #     plt.scatter(list_sigma, y_numbercellsplane, color = 'xkcd:green', alpha = 0.7)

        y_numbercellsplane_average = runner.AverageFunction(y_numbercellsplane_average,y_numbercellsplane,TotalTestNb)
        y_areaplane_average = runner.AverageFunction(y_areaplane_average, y_areaplane,TotalTestNb)

    fig, ax1 = plt.subplots()
    ax1.set_xlabel(r'$\sigma$', fontsize = 32)
    ax1.set_xscale('log')
    ax1.tick_params(axis='x', labelsize = 18)
    color1 = 'tab:green'
    ax1.set_ylabel('Number of Cells in Plane', fontsize = 28, color = color1)
    ax1.plot(list_sigma, y_numbercellsplane_average, color = color1, linestyle = '-', marker = '.', linewidth = 3.0, markersize = 12.0)
    ax1.tick_params(axis='y', labelcolor=color1, labelsize = 18)
    ax2 = ax1.twinx()
    color2 = 'tab:red'
    ax2.set_ylabel('Area of Plane', fontsize = 28, color = color2)
    ax2.plot(list_sigma, y_areaplane_average, color = color2, linestyle = '-', marker = '.', linewidth = 3.0, markersize = 12.0)
    ax2.tick_params(axis='y', labelcolor=color2, labelsize = 18)
    fig.tight_layout()

if GraphNbCellsAreaPlane1DOmegaa:
    y_numbercellsplane_average = np.zeros(9)
    y_areaplane_average = np.zeros(9)
    for l in range(1,TotalTestNb+1,1):
        y_numbercellsplane = []
        y_areaplane = []
        for k in range(2,11,1):
            file_nodescoordinates = file_path + str(l) + str(TestCounter1DSigma) + str(TestCounter1DOmegap) + str(1) + "/results_from_time_0/" + "results.viznodes"
            y_numbercellsplane.append(runner.NbCellsPlane(file_nodescoordinates, plane_normal_vector))
            y_areaplane.append(runner.AreaCellsPlane(file_nodescoordinates, plane_normal_vector))
        
        # if(l == 9):
        #     plt.scatter(list_sigma, y_numbercellsplane, color = 'xkcd:green', alpha = 0.7, label = 'Cell Density in Plane for one Realisation')
        # else:
        #     plt.scatter(list_sigma, y_numbercellsplane, color = 'xkcd:green', alpha = 0.7)

        y_numbercellsplane_average = runner.AverageFunction(y_numbercellsplane_average,y_numbercellsplane,TotalTestNb)
        y_areaplane_average = runner.AverageFunction(y_areaplane_average, y_areaplane,TotalTestNb)

    fig, ax1 = plt.subplots()
    ax1.set_xlabel(r'$\omega_a$', fontsize = 32)
    ax1.set_xscale('log')
    ax1.tick_params(axis='x', labelsize = 18)
    color1 = 'tab:green'
    ax1.set_ylabel('Number of Cells in Plane', fontsize = 28, color = color1)
    ax1.plot(list_omegaa, y_numbercellsplane_average, color = color1, linestyle = '-', marker = '.', linewidth = 3.0, markersize = 12.0)
    ax1.tick_params(axis='y', labelcolor=color1, labelsize = 18)
    ax2 = ax1.twinx()
    color2 = 'tab:red'
    ax2.set_ylabel('Area of Plane', fontsize = 28, color = color2)
    ax2.plot(list_omegaa, y_areaplane_average, color = color2, linestyle = '--', marker = '.', linewidth = 3.0, markersize = 8.0)
    ax2.tick_params(axis='y', labelcolor=color2, labelsize = 18)
    fig.tight_layout()

## 2D GRAPHS ##

# Initialisation of the lists
x_sigma_2dgraph = []
x_omegap_2dgraph = []
x_omegaa_2dgraph = []
y_totalnbcells_sigma_omegap = []
y_area_sigma_omegap = []
y_volume_sigma_omegap = []
y_tortuosity_sigma_omegap = []
list_color_area_sigma_omegap = []
list_color_volume_sigma_omegap = []
list_color_tortuosity_sigma_omegap = []
y_totalnbcells_sigma_omegaa = []
y_area_sigma_omegaa = []
y_volume_sigma_omegaa = []
y_tortuosity_sigma_omegaa = []
list_color_area_sigma_omegaa = []
list_color_volume_sigma_omegaa = []
list_color_tortuosity_sigma_omegaa = []
y_totalnbcells_omegap_omegaa = []
y_area_omegap_omegaa = []
y_volume_omegap_omegaa = []
y_tortuosity_omegap_omegaa = []
list_color_area_omegap_omegaa = []
list_color_volume_omegap_omegaa = []
list_color_tortuosity_omegap_omegaa = []

if Graph2DSigmaOmegap:
    for k in range(2,11):
        for j in range(2, 11):
            x_sigma_2dgraph.append(10**(-k))
            x_omegap_2dgraph.append(10**(-j))

            if MultipleRandomSeedTest:
                for l in range(1,11,1):
                    # initialisation 
                    value_totalnbcells_multiple_seeds = []
                    value_area_multiple_seeds = []
                    value_volume_multiple_seeds = []
                    value_tortuosity_multiple_seeds = []

                    # we read the files 
                    main_file_path = file_path + str(l) + str(k) + str(j) + str(TestCounter2DOmegaa) + "/results_from_time_0/"

                    file_cellmutation = main_file_path + "results.vizmutationstates"
                    file_nodescoordinates = main_file_path + "results.viznodes"
                    file_tortuosity = main_file_path + "results.viztortuosity"

                    value_totalnbcells_seed = runner.TotalNumberCells(file_cellmutation)
                    value_area_seed = runner.AverageArea(file_nodescoordinates, file_cellmutation)
                    value_volume_seed = runner.AverageVolume(file_nodescoordinates, file_cellmutation)
                    value_tortuosity_seed = runner.Tortuosity(file_tortuosity)

                    value_totalnbcells_multiple_seeds.append(value_totalnbcells_seed)
                    value_area_multiple_seeds.append(value_area_seed)
                    value_volume_multiple_seeds.append(value_volume_seed)
                    value_tortuosity_multiple_seeds.append(value_tortuosity_seed)

                value_totalnbcells = stats.mean(value_totalnbcells_multiple_seeds)
                value_area = stats.mean(value_area_multiple_seeds)
                value_volume = stats.mean(value_volume_multiple_seeds)
                value_tortuosity = stats.mean(value_tortuosity_multiple_seeds)

                y_totalnbcells_sigma_omegap.append(value_totalnbcells)
                y_area_sigma_omegap.append(value_area)
                y_volume_sigma_omegap.append(value_volume)
                y_tortuosity_sigma_omegap.append(value_tortuosity)

            else:
                # we read the files 
                main_file_path = file_path + str(k) + str(j) + str(TestCounter2DOmegaa) + "/results_from_time_0/"

                file_cellmutation = main_file_path + "results.vizmutationstates"
                file_nodescoordinates = main_file_path + "results.viznodes"
                file_tortuosity = main_file_path + "results.viztortuosity"

                value_totalnbcells = runner.TotalNumberCells(file_cellmutation)
                value_area = runner.AverageArea(file_nodescoordinates, file_cellmutation)
                value_volume = runner.AverageVolume(file_nodescoordinates, file_cellmutation)
                value_tortuosity = runner.Tortuosity(file_tortuosity)

                y_totalnbcells_sigma_omegap.append(value_totalnbcells)
                y_area_sigma_omegap.append(value_area)
                y_volume_sigma_omegap.append(value_volume)
                y_tortuosity_sigma_omegap.append(value_tortuosity)

    for value_area in y_area_sigma_omegap:
        if(value_area < 126):
            list_color_area_sigma_omegap.append('xkcd:lightgreen')
        elif(126 <= value_area < 160):
            list_color_area_sigma_omegap.append('xkcd:green')
        elif(160 <= value_area < 190):
            list_color_area_sigma_omegap.append('xkcd:lightblue')
        elif(190 <= value_area < 210):
            list_color_area_sigma_omegap.append('xkcd:blue')
        elif(210 <= value_area < 340):
            list_color_area_sigma_omegap.append('xkcd:orange')
        else:
            list_color_area_sigma_omegap.append('xkcd:red')

    for value_volume in y_volume_sigma_omegap:
        if(value_volume < 200):
            list_color_volume_sigma_omegap.append('xkcd:lightgreen')
        elif(200 <= value_volume < 210):
            list_color_volume_sigma_omegap.append('xkcd:green')
        elif(360 <= value_volume < 420):
            list_color_volume_sigma_omegap.append('xkcd:lightblue')
        elif(450 <= value_volume < 575):
            list_color_volume_sigma_omegap.append('xkcd:blue')
        elif(790 <= value_volume < 1000):
            list_color_volume_sigma_omegap.append('xkcd:orange')
        else:
            list_color_volume_sigma_omegap.append('xkcd:red')

    for value_tortuosity in y_tortuosity_sigma_omegap:
        if(value_tortuosity < 1.05):
            list_color_tortuosity_sigma_omegap.append('xkcd:lightgreen')
        elif(1.05 <= value_tortuosity < 1.2):
            list_color_tortuosity_sigma_omegap.append('xkcd:green')
        elif(1.2 <= value_tortuosity < 1.4):
            list_color_tortuosity_sigma_omegap.append('xkcd:lightblue')
        elif(1.4 <= value_tortuosity < 1.6):
            list_color_tortuosity_sigma_omegap.append('xkcd:blue')
        elif(1.6 <= value_tortuosity < 2.0):
            list_color_tortuosity_sigma_omegap.append('xkcd:orange')
        else:
            list_color_tortuosity_sigma_omegap.append('xkcd:red')

    if TotalNbCellsGraph2DSigmaOmegap:
        fig = plt.figure(figsize=(6.4,4.8), dpi = 300)
        print(y_totalnbcells_sigma_omegap)
        plt.scatter(np.log10(x_sigma_2dgraph), np.log10(x_omegap_2dgraph), s=[0.3*item for item in y_totalnbcells_sigma_omegap])
        plt.xlabel(r'$\sigma$', fontsize = 14)
        plt.ylabel(r'$\omega_p$', fontsize = 14)
        plt.xticks([-10, -8, -6, -4, -2], labels = ['1e-10', '1e-8', '1e-6', '1e-4', '1e-2'])
        plt.yticks([-10, -8, -6, -4, -2], labels = ['1e-10', '1e-8', '1e-6', '1e-4', '1e-2'])
        plt.title('Total Number of Cells', loc = 'left', fontsize = 14)

    if AreaGraph2DSigmaOmegap:
        fig = plt.figure(figsize=(6.4,4.8), dpi = 300)
        plt.scatter(np.log10(x_sigma_2dgraph), np.log10(x_omegap_2dgraph), s=[0.3*item for item in y_area_sigma_omegap], c=list_color_area_sigma_omegap)
        plt.xlabel(r'$\sigma$', fontsize = 14)
        plt.ylabel(r'$\omega_p$', fontsize = 14)
        plt.xticks([-10, -8, -6, -4, -2], labels = ['1e-10', '1e-8', '1e-6', '1e-4', '1e-2'])
        plt.yticks([-10, -8, -6, -4, -2], labels = ['1e-10', '1e-8', '1e-6', '1e-4', '1e-2'])
        plt.title('Area', loc = 'left', fontsize = 14)

    if VolumeGraph2DSigmaOmegap:
        fig = plt.figure(figsize=(6.4,4.8), dpi = 300)
        plt.scatter(np.log10(x_sigma_2dgraph), np.log10(x_omegap_2dgraph), s=[0.3*item for item in y_volume_sigma_omegap], c=list_color_volume_sigma_omegap)
        plt.xlabel(r'$\sigma$', fontsize = 14)
        plt.ylabel(r'$\omega_p$', fontsize = 14)
        plt.xticks([-10, -8, -6, -4, -2], labels = ['1e-10', '1e-8', '1e-6', '1e-4', '1e-2'])
        plt.yticks([-10, -8, -6, -4, -2], labels = ['1e-10', '1e-8', '1e-6', '1e-4', '1e-2'])
        plt.title('Volume', loc = 'left', fontsize = 14)

    if TortuosityGraph2DSigmaOmegap:
        fig = plt.figure(figsize=(6.4,4.8), dpi = 300)
        plt.scatter(np.log10(x_sigma_2dgraph), np.log10(x_omegap_2dgraph), s=[20*item for item in y_tortuosity_sigma_omegap], c=list_color_tortuosity_sigma_omegap)
        plt.xlabel(r'$\sigma$', fontsize = 14)
        plt.ylabel(r'$\omega_p$', fontsize = 14)
        plt.xticks([-10, -8, -6, -4, -2], labels = ['1e-10', '1e-8', '1e-6', '1e-4', '1e-2'])
        plt.yticks([-10, -8, -6, -4, -2], labels = ['1e-10', '1e-8', '1e-6', '1e-4', '1e-2'])
        plt.title('Tortuosity', loc = 'left', fontsize = 14)

if Graph2DSigmaOmegaa:
    for k in range(2,11):
        for j in range(0,9):
            x_sigma_2dgraph.append(10**(-k))
            x_omegaa_2dgraph.append(10**(-j))

            if MultipleRandomSeedTest:
                for l in range(1,11,1):
                    # initialisation 
                    value_totalnbcells_multiple_seeds = []
                    value_area_multiple_seeds = []
                    value_volume_multiple_seeds = []
                    value_tortuosity_multiple_seeds = []

                    # we read the files 
                    main_file_path = file_path + str(l) + str(k) + str(TestCounter2DOmegap) + str(j) + "/results_from_time_0/"

                    file_cellmutation = main_file_path + "results.vizmutationstates"
                    file_nodescoordinates = main_file_path + "results.viznodes"
                    file_tortuosity = main_file_path + "results.viztortuosity"

                    value_totalnbcells_seed = runner.TotalNumberCells(file_cellmutation)
                    value_area_seed = runner.AverageArea(file_nodescoordinates, file_cellmutation)
                    value_volume_seed = runner.AverageVolume(file_nodescoordinates, file_cellmutation)
                    value_tortuosity_seed = runner.Tortuosity(file_tortuosity)

                    value_totalnbcells_multiple_seeds.append(value_totalnbcells_seed)
                    value_area_multiple_seeds.append(value_area_seed)
                    value_volume_multiple_seeds.append(value_volume_seed)
                    value_tortuosity_multiple_seeds.append(value_tortuosity_seed)

                value_totalnbcells = stats.mean(value_totalnbcells_multiple_seeds)    
                value_area = stats.mean(value_area_multiple_seeds)
                value_volume = stats.mean(value_volume_multiple_seeds)
                value_tortuosity = stats.mean(value_tortuosity_multiple_seeds)

                y_totalnbcells_sigma_omegaa.append(value_totalnbcells)
                y_area_sigma_omegaa.append(value_area)
                y_volume_sigma_omegaa.append(value_volume)
                y_tortuosity_sigma_omegaa.append(value_tortuosity)

            else:
                # we read the files 
                main_file_path = file_path + str(k) + str(TestCounter2DOmegap) + str(j) + "/results_from_time_0/"

                file_cellmutation = main_file_path + "results.vizmutationstates"
                file_nodescoordinates = main_file_path + "results.viznodes"
                file_tortuosity = main_file_path + "results.viztortuosity"

                value_totalnbcells = runner.TotalNumberCells(file_cellmutation)
                value_area = runner.AverageArea(file_nodescoordinates, file_cellmutation)
                value_volume = runner.AverageVolume(file_nodescoordinates, file_cellmutation)
                value_tortuosity = runner.Tortuosity(file_tortuosity)

                y_totalnbcells_sigma_omegaa.append(value_totalnbcells)
                y_area_sigma_omegaa.append(value_area)
                y_volume_sigma_omegaa.append(value_volume)
                y_tortuosity_sigma_omegaa.append(value_tortuosity)

    for value_area in y_area_sigma_omegaa:
        if(value_area < 126):
            list_color_area_sigma_omegaa.append('xkcd:lightgreen')
        elif(126 <= value_area < 160):
            list_color_area_sigma_omegaa.append('xkcd:green')
        elif(160 <= value_area < 190):
            list_color_area_sigma_omegaa.append('xkcd:lightblue')
        elif(190 <= value_area < 210):
            list_color_area_sigma_omegaa.append('xkcd:blue')
        elif(210 <= value_area < 340):
            list_color_area_sigma_omegaa.append('xkcd:orange')
        else:
            list_color_area_sigma_omegaa.append('xkcd:red')

    for value_volume in y_volume_sigma_omegaa:
        if(value_volume < 200):
            list_color_volume_sigma_omegaa.append('xkcd:lightgreen')
        elif(200 <= value_volume < 210):
            list_color_volume_sigma_omegaa.append('xkcd:green')
        elif(360 <= value_volume < 420):
            list_color_volume_sigma_omegaa.append('xkcd:lightblue')
        elif(450 <= value_volume < 575):
            list_color_volume_sigma_omegaa.append('xkcd:blue')
        elif(790 <= value_volume < 1000):
            list_color_volume_sigma_omegaa.append('xkcd:orange')
        else:
            list_color_volume_sigma_omegaa.append('xkcd:red')

    for value_tortuosity in y_tortuosity_sigma_omegaa:
        if(value_tortuosity < 1.2):
            list_color_tortuosity_sigma_omegaa.append('xkcd:lightgreen')
        elif(1.2 <= value_tortuosity < 1.4):
            list_color_tortuosity_sigma_omegaa.append('xkcd:green')
        elif(1.4 <= value_tortuosity < 1.6):
            list_color_tortuosity_sigma_omegaa.append('xkcd:lightblue')
        elif(1.6 <= value_tortuosity < 1.8):
            list_color_tortuosity_sigma_omegaa.append('xkcd:blue')
        elif(1.8 <= value_tortuosity < 2.0):
            list_color_tortuosity_sigma_omegaa.append('xkcd:orange')
        else:
            list_color_tortuosity_sigma_omegaa.append('xkcd:red')

    if TotalNbCellsGraph2DSigmaOmegaa:
        fig = plt.figure(figsize=(6.4,4.8), dpi = 300)
        print(y_totalnbcells_sigma_omegaa)
        plt.scatter(np.log10(x_sigma_2dgraph), np.log10(x_omegaa_2dgraph), s=[0.3*item for item in y_totalnbcells_sigma_omegaa])
        plt.xlabel(r'$\sigma$', fontsize = 14)
        plt.ylabel(r'$\omega_a$', fontsize = 14)
        plt.xticks([-10, -8, -6, -4, -2], labels = ['1e-10', '1e-8', '1e-6', '1e-4', '1e-2'])
        plt.yticks([-6, -4, -2, 0], labels = ['1e-6', '1e-4', '1e-2', '1'])
        plt.title('Total Number of Cells', loc = 'left', fontsize = 14)

    if AreaGraph2DSigmaOmegaa:
        fig = plt.figure(figsize=(6.4,4.8), dpi = 300)
        plt.scatter(np.log10(x_sigma_2dgraph), np.log10(x_omegaa_2dgraph), s=[0.3*item for item in y_area_sigma_omegaa], c=list_color_area_sigma_omegaa)
        plt.xlabel(r'$\sigma$', fontsize = 14)
        plt.ylabel(r'$\omega_a$', fontsize = 14)
        plt.xticks([-10, -8, -6, -4, -2], labels = ['1e-10', '1e-8', '1e-6', '1e-4', '1e-2'])
        plt.yticks([-6, -4, -2, 0], labels = ['1e-6', '1e-4', '1e-2', '1'])
        plt.title('Area', loc = 'left', fontsize = 14)

    if VolumeGraph2DSigmaOmegaa:
        fig = plt.figure(figsize=(6.4,4.8), dpi = 300)
        plt.scatter(np.log10(x_sigma_2dgraph), np.log10(x_omegaa_2dgraph), s=[0.3*item for item in y_volume_sigma_omegaa], c=list_color_volume_sigma_omegaa)
        plt.xlabel(r'$\sigma$', fontsize = 14)
        plt.ylabel(r'$\omega_a$', fontsize = 14)
        plt.xticks([-10, -8, -6, -4, -2], labels = ['1e-10', '1e-8', '1e-6', '1e-4', '1e-2'])
        plt.yticks([-6, -4, -2, 0], labels = ['1e-6', '1e-4', '1e-2', '1'])
        plt.title('Volume', loc = 'left', fontsize = 14)

    if TortuosityGraph2DSigmaOmegaa:
        fig = plt.figure(figsize=(6.4,4.8), dpi = 300)
        plt.scatter(np.log10(x_sigma_2dgraph), np.log10(x_omegaa_2dgraph), s=[20*item for item in y_tortuosity_sigma_omegaa], c=list_color_tortuosity_sigma_omegaa)
        plt.xlabel(r'$\sigma$', fontsize = 14)
        plt.ylabel(r'$\omega_a$', fontsize = 14)
        plt.xticks([-10, -8, -6, -4, -2], labels = ['1e-10', '1e-8', '1e-6', '1e-4', '1e-2'])
        plt.yticks([-6, -4, -2, 0], labels = ['1e-6', '1e-4', '1e-2', '1'])
        plt.title('Tortuosity', loc = 'left', fontsize = 14)

if Graph2DOmegapOmegaa:
    for k in range(2,11):
        for j in range(0,9):
            x_omegap_2dgraph.append(10**(-k))
            x_omegaa_2dgraph.append(10**(-j))

            if MultipleRandomSeedTest:
                for l in range(1,11,1):
                    # initialisation 
                    value_totalnbcells_multiple_seeds = []
                    value_area_multiple_seeds = []
                    value_volume_multiple_seeds = []
                    value_tortuosity_multiple_seeds = []

                    # we read the files 
                    main_file_path = file_path + str(l) + str(TestCounter2DSigma) + str(k) + str(j) + "/results_from_time_0/"

                    file_cellmutation = main_file_path + "results.vizmutationstates"
                    file_nodescoordinates = main_file_path + "results.viznodes"
                    file_tortuosity = main_file_path + "results.viztortuosity"

                    value_totalnbcells_seed = runner.TotalNumberCells(file_cellmutation)
                    value_area_seed = runner.AverageArea(file_nodescoordinates, file_cellmutation)
                    value_volume_seed = runner.AverageVolume(file_nodescoordinates, file_cellmutation)
                    value_tortuosity_seed = runner.Tortuosity(file_tortuosity)

                    value_totalnbcells_multiple_seeds.append(value_totalnbcells_seed)
                    value_area_multiple_seeds.append(value_area_seed)
                    value_volume_multiple_seeds.append(value_volume_seed)
                    value_tortuosity_multiple_seeds.append(value_tortuosity_seed)

                value_totalnbcells = stats.mean(value_totalnbcells_multiple_seeds)
                value_area = stats.mean(value_area_multiple_seeds)
                value_volume = stats.mean(value_volume_multiple_seeds)
                value_tortuosity = stats.mean(value_tortuosity_multiple_seeds)

                y_totalnbcells_omegap_omegaa.append(value_totalnbcells)
                y_area_omegap_omegaa.append(value_area)
                y_volume_omegap_omegaa.append(value_volume)
                y_tortuosity_omegap_omegaa.append(value_tortuosity)

            else:
                # we read the files 
                main_file_path = file_path + str(TestCounter2DSigma) + str(k) + str(j) + "/results_from_time_0/"

                file_cellmutation = main_file_path + "results.vizmutationstates"
                file_nodescoordinates = main_file_path + "results.viznodes"
                file_tortuosity = main_file_path + "results.viztortuosity"

                value_totalnbcells = runner.TotalNumberCells(file_cellmutation)
                value_area = runner.AverageArea(file_nodescoordinates, file_cellmutation)
                value_volume = runner.AverageVolume(file_nodescoordinates, file_cellmutation)
                value_tortuosity = runner.Tortuosity(file_tortuosity)

                y_totalnbcells_omegap_omegaa.append(value_totalnbcells)
                y_area_omegap_omegaa.append(value_area)
                y_volume_omegap_omegaa.append(value_volume)
                y_tortuosity_omegap_omegaa.append(value_tortuosity)

    for value_area in y_area_omegap_omegaa:
        if(value_area < 126):
            list_color_area_omegap_omegaa.append('xkcd:lightgreen')
        elif(126 <= value_area < 160):
            list_color_area_omegap_omegaa.append('xkcd:green')
        elif(160 <= value_area < 190):
            list_color_area_omegap_omegaa.append('xkcd:lightblue')
        elif(190 <= value_area < 210):
            list_color_area_omegap_omegaa.append('xkcd:blue')
        elif(210 <= value_area < 340):
            list_color_area_omegap_omegaa.append('xkcd:orange')
        else:
            list_color_area_omegap_omegaa.append('xkcd:red')

    for value_volume in y_volume_omegap_omegaa:
        if(value_volume < 200):
            list_color_volume_omegap_omegaa.append('xkcd:lightgreen')
        elif(200 <= value_volume < 210):
            list_color_volume_omegap_omegaa.append('xkcd:green')
        elif(360 <= value_volume < 420):
            list_color_volume_omegap_omegaa.append('xkcd:lightblue')
        elif(450 <= value_volume < 575):
            list_color_volume_omegap_omegaa.append('xkcd:blue')
        elif(790 <= value_volume < 1000):
            list_color_volume_omegap_omegaa.append('xkcd:orange')
        else:
            list_color_volume_omegap_omegaa.append('xkcd:red')

    for value_tortuosity in y_tortuosity_omegap_omegaa:
        if(value_tortuosity < 1.2):
            list_color_tortuosity_omegap_omegaa.append('xkcd:lightgreen')
        elif(1.2 <= value_tortuosity < 1.4):
            list_color_tortuosity_omegap_omegaa.append('xkcd:green')
        elif(1.4 <= value_tortuosity < 1.6):
            list_color_tortuosity_omegap_omegaa.append('xkcd:lightblue')
        elif(1.6 <= value_tortuosity < 1.8):
            list_color_tortuosity_omegap_omegaa.append('xkcd:blue')
        elif(1.8 <= value_tortuosity < 2.0):
            list_color_tortuosity_omegap_omegaa.append('xkcd:orange')
        else:
            list_color_tortuosity_omegap_omegaa.append('xkcd:red')

    if TotalNbCellsGraph2DOmegapOmegaa:
        fig = plt.figure(figsize=(6.4,4.8), dpi = 300)
        print(y_totalnbcells_omegap_omegaa)
        plt.scatter(np.log10(x_omegap_2dgraph), np.log10(x_omegaa_2dgraph), s=[0.3*item for item in y_totalnbcells_omegap_omegaa])
        plt.xlabel(r'$\omega_p$', fontsize = 14)
        plt.ylabel(r'$\omega_a$', fontsize = 14)
        plt.xticks([-10, -8, -6, -4, -2], labels = ['1e-10', '1e-8', '1e-6', '1e-4', '1e-2'])
        plt.yticks([-6, -4, -2, 0], labels = ['1e-6', '1e-4', '1e-2', '1'])        
        plt.title('Total Number of Cells', loc = 'left', fontsize = 14)

    if AreaGraph2DOmegapOmegaa:
        fig = plt.figure(figsize=(6.4,4.8), dpi = 300)
        plt.scatter(np.log10(x_omegap_2dgraph), np.log10(x_omegaa_2dgraph), s=[0.3*item for item in y_area_omegap_omegaa], c=list_color_area_omegap_omegaa)
        plt.xlabel(r'$\omega_p$', fontsize = 14)
        plt.ylabel(r'$\omega_a$', fontsize = 14)
        plt.xticks([-10, -8, -6, -4, -2], labels = ['1e-10', '1e-8', '1e-6', '1e-4', '1e-2'])
        plt.yticks([-6, -4, -2, 0], labels = ['1e-6', '1e-4', '1e-2', '1'])
        plt.title('Area', loc = 'left', fontsize = 14)

    if VolumeGraph2DOmegapOmegaa:
        print(len(x_omegap_2dgraph))
        print(len(x_omegaa_2dgraph))
        fig = plt.figure(figsize=(6.4,4.8), dpi = 300)
        plt.scatter(np.log10(x_omegap_2dgraph), np.log10(x_omegaa_2dgraph), s=[0.3*item for item in y_volume_omegap_omegaa], c=list_color_volume_omegap_omegaa)
        plt.xlabel(r'$\omega_p$', fontsize = 14)
        plt.ylabel(r'$\omega_a$', fontsize = 14)
        plt.xticks([-10, -8, -6, -4, -2], labels = ['1e-10', '1e-8', '1e-6', '1e-4', '1e-2'])
        plt.yticks([-6, -4, -2, 0], labels = ['1e-6', '1e-4', '1e-2', '1'])
        plt.title('Volume', loc = 'left', fontsize = 14)

    if TortuosityGraph2DOmegapOmegaa:
        fig = plt.figure(figsize=(6.4,4.8), dpi = 300)
        plt.scatter(np.log10(x_omegap_2dgraph), np.log10(x_omegaa_2dgraph), s=[20*item for item in y_tortuosity_omegap_omegaa], c=list_color_tortuosity_omegap_omegaa)
        plt.xlabel(r'$\omega_p$', fontsize = 14)
        plt.ylabel(r'$\omega_a$', fontsize = 14)
        plt.xticks([-10, -8, -6, -4, -2], labels = ['1e-10', '1e-8', '1e-6', '1e-4', '1e-2'])
        plt.yticks([-6, -4, -2, 0], labels = ['1e-6', '1e-4', '1e-2', '1'])
        plt.title('Tortuosity', loc = 'left', fontsize = 14)

## PLOT ##
# plt.legend(loc = 'upper left', fontsize = 24, scatterpoints = 1)
# plt.legend(loc = 'upper right', fontsize = 24, ncol = 1)
# plt.legend(loc = 'lower left', fontsize = 18)
plt.show()