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

file_path = '/Users/coli171/Chaste/Output/TestParametersSensitivity/TestAreaVolumeTortuosity/Test'

list_sigma = np.array([1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10])
list_omegap = np.array([1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10])
list_omegaa = np.array([1.0, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8])
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

for i in range(0,8,2):
    for j in range(2,11,2):
        for k in range(2,11,2):
            x_omegaa.append(10**(-i))
            x_omegap.append(10**(-j))
            x_sigma.append(10**(-k))

            # we read the files 
            main_file_path = file_path + str(k) + str(j) + str(i) + "/results_from_time_0/"

            file_cellmutation = main_file_path + "results.vizmutationstates"
            file_nodescoordinates = main_file_path + "results.viznodes"
            file_branchesnumber = main_file_path + "results.vizbranchesnumber"
            file_consecutivebranches = main_file_path + "results.vizconsecutivebranches"
            file_birthtimecells = main_file_path + "results.vizbirthtimecell"
            file_tortuosity = main_file_path + "results.viztortuosity"

            # we collect the results
            list_consecutivebranches = ConsecutiveBranches(file_consecutivebranches)

            y_totalnumbercells.append(TotalNumberCells(file_cellmutation))
            y_totalnumbertipcells.append(TotalNumberTipCells(file_cellmutation))
            y_totalnumberstalkcells.append(TotalNumberStalkCells(file_cellmutation))
            y_totalnumberbranches.append(TotalNumberBranches(file_cellmutation))

            value_area = AverageArea(file_nodescoordinates, file_cellmutation)
            value_volume = AverageVolume(file_nodescoordinates, file_cellmutation)
            value_tortuosity = Tortuosity(file_tortuosity)

            y_area.append(value_area)
            y_volume.append(value_volume)
            y_tortuosity.append(value_tortuosity)

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

            if(value_volume < 70):
                list_color_volume.append('xkcd:green')
            elif(70 <= value_volume < 190):
                list_color_volume.append('xkcd:lightblue')
            elif(190 <= value_volume < 219):
                list_color_volume.append('xkcd:blue')
            elif(219 <= value_volume < 300):
                list_color_volume.append('xkcd:orange')
            else:
                list_color_volume.append('xkcd:red')

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

# List of the colors for colormap
# list_color = ['xkcd:lightgreen', 'xkcd:white', 'xkcd:green', 'xkcd:white', 'xkcd:lightblue', 'xkcd:white', 'xkcd:blue', 'xkcd:white', 'xkcd:orange', 'xkcd:white', 'xkcd:red']
# list_color = ['xkcd:green', 'xkcd:white', 'xkcd:lightblue', 'xkcd:white', 'xkcd:blue', 'xkcd:white', 'xkcd:orange', 'xkcd:white', 'xkcd:red'] # Volume color bar 

## COLORMAPS ##
# colormap for tortuosity
# fig, ax = plt.subplots(figsize=(1,6), layout ='constrained', dpi=300)
# bounds_tortuosity = [1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2]
# cmap_tortuosity = mpl.colors.ListedColormap(list_color)
# norm_tortuosity = mpl.colors.BoundaryNorm(bounds_tortuosity, cmap_tortuosity.N)
# fig.colorbar(mpl.cm.ScalarMappable(cmap=cmap_tortuosity, norm=norm_tortuosity), cax=ax, orientation = 'vertical', spacing = 'proportional')
# ax.tick_params(labelsize=10)

# colormap for area
# fig, ax = plt.subplots(figsize=(1,6), layout ='constrained', dpi=300)
# bounds_area = [100, 120, 160, 190, 210, 340, 410]
# cmap_area = mpl.colors.ListedColormap(list_color)
# norm_area = mpl.colors.BoundaryNorm(bounds_area, cmap_area.N)
# fig.colorbar(mpl.cm.ScalarMappable(cmap=cmap_area, norm=norm_area), cax=ax, orientation = 'vertical', spacing = 'proportional')
# ax.tick_params(labelsize=10)

# colormap for volume
# fig, ax = plt.subplots(figsize=(1,6), layout ='constrained', dpi=300)
# bounds_volume = [50, 60, 73, 90, 197, 202, 223, 229, 358, 360]
# cmap_volume = mpl.colors.ListedColormap(list_color)
# norm_volume = mpl.colors.BoundaryNorm(bounds_volume, cmap_volume.N)
# fig.colorbar(mpl.cm.ScalarMappable(cmap=cmap_volume, norm=norm_volume), cax=ax, orientation = 'vertical', spacing = 'proportional')
# ax.tick_params(labelsize=9)

## TORTUOSITY ##
# fig = plt.figure(figsize=(5,4), dpi = 300)
# ax = fig.add_subplot(111, projection='3d')
# ax.scatter(np.log10(x_sigma), np.log10(x_omegap), np.log10(x_omegaa), s=[20*item for item in y_tortuosity], c=list_color_tortuosity)
# ax.set_xlabel(r'$\sigma$', fontsize = 14)
# ax.set_ylabel(r'$\omega_p$', fontsize = 14)
# ax.set_zlabel(r'$\omega_a$', fontsize = 14)
# ax.set_xticks([-10, -8, -6, -4, -2], labels = ['1e-10', '1e-8', '1e-6', '1e-4', '1e-2'])
# ax.set_yticks([-10, -8, -6, -4, -2], labels = ['1e-10', '1e-8', '1e-6', '1e-4', '1e-2'])
# ax.set_zticks([-6, -4, -2], labels = ['1e-6', '1e-4', '1e-2'])
# ax.set_title('Tortuosity', loc = 'left', fontsize = 14)

## AREA ##
# fig = plt.figure(figsize=(5,4), dpi = 300)
# ax = fig.add_subplot(111, projection='3d')
# ax.scatter(np.log10(x_sigma), np.log10(x_omegap), np.log10(x_omegaa), s=[0.5*item for item in y_area], c=list_color_area)
# ax.set_xlabel(r'$\sigma$', fontsize = 14)
# ax.set_ylabel(r'$\omega_p$', fontsize = 14)
# ax.set_zlabel(r'$\omega_a$', fontsize = 14)
# ax.set_xticks([-10, -8, -6, -4, -2], labels = ['1e-10', '1e-8', '1e-6', '1e-4', '1e-2'])
# ax.set_yticks([-10, -8, -6, -4, -2], labels = ['1e-10', '1e-8', '1e-6', '1e-4', '1e-2'])
# ax.set_zticks([-6, -4, -2], labels = ['1e-6', '1e-4', '1e-2'])
# ax.set_title('Area', loc = 'left', fontsize = 14)

## VOLUME ##
# fig = plt.figure(figsize=(5,4), dpi = 300)
# ax = fig.add_subplot(111, projection='3d')
# ax.scatter(np.log10(x_sigma), np.log10(x_omegap), np.log10(x_omegaa), s=[0.5*item for item in y_volume], c=list_color_volume)
# ax.set_xlabel(r'$\sigma$', fontsize = 14)
# ax.set_ylabel(r'$\omega_p$', fontsize = 14)
# ax.set_zlabel(r'$\omega_a$', fontsize = 14)
# ax.set_xticks([-10, -8, -6, -4, -2], labels = ['1e-10', '1e-8', '1e-6', '1e-4', '1e-2'])
# ax.set_yticks([-10, -8, -6, -4, -2], labels = ['1e-10', '1e-8', '1e-6', '1e-4', '1e-2'])
# ax.set_zticks([-6, -4, -2], labels = ['1e-6', '1e-4', '1e-2'])
# ax.set_title('Volume', loc = 'left', fontsize = 14)

# x_sigma = []
# x_omegap = []
# x_omegaa = []
# y_area_sigma = []
# y_volume_sigma = []
# y_area_omegap = []
# y_volume_omegap = []
# y_area_omegaa = []
# y_volume_omegaa = []

## VARIATION OF ONLY ONE PARAMETER ##
# # sigma 
# for k in range(2,11):
#     x_sigma.append(10**(-k))

#     # we read the files 
#     main_file_path = file_path + str(k) + str(4) + str(2) + "/results_from_time_0/"

#     file_cellmutation = main_file_path + "results.vizmutationstates"
#     file_nodescoordinates = main_file_path + "results.viznodes"

#     value_area = AverageArea(file_nodescoordinates, file_cellmutation)
#     value_volume = AverageVolume(file_nodescoordinates, file_cellmutation)

#     y_area_sigma.append(value_area)
#     y_volume_sigma.append(value_volume)

# # omega_p
# for k in range(2,11):
#     x_omegap.append(10**(-k))

#     # we read the files 
#     main_file_path = file_path + str(8) + str(k) + str(4) + "/results_from_time_0/"

#     file_cellmutation = main_file_path + "results.vizmutationstates"
#     file_nodescoordinates = main_file_path + "results.viznodes"

#     value_area = AverageArea(file_nodescoordinates, file_cellmutation)
#     value_volume = AverageVolume(file_nodescoordinates, file_cellmutation)

#     y_area_omegap.append(value_area)
#     y_volume_omegap.append(value_volume)

# # omega_a
# for k in range(0,9):
#     x_omegaa.append(10**(-k))

#     # we read the files 
#     main_file_path = file_path + str(2) + str(2) + str(k) + "/results_from_time_0/"

#     file_cellmutation = main_file_path + "results.vizmutationstates"
#     file_nodescoordinates = main_file_path + "results.viznodes"

#     value_area = AverageArea(file_nodescoordinates, file_cellmutation)
#     value_volume = AverageVolume(file_nodescoordinates, file_cellmutation)

#     y_area_omegaa.append(value_area)
#     y_volume_omegaa.append(value_volume)

# Area (red) and Volume (green) for sigma 
# plt.plot(x_sigma, y_area_sigma, 'rx-', label = 'Area')
# plt.plot(x_sigma, y_volume_sigma, 'gx-', label = 'Volume')
# plt.xscale('log')
# plt.xlabel(r'$\sigma$', fontsize = 18)
# plt.xticks(fontsize = 18)
# plt.yticks(fontsize = 18)

# Area (red) and Volume (green) for omegap
# plt.plot(x_omegap, y_area_omegap, 'rx-', label = 'Area')
# plt.plot(x_omegap, y_volume_omegap, 'gx-', label = 'Volume')
# plt.xscale('log')
# plt.xlabel(r'$\omega_p$', fontsize = 18)
# plt.xticks(fontsize = 18)
# plt.yticks(fontsize = 18)

# Area (red) and Volume (green) for omegaa 
# plt.plot(x_omegaa, y_area_omegaa, 'rx-', label = 'Area')
# plt.plot(x_omegaa, y_volume_omegaa, 'gx-', label = 'Volume')
# plt.xscale('log')
# plt.xlabel(r'$\omega_a$', fontsize = 18)
# plt.xticks(fontsize = 18)
# plt.yticks(fontsize = 18)


## VARIATION OF TWO PARAMETERS ##
# x_sigma = []
# x_omegap = []
# x_omegaa = []
# y_area_sigma_omegap = []
# y_volume_sigma_omegap = []
# y_tortuosity_sigma_omegap = []
# list_color_area_sigma_omegap = []
# list_color_volume_sigma_omegap = []
# list_color_tortuosity_sigma_omegap = []
# y_area_sigma_omegaa = []
# y_volume_sigma_omegaa = []
# y_tortuosity_sigma_omegaa = []
# list_color_area_sigma_omegaa = []
# list_color_volume_sigma_omegaa = []
# list_color_tortuosity_sigma_omegaa = []
# y_area_omegap_omegaa = []
# y_volume_omegap_omegaa = []
# y_tortuosity_omegap_omegaa = []
# list_color_area_omegap_omegaa = []
# list_color_volume_omegap_omegaa = []
# list_color_tortuosity_omegap_omegaa = []


# # Area for (sigma, omegap) / omegaa = 1e-2
# # (sigma, omegap) 
# for k in range(2,11):
#     for j in range(2, 11):
#         x_sigma.append(10**(-k))
#         x_omegap.append(10**(-j))

#         # we read the files 
#         main_file_path = file_path + str(k) + str(j) + str(2) + "/results_from_time_0/"

#         file_cellmutation = main_file_path + "results.vizmutationstates"
#         file_nodescoordinates = main_file_path + "results.viznodes"
#         file_tortuosity = main_file_path + "results.viztortuosity"

#         value_area = AverageArea(file_nodescoordinates, file_cellmutation)
#         value_volume = AverageVolume(file_nodescoordinates, file_cellmutation)
#         value_tortuosity = Tortuosity(file_tortuosity)

#         y_area_sigma_omegap.append(value_area)
#         y_volume_sigma_omegap.append(value_volume)
#         y_tortuosity_sigma_omegap.append(value_tortuosity)

#         if(value_area < 126):
#             list_color_area_sigma_omegap.append('xkcd:lightgreen')
#         elif(126 <= value_area < 160):
#             list_color_area_sigma_omegap.append('xkcd:green')
#         elif(160 <= value_area < 190):
#             list_color_area_sigma_omegap.append('xkcd:lightblue')
#         elif(190 <= value_area < 210):
#             list_color_area_sigma_omegap.append('xkcd:blue')
#         elif(210 <= value_area < 340):
#             list_color_area_sigma_omegap.append('xkcd:orange')
#         else:
#             list_color_area_sigma_omegap.append('xkcd:red')

#         if(value_volume < 70):
#             list_color_volume_sigma_omegap.append('xkcd:lightgreen')
#         elif(70 <= value_volume < 80):
#             list_color_volume_sigma_omegap.append('xkcd:green')
#         elif(80 <= value_volume < 190):
#             list_color_volume_sigma_omegap.append('xkcd:lightblue')
#         elif(190 <= value_volume < 219):
#             list_color_volume_sigma_omegap.append('xkcd:blue')
#         elif(219 <= value_volume < 300):
#             list_color_volume_sigma_omegap.append('xkcd:orange')
#         else:
#             list_color_volume_sigma_omegap.append('xkcd:red')

#         if(value_tortuosity < 1.05):
#             list_color_tortuosity_sigma_omegap.append('xkcd:lightgreen')
#         elif(1.05 <= value_tortuosity < 1.2):
#             list_color_tortuosity_sigma_omegap.append('xkcd:green')
#         elif(1.2 <= value_tortuosity < 1.4):
#             list_color_tortuosity_sigma_omegap.append('xkcd:lightblue')
#         elif(1.4 <= value_tortuosity < 1.6):
#             list_color_tortuosity_sigma_omegap.append('xkcd:blue')
#         elif(1.6 <= value_tortuosity < 2.0):
#             list_color_tortuosity_sigma_omegap.append('xkcd:orange')
#         else:
#             list_color_tortuosity_sigma_omegap.append('xkcd:red')

# fig = plt.figure(figsize=(6.4,4.8), dpi = 300)
# plt.scatter(np.log10(x_sigma), np.log10(x_omegap), s=[0.3*item for item in y_area_sigma_omegap], c=list_color_area_sigma_omegap)
# plt.xlabel(r'$\sigma$', fontsize = 14)
# plt.ylabel(r'$\omega_p$', fontsize = 14)
# plt.xticks([-10, -8, -6, -4, -2], labels = ['1e-10', '1e-8', '1e-6', '1e-4', '1e-2'])
# plt.yticks([-10, -8, -6, -4, -2], labels = ['1e-10', '1e-8', '1e-6', '1e-4', '1e-2'])
# plt.title('Area', loc = 'left', fontsize = 14)

# fig = plt.figure(figsize=(6.4,4.8), dpi = 300)
# plt.scatter(np.log10(x_sigma), np.log10(x_omegap), s=[0.3*item for item in y_volume_sigma_omegap], c=list_color_volume_sigma_omegap)
# plt.xlabel(r'$\sigma$', fontsize = 14)
# plt.ylabel(r'$\omega_p$', fontsize = 14)
# plt.xticks([-10, -8, -6, -4, -2], labels = ['1e-10', '1e-8', '1e-6', '1e-4', '1e-2'])
# plt.yticks([-10, -8, -6, -4, -2], labels = ['1e-10', '1e-8', '1e-6', '1e-4', '1e-2'])
# plt.title('Volume', loc = 'left', fontsize = 14)

# fig = plt.figure(figsize=(6.4,4.8), dpi = 300)
# plt.scatter(np.log10(x_sigma), np.log10(x_omegap), s=[20*item for item in y_tortuosity_sigma_omegap], c=list_color_tortuosity_sigma_omegap)
# plt.xlabel(r'$\sigma$', fontsize = 14)
# plt.ylabel(r'$\omega_p$', fontsize = 14)
# plt.xticks([-10, -8, -6, -4, -2], labels = ['1e-10', '1e-8', '1e-6', '1e-4', '1e-2'])
# plt.yticks([-10, -8, -6, -4, -2], labels = ['1e-10', '1e-8', '1e-6', '1e-4', '1e-2'])
# plt.title('Tortuosity', loc = 'left', fontsize = 14)

# Area for (sigma, omegaa) / omegap = 1e-2
# (sigma, omega_a)
# for k in range(2,11):
#     for j in range(0,9):
#         x_sigma.append(10**(-k))
#         x_omegaa.append(10**(-j))

#         # we read the files 
#         main_file_path = file_path + str(k) + str(2) + str(j) + "/results_from_time_0/"

#         file_cellmutation = main_file_path + "results.vizmutationstates"
#         file_nodescoordinates = main_file_path + "results.viznodes"
#         file_tortuosity = main_file_path + "results.viztortuosity"

#         value_area = AverageArea(file_nodescoordinates, file_cellmutation)
#         value_volume = AverageVolume(file_nodescoordinates, file_cellmutation)
#         value_tortuosity = Tortuosity(file_tortuosity)

#         y_area_sigma_omegaa.append(value_area)
#         y_volume_sigma_omegaa.append(value_volume)
#         y_tortuosity_sigma_omegaa.append(value_tortuosity)

#         if(value_area < 126):
#             list_color_area_sigma_omegaa.append('xkcd:lightgreen')
#         elif(126 <= value_area < 160):
#             list_color_area_sigma_omegaa.append('xkcd:green')
#         elif(160 <= value_area < 190):
#             list_color_area_sigma_omegaa.append('xkcd:lightblue')
#         elif(190 <= value_area < 210):
#             list_color_area_sigma_omegaa.append('xkcd:blue')
#         elif(210 <= value_area < 340):
#             list_color_area_sigma_omegaa.append('xkcd:orange')
#         else:
#             list_color_area_sigma_omegaa.append('xkcd:red')

#         if(value_volume < 70):
#             list_color_volume_sigma_omegaa.append('xkcd:green')
#         elif(70 <= value_volume < 190):
#             list_color_volume_sigma_omegaa.append('xkcd:lightblue')
#         elif(190 <= value_volume < 219):
#             list_color_volume_sigma_omegaa.append('xkcd:blue')
#         elif(219 <= value_volume < 300):
#             list_color_volume_sigma_omegaa.append('xkcd:orange')
#         else:
#             list_color_volume_sigma_omegaa.append('xkcd:red')

#         if(value_tortuosity < 1.2):
#             list_color_tortuosity_sigma_omegaa.append('xkcd:lightgreen')
#         elif(1.2 <= value_tortuosity < 1.4):
#             list_color_tortuosity_sigma_omegaa.append('xkcd:green')
#         elif(1.4 <= value_tortuosity < 1.6):
#             list_color_tortuosity_sigma_omegaa.append('xkcd:lightblue')
#         elif(1.6 <= value_tortuosity < 1.8):
#             list_color_tortuosity_sigma_omegaa.append('xkcd:blue')
#         elif(1.8 <= value_tortuosity < 2.0):
#             list_color_tortuosity_sigma_omegaa.append('xkcd:orange')
#         else:
#             list_color_tortuosity_sigma_omegaa.append('xkcd:red')

# fig = plt.figure(figsize=(6.4,4.8), dpi = 300)
# plt.scatter(np.log10(x_sigma), np.log10(x_omegaa), s=[0.3*item for item in y_area_sigma_omegaa], c=list_color_area_sigma_omegaa)
# plt.xlabel(r'$\sigma$', fontsize = 14)
# plt.ylabel(r'$\omega_a$', fontsize = 14)
# plt.xticks([-10, -8, -6, -4, -2], labels = ['1e-10', '1e-8', '1e-6', '1e-4', '1e-2'])
# plt.yticks([-6, -4, -2, 0], labels = ['1e-6', '1e-4', '1e-2', '1'])
# plt.title('Area', loc = 'left', fontsize = 14)

# fig = plt.figure(figsize=(6.4,4.8), dpi = 300)
# plt.scatter(np.log10(x_sigma), np.log10(x_omegaa), s=[0.3*item for item in y_volume_sigma_omegaa], c=list_color_volume_sigma_omegaa)
# plt.xlabel(r'$\sigma$', fontsize = 14)
# plt.ylabel(r'$\omega_a$', fontsize = 14)
# plt.xticks([-10, -8, -6, -4, -2], labels = ['1e-10', '1e-8', '1e-6', '1e-4', '1e-2'])
# plt.yticks([-6, -4, -2, 0], labels = ['1e-6', '1e-4', '1e-2', '1'])
# plt.title('Volume', loc = 'left', fontsize = 14)

# fig = plt.figure(figsize=(6.4,4.8), dpi = 300)
# plt.scatter(np.log10(x_sigma), np.log10(x_omegaa), s=[20*item for item in y_tortuosity_sigma_omegaa], c=list_color_tortuosity_sigma_omegaa)
# plt.xlabel(r'$\sigma$', fontsize = 14)
# plt.ylabel(r'$\omega_a$', fontsize = 14)
# plt.xticks([-10, -8, -6, -4, -2], labels = ['1e-10', '1e-8', '1e-6', '1e-4', '1e-2'])
# plt.yticks([-6, -4, -2, 0], labels = ['1e-6', '1e-4', '1e-2', '1'])
# plt.title('Tortuosity', loc = 'left', fontsize = 14)

# # Area for (omegap, omegaa) / sigma = 1e-2
# # (omega_p, omega_a)
# for k in range(2,11):
#     for j in range(0,9):
#         x_omegap.append(10**(-k))
#         x_omegaa.append(10**(-j))

#         # we read the files 
#         main_file_path = file_path + str(2) + str(k) + str(j) + "/results_from_time_0/"

#         file_cellmutation = main_file_path + "results.vizmutationstates"
#         file_nodescoordinates = main_file_path + "results.viznodes"
#         file_tortuosity = main_file_path + "results.viztortuosity"

#         value_area = AverageArea(file_nodescoordinates, file_cellmutation)
#         value_volume = AverageVolume(file_nodescoordinates, file_cellmutation)
#         value_tortuosity = Tortuosity(file_tortuosity)

#         y_area_omegap_omegaa.append(value_area)
#         y_volume_omegap_omegaa.append(value_volume)
#         y_tortuosity_omegap_omegaa.append(value_tortuosity)

#         if(value_area < 126):
#             list_color_area_omegap_omegaa.append('xkcd:lightgreen')
#         elif(126 <= value_area < 160):
#             list_color_area_omegap_omegaa.append('xkcd:green')
#         elif(160 <= value_area < 190):
#             list_color_area_omegap_omegaa.append('xkcd:lightblue')
#         elif(190 <= value_area < 210):
#             list_color_area_omegap_omegaa.append('xkcd:blue')
#         elif(210 <= value_area < 340):
#             list_color_area_omegap_omegaa.append('xkcd:orange')
#         else:
#             list_color_area_omegap_omegaa.append('xkcd:red')

#         if(value_volume < 70):
#             list_color_volume_omegap_omegaa.append('xkcd:green')
#         elif(70 <= value_volume < 190):
#             list_color_volume_omegap_omegaa.append('xkcd:lightblue')
#         elif(190 <= value_volume < 219):
#             list_color_volume_omegap_omegaa.append('xkcd:blue')
#         elif(219 <= value_volume < 300):
#             list_color_volume_omegap_omegaa.append('xkcd:orange')
#         else:
#             list_color_volume_omegap_omegaa.append('xkcd:red')

#         if(value_tortuosity < 1.2):
#             list_color_tortuosity_omegap_omegaa.append('xkcd:lightgreen')
#         elif(1.2 <= value_tortuosity < 1.4):
#             list_color_tortuosity_omegap_omegaa.append('xkcd:green')
#         elif(1.4 <= value_tortuosity < 1.6):
#             list_color_tortuosity_omegap_omegaa.append('xkcd:lightblue')
#         elif(1.6 <= value_tortuosity < 1.8):
#             list_color_tortuosity_omegap_omegaa.append('xkcd:blue')
#         elif(1.8 <= value_tortuosity < 2.0):
#             list_color_tortuosity_omegap_omegaa.append('xkcd:orange')
#         else:
#             list_color_tortuosity_omegap_omegaa.append('xkcd:red')

# fig = plt.figure(figsize=(6.4,4.8), dpi = 300)
# plt.scatter(np.log10(x_omegap), np.log10(x_omegaa), s=[0.3*item for item in y_area_omegap_omegaa], c=list_color_area_omegap_omegaa)
# plt.xlabel(r'$\omega_p$', fontsize = 14)
# plt.ylabel(r'$\omega_a$', fontsize = 14)
# plt.xticks([-10, -8, -6, -4, -2], labels = ['1e-10', '1e-8', '1e-6', '1e-4', '1e-2'])
# plt.yticks([-6, -4, -2, 0], labels = ['1e-6', '1e-4', '1e-2', '1'])
# plt.title('Area', loc = 'left', fontsize = 14)

# fig = plt.figure(figsize=(6.4,4.8), dpi = 300)
# plt.scatter(np.log10(x_omegap), np.log10(x_omegaa), s=[0.3*item for item in y_volume_omegap_omegaa], c=list_color_volume_omegap_omegaa)
# plt.xlabel(r'$\omega_p$', fontsize = 14)
# plt.ylabel(r'$\omega_a$', fontsize = 14)
# plt.xticks([-10, -8, -6, -4, -2], labels = ['1e-10', '1e-8', '1e-6', '1e-4', '1e-2'])
# plt.yticks([-6, -4, -2, 0], labels = ['1e-6', '1e-4', '1e-2', '1'])
# plt.title('Volume', loc = 'left', fontsize = 14)

# fig = plt.figure(figsize=(6.4,4.8), dpi = 300)
# plt.scatter(np.log10(x_omegap), np.log10(x_omegaa), s=[20*item for item in y_tortuosity_omegap_omegaa], c=list_color_tortuosity_omegap_omegaa)
# plt.xlabel(r'$\omega_p$', fontsize = 14)
# plt.ylabel(r'$\omega_a$', fontsize = 14)
# plt.xticks([-10, -8, -6, -4, -2], labels = ['1e-10', '1e-8', '1e-6', '1e-4', '1e-2'])
# plt.yticks([-6, -4, -2, 0], labels = ['1e-6', '1e-4', '1e-2', '1'])
# plt.title('Tortuosity', loc = 'left', fontsize = 14)

## PLOT ##
# plt.legend(loc = 'upper left', fontsize = 24, scatterpoints = 1)
plt.legend(loc = 'upper right', fontsize = 24, ncol = 1)
# plt.legend(loc = 'lower left', fontsize = 18)
plt.show()