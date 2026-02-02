import numpy as np
import matplotlib.pyplot as plt
import statistics as stats
import polyscope as ps

from DataExtractionFetoFlow import DataExtractionFetoFlow
runner = DataExtractionFetoFlow

## Settings ##
hpc = False
local = True

InitialisationFiles = True
Visualisation = False

QuantitativeAnalysis = True
AverageFlow = False
FlowDistributionOneRealisation = False
FlowInsideLesion = True
FlowInsideLesionTwoGraphs = False

dim = 2
ref_point = 40

## Parameters ##
# at the inlet, we have an arterioles (pressure between 30 and 60 mmHg) and at the outlet, we have capillaries (average of 20 mmHg) -> Al-Nuaimi et al.
# pressure in Pa, 50 mmHg = 6666.1 Pa and 20 mmHg = 2666.45 Pa
# pressure difference between inlet and outlet of 30 mmHg -> Bazmara et al. but pressure difference of 9mmHg -> Chaplain et al.
inlet_pressure, outlet_pressure = 6666.1, 1333.22*4 
umbilical_artery_radius, decay_factor = 5e-6, 1.38 # capillaries have a radius of approx. 5 micrometers, 5e-6
viscosity_type = 'constant'
SmallSystem = False


## Quantitative Analysis ## 
if QuantitativeAnalysis:
    if AverageFlow:
        fig, ax = plt.subplots(figsize = (12,8), dpi = 300, layout='constrained')

        sourceterm = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
        flow_analyticalapproxpde_average = np.zeros(10)
        flow_constant_average = np.zeros(10)

        # Constant 
        for SourceNb in range(1, 7):
            flow_analyticalapproxpde_list = []
            flow_constant_list = []
            for SeedNb in [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]: 
                if hpc:
                    main_pathway_constant = "/hpc/coli171/Results/PaperAngiogenesisModel2025/PaperModel2025Analysis2D/" + "CoupledModel2DConstant/CoupledModel2DConstantSeed" 
                    main_pathway_analyticalapproxpde = "/hpc/coli171/Results/PaperAngiogenesisModel2025/PaperModel2025Analysis2D/" + "CoupledModel2DAnalyticalApproxPde/CoupledModel2DAnalyticalApproxPdeSeed"
                    main_pathway_constant_seedsource = main_pathway_constant + str(SeedNb) + "Source" + str(SourceNb) + "/results_from_time_0/"
                    main_pathway_analyticalapproxpde_seedsource = main_pathway_analyticalapproxpde + str(SeedNb) + "Source" + str(SourceNb) + "/results_from_time_0/"
                elif local:
                    main_pathway_constant = "/Users/coli171/Desktop/test/Constant/Seed" 
                    main_pathway_analyticalapproxpde = "/Users/coli171/Desktop/test/AnalyticalApproxPde/Seed" 
                    main_pathway_constant_seedsource = main_pathway_constant + str(SeedNb) + "Source" + str(SourceNb) + "/"
                    main_pathway_analyticalapproxpde_seedsource = main_pathway_analyticalapproxpde + str(SeedNb) + "Source" + str(SourceNb) + "/"

                # nodes and elements calculations 
                nodes_analyticalapproxpde, edges_analyticalapproxpde = runner.nodes_elements_calculation(main_pathway_analyticalapproxpde, InitialisationFiles, hpc, local, SeedNb, SourceNb, dim)
                nodes_constant, edges_constant = runner.nodes_elements_calculation(main_pathway_constant, InitialisationFiles, hpc, local, SeedNb, SourceNb, dim)

                # calculate flow and pressure for one random seed and one source term
                nodes_analyticalapproxpde_ps, edges_analyticalapproxpde_ps, pressure_analyticalapproxpde_ps, flow_analyticalapproxpde_ps = runner.flow_pressure_calculation(main_pathway_analyticalapproxpde, nodes_analyticalapproxpde, edges_analyticalapproxpde, inlet_pressure, outlet_pressure, umbilical_artery_radius, decay_factor, viscosity_type, SmallSystem, hpc, local, SeedNb, SourceNb, dim)
                nodes_constant_ps, edges_constant_ps, pressure_constant_ps, flow_constant_ps = runner.flow_pressure_calculation(main_pathway_constant, nodes_constant, edges_constant, inlet_pressure, outlet_pressure, umbilical_artery_radius, decay_factor, viscosity_type, SmallSystem, hpc, local, SeedNb, SourceNb, dim)

                flow_analyticalapproxpde = stats.mean(flow_analyticalapproxpde_ps)
                flow_constant = stats.mean(flow_constant_ps)

                # get the average flow for this realisation 
                flow_analyticalapproxpde_list.append(flow_analyticalapproxpde)
                flow_constant_list.append(flow_constant)

            # scatter plots 
            ax.scatter([sourceterm[SourceNb-1] for i in range(10)], flow_analyticalapproxpde_list, alpha=0.75, color='#757575', marker='.', s = 180)
            ax.scatter([sourceterm[SourceNb-1] for i in range(10)], flow_constant_list, alpha=0.75, color='#c90000', marker='D', s = 35)
            
            # add it to the average for this source term
            flow_analyticalapproxpde_average[SourceNb-1] = stats.mean(flow_analyticalapproxpde_list)
            flow_constant_average[SourceNb-1] = stats.mean(flow_constant_list)

        # Steady-state
        ax.plot(sourceterm, flow_constant_average, marker = 'D', markersize = 6, label = 'ECM Hypothesis', color='#c90000')
        ax.plot(sourceterm, flow_analyticalapproxpde_average, label = 'Lesion Hypothesis', marker = '.', markersize = 15.0, color='#757575')

        ax.legend(loc='upper right', fontsize = 12)
        ax.set_xlabel(r'$c_{max}$', fontsize = 15)
        ax.set_ylabel('Average Flow', fontsize = 13)
        ax.set_xticks(sourceterm)
        ax.tick_params(axis = 'both', labelsize = 12)
        plt.show()

    if FlowDistributionOneRealisation:
        SeedNb = 10
        SourceNb = 10
        hpc = False
        local = True
        dim = 2

        # constant 
        main_pathway_constant = "/Users/coli171/Desktop/test/Constant/Seed" 
        nodes_constant, edges_constant = runner.nodes_elements_calculation(main_pathway_constant, InitialisationFiles, hpc, local, SeedNb, SourceNb, dim)
        indices_vesseltips_constant = runner.find_vesseltips(nodes_constant, edges_constant)
        nodes_constant_ps, edges_constant_ps, pressure_constant, flow_constant = runner.flow_pressure_calculation(main_pathway_constant, nodes_constant, edges_constant, inlet_pressure, outlet_pressure, umbilical_artery_radius, decay_factor, viscosity_type, SmallSystem, hpc, local, SeedNb, SourceNb, dim)
        flow_constant_vesseltips = [flow_constant[i] for i in indices_vesseltips_constant]
        bins_constant = np.histogram_bin_edges(flow_constant_vesseltips, bins=15)

        # steady-state
        main_pathway_analyticalapproxpde = "/Users/coli171/Desktop/test/AnalyticalApproxPde/Seed" 
        nodes_analyticalapproxpde, edges_analyticalapproxpde = runner.nodes_elements_calculation(main_pathway_analyticalapproxpde, InitialisationFiles, hpc, local, SeedNb, SourceNb, dim)
        indices_vesseltips_analyticalapproxpde = runner.find_vesseltips(nodes_analyticalapproxpde, edges_analyticalapproxpde)
        nodes_analyticalapproxpde_ps, edges_analyticalapproxpde_ps, pressure_analyticalapproxpde, flow_analyticalapproxpde = runner.flow_pressure_calculation(main_pathway_analyticalapproxpde, nodes_analyticalapproxpde, edges_analyticalapproxpde, inlet_pressure, outlet_pressure, umbilical_artery_radius, decay_factor, viscosity_type, SmallSystem, hpc, local, SeedNb, SourceNb, dim)
        flow_analyticalapproxpde_vesseltips = [flow_analyticalapproxpde[i] for i in indices_vesseltips_analyticalapproxpde]
        bins_analyticalapproxpde = np.histogram_bin_edges(flow_analyticalapproxpde_vesseltips, bins=15)

        fig, ax = plt.subplots(figsize=(7, 5), dpi=150)

        plt.hist(flow_constant_vesseltips, bins=bins_constant, histtype='step', linewidth=2, label='ECM Hypothesis', color='#008080')
        plt.hist(flow_analyticalapproxpde_vesseltips,  bins=bins_analyticalapproxpde, histtype='step', linewidth=2, label='Lesion Hypothesis', color='#C00000')

        ax.legend(loc='upper right', fontsize = 12)
        ax.set_ylabel('Distribution of Flow at Vessel Tips', fontsize = 13)
        ax.set_xlabel('Blood Flow')
        ax.tick_params(axis = 'both', labelsize = 12)
        plt.show()

    if FlowInsideLesion:
        fig, ax = plt.subplots(dpi = 300, layout='constrained')

        sourceterm = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
        flow_analyticalapproxpde_average = np.zeros(10)
        flow_constant_average = np.zeros(10)
        flowinlet_analyticalapproxpde_average = np.zeros(10)
        flowinlet_constant_average = np.zeros(10)

        # Constant 
        for SourceNb in range(1, 11):
            flow_analyticalapproxpde_list = []
            flow_constant_list = []
            flowinlet_analyticalapproxpde_list = []
            flowinlet_constant_list = []
            for SeedNb in [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]: 
                if hpc:
                    main_pathway_constant = "/hpc/coli171/Results/PaperAngiogenesisModel2025/PaperModel2025Analysis2D/" + "CoupledModel2DConstant/CoupledModel2DConstantSeed" 
                    main_pathway_analyticalapproxpde = "/hpc/coli171/Results/PaperAngiogenesisModel2025/PaperModel2025Analysis2D/" + "CoupledModel2DAnalyticalApproxPde/CoupledModel2DAnalyticalApproxPdeSeed"
                    main_pathway_constant_seedsource = main_pathway_constant + str(SeedNb) + "Source" + str(SourceNb) + "/results_from_time_0/"
                    main_pathway_analyticalapproxpde_seedsource = main_pathway_analyticalapproxpde + str(SeedNb) + "Source" + str(SourceNb) + "/results_from_time_0/"
                elif local:
                    main_pathway = "/Users/coli171/Desktop/test/"
                    main_pathway_constant = main_pathway + "Constant/Seed" 
                    main_pathway_analyticalapproxpde = main_pathway + "AnalyticalApproxPde/Seed" 
                    main_pathway_constant_seedsource = main_pathway_constant + str(SeedNb) + "Source" + str(SourceNb) + "/"
                    main_pathway_analyticalapproxpde_seedsource = main_pathway_analyticalapproxpde + str(SeedNb) + "Source" + str(SourceNb) + "/"

                # nodes and elements calculations 
                nodes_analyticalapproxpde, edges_analyticalapproxpde = runner.nodes_elements_calculation(main_pathway_analyticalapproxpde, InitialisationFiles, hpc, local, SeedNb, SourceNb, dim)
                nodes_constant, edges_constant = runner.nodes_elements_calculation(main_pathway_constant, InitialisationFiles, hpc, local, SeedNb, SourceNb, dim)

                indices_vesseltips_constant = runner.find_vesseltips_insidelesion(nodes_constant, edges_constant, ref_point)
                print(len(indices_vesseltips_constant))
                indices_vesseltips_analyticalapproxpde = runner.find_vesseltips_insidelesion(nodes_analyticalapproxpde, edges_analyticalapproxpde, ref_point)
                print(len(indices_vesseltips_analyticalapproxpde))

                # calculate flow and pressure for one random seed and one source term
                nodes_analyticalapproxpde_ps, edges_analyticalapproxpde_ps, pressure_analyticalapproxpde_ps, flow_analyticalapproxpde_ps = runner.flow_pressure_calculation(main_pathway_analyticalapproxpde, nodes_analyticalapproxpde, edges_analyticalapproxpde, inlet_pressure, outlet_pressure, umbilical_artery_radius, decay_factor, viscosity_type, SmallSystem, hpc, local, SeedNb, SourceNb, dim)
                nodes_constant_ps, edges_constant_ps, pressure_constant_ps, flow_constant_ps = runner.flow_pressure_calculation(main_pathway_constant, nodes_constant, edges_constant, inlet_pressure, outlet_pressure, umbilical_artery_radius, decay_factor, viscosity_type, SmallSystem, hpc, local, SeedNb, SourceNb, dim)

                flow_constant_vesseltips = [flow_constant_ps[i-1] for i in indices_vesseltips_constant]
                flow_analyticalapproxpde_vesseltips = [flow_analyticalapproxpde_ps[i-1] for i in indices_vesseltips_analyticalapproxpde]

                flow_analyticalapproxpde = 0
                if len(flow_analyticalapproxpde_vesseltips)!=0:
                    flow_analyticalapproxpde = stats.mean(flow_analyticalapproxpde_vesseltips)
                
                flow_constant = 0
                if len(flow_constant_vesseltips)!=0:
                    flow_constant = stats.mean(flow_constant_vesseltips)

                # get the average flow for this realisation 
                flowinlet_analyticalapproxpde_list.append(flow_analyticalapproxpde_ps[0])
                flow_analyticalapproxpde_list.append(flow_analyticalapproxpde)
                flow_constant_list.append(flow_constant)
                flowinlet_constant_list.append(flow_constant_ps[0])

            # scatter plots 
            ax.scatter([sourceterm[SourceNb-1] for i in range(10)], flow_analyticalapproxpde_list, alpha=0.75, color='#757575', marker='.', s = 90)
            ax.scatter([sourceterm[SourceNb-1] for i in range(10)], flow_constant_list, alpha=0.75, color='#c90000', marker='D', s = 17.5)
            ax.scatter([sourceterm[SourceNb-1] for i in range(10)], flowinlet_analyticalapproxpde_list, marker = '.', color='xkcd:teal', alpha=0.75)
            ax.scatter([sourceterm[SourceNb-1] for i in range(10)], flowinlet_constant_list, marker = '.', color='xkcd:green', alpha=0.75)
            
            # add it to the average for this source term
            flow_analyticalapproxpde_average[SourceNb-1] = stats.mean(flow_analyticalapproxpde_list)
            flow_constant_average[SourceNb-1] = stats.mean(flow_constant_list)
            flowinlet_analyticalapproxpde_average[SourceNb-1] = stats.mean(flowinlet_analyticalapproxpde_list)
            flowinlet_constant_average[SourceNb-1] = stats.mean(flowinlet_constant_list)

        # Steady-state
        ax.plot(sourceterm, flow_constant_average, marker = 'D', markersize = 3, label = 'ECM Hypothesis', color='#c90000')
        ax.plot(sourceterm, flow_analyticalapproxpde_average, label = 'Lesion Hypothesis', marker = '.', markersize = 7.5, color='#757575')
        ax.plot(sourceterm, flowinlet_analyticalapproxpde_average, label = 'Inlet Flow Lesion Hyp', marker = '.', markersize = 7.5, color='xkcd:teal')
        ax.plot(sourceterm, flowinlet_constant_average, label = 'Inlet Flow ECM Hyp', marker = 'D', markersize = 7.5, color='xkcd:green')

        ax.set_title('Comparison Outlet Flow / Inlet Flow')
        ax.legend(loc='upper right', fontsize = 12)
        ax.set_xlabel(r'$c_{max}$', fontsize = 15)
        ax.set_ylabel('Flow', fontsize = 13)
        #ax.set_yscale('log')
        ax.set_xticks(sourceterm)
        ax.tick_params(axis = 'both', labelsize = 12)
        plt.show()
        #plt.savefig("/Users/coli171/Library/CloudStorage/OneDrive-TheUniversityofAuckland/Images/ANZIAM2026" + "ProportionFlowInsideLesion.png")

    if FlowInsideLesionTwoGraphs:
        fig, ax = plt.subplots(1, 2, dpi = 300, layout='constrained', sharey = True)

        sourceterm = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
        flow_analyticalapproxpde_average = np.zeros(10)
        flow_constant_average = np.zeros(10)
        flowinlet_analyticalapproxpde_average = np.zeros(10)
        flowinlet_constant_average = np.zeros(10)

        # Constant 
        for SourceNb in range(1, 7):
            flow_analyticalapproxpde_list = []
            flow_constant_list = []
            flowinlet_analyticalapproxpde_list = []
            flowinlet_constant_list = []
            for SeedNb in [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]: 
                if hpc:
                    main_pathway_constant = "/hpc/coli171/Results/PaperAngiogenesisModel2025/PaperModel2025Analysis2D/" + "CoupledModel2DConstant/CoupledModel2DConstantSeed" 
                    main_pathway_analyticalapproxpde = "/hpc/coli171/Results/PaperAngiogenesisModel2025/PaperModel2025Analysis2D/" + "CoupledModel2DAnalyticalApproxPde/CoupledModel2DAnalyticalApproxPdeSeed"
                    main_pathway_constant_seedsource = main_pathway_constant + str(SeedNb) + "Source" + str(SourceNb) + "/results_from_time_0/"
                    main_pathway_analyticalapproxpde_seedsource = main_pathway_analyticalapproxpde + str(SeedNb) + "Source" + str(SourceNb) + "/results_from_time_0/"
                elif local:
                    main_pathway_constant = "/Users/coli171/Desktop/test/Constant/Seed" 
                    main_pathway_analyticalapproxpde = "/Users/coli171/Desktop/test/AnalyticalApproxPde/Seed" 
                    main_pathway_constant_seedsource = main_pathway_constant + str(SeedNb) + "Source" + str(SourceNb) + "/"
                    main_pathway_analyticalapproxpde_seedsource = main_pathway_analyticalapproxpde + str(SeedNb) + "Source" + str(SourceNb) + "/"

                # nodes and elements calculations 
                nodes_analyticalapproxpde, edges_analyticalapproxpde = runner.nodes_elements_calculation(main_pathway_analyticalapproxpde, InitialisationFiles, hpc, local, SeedNb, SourceNb, dim)
                nodes_constant, edges_constant = runner.nodes_elements_calculation(main_pathway_constant, InitialisationFiles, hpc, local, SeedNb, SourceNb, dim)

                indices_vesseltips_constant = runner.find_vesseltips_insidelesion(nodes_constant, edges_constant, ref_point)
                indices_vesseltips_analyticalapproxpde = runner.find_vesseltips_insidelesion(nodes_analyticalapproxpde, edges_analyticalapproxpde, ref_point)

                # calculate flow and pressure for one random seed and one source term
                nodes_analyticalapproxpde_ps, edges_analyticalapproxpde_ps, pressure_analyticalapproxpde_ps, flow_analyticalapproxpde_ps = runner.flow_pressure_calculation(main_pathway_analyticalapproxpde, nodes_analyticalapproxpde, edges_analyticalapproxpde, inlet_pressure, outlet_pressure, umbilical_artery_radius, decay_factor, viscosity_type, SmallSystem, hpc, local, SeedNb, SourceNb, dim)
                nodes_constant_ps, edges_constant_ps, pressure_constant_ps, flow_constant_ps = runner.flow_pressure_calculation(main_pathway_constant, nodes_constant, edges_constant, inlet_pressure, outlet_pressure, umbilical_artery_radius, decay_factor, viscosity_type, SmallSystem, hpc, local, SeedNb, SourceNb, dim)

                flow_constant_vesseltips = [flow_constant_ps[i-1] for i in indices_vesseltips_constant]
                flow_analyticalapproxpde_vesseltips = [flow_analyticalapproxpde_ps[i-1] for i in indices_vesseltips_analyticalapproxpde]

                flow_analyticalapproxpde = 0
                if len(flow_analyticalapproxpde_vesseltips)!=0:
                    flow_analyticalapproxpde = stats.mean(flow_analyticalapproxpde_vesseltips)
                
                flow_constant = 0
                if len(flow_constant_vesseltips)!=0:
                    flow_constant = stats.mean(flow_constant_vesseltips)

                # get the average flow for this realisation 
                flowinlet_analyticalapproxpde_list.append(flow_analyticalapproxpde_ps[0])
                flow_analyticalapproxpde_list.append(flow_analyticalapproxpde)
                flow_constant_list.append(flow_constant)
                flowinlet_constant_list.append(flow_constant_ps[0])

            # scatter plots 
            # ax[1].scatter([sourceterm[SourceNb-1] for i in range(10)], flow_analyticalapproxpde_list, alpha=0.75, color='#757575', marker='.', s = 90)
            # ax[0].scatter([sourceterm[SourceNb-1] for i in range(10)], flow_constant_list, alpha=0.75, color='#c90000', marker='D', s = 17.5)
            
            # add it to the average for this source term
            flow_analyticalapproxpde_average[SourceNb-1] = stats.mean(flow_analyticalapproxpde_list)
            flow_constant_average[SourceNb-1] = stats.mean(flow_constant_list)
            flowinlet_analyticalapproxpde_average[SourceNb-1] = stats.mean(flowinlet_analyticalapproxpde_list)
            flowinlet_constant_average[SourceNb-1] = stats.mean(flowinlet_constant_list)

        # Steady-state
        ax[0].plot(sourceterm, flow_constant_average, marker = 'D', markersize = 3, label = 'Outlet Flow Inside Lesion', color='#c90000')
        ax[1].plot(sourceterm, flow_analyticalapproxpde_average, label = 'Outlet Flow Inside Lesion', marker = '.', markersize = 7.5, color='#757575')
        ax[0].plot(sourceterm, flowinlet_constant_average, marker = 'D', markersize = 3, label = 'Inlet Flow', color='xkcd:blue')
        ax[1].plot(sourceterm, flowinlet_analyticalapproxpde_average, label = 'Inlet Flow', marker = '.', markersize = 7.5, color='xkcd:green')

        ax[0].legend(loc='upper right', fontsize = 12)
        ax[1].legend(loc='upper right', fontsize = 12)
        ax[0].set_xlabel(r'$c_{max}$', fontsize = 15)
        ax[1].set_xlabel(r'$c_{max}$', fontsize = 15)
        ax[0].set_ylabel('Flow', fontsize = 13)
        ax[0].set_xticks(sourceterm)
        ax[1].set_xticks(sourceterm)
        ax[0].tick_params(axis = 'both', labelsize = 12)
        ax[1].tick_params(axis = 'both', labelsize = 12)
        ax[1].set_title('Lesion Hypothesis', fontsize = 12)
        ax[0].set_title('ECM Hypothesis', fontsize = 12)
        plt.show()
        #plt.savefig("/Users/coli171/Library/CloudStorage/OneDrive-TheUniversityofAuckland/Images/ANZIAM2026" + "ProportionFlowInsideLesion.png")


## Visualisation ##
if Visualisation:
    # enter pathway, seed and source term to visualise
    SeedNb = 10
    SourceNb = 7
    SteadyState = False
    Constant = True 

    if hpc:
        if SteadyState:
            main_pathway = "/hpc/coli171/Results/PaperAngiogenesisModel2025/PaperModel2025Analysis2D/" + "CoupledModel2DAnalyticalApproxPde/CoupledModel2DAnalyticalApproxPdeSeed"
        elif Constant:
            main_pathway = "/hpc/coli171/Results/PaperAngiogenesisModel2025/PaperModel2025Analysis2D/" + "CoupledModel2DConstant/CoupledModel2DConstantSeed" 
        main_pathway_seedsource = main_pathway + str(SeedNb) + "Source" + str(SourceNb) + "/results_from_time_0/"
    elif local:
        if SteadyState:
            main_pathway = "/Users/coli171/Desktop/test/AnalyticalApproxPde/Seed" 
        elif Constant:
            main_pathway = "/Users/coli171/Desktop/test/Constant/Seed" 
        main_pathway_seedsource = main_pathway + str(SeedNb) + "Source" + str(SourceNb) + "/"
    dim = 2

    if hpc:
        ps.init("openGL3_egl")
    elif local:
        ps.init()

    # nodes and elements calculations 
    nodes, edges = runner.nodes_elements_calculation(main_pathway, InitialisationFiles, hpc, local, SeedNb, SourceNb, dim)

    # calculate flow and pressure for one random seed and one source term
    nodes_ps, edges_ps, pressure_ps, flow_ps = runner.flow_pressure_calculation(main_pathway, nodes, edges, inlet_pressure, outlet_pressure, umbilical_artery_radius, decay_factor, viscosity_type, SmallSystem, hpc, local, SeedNb, SourceNb, dim)

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