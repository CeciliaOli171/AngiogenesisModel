# Setup details
import os, sys
sys.path.insert(1, "/".join(os.path.realpath(__file__).split("/")[0:-2]))

from ChasteMultiProcRunner import ChasteMultiProcRunner
import numpy as np 
import cmd 
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-valuerandomseed")
parser.add_argument("-sourceterm")
args = parser.parse_args()
path_executable_constantvegf = "Chaste/Build/projects/AngiogenesisModel/test/TestAngiogenesisModelWithVegfConcentrationConstant"


input_list_constantvegf = []
list_random_seed = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
list_source_term = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]


# Constant
for random_seed in list_random_seed:
    for source_term in list_source_term:
        RandomSeedNb = int(random_seed)
        SourceTermNb = int(float(source_term)*10)
        input_list_constantvegf.append({"vegfdudtcoeff": 0.1, "vegfdiffusioncoeff": 0.44e4, "vegfdecaycoeff" : 1.0, "vegfcreationcoeff": 0.0, "vegfconsumptioncoeff": 1.0, "vegfinitialvalue": 0.1, "vegfboundaryvalue": source_term, "vegfconstantbackground": source_term, "S": 1e-5, "sigma": 1e-1, "chi": 1e-4, "omegap": 1e-1, "omegaa": 1e-5, "maxsproutingrate": 0.08, "psproutfunctiontestnb": 1, "cmax": 0.8, "cmin": 0.3, "pmax": 0.98, "pmin": 0.5, "anastomosislength": 0.75, "thresholdlength": 5.0, "time": 2000, "seed": random_seed, "output_directory": "PaperAngiogenesisModel2025/PaperModel2025Analysis2D/CoupledModel2DPde/CoupledModel2DPdeSeed" + str(RandomSeedNb) + "Source" + str(SourceTermNb), "output_directory_analyticalapproxvegf": "PaperAngiogenesisModel2025/PaperModel2025Analysis2D/CoupledModel2DAnalyticalApproxPde/CoupledModel2DAnalyticalApproxPdeSeed" + str(RandomSeedNb) + "Source" + str(SourceTermNb), "output_directory_vegfconstant": "PaperAngiogenesisModel2025/PaperModel2025Analysis2D/CoupledModel2DConstant/CoupledModel2DConstantSeed" + str(RandomSeedNb) + "Source" + str(SourceTermNb)})


# Setup the multiproc class runner
runner_constantvegf = ChasteMultiProcRunner(f_executable=path_executable_constantvegf)
# Run the simulations and print output to log (print error to stdout if exception thrown)
try:
    output = runner_constantvegf.RunMultipleChasteInstancesFromList(input_list_constantvegf)
    output = f"\n".join(output)
except Exception as e:
    output = e
    print("\nERROR (also copied to log): " + str(e))