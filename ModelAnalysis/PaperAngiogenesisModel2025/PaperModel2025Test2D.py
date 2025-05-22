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
path_executable = "build/projects/AngiogenesisModel/test/TestAngiogenesisModelWithVegfConcentrationPdes"

input_list = []
value_sourceterm = args.sourceterm
RandomSeedNb = int(args.valuerandomseed)
SourceTermNb = int(float(args.sourceterm)*10)

# Analytical Approx Pde: hill function
input_list.append({"vegfdudtcoeff": 0.1, "vegfdiffusioncoeff": 1e4, "vegfdecaycoeff" : 1.0, "vegfcreationcoeff": 0.0, "vegfconsumptioncoeff": 1.0, "vegfinitialvalue": 0.1, "vegfboundaryvalue": args.sourceterm, "vegfconstantbackground": 0.1, "sigma": 1e-3, "chi": 1e-2, "omegap": 1e-3, "omegaa": 1e-1, "maxsproutingrate": 0.02, "psproutfunctiontestnb": 1, "time": 3000, "seed": args.valuerandomseed, "output_directory": "PaperAngiogenesisModel2025/PaperModel2025Analysis2D/CoupledModel2DPde/CoupledModel2DPdeSeed" + str(RandomSeedNb) + "Source" + str(SourceTermNb), "output_directory_analyticalapproxvegf": "PaperAngiogenesisModel2025/PaperModel2025Analysis2D/CoupledModel2DAnalyticalApproxPde/CoupledModel2DAnalyticalApproxPdeSeed" + str(RandomSeedNb) + "Source" + str(SourceTermNb), "output_directory_vegfconstant": "PaperAngiogenesisModel2025/PaperModel2025Analysis2D/CoupledModel2DConstant/CoupledModel2DConstantSeed" + str(RandomSeedNb) + "Source" + str(SourceTermNb)})

# Analytical Approx Pde: linear function
# input_list.append({"vegfdudtcoeff": 0.1, "vegfdiffusioncoeff": 200, "vegfdecaycoeff" : 10.0, "vegfcreationcoeff": 8.3, "vegfconsumptioncoeff": 1.0, "vegfinitialvalue": 0.1, "vegfboundaryvalue": args.sourceterm, "vegfconstantbackground": 0.1, "sigma": 1e-5, "chi": 1e-2, "omegap": 1e-3, "omegaa": 1e-1, "maxsproutingrate": 0.02, "psproutfunctiontestnb": 0, "time": 1000, "seed": args.valuerandomseed, "output_directory": "PaperModel2025Analysis2D/CoupledModel2DPdeSeed" + str(RandomSeedNb) + "Source" + str(SourceTermNb), "output_directory_analyticalapproxvegf": "PaperAngiogenesisModel2025/PaperModel2025Analysis2D/CoupledModel2DAnalyticalApproxPdeLinearFunction/CoupledModel2DAnalyticalApproxPdeSeed" + str(RandomSeedNb) + "Source" + str(SourceTermNb), "output_directory_vegfconstant": "PaperAngiogenesisModel2025/PaperModel2025Analysis2D/CoupledModel2DConstantLinearFunction/CoupledModel2DConstantSeed" + str(RandomSeedNb) + "Source" + str(SourceTermNb)})

# Constant: problem for random seed 3
# input_list.append({"vegfdudtcoeff": 0.1, "vegfdiffusioncoeff": 200, "vegfdecaycoeff" : 10.0, "vegfcreationcoeff": 8.3, "vegfconsumptioncoeff": 1.0, "vegfinitialvalue": 0.1, "vegfboundaryvalue": 0.1, "vegfconstantbackground": args.sourceterm, "sigma": 1e-5, "chi": 1e-2, "omegap": 1e-3, "omegaa": 1e-1, "maxsproutingrate": 0.02, "psproutfunctiontestnb": 1, "time": 1000, "seed": args.valuerandomseed, "output_directory": "PaperModel2025Analysis2D/CoupledModel2DPdeSeed" + str(RandomSeedNb) + "Source" + str(SourceTermNb), "output_directory_analyticalapproxvegf": "PaperAngiogenesisModel2025/PaperModel2025Analysis2D/CoupledModel2DAnalyticalApproxPde/CoupledModel2DAnalyticalApproxPdeSeed" + str(RandomSeedNb) + "Source" + str(SourceTermNb), "output_directory_vegfconstant": "PaperAngiogenesisModel2025/PaperModel2025Analysis2D/CoupledModel2DConstant/CoupledModel2DConstantSeed" + str(RandomSeedNb) + "Source" + str(SourceTermNb)})

# Setup the multiproc class runner
runner = ChasteMultiProcRunner(f_executable=path_executable)
# Run the simulations and print output to log (print error to stdout if exception thrown)
try:
    output = runner.RunMultipleChasteInstancesFromList(input_list)
    output = f"\n".join(output)
except Exception as e:
    output = e
    print("\nERROR (also copied to log): " + str(e))