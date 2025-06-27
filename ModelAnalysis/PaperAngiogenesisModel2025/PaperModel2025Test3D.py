# Setup details: to modify when 3D tests
from ChasteMultiProcRunner import ChasteMultiProcRunner
import numpy as np 
import cmd 
import os 
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-valuerandomseed")
parser.add_argument("-sourceterm")
args = parser.parse_args()
path_executable = "build/projects/AngiogenesisModel/test/TestAngiogenesisModelWithVegfConcentrationPde"

input_list = []
value_sourceterm = args.sourceterm
RandomSeedNb = int(args.valuerandomseed)
SourceTermNb = int(float(args.sourceterm)*10)
input_list.append({"vegfdudtcoeff": 0.1, "vegfdiffusioncoeff": 200, "vegfdecaycoeff" : 10.0, "vegfcreationcoeff": 8.3, "vegfconsumptioncoeff": 1.0, "vegfinitialvalue": 0.1, "vegfboundaryvalue": args.sourceterm, "vegfconstantbackground": 0.1, "sigma": 1e-5, "chi": 1e-2, "omegap": 1e-3, "omegaa": 1e-1, "maxsproutingrate": 0.2, "psproutfunctiontestnb": 1, "time":500, "seed": args.valuerandomseed, "output_directory": "PaperModel2025Analysis3D/CoupledModel3DPdeSeed" + str(RandomSeedNb) + "Source" + str(SourceTermNb), "output_directory_analyticalapproxvegf": "PaperAngiogenesisModel2025/PaperModel2025Analysis3D/CoupledModel3DAnalyticalApproxPde/CoupledModel3DAnalyticalApproxPdeSeed" + str(RandomSeedNb) + "Source" + str(SourceTermNb), "output_directory_vegfconstant": "PaperAngiogenesisModel2025/PaperModel2025Analysis3D/CoupledModel3DConstant/CoupledModel3DConstantSeed" + str(RandomSeedNb) + "Source" + str(SourceTermNb)})

# Setup the multiproc class runner
runner = ChasteMultiProcRunner(f_executable=path_executable)
# Run the simulations and print output to log (print error to stdout if exception thrown)
try:
    output = runner.RunMultipleChasteInstancesFromList(input_list)
    output = f"\n".join(output)
except Exception as e:
    output = e
    print("\nERROR (also copied to log): " + str(e))