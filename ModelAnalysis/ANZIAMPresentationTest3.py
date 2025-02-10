# Setup details
from ChasteMultiProcRunner import ChasteMultiProcRunner
import numpy as np 
import cmd 
import os 
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-valuerandomseed")
parser.add_argument("-sourceterm")
args = parser.parse_args()
path_executable = "build/projects/AngiogenesisModel/test/TestAngiogenesisModelWithMolecularConcentrationsPdes"

input_list = []
value_sourceterm = args.sourceterm
RandomSeedNb = int(args.valuerandomseed)
SourceTermNb = int(float(args.sourceterm)*10)
input_list.append({"vegfdudtcoeff": 0.1, "vegfdiffusioncoeff": 200, "vegfdecaycoeff" : 10.0, "vegfcreationcoeff": 8.3, "vegfconsumptioncoeff": 1.0, "vegfinitialvalue": 0.1, "vegfboundaryvalue": args.sourceterm, "output_directory": "ModelWithPdes2DTestAnalyticalSeedSource/Test3/ModelWithPdes2DTest3PdeSeed" + str(RandomSeedNb) + "Source" + str(SourceTermNb), "output_directory_analyticalvegf": "ModelWithPdes2DTestAnalyticalSeedSource/Test3/ModelWithPdes2DTest3AnalyticalSeed" + str(RandomSeedNb) + "Source" + str(SourceTermNb), "maxsproutingrate": 0.2, "sigma": 1e-5, "chi": 1e-2, "omegap": 1e-3, "omegaa": 1e-1, "seed": args.valuerandomseed, "testnb": 3})

# Setup the multiproc class runner
runner = ChasteMultiProcRunner(f_executable=path_executable)
# Run the simulations and print output to log (print error to stdout if exception thrown)
try:
    output = runner.RunMultipleChasteInstancesFromList(input_list)
    output = f"\n".join(output)
except Exception as e:
    output = e
    print("\nERROR (also copied to log): " + str(e))