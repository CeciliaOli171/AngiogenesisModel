# Setup details
from ChasteMultiProcRunner import ChasteMultiProcRunner
import numpy as np 
import cmd 
import os 
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-valuerandomseed")
parser.add_argument("-countersigma")
parser.add_argument("-counteromegap")
args = parser.parse_args()
path_executable = "build/projects/AngiogenesisModel/test/TestAngiogenesisModel"
#log_file = "~/chaste.log"

#def MultipleTestParametersSensitivity(args.value_sigma, args.value_omegap, args.value_random_seed, RandomSeedNb):
input_list = []
listvalues_omegaa = [1.0, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8]  
value_sigma = 10**(-int(args.countersigma))
value_omegap = 10**(-int(args.counteromegap))
counter_omegaa = 0
RandomSeedNb = int(int(args.valuerandomseed)/10)
for value_omegaa in listvalues_omegaa:
    input_list.append({"output_directory": "TestParametersSensitivity/TestAreaVolumeTortuosity/Test" + str(RandomSeedNb) + str(args.countersigma) + str(args.counteromegap) + str(counter_omegaa), "sigma": value_sigma, "chi": 0.0, "omegap": value_omegap, "omegaa": value_omegaa, "Psprout": 0.1, "time": 300.0, "seed": args.valuerandomseed})
    counter_omegaa += 1

# Setup the multiproc class runner
runner = ChasteMultiProcRunner(f_executable=path_executable)
# Run the simulations and print output to log (print error to stdout if exception thrown)
try:
    output = runner.RunMultipleChasteInstancesFromList(input_list)
    output = f"\n".join(output)
except Exception as e:
    output = e
    print("\nERROR (also copied to log): " + str(e))
#finally:
#    with open(log_file, 'w') as f:
#        print(output,file=f)
