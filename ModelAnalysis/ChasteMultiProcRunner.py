# This is a class to run multiple instances of a Chaste simulation with different parameters
# Note 1: You NEED to include a flag for the output directory that is read into the Chaste code
#         Otherwise all the results will write over each other into the same directory
# Note 2: Change the paths in the init for what you use
#

from spython.main import Client # Singularity client
import multiprocess as mp # Parallel processing
from re import sub
from pathlib import Path
from subprocess import run # For running docker
import platform


class ChasteMultiProcRunner:
    """A class to deal with running several instances of chaste at once"""
    def __init__(self,  
                f_executable,
                num_processors = 10,
                verbose = False):
        self.f_executable = f_executable
        self.n_procs = num_processors
        self.verbose = verbose
        self.use_docker = (platform.system()=="Darwin") # Mac use docker, linux use singularity
        # Setup our defaults based on our platform
        self.local_build_dir = "/Users/coli171/Chaste/Build/"
        if (self.use_docker):
            self.local_output_dir = "/Users/coli171/Chaste/Output/"
            self.container_path = "chaste/base"
        else:
            self.local_output_dir ="/hpc/coli171/Results"
            self.container_path = "/people/coli171/Containers/base_latest.sif"
            

    def RunChasteInstance(self, setup_params):
        # Convert the dictionary to an array, checking that there is a '-' at the front of the key
        param_array = [ x for key,value in setup_params.items() for x in [sub('(^\w)','-\\1',key), str(value)]]

        # Add to the executable
        command_array = [self.f_executable, *param_array]
        options = ["--bind", f"{self.local_output_dir}:/home/chaste/testoutput"]
        if (self.verbose):
            print("Running command:")
            print(command_array)
            print(options)

        # add a file in the output directory where to store the parameters 
        f = open(self.local_output_dir + "parameters.txt", "w")
        for item in param_array :
            f.write(item)
        f.close()

        # Run the instance
        return_code = 0
        if (self.use_docker):
            docker_command_array = ["docker","run", "--rm",
                                    "-v", self.local_build_dir + ":/home/chaste/build",
                                    "-v", self.local_output_dir + ":/home/chaste/testoutput",
                                    self.container_path]
            command_array = [*docker_command_array, *command_array]
            run_result = run(command_array, capture_output=True)
            return_code = run_result.returncode
            output = run_result.stdout.decode()
        else: 
            singularity_instance = Client.instance(self.container_path,options=options)
            run_result = Client.execute(singularity_instance,command_array,return_result=True)
            singularity_instance.stop()
            return_code = run_result['return_code']
            output = run_result['message']

        # Deal with any run errors
        if (return_code != 0):
            msg = f"Non-zero return code ({return_code}) in chaste runner.\n\nCommand run: {command_array}\n\nExecutable output:\n{output}\n" 
            raise Exception(msg)

        # Return the simulation output string
        return output
    

    def RunMultipleChasteInstancesFromList(self,sim_input_list):
        # Not in parallel 
        # output = map(self.RunChasteInstance, sim_input_list)
        # return(output)

        # Run the instances in parallel
        pool = mp.Pool(self.n_procs)
        output = pool.map(self.RunChasteInstance, sim_input_list)
        return(output)