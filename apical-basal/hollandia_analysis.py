import os
import time

# function that writes a jobscript file for PBS
# nodeName can be give to specify a specific node to run on
def write_pbs_jobfile(jobName, walltime, programName, nodeName=None):
    pbsFileName = jobName+".pbs"

    file = open(pbsFileName, "w")
    jobString = """#!/bin/bash
#PBS -N %s
#PBS -l walltime=%s
#PBS -l select=1:ncpus=1
#PBS -l mem=4gb
#PBS -V

cd $PBS_O_WORKDIR

mkdir -p $TMPDIR
pwd
python36 %s > output.txt 2> error.txt""" % (jobName, walltime, programName)

    if (nodeName):
        jobString = """#!/bin/bash
        #PBS -N %s
        #PBS -l walltime=%s
        #PBS -l nodes=%s:ppn=1
        #PBS -V
        
        cd $PBS_O_WORKDIR
        mkdir -p $TMPDIR
        pwd
        python36 %s > output.txt 2> error.txt
        """ % (jobName, walltime, nodeName, programName)

    file.write(jobString)

    file.close()

    return pbsFileName

# do something (this allows this file to be included somewhere without executing the code below)    
if __name__ == '__main__':
    """
    Small example of a job submission python script 
    """

    programName = "analysing.py"
    simDirName = "sim_"

    # maximum time the simulation is allowed to take:
    walltime = "8:00:00"
    
    startDir = os.getcwd()
    
    result_folders = []
    for file in os.listdir():
        if 'sim_' in file:
            result_folders.append(file)
    
    for resultsDir in result_folders:
        os.system("cp analysing.py %s" %resultsDir)
        
        os.chdir(resultsDir)
        
        print(os.getcwd())
        
        jobName = "analysing_" + resultsDir
        
        pbsFileName = write_pbs_jobfile(jobName, walltime, programName)
        
        os.system("qsub %s" %pbsFileName)
        os.chdir(startDir)
        
        time.sleep(0.1)
