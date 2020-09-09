import os # use os.system("command") to execute a shell command
import time
import pickle

# function that writes a jobscript file for PBS
# nodeName can be give to specify a specific node to run on
def write_pbs_jobfile(jobName, walltime, programName, nodeName=None):
    pbsFileName = jobName+".pbs"

    file = open(pbsFileName, "w")
    jobString = """#!/bin/bash
#PBS -N %s
#PBS -l walltime=%s
#PBS -l select=1:ncpus=1
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

if __name__ == '__main__':
    """
    Job submission script for doing simulations 
    """

    programName = "main.py"
    
    t = time.localtime()
    simDirName = "sim_%d-%d-%d-%d_"%(t.tm_mday, t.tm_mon, t.tm_year, t.tm_hour)

    # maximum time the simulation is allowed to take:
    walltime = "30:00:00"

    # number of copies per parameter set
    nsims = 3

    # vary some parameter
    varyParameter = [.5 , .1, .05, .01]
    
    # create a set of datapoints running a simulation for each datapoint
    startDir = os.getcwd()
    
    for sim in range(nsims):
        for paraM in varyParameter:
            simDir = startDir + "/" + simDirName + str(sim) + "_" + str(paraM) + "/" 
    
            pars = {'seconds' : 100000, 'stepsize' : .01, 'tension' : 0,\
                    'basal_tension' : .25, 'apical_tension' : .25,\
                    'lateral_tension' : .1, 'compression' : 4, 'growthrate' : .001,\
                    'init_volume' : 2.6, 'pressure' : 0, 'transition_boundary' : .05,\
                    'new_edge_length' : .05, 'paneth_dist' : 2}
            pars['transition_boundary'] = paraM
            pars['new_edge_length'] = paraM
    
            if os.path.isfile(programName):
                if not os.path.exists(simDir):
                    os.makedirs(simDir)
    
                os.system("cp main.py parameters.py transitions.py vertex_objects.py ?_rings.txt %s" %(simDir))
    
                os.chdir(simDir)
                
                pkl = open('parameters.pickle','wb')
                pickle.dump(pars, pkl)
                pkl.close()
                
                print(os.getcwd())
                # write_config(config)
    
                jobName = "vertex_" + str(sim) + "_" + str(paraM)
    
                pbsFileName = write_pbs_jobfile(jobName, walltime, programName)
    
                os.system("qsub %s" %pbsFileName)
                os.chdir(startDir)
                
                time.sleep(0.1)
    
            else:
                print("The program named: " + programName + " could not be found")
