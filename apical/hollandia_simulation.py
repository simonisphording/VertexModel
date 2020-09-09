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
        python36 %s > output.txt
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
    simDirName = "sim_%d-%d-%d_"%(t.tm_mday, t.tm_mon, t.tm_year)

    # maximum time the simulation is allowed to take:
    walltime = "8:00:00"

    # number of copies per parameter set
    nsims = 5 

    # vary some parameter
    varyParameter = [1, 1.2, 1.4, 1.6, 1.8, 2]
    
    # create a set of datapoints running a simulation for each datapoint
    startDir = os.getcwd()
    
    for sim in range(nsims):
        for paraM in varyParameter:
            simDir = startDir + "/" + simDirName + str(sim) + "_" + str(paraM) + "/" 
            
            # Setting parameters, of which one is varied per simulation
            pars = {'seconds' : 1000, 'stepsize' : .05, 'radius' : 6,\
                    'center' : [0., 0., 0.], 'tension' : .08, 'contract' : .1,\
                    'init_surface' : 1, 'elastic' : 4, 'gradient' : 1.5,\
                    'transition_boundary' : .05, 'new_edge_length' : .1, 'paneth_dist' : 2}
            pars['gradient'] = paraM
    
            if os.path.isfile(programName):
                if not os.path.exists(simDir):
                    os.makedirs(simDir)
    
                os.system("cp main.py parameters.py transitions.py vertex_objects.py voronoi.py %s" %(simDir))
    
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