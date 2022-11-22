import os, sys, json
import pandas as pd # for manipulating dataframes
from functools import lru_cache
import math, time, tqdm
import parsl
from parsl.app.app import python_app

from parsl.config import Config
from parsl.providers import CobaltProvider
from parsl.launchers import MpiExecLauncher
from parsl.executors import HighThroughputExecutor

print(f'Parsl version: {parsl.__version__}')

# I defined these locations otherwise parsl will use the current directory to output the run informations and log messages
rundir = '/projects/covid-ct/imlab/users/temi/projects/TFXcan/enformer-predict/runinfo'
workingdir = '/projects/covid-ct/imlab/users/temi/projects/TFXcan/enformer-predict'

# I want to put the cobalt directives 
sch_options = ['#COBALT --attrs filesystems=home,theta-fs0,grand,eagle:enable_ssh=1',
                '#COBALT --jobname=enformer-predict-personalized',
                '#COBALT -o /projects/covid-ct/imlab/users/temi/projects/TFXcan/enformer-predict/cobalt-log/enformer-predict-personalized.out',
                '#COBALT -e /projects/covid-ct/imlab/users/temi/projects/TFXcan/enformer-predict/cobalt-log/enformer-predict-personalized.err',
                '#COBALT --debuglog /projects/covid-ct/imlab/users/temi/projects/TFXcan/enformer-predict/cobalt-log/enformer-predict-personalized.cobalt'
]

sch_options = '\n'.join(sch_options)

#just in case there is a config loaded already
cobalt_htex = Config(
    executors=[
        HighThroughputExecutor(
            label='htex-Cobalt',
            max_workers=8,
            available_accelerators=8,
            worker_debug=True,
            cores_per_worker=1, # how many cores per worker #nodes_per_block, 2 is usually enough or 1.
            working_dir=workingdir,
            provider=CobaltProvider(
                queue='full-node',
                account='covid-ct',
                launcher=MpiExecLauncher(),
                walltime='00:10:00',
                nodes_per_block=1, # number of full-nodes - 3 will launch 3 full nodes at a time for one instance for each `cores_per_worker`
                min_blocks=1,
                #max_blocks=3,
                worker_init='source ~/.bashrc; conda activate dl-tools; which python',
                cmd_timeout=120,
                scheduler_options=sch_options
            ),
        )   
    ],
    run_dir=rundir,
    retries=6
)

parsl.load(cobalt_htex)

# function

# @python_app
# def sumtwonums(a, b):
#     return(a+b)

@python_app
def tf_version(i):
    import tensorflow as tf
    return(f'[INFO {i}] Tensorflow found {tf.config.list_physical_devices()}')

mpath = os.path.abspath(os.path.dirname(__file__))
print(mpath)
#sys.path.append(mpath)

futures = [tf_version(i) for i in range(0, 20)]
print(futures)
execution = [o.result() for o in futures]
print(execution)

#saving
for i, value in enumerate(execution):
    with open(f'{mpath}/output/file_{i}.txt', 'w') as wf:
        wf.writelines(value)