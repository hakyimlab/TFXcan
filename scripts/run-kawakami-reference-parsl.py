import parsl
import os
import csv
from parsl.app.app import python_app, bash_app
# from parsl.configs.local_threads import config

from parsl.config import Config
from parsl.providers import CobaltProvider
from parsl.launchers import MpiExecLauncher
from parsl.launchers import AprunLauncher
from parsl.executors import HighThroughputExecutor

import argparse

# Set up parsl config
nodes_per_block = 4

config = Config(
    executors=[
        HighThroughputExecutor(
            label='theta_local_gpu',
            max_workers=8,
            available_accelerators=8,
            provider=CobaltProvider(
                queue='full-node',
                account='covid-ct',
                launcher=MpiExecLauncher(),
                walltime='00:360:00',
                nodes_per_block=nodes_per_block,
                init_blocks=1,
                min_blocks=1,
                max_blocks=1,
                # string to prepend to #COBALT blocks in the submit
                # script to the scheduler eg: '#COBALT -t 50'
                scheduler_options='#COBALT --attrs enable_ssh=1',
                # Command to be run before starting a worker, such as:
                # 'module load Anaconda; source activate parsl_env'.
                worker_init='module load conda;conda activate; which python',
                cmd_timeout=120
            ),
        )   
    ],
)

print(parsl.__version__)
parsl.load(config)