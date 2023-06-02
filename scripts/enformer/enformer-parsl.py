#


import parsl
import os
from parsl.app.app import python_app, bash_app
from parsl.configs.local_threads import config

#parsl.set_stream_logger() # <-- log everything to stdout

print(parsl.__version__)

config = Config(
    executors=[
        HighThroughputExecutor(
            label='theta_local_gpu',
            max_workers=8,
            available_accelerators=8,
            provider=CobaltProvider(
                queue=args.queue,
                account='covid-ct',
                launcher=MpiExecLauncher(),
                walltime=args.walltime,
                nodes_per_block=nodes_per_block,
                init_blocks=1,
                min_blocks=1,
                max_blocks=1,
                # string to prepend to #COBALT blocks in the submit
                # script to the scheduler eg: '#COBALT -t 50'
                scheduler_options='#COBALT --attrs enable_ssh=1',
                # Command to be run before starting a worker, such as:
                # 'module load Anaconda; source activate parsl_env'.
                worker_init='module load conda;conda activate /projects/covid-ct/imlab/users/saideep/enformerxcan_naive/envs; which python',
                cmd_timeout=120
            ),
        )   
    ],
)

parsl.load(config)