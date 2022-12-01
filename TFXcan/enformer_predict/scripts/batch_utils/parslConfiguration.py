def localParslConfig():
    from tqdm import tqdm
    from parsl import python_app
    import parsl
    import json

    # Make a config that runs on two nodes
    from parsl.executors import HighThroughputExecutor
    from parsl.providers import LocalProvider
    from parsl.config import Config
    from parsl.channels import LocalChannel
    from parsl.launchers import MpiExecLauncher

    rundir = '/projects/covid-ct/imlab/users/temi/projects/TFXcan/enformer_predict/runinfo'
    workingdir = '/projects/covid-ct/imlab/users/temi/projects/TFXcan/enformer_predict'

    #parsl.clear()

    local_htex = Config(
        executors=[
            HighThroughputExecutor(
                label="htex_Local",
                max_workers=8, # vs max_workers
                available_accelerators=8,
                worker_debug=True,
                cores_per_worker=1, # how many cores per worker #nodes_per_block, 2 is usually enough or 1.
                working_dir=workingdir,
                provider=LocalProvider(
                    channel=LocalChannel(),
                    init_blocks=1,
                    max_blocks=1,
                    launcher=MpiExecLauncher(),
                    worker_init='source ~/.bashrc; conda activate dl-tools; which python; export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/temi/miniconda3/envs/dl-tools/lib'
                ),
            )
        ],
        strategy=None,
        run_dir=rundir
    )

    return(local_htex)


def htParslConfig(workingdir=None):

    import parsl
    from parsl.config import Config
    from parsl.providers import CobaltProvider
    from parsl.launchers import MpiExecLauncher
    from parsl.executors import HighThroughputExecutor

    print(f'Parsl version: {parsl.__version__}')

    # I defined these locations otherwise parsl will use the current directory to output the run informations and log messages
    if workingdir is None:
        workingdir = '/projects/covid-ct/imlab/users/temi/projects/TFXcan/enformer_predict'
        rundir = f'{workingdir}/runinfo'
    
    # I want to put the cobalt directives 
    sch_options = ['#COBALT --attrs filesystems=home,theta-fs0,grand,eagle:enable_ssh=1',
                    '#COBALT --jobname=enformer-predict-personalized',
                    '#COBALT -o /projects/covid-ct/imlab/users/temi/projects/TFXcan/enformer_predict/cobalt-log/enformer-predict-personalized.out',
                    '#COBALT -e /projects/covid-ct/imlab/users/temi/projects/TFXcan/enformer_predict/cobalt-log/enformer-predict-personalized.err',
                    '#COBALT --debuglog /projects/covid-ct/imlab/users/temi/projects/TFXcan/enformer_predict/cobalt-log/enformer-predict-personalized.cobalt'
    ]

    sch_options = '\n'.join(sch_options)

    # worker init
    workerinit = 'source ~/.bashrc; conda activate dl-tools; which python; export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/temi/miniconda3/envs/dl-tools/lib'

    #just in case there is a config loaded already
    cobalt_htex = Config(
        executors=[
            HighThroughputExecutor(
                label='htex-Cobalt',
                max_workers=8, # vs max_workers
                available_accelerators=8,
                worker_debug=True,
                cores_per_worker=1, # how many cores per worker #nodes_per_block, 2 is usually enough or 1.
                working_dir=workingdir,
                provider=CobaltProvider(
                    queue='full-node',
                    account='covid-ct',
                    launcher=MpiExecLauncher(),
                    walltime='00:20:00',
                    nodes_per_block=2, # number of full-nodes - 3 will launch 3 full nodes at a t   ime for one instance for each `cores_per_worker`
                    min_blocks=1,
                    #max_blocks=2,
                    worker_init=workerinit,
                    cmd_timeout=120,
                    scheduler_options=sch_options
                ),
            )   
        ],
        run_dir=rundir,
        retries=6
    )

    return(cobalt_htex)


def polaris_htParslConfig(workingdir=None):

    import parsl
    from parsl.config import Config
    from parsl.executors import HighThroughputExecutor
    from parsl.launchers import MpiExecLauncher
    from parsl.addresses import address_by_hostname
    from parsl.providers import PBSProProvider

    print(f'Parsl version: {parsl.__version__}')

    # I defined these locations otherwise parsl will use the current directory to output the run informations and log messages
    if workingdir is None:
        workingdir = '/grand/projects/covid-ct/imlab/users/temi/projects/TFXcan/enformer_predict'
        rundir = f'{workingdir}/runinfo'
    
    # # I want to put the cobalt directives 
    # sch_options = ['#COBALT --attrs filesystems=home,theta-fs0,grand,eagle:enable_ssh=1',
    #                 '#COBALT --jobname=enformer-predict-personalized',
    #                 '#COBALT -o /grand/projects/covid-ct/imlab/users/temi/projects/TFXcan/enformer_predict/cobalt-log/enformer-predict-personalized.out',
    #                 '#COBALT -e /grand/projects/covid-ct/imlab/users/temi/projects/TFXcan/enformer_predict/cobalt-log/enformer-predict-personalized.err',
    #                 '#COBALT --debuglog /grand/projects/covid-ct/imlab/users/temi/projects/TFXcan/enformer_predict/cobalt-log/enformer-predict-personalized.cobalt'
    # ]

    # sch_options = '\n'.join(sch_options)

    # worker init
    #workerinit = 'source ~/.bashrc; conda activate dl-tools; which python; export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/temi/miniconda3/envs/dl-tools/lib'


    user_opts = {
        'polaris': {
            # Node setup: activate necessary conda environment and such.
            'worker_init': 'source ~/.bashrc; conda activate dl-tools; which python; export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/temi/miniconda3/envs/dl-tools/lib',
            'scheduler_options': '#PBS -l filesystems=home:grand:eagle',
            # ALCF allocation to use
            'account': 'covid-ct',
        }
    }

    config = Config(
        executors=[
            HighThroughputExecutor(
                available_accelerators=4,  # Pin each worker to a different GPU
                            max_workers=4,
                #strategy=SimpleStrategy(max_idletime=300),
                address=address_by_hostname(),
                provider=PBSProProvider(
                    launcher=MpiExecLauncher(
                        bind_cmd="--cpu-bind", overrides="--depth=64 --ppn 1"
                    ),  # Ensures 1 manger per node, work on all 64 cores
                    account=user_opts['polaris']['account'],
                    queue='preemptable',
                    cpus_per_node=32,
                    select_options='ngpus=4',
                    worker_init=user_opts['polaris']['worker_init'],
                    scheduler_options=user_opts['polaris']['scheduler_options'],
                    walltime='00:20:00',
                    nodes_per_block=1,
                    init_blocks=0,
                    min_blocks=0,
                    max_blocks=2,
                ),
            )
        ],
        run_dir=rundir
    )

    return config