

def polaris_htParslConfig(params):

    import parsl
    from parsl.config import Config
    from parsl.executors import HighThroughputExecutor
    from parsl.launchers import MpiExecLauncher
    from parsl.addresses import address_by_hostname
    from parsl.providers import PBSProProvider

    print(f'Parsl version: {parsl.__version__}')
    # I defined these locations otherwise parsl will use the current directory to output the run informations and log messages
    workingdir = params['working_dir']
    rundir = f'{workingdir}/runinfo'
    job_name = params['job_name']
    
    # I want to put the cobalt directives 
    sch_options = ['#PBS -l filesystems=home:grand:eagle',
                    f'#PBS -N {job_name}',
                    f'#PBS -k doe'
    ]

    sch_options = '\n'.join(sch_options)

    user_opts = {
        'polaris': {
            # Node setup: activate necessary conda environment and such.
            'worker_init': 'source /home/temi/.bashrc; conda activate dl-tools; which python; export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/temi/miniconda3/envs/dl-tools/lib; echo Running on host `hostname`; echo Running on nodes `cat $PBS_NODEFILE`',
            'scheduler_options': f'{sch_options}',
            # ALCF allocation to use
            'account': 'covid-ct',
        }
    }

    pbs_htex = Config(
        executors=[
            HighThroughputExecutor(
                label='htex_PBS',
                available_accelerators=4,  # Pin each worker to a different GPU
                max_workers=4,
                cpu_affinity = 'alternating',
                address=address_by_hostname(),
                provider=PBSProProvider(
                    launcher=MpiExecLauncher(
                        bind_cmd="--cpu-bind", overrides="--depth=64 --ppn 1"
                    ),  # Ensures 1 manger per node, work on all 64 cores
                    account=user_opts['polaris']['account'],
                    queue=params['queue'], #preemptable',
                    cpus_per_node=64,
                    select_options='ngpus=4',
                    worker_init=user_opts['polaris']['worker_init'],
                    scheduler_options=user_opts['polaris']['scheduler_options'],
                    walltime=params['walltime'],
                    nodes_per_block=params['num_of_full_nodes'],
                    init_blocks=1,
                    min_blocks=params['min_num_blocks'],
                    max_blocks=params['max_num_blocks'],
                ),
            )
        ],
        run_dir=rundir,
        retries=6
    )
    return pbs_htex

def polaris_localParslConfig(params):

    import parsl
    # Make a config that runs on two nodes
    from parsl.executors import HighThroughputExecutor
    from parsl.providers import LocalProvider
    from parsl.config import Config
    from parsl.channels import LocalChannel
    from parsl.launchers import MpiExecLauncher
    from parsl.addresses import address_by_hostname

    import os

    print(f'Parsl version: {parsl.__version__}')
    # I defined these locations otherwise parsl will use the current directory to output the run informations and log messages
    workingdir = params['working_dir']
    rundir = os.path.join(workingdir, 'runinfo')
    job_name = params['job_name']

    local_htex = Config(
        executors=[
            HighThroughputExecutor(
                label="htex_Local",
                max_workers=4, # vs max_workers
                available_accelerators=4,
                worker_debug=True,
                cores_per_worker=32, # how many cores per worker #nodes_per_block, 2 is usually enough or 1.
                working_dir=workingdir,
                address=address_by_hostname(),
                provider=LocalProvider(
                    channel=LocalChannel(),
                    init_blocks=1,
                    nodes_per_block=params['num_of_full_nodes'],
                    min_blocks=params['min_num_blocks'],
                    max_blocks=params['max_num_blocks'],
                    launcher=MpiExecLauncher(bind_cmd="--cpu-bind", overrides="--depth=64 --ppn 1"),
                    worker_init='source ~/.bashrc; conda activate dl-tools; which python; export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/temi/miniconda3/envs/dl-tools/lib'
                ),
            )
        ],
        strategy=None,
        run_dir=rundir
    )

    return(local_htex)