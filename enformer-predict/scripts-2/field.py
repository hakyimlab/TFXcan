g = tf.config.experimental.list_physical_devices('GPU')

tf.config.experimental.get_memory_info(g[0][0])

tf.config.experimental.get_memory_info('GPU:0')


import subprocess


def get_gpu_memory():
    import subprocess
    command = "nvidia-smi --query-gpu=memory.free,memory.used --format=csv"
    memory_info = subprocess.check_output(command.split()).decode('ascii').split('\n')[1].split(',')
    memory_values = [int(x.strip().split(' ')[0]) for i, x in enumerate(memory_info)]
    return memory_values

get_gpu_memory()



command = "nvidia-smi -q"
subprocess.check_output(command.split())