g = tf.config.experimental.list_physical_devices('GPU')

tf.config.experimental.get_memory_info(g[0][0])

tf.config.experimental.get_memory_info('GPU:0')


import subprocess as sp
import os

def get_gpu_memory():
    command = "nvidia-smi --query-gpu=memory.free --format=csv"
    memory_free_info = sp.check_output(command.split()).decode('ascii').split('\n')[:-1][1:]
    memory_free_values = [int(x.split()[0]) for i, x in enumerate(memory_free_info)]
    return memory_free_values

get_gpu_memory()



command = "nvidia-smi -q"
sp.check_output(command.split())