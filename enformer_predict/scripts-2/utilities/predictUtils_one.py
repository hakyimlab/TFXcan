# import parsl
# from parsl.app.app import python_app

import os, sys
from functools import lru_cache

whereis_script = os.path.dirname(__file__) #os.path.dirname(sys.argv[0]) # or os.path.dirname(__file__)
global module_path # this variable should be a global one and will be used by the functions defined below
module_path = os.path.abspath(whereis_script)

