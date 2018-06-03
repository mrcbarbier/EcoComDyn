from libc.stdlib cimport rand, srand, RAND_MAX
from libc.math cimport log,exp,sqrt
import random as rnd
import numpy as np
cimport numpy as np

srand(np.random.randint(1,RAND_MAX))


