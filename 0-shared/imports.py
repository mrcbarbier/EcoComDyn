# -*- coding: utf-8 -*-
import numpy as np
from datatools import Path
import scipy.stats as scipystats
from copy import deepcopy
from itertools import product as iproduct


import os, inspect
REFPATH=Path(os.path.realpath(os.path.abspath(os.path.split(inspect.getfile(
    inspect.currentframe() ))[0])))