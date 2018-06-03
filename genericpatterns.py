from datatools import *




### FIGURE 5
ADD_FOLDER_TO_PATH('/1-assembly/cavity_interp/')
from network import ordering,show_ordering
path = Path('RESULTS/ordering')

for mode in ['competition', 'mutualism', 'predation']:
    ordering(path=Path(path),funcresp=False,
             mode=mode, struct='all',
             rerun=1, measure=1, nsys=100, sysmin=0, show=0)

show_ordering(path=path)

