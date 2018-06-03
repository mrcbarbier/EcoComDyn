from datatools import *
ADD_FOLDER_TO_PATH('/1-assembly/cavity_interp/')


### FIGURE 2 -- RESOURCE COMPETITION
from resource_comp import *
path = Path('RESULTS/resource_comp')

resource_gen(path=path,
            mode='RX',
            rerun=1,measure=1,nsys=70,show=1,imit=0)


### FIGURE 3 -- PARAMETER SPACE
from run import *
path = Path('RESULTS/intersection')
intersection(path=path,rerun=1,nsys=100,measure=1,show=1)
compare_prms_models(
    [
        (model_resource_generalists,
         {'legend':'Resource', 'vars': ('R', 'X'), 'S': 20, 'R': np.logspace(1, 4, 15), 'X': np.linspace(.3, .9, 15), 'sigma_rho': .25}),
        (model_trophic_satur,
         {'legend':'Predation','vars': ('m0', 'e'), 'm0': np.linspace(0.1, 25, 10), 'e': np.linspace(0, 1, 10)[1:-1]}),
    ],
    mode='hull',
)

### FIGURE 4 -- NETWORK METRICS


from network import ordering,show_ordering
path = Path('RESULTS/ordering')

for mode in ['competition', 'mutualism', 'predation']:
    ordering(path=Path(path),funcresp=False,
             mode=mode, struct='all',
             rerun=1, measure=1, nsys=100, sysmin=0, show=0)

show_ordering(path=path)


### FIGURE 5 -- BIPARTITE COMMUNITY
from mix import *
path=Path('RESULTS/mixture')
mixture(path=Path(path),
            mode=['dft'][0],
            rerun=0,measure=0,nsys=30,show=1)