from datatools import *
ADD_FOLDER_TO_PATH('/1-assembly/cavity_interp/')


RERUN='rerun' in sys.argv


### FIGURE 2 -- RESOURCE COMPETITION
print 'FIGURE 2'
from resource_comp import *
path = Path('RESULTS/resource_comp')

resource_gen(path=path,
            mode='RX',
            rerun=RERUN,measure=1,nsys=40,show=1,imit=0,hold=1)
plt.suptitle('Figure 2')
plt.show()
### FIGURE 3 -- PARAMETER SPACE
print 'FIGURE 3'

from run import *
path = Path('RESULTS/intersection')
intersection(path=path,rerun=RERUN,nsys=100,measure=0,show=1,hold=1)
plt.title('Figure 3 inset')

compare_prms_models(
    [
        (model_resource_generalists,
         {'legend':'Resource', 'vars': ('R', 'X'), 'S': 20, 'R': np.logspace(1, 4, 15), 'X': np.linspace(.3, .9, 15), 'sigma_rho': .25}),
        (model_trophic_satur,
         {'legend':'Predation','vars': ('m0', 'e'), 'm0': np.linspace(0.1, 25, 10), 'e': np.linspace(0, 1, 10)[1:-1]}),
    ],
    mode='hull',hold=1,
)
plt.suptitle('Figure 3')

### FIGURE 4 -- NETWORK METRICS
print 'FIGURE 4'
from network import ordering,show_ordering
path = Path('RESULTS/ordering')

plt.figure()
for mode in ['competition', 'mutualism', 'predation']:
    ordering(path=Path(path),funcresp=False,
             mode=mode, struct='all',
             rerun=RERUN, measure=0, nsys=100, sysmin=0, show=0)

show_ordering(path=path,hold=1)
plt.suptitle('Figure 4')

### FIGURE 5 -- BIPARTITE COMMUNITY
print 'FIGURE 5'

from mix import *
path=Path('RESULTS/mixture')
if 0:
    mixture(path=Path(path),
            mode=['dft'][0],
            rerun=RERUN,measure=0,nsys=30,show=1,hold=1)

    plt.suptitle('Figure 5')


plt.show()

