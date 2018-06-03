# -*- coding: utf-8 -*-
## Ranked species
from interp import *
from show import *


dft_dyn={
    'nlv':{'type':'simplelv',
        'variables':[('n','com'),],
        },
    }

dft_prm={
        'n':{
    'type':'variable',
    'axes': ('com',),
    'death': 10.**-8,
    'shape':(100,),
    'mean':1.,
    'std':1.,
    'sign':1,
    },

        'niches':{
    'type':'constant',
    'role':'niche',
    'axes': [('n','com')],
    'min':0,
    'max':1,
    'distribution':'uniform',
    },



        'growth':{
    'type':'matrix',
    'variables':['n'],
    'mean':1.,
    'meanslope':0,
    'std':0,
    'stdslope':0,
    'dynamics':'nlv',
    'role':'growth',
    'mode':'exponential',
    'structure':'ranked',
    'requires':['niches'],
    },

        'selfint':{
    'type':'matrix',
    'variables':['n'],
    'dynamics':'nlv',
    'mean':1,
    'std':0,
    'role':'diagonal',
    'structure':'random',
    'requires':['niches'],
    },

        'community':{
    'type':'matrix',
    'variables':[('n','com'),('n','com') ],
    'role':'interactions',
    'mean':-.05,
    'meanslope':0,
    'std':0.02,
    'stdslope':0,
    'symmetry':-.5,
    'dynamics':'nlv',
    'mode':'exponential',
    'structure':'ranked',
    'requires':['niches'],
    },

    }

def ranked(path=Path('ranked'),mode=['dft'][0],
        rerun=0,measure=0,nsys=3,show=0):
    #Ordered competitive community (Tilman's colonization-competition tradeoff)

    if mode != 'dft':
        path=Path(path.strip()+'_{}'.format(mode))
    #Niche axis competition (MacArthur Levins)
    try:
        if not os.listdir(path):
            rerun=1
    except:
        rerun =1
        #If not already run before

    kwargs={}
    prm=deepcopy(dft_prm)
    kwargs['n_shape']=(100,)
    prm['community']['symmetry']=.5

    sigmas=np.linspace(0.,1.,11)
    modes=['exponential', 'linear']
    if mode=='dft':
        prm['growth']['std']=.25
        axes=[
                ( 'community_meanslope', sigmas ),
                #('community_stdslope', [0.,0.01]  ) ,
                (('community_mode','growth_mode'),np.array([modes,modes]).T ),
                ( 'sys',range(nsys) )
            ]
    if mode=='K':
        prm['growth']['std']=0
        axes=[
                ( 'growth_meanslope', sigmas ),
                (('community_mode','growth_mode'),np.array([modes,modes]).T ),
                ( 'sys',range(nsys) )
            ]


    if rerun:
        loop(path=path,axes=axes,check_key=0,print_msg=1,tmax=500000,tsample=10000.,
            parameters=prm,dynamics=dft_dyn,
            **kwargs )

    if rerun or measure:
        rebuild_filelist(path)
        make_measures(path,
            #read_filter=('community_meanslope','>',.0),
            #read_filter=('community_mode','==','linear'),
            split_arrays={'contains':'n_groups','add':group_props_list+('n_abundance',) },

            use_measures=['usual',
            'effective',
            'matrix',
            #measure_cavity,
            #measure_ordered,
            #measure_ordered_cavity,
            measure_groups,
            measure_groups_cavity,
            #measure_continuum,
            #measure_continuum_cavity,

            ],death=10**-10)


    if not show:
        return

    measures=open_measures(path)
    axes=[ax[0] for ax in axes][:1]
    axes=[ax if not hasattr(ax,'__iter__') else ax[0] for ax in axes]
    split_by=['community_mode']

    dico={
        'community_meanslope':r'$c_\mu$ (linear) or $\log c_\mu$ (exponential)','community_mode':'Type'
        }


    if 0:
        plot_data(measures=measures, axes=axes,
            split_by=split_by,
            hold=1,dictionary=dico,
            values=[('mean_rank',{'style':'plot'}),('n_mean_rank',{'style':'scatter'}) ]
            )
    plt.ylabel('')
    plt.title('Mean niche position')

    #plot_main(measures=measures,axes=axes,split_by=split_by,dictionary=dico)

    plot_groups(measures=measures,axes=axes,split_by=split_by)

    plt.show()



path='interp_ranked'

if 1:
    ranked(path=Path(path),mode=['dft','K'][0],
        rerun=0,measure=1,nsys=1,show=1)

