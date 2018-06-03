# -*- coding: utf-8 -*-

## McArthur resource competition
from show import *


dft_dyn={
    'competition':{'type':'macarthur',
        'variables':[('n','com'),('r','com')],
    },

}

dft_prm={
        'n':{
    'type':'variable',
    'axes': ('com',),
    'death': 10.**-8,
    'shape':(80,),
    'mean':1.,
    'std':1,
    'sign':1,
    },
        'r':{
    'type':'constant',
    'axes': ('com',),
    'shape':(80,),
    'mean':1.,
    'std':.2,
    'sign':1,
    },

        'consumption':{
    'type':'matrix',
    'variables':[('n','com'),('r','com') ],
    'mean':1./80, 'std':.3,
    'sign':1,
    'dynamics':'competition',
    'role':'consumption',
    },

        'mortality':{
    'type':'matrix',
    'variables':['n'],
    'dynamics':'competition',
    'role':'mortality',
    'mean':1., 'std':0,
    },
    }

def test_simple(profile=0):
    #TEST RESOURCE COMPETITION

    path=Path('tests/resource')
    tmax=200
    tsample=1.
    m=Model(parameters=dft_prm,dynamics=dft_dyn)
    if not profile:
        m.evol(tmax=tmax,print_msg=1,tsample=tsample)
    else:
        profiler.runcall(m.evol,tmax=tmax,tsample=tsample,print_msg=1)

        prolog('tmp.dat')

    m.save(path)
    m=Model.load(path,initialize=False)

    analyze_trajectory(m,hold=1,log='')

    if 0:
        m.evol(tmax=tmax,print_msg=1,tsample=tsample)
        analyze_trajectory(m,hold=1,log='y')
    plt.show()


def resource_loop(rerun=0,measure=0,param='equal'):

    path=Path('resource_loop')
    try:
        if not os.listdir(path):
            rerun=1
    except:
        rerun =1
    if rerun:
        #If not already run before

        nsys=10
        dft_prm['n']['shape']=(50,)
        sigmas=np.linspace(.2,2.2,11)
        specs=[(N,) for N in np.linspace(50,200,4).astype("int")]

        axes=[
                #('r_shape',specs),
                ('r_std', sigmas ),
                ('consumption_std', sigmas) ,
                ('sys',range(nsys))]

        attrs={'parameters':dft_prm,'dynamics':dft_dyn}


        loop(path=path,axes=axes,check_key=0,print_msg=1,tsample=100.,
            **attrs )
    if rerun or measure:
        make_measures(path)
        rebuild_filelist(path)

    #plot_function(path,axes=['community_std','community_symmetry'],
        #function= cavity_output ,
        #hold=1,
        #values='phi',
        #)

    plot_data(path,axes=[
            #'r_shape',
            'r_std',
            'consumption_std',
            ],
        values=['n_%alive'],
        #split_by=[]
        hold=1,newfig=1,
        )
    plt.show()

#test_simple(1)
resource_loop(0)