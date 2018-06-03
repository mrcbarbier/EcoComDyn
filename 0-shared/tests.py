# -*- coding: utf-8 -*-

from models import DynamicModel,dftparams
from analysis import *
from loop import loop
from datatools import profiler,prolog

def test_noise():
    #TEST NOISE GENERATION AND PROPERTIES
    pass


def test_model(profile=1):
    #TEST MODEL RUN WITH GIVEN PARAMETERS, SAVING TO A FOLDER
    #LOADING AND PLOTTING TRAJECTORY

    path=Path('tests/test_model')
    tmax=5
    prm={}
    prm.update(dftparams)
    prm['growth_std']=.5
    prm['capacity_std']=.2

    m=DynamicModel(**prm)
    #c=m.data['community'].matrix
    #c=c[c!=0]
    #print np.mean(c),np.std(c)
    #hist(m.data['community'].matrix.ravel() )
    if not profile:
        m.evol(tmax=tmax,print_msg=1,tsample=.1,**prm)
    else:
        profiler.runcall(m.evol,tmax=tmax,tsample=.1,print_msg=1,**prm)

        prolog('tests/tmp.dat')

    m.save(path)
    #run_movie(path,destination=path+Path('Movie') )
    m=DynamicModel.load(path,initialize=False)

    analyze_trajectory(m,hold=1,log='')
    plt.figure()

    if 0:
        prm['growth_dynamics']='nlv'
        prm['capacity_dynamics']='nlv'
        prm['community_dynamics']='nlv'

    m.initialize()
    m.evol(reseed=m.variables.keys() ,tmax=tmax,tsample=.1,print_msg=1,**prm)
    analyze_trajectory(m,hold=0,log='')

def test_var(profile=0):
    #TEST MULTIPLE VARIABLES

    path=Path('tests/test_var')
    tmax=600
    tsample=2

    dyn={
    'blog':{'type':'logistic',
        'variables':[('b','com'),],
        },
    'nlin':{'type':'linear',
        'variables':[('n','com'),],
        },
    'bncom':{'type':'community',
        'variables':[('b','com'),('n','com')],
        },
    'nbcom':{'type':'community',
        'variables':[('n','com'),('b','com')],
        },

    }

    prm={
        'n':{
    'type':'variable',
    'axes': ('com',),
    'death': 10.**-7,
    'shape':(2,),
    'mean':.1,
    'std':.1,
    'sign':True,
    },
        'b':{
    'type':'variable',
    'axes': ('com',),
    'death': 10.**-7,
    'shape':(2,),
    'mean':1.,
    'std':.5,
    'sign':True,
    },

        'growth':{
    'type':'matrix',
    'variables':['b'],
    'dynamics':'blog',
    'role':'growth',
    'mean':1., 'std':0,
    },
        'capacity':{
    'type':'matrix',
    'variables':['b'],
    'dynamics':'blog',
    'role':'capacity',
    'mean':1., 'std':0,
    },
        'mortality':{
    'type':'matrix',
    'variables':['n'],
    'dynamics':'nlin',
    'role':'diagonal',
    'mean':-.1, 'std':0,
    },
        'bn':{
    'type':'matrix',
    'variables':[('b','com'),('n','com') ],
    'dynamics':'bncom',
    'role':'interactions',
    'mean':.5, 'std':.3,
    'sign':-1,
    },
        'epsilon':{
    'type':'matrix',
    'variables':[('n','com'),('b','com') ],
    'mean':.7, 'std':0.05,
    'sign':True,
    },
        'nb':{
    'type':'matrix',
    'variables':[('n','com'),('b','com') ],
    'dynamics':'nbcom',
    'role':'interactions',
    'formula':'-bn * epsilon'
    },

        }

    m=DynamicModel(parameters=prm,dynamics=dyn)
    print m.data['nb']
    print m.data['bn']

    if not profile:
        m.evol(tmax=tmax,print_msg=1,tsample=tsample)
    else:
        profiler.runcall(m.evol,tmax=tmax,tsample=tsample,print_msg=1)

        prolog('tmp.dat')

    m.save(path)

    analyze_trajectory(m,hold=1,log='')
    plt.show()


def test_loop(rerun=0):
    #TEST LOOPING OVER PARAMETERS, SAVING TO FOLDER HIERARCHY
    #LOADING AND PLOTTING AGGREGATE MEASURES

    path=Path('tests/test_loop')
    try:
        if not os.listdir(path):
            rerun=1
    except:
        rerun =1
    if rerun:
        #If not already run before
        prm={}
        prm.update(dftparams)

        axes=[
            ('community_std', [.1,1,10.] ),
            ('growth_std',[0,1.])
             ]

        loop(path=path,axes=axes,check_key=1,print_msg=1,tsample=.1,**prm )
    make_measures(path)
    rebuild_filelist(path)
    plot_data(path,axes=['community_std','growth_std'],values=['n_biomass'],style='3d',
        #split_by=[]
        )

#test_model()
test_var()
#test_loop()