# -*- coding: utf-8 -*-

#MIXED COMMUNITY: COMPETITION AND MUTUALISM

from run import *



dft_prm={
        'n':{
    'type':'variable',
    'axes': ('com',),
    'death': 10.**-14,
    'shape':(80,),
    'mean':1.,
    'std':1,
    'sign':1,
    },
        'growth':{
    'type':'matrix',
    'variables':['n'],
    'requires':['niches'],
    'dynamics':'nlv',
    'role':'growth',
    'structure':'mixture',
    'mean':1., 'std':0.4,
    'sign':1,
    },
        'selfint':{
    'type':'matrix',
    'variables':['n'],
    'dynamics':'nlv',
    'role':'diagonal',
    'mean':1, 'std':0.,
    'save':1,
    'sign':1,
    },


        'niches':{
    'type':'constant',
    'role':'niche',
    'axes': [('n','com')],
    'distribution':'randint',
    'min':0,
    'max':1,
    'save':True,
    },

        'edges':{
    'type':'matrix',
    'variables':[('n','com'),('n','com') ],
    'mode':'network',
    'connectivity':.4,
    'role':'edges',
    'structure':'random',
    },



        'community':{
    'type':'matrix',
    'variables':[('n','com'),('n','com') ],
    'role':'interactions',
    'requires':['niches','edges'],
    'structure':'mixture',
    'diagonal':0,
    'mean':.01,
    'stdrel':.3,
    'symmetry':-1.,
    'dynamics':'nlv',
    'edgefrom':'edges',
    },

}

dft_dyn={
    'nlv':{'type':'simplelv',
        'variables':[('n','com'),],
        },
    }


def measure_mixture_simple(model,measure,**kwargs):
    measure_groups(model,measure,**kwargs)
    #print measure
    for z in ('sigma','mu'):
        measure['n_groups_{}'.format(z) ][:]=np.mean(measure['n_groups_{}'.format(z) ])
    res=group_cavity_output(measure,**kwargs)
    measure.update(res )


def measure_mixture(model,measure,**kwargs):
    measure_groups(model,measure,**kwargs)
    res=group_cavity_output(measure,**kwargs)
    measure.update(res )

def mixture(path=Path('mixture'),mode=['dft'][0],
        rerun=0,measure=0,nsys=3,show=0,imit=0,hold=0):

    #MIXTURE MODEL EXAMPLE FOR PAPER

    if mode != 'dft':
        path=Path(path.strip()+'_{}'.format(mode))
    try:
        if not os.listdir(path):
            rerun=1
    except:
        rerun =1

    kwargs={}
    prm=deepcopy(dft_prm)
    prm['growth']['sign']=0

    if mode=='dft':
        weights=[ (0.5,0.5),(0.7,0.3),(0.9,0.1)   ]
        axes=[
                ( 'community_ordered',np.linspace(0,1,11) ),
                ( 'niches_weights',weights[::-1] ),
                ( 'sys',range(nsys) )
            ]

    if rerun:
        loop(path=path,axes=axes,check_key=0,print_msg=1,tmax=500000,ftol=0.0002,tsample=10000.,
            parameters=prm,dynamics=dft_dyn,converge=1,
            **kwargs )

    if rerun or measure:
        #rebuild_filelist(path)
        make_measures(path,
        #split_arrays={'contains':'n_groups','add':group_props_list+('n_abundance',) },
            use_measures=
        [
            [
                measure_mixture,
                #measure_mixture_simple,
            'usual','matrix'
                ],
            ['usual',
            'effective',
            measure_cavity,
            'matrix',
            'degree',
            measure_cavity_compare,
            #'finalprm',
            #'removal','testcavity'
            ]
        ][0],
            #read_filter=('community_ordered','>=',.95),
            death=10**-10)

    if not show:
        return

    axes=[ax[0] for ax in axes if ax[0]!='sys']
    axes=[ax if not hasattr(ax,'__iter__') else ax[0] for ax in axes]
    split_by=axes[1:]
    if not split_by:
        split_by=None
    axes=axes[:1]

    measures=open_measures(path)
    if 1:
        measures_kstruct=open_measures(path+'kstruct_')


        measures=measures[measures['niches_weights']==(0.5,0.5)]
        measures_kstruct=measures_kstruct[measures_kstruct['niches_weights']==(0.5,0.5)]

    dico={'community_edgemean':r'$s$','n_shape':r'$S$','community_nestedness':'nestedness',
        'community_partition':'partition','community_edgestruct':'network','niches_weights':'Facultative',
        'rewiring':'Rewiring','community_ordered':'Ordering'
        }
    #plot_main(measures=aggregate_groups(measures),
        #axes=axes,split_by=split_by,dictionary=dico,aggregator='mean',log='')


    #plot_modelprm(measures=measures,axes=axes,split_by=split_by,dictionary=dico)
    if 0:
        plot_data(measures=measures_nostruct, axes=axes,
            split_by=split_by,
            hold=1,
            values=['cavity_compare']
            )


    if 1:
        #axes=['rewiring']
        plt.figure()
        panel=MutableInt(0)
        mlist=[aggregate_groups(measures_kstruct),aggregate_groups(measures)]
        for m in mlist:
            auto_subplot(plt,len(mlist),panel)
            plot_data(measures=m, axes=axes,
                #split_by=split_by,
                 dictionary=dico,
                newfig=0,
                hold=1,
                values=[
                    ('n_%alive',{'style':'plot','marker':'^','linestyle':'None'})  ],
                )
            plot_data(measures=m, axes=axes,
                #split_by=split_by,
                 dictionary=dico,
                hold=1,
                newfig=0,
                ylabel=ifelse(panel.val==1,'Fraction of survivors',''),
                values=[('phi',{'style':'plot',
                        #'linestyle':'--'
                        } ),])
    if not hold:
        plt.show()


def groups(path=Path('groups'),mode=['dft'][0],
        rerun=0,measure=0,nsys=3,show=0,imit=0):

    if mode != 'dft':
        path=Path(path.strip()+'_{}'.format(mode))
    try:
        if not os.listdir(path):
            rerun=1
    except:
        rerun =1

    kwargs={}
    prm=deepcopy(dft_prm)


    prm['niches']['min']=1 ##CRUCIAL SO THAT niches_max GIVES NUMBER OF GROUPS

    prm['rmean']={
    'type':'constant',
    'role':'rmean',
    'shape': ['niches_max'],
    'mean':1,
    'std':.5,
    }
    prm['rstdrel']={
    'type':'constant',
    'role':'rstdrel',
    'shape': ['niches_max'],
    'mean':0.25,
    'stdrel':0.25,
    }
    prm['Amean']={
    'type':'constant',
    'role':'Amean',
    'shape': ['niches_max','niches_max'],
    'mean':-0.1,
    'stdrel':.25,
    'diagonal':0,
    }
    prm['Astdrel']={
    'type':'constant',
    'role':'Astdrel',
    'shape': ['niches_max','niches_max'],
    'mean':.2,
    'std':.1,
    'sign':1,
    }
    prm['Asym']={
    'type':'constant',
    'role':'Asym',
    'shape': ['niches_max','niches_max'],
    'distribution':'uniform',
    'min':-1,
    'max':1,
    'symmetry':1,
    }

    prm['growth']['sign']=0
    prm['growth']['stdrel']=0
    prm['growth']['requires']+=['rmean','rstdrel']
    prm['community']['requires']+=['Amean','Astdrel','Asym']
    del prm['community']['symmetry']


    if mode=='dft':
        ng=range(1,6)
        axes=[
                ('niches_max', ng  ),
                ( 'sys',range(nsys) )
            ]
    elif mode=='order':
        order=np.linspace(0,1,11)
        prm['niches']['max']=2
        prm['Amean']['stdrel']=prm['Amean']['std']=0
        prm['Astdrel']['stdrel']=prm['Astdrel']['std']=0
        prm['Asym']['stdrel']=prm['Asym']['std']=0
        prm['rmean']['stdrel']=prm['rmean']['std']=0
        prm['rstdrel']['stdrel']=prm['rstdrel']['std']=0
        prm['n']['shape']=(120,)
        axes=[
                ('community_ordered', order  ),
                ( 'sys',range(nsys) )
            ]

    if rerun:
        loop(path=path,axes=axes,check_key=0,print_msg=1,tmax=500000,ftol=0.0002,tsample=10000.,
            parameters=prm,dynamics=dft_dyn,converge=1,
            **kwargs )




    def measgp(model,measure,*args,**kwargs):
        #kwargs.setdefault('groupby','niche')
        niches=model.find_data(role='niche').matrix
        names=sorted(set(niches))
        bioms=[np.sum(model.results['n'][-1][ niches==n ] ) for n in sorted(set(niches)) ]
        names={n:b  for n,b in zip(names,np.argsort(bioms)) }
        #from datatools import code_debugger
        # code_debugger()
        orderedniches=np.array([names[n] for n in niches]  )


        kwargs.setdefault('groupby', orderedniches  )
        measure_groups(model,measure,*args,**kwargs)
        measure_groups_cavity(model,measure,*args,**kwargs)
        from networkx.algorithms import community
        Aij=model.find_data(role='interactions').matrix.copy()

        ords=[]
        for i in range(5):
            parts = community.kernighan_lin_bisection(nx.Graph(Aij))
            pbioms=[np.sum(model.results['n'][-1][ list(part) ] ) for part in parts ]
            pnames=np.argsort(pbioms)
            orderedparts=np.zeros(niches.shape)
            for i,p in enumerate(parts):
                orderedparts[list(p)]=pnames[i]
            ords.append(orderedparts)

        sbm=SBM(Aij)
        measure['partition_compare']=np.mean([ max( np.mean(ord==orderedniches), np.mean(1-ord==orderedniches)) for ord in ords]  )
        measure['SBM_compare']=max(np.mean(sbm==orderedniches),np.mean(1-sbm==orderedniches))
        # code_debugger()

    if rerun or measure:
        rebuild_filelist(path)
        make_measures(path,
        split_arrays={'contains':'groups','!contains':['ranks','selfint'] },
            use_measures=
            ['usual',
            'effective',
            'matrix',
            'degree',
             measure_cavity,
             measgp,
             measure_cavity_compare,
            ],
            # read_filter=('sys','<=',1),
            death=10**-10)

    if not show:
        return

    axes=[ax[0] for ax in axes if ax[0]!='sys']
    axes=[ax if not hasattr(ax,'__iter__') else ax[0] for ax in axes]
    split_by=axes[1:]
    if not split_by:
        split_by=None
    axes=axes[:1]


    measures=open_measures(path)

    #tmp=measures['N1']/measures['n_groups_biomass*']
    #scatter(measures['Phi'],np.clip(tmp,0,10),log='' )


    #plot_main(measures=aggregate_groups(measures),axes=axes,split_by=split_by)

    if 0:
        for idx,r in measures.sort_values('cavity_compare',ascending=False).iterrows():
            print r['cavity_compare']
            print r["path"]
            print r['n_groups_S']
            print r['n_groups_biomass*'],r["N1"]
            print r['n_groups_biomass_std'],r["stdN"]
            print r['n_groups_mu']
            print r['n_groups_sigma']
            print r['n_groups_gamma']
            try:
                t=input()
                if t=='q':
                    break
            except:
                pass

    dico={'community_edgemean':r'$s$','n_shape':r'$S$','community_nestedness':'nestedness',
        'community_partition':'partition','community_edgestruct':'network','niches_weights':'Facultative',
        'rewiring':'Rewiring','community_ordered':'Ordering'
        }


    m=measures.groupby('community_ordered').mean()
    plot(m.index,np.abs(m['partition_compare' ]-0.5)*2,hold=1,log='',marker='o')
    plot(m.index,np.abs(m['SBM_compare' ]-0.5)*2,hold=1,log='',marker='o')
    # plt.ylim(ymin=0)
    plot(m.index,m['cavity_compare' ]*2,title='Bipartite graph',xlabel='Order parameter',hold=1,marker='o')
    plt.legend(['Partitioning','SBM','Dynamics'])
    plt.show()
    code_debugger()
    plot_groups(measures=measures,axes=axes,split_by=split_by)


    plt.show()

if __name__=='__main__':
    if 0:
        path=Path('interp_mixture')
        mixture(path=Path(path),
            mode=['dft'][0],
            rerun=0,measure=0,nsys=30,show=1)
    if 1:
        path=Path('interp_groups')
        groups(path=Path(path),
            mode=['dft','order'][1],
            rerun='rerun' in sys.argv,measure='measure' in sys.argv,nsys=10,show=1)
