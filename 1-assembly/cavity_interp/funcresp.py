# -*- coding: utf-8 -*-

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

        'threshold':{
    'type':'matrix',
    'variables':['n'],
    'dynamics':'nlv',
    'role':'threshold',
    'mean':1., 'std':0.,
    'sign':1,
    },
        'growth':{
    'type':'matrix',
    'variables':['n'],
    'dynamics':'nlv',
    'role':'growth',
    'mean':1., 'std':0.,
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


        'community':{
    'type':'matrix',
    'variables':[('n','com'),('n','com') ],
    'role':'interactions',
    'structure':'random',
    'diagonal':0,
    'mean':.5,
    'std':.1,
    'symmetry':-1.,
    'dynamics':'nlv',
    'diagonal':0,
    },

}

dft_dyn={
    'nlv':{'type':'simplefr',
        'variables':[('n','com'),],
        },
    }


def measure_funcresp(model,measure,**kwargs):

    A=model.find_data(role='interactions')[0]
    D=model.get_labeled_data('selfint')
    r=model.get_labeled_data('growth')
    K=r/D
    Nf=model.results['n'][-1]
    Nc=measure['threshold_mean']
    ref=np.abs(np.mean(A.matrix[np.abs(A.matrix)>10**-16]))
    flux=np.abs(np.dot(A.matrix,Nf))
    measure['n_saturation_fraction']=np.mean( (flux>Nc * ref).astype('float')  )  #np.mean(1./(1+ref*Nc/flux) )

    mat=A.matrix/(1+ np.abs(np.dot(A.matrix,Nf))/Nc/ref )

    #print 'note:',
    measure['n_saturation_flux']=np.abs(np.mean(np.sum(mat,axis=1)))
    return

def mutualism(path=Path('interp_mutualism'),mode=['dft'][0],
        rerun=0,measure=0,nsys=3,show=0,imit=0):
    #Mutualism with functional response

    if mode != 'dft':
        path=Path(path.strip()+'_{}'.format(mode))
    try:
        if not os.listdir(path):
            rerun=1
    except:
        rerun =1

    kwargs={}
    prm=deepcopy(dft_prm)
    prm['n']['shape']=(100,)
    prm['n']['death']= 10.**-15
    prm['n']['mean']=3.
    prm['growth']['mean']=1
    prm['growth']['std']=0.5
    prm['growth']['sign']=None
    prm['community']['structure']='random'
    prm['community']['mean']=mean=.02
    prm['community']['stdrel']=.2
    prm['community']['sign']=1
    #prm['community']['connectivity']=.2
    prm['threshold']['mean']=5.



    edge=np.linspace(2.,20.,3)
    specs=[(N,) for N in np.logspace(np.log10(20),np.log10(200),10).astype("int")]
    if mode=='dft':
        axes=[
            #('n_shape',specs),
                #( 'community_edgemean', edge ),
                #( 'community_modularity', np.linspace(0,1,3) ),
                #( 'community_nestedness', np.linspace(0,1,3) ),
                ( 'threshold_mean',np.logspace(0,2.5,7) ),
                ( 'sys',range(nsys) )
            ]

    elif mode=='nested':
        prm['community']['structure']='nested'
        axes=[
                ( 'threshold_mean',np.logspace(0,2.5,15) ),
                ( 'community_nestedness', np.linspace(0,1,3) ),
                ( 'sys',range(nsys) )
            ]

    elif mode=='bipartite':
        prm['community']['structure']='multipartite'
        part=np.array([0.,.5,.75,1])
        axes=[
                ( 'threshold_mean',np.logspace(0,2.5,15) ),
                ( 'community_partition',part ),
                #( ('community_partition','community_mean'), zip(part,mean/(1-part/2.) ) ),
                ( 'sys',range(nsys) )
            ]

    if rerun:
        loop(path=path,axes=axes,check_key=0,print_msg=1,tmax=500000,ftol=0.0002,
            tsample=10000.,converge=1,
            parameters=prm,dynamics=dft_dyn,
            **kwargs )

    if rerun or measure:
        rebuild_filelist(path)
        make_measures(path,use_measures=['usual',
            'effective',
            measure_funcresp_cavity,
            'matrix',
            measure_funcresp,
            #'finalprm',
            #'removal','testcavity'
            ],death=10**-10)

    if not show:
        return

    axes=[ax[0] for ax in axes if ax[0]!='sys']
    split_by=axes[1:]
    if not split_by:
        split_by=None
    axes=axes[:1]

    measures=open_measures(path)
    if imit:
        pathimit=Path(str(path).strip('\/')+'_imit')
        try:
            measures_imit=open_measures(pathimit)
        except Exception as e:
            print e
            measures_imit=None


    dico={'community_edgemean':r'$s$','n_shape':r'$S$','community_nestedness':'nestedness',
        'community_partition':'partition',
        'threshold_mean':r'$N_c$'}
    plot_main(measures=measures,axes=axes,split_by=split_by,dictionary=dico,aggregator='mean',log='yx')

    #plot_modelprm(measures=measures,axes=axes,split_by=split_by,dictionary=dico)


    plot_data(measures=measures, axes=axes,
        split_by=split_by, hold=1,
        values=[
             ('saturation_fraction',{'style':'plot'}),
             ('n_saturation_fraction',{'style':'scatter'}),]
        )
    #plot_cavity(measures=measures, axes=axes,split_by=split_by)
    #plot_removal(measures=measures, axes=axes,split_by=split_by)

    plt.show()


def competition(path=Path('funcresp_competition'),mode=['dft'][0],
        rerun=0,measure=0,nsys=3,show=0,imit=0):
    #Competition with functional response

    if mode != 'dft':
        path=Path(path.strip()+'_{}'.format(mode))
    try:
        if not os.listdir(path):
            rerun=1
    except:
        rerun =1

    kwargs={}
    prm=deepcopy(dft_prm)
    prm['n']['shape']=(100,)
    prm['n']['death']= 10.**-15
    prm['growth']['mean']=1
    prm['growth']['std']=0.2
    prm['growth']['sign']=1
    prm['community']['symmetry']=.5
    #prm['community']['connectivity']=.2
    prm['threshold']['mean']=5.



    #edge=np.linspace(2.,20.,3)
    #specs=[(N,) for N in np.logspace(np.log10(20),np.log10(200),10).astype("int")]
    if mode=='dft':
        prm['community']['structure']='bunin'
        prm['community']['mean']=10.
        prm['community']['std']=.5
        axes=[
                #( 'community_edgemean', edge ),
                #( 'community_nestedness', np.linspace(0,1,3) ),
                ( 'threshold_mean',np.logspace(2,3,11) ),
                ( 'sys',range(nsys) )
            ]

    elif mode=='nested':
        prm['community']['structure']='nested'
        prm['community']['mean']=mean=-.05
        prm['community']['stdrel']=.3
        axes=[
                ( 'threshold_mean',np.logspace(0,2.5,15)[12:] ),
                ( 'community_nestedness', np.linspace(0,1,3)[::-1] ),
                ( 'sys',range(nsys) )
            ]
    elif mode=='bipartite':

        prm['community']['structure']='multipartite'
        prm['community']['mean']=mean=-.05
        prm['community']['stdrel']=.3

        part=np.array([0.,.5,.75,1])
        axes=[
                ( 'threshold_mean',np.logspace(1,3,11) ),
                ( 'community_partition', part ),
                #( ('community_partition','community_mean'), zip(part,mean/(1-part/2.) ) ),
                ( 'sys',range(nsys) )
            ]


    if rerun:
        loop(path=path,axes=axes,check_key=0,print_msg=1,tmax=500000,ftol=0.0002,tsample=10000.,
            parameters=prm,dynamics=dft_dyn,converge=1,
            **kwargs )

    if rerun or measure:
        rebuild_filelist(path)
        make_measures(path,
            read_filter=['sys','<=',0],
            use_measures=['usual',
            'effective',
            #measure_cavity,
            measure_funcresp_cavity,
            'matrix',
            measure_funcresp,
            #'finalprm',
            #'removal','testcavity'
            ],death=10**-10)

    if not show:
        return

    axes=[ax[0] for ax in axes if ax[0]!='sys']
    axes=[ax if not hasattr(ax,'__iter__') else ax[0] for ax in axes]
    split_by=axes[1:]
    if not split_by:
        split_by=None
    axes=axes[:1]

    measures=open_measures(path)

    #if mode=='nested':
        #measures=measures[measures['community_nestedness']>0]
    if imit:
        pathimit=Path(str(path).strip('\/')+'_imit')
        try:
            measures_imit=open_measures(pathimit)
        except Exception as e:
            print e
            measures_imit=None


    dico={'community_edgemean':r'$s$','n_shape':r'$S$','community_nestedness':'nestedness',
        'community_partition':'partition',
        'threshold_mean':r'$N_c$'}
    plot_main(measures=measures,axes=axes,split_by=split_by,dictionary=dico,aggregator='mean',log='x',
        subplots=True)

    plot_modelprm(measures=measures,axes=axes,split_by=split_by,dictionary=dico)

    #plt.show()

    plot_data(measures=measures, axes=axes,
        split_by=split_by, hold=1,
        values=[
             ('saturation_fraction',{'style':'plot','log':'x'}),
             ('n_saturation_fraction',{'style':'scatter'}),]
        )

    if 0:
        plot_data(measures=measures, axes=axes,
            split_by=split_by, hold=1,
            values=[
                 ('saturation_flux',{'style':'plot'}),
                 ('n_saturation_flux',{'style':'scatter'}),]
            )

    #plot_cavity(measures=measures, axes=axes,split_by=split_by)
    #plot_removal(measures=measures, axes=axes,split_by=split_by)

    plt.show()

if __name__=='__main__':
    if 0:
        path=Path('interp_mutualism')
        #path=Path('interp_mutualism_edge')
        mutualism(path=Path(path),
            mode=['dft','nested','bipartite'][0],
            rerun=0,measure=1,nsys=30,show=1,imit=False)

    if 1:
        path=Path('funcresp_comp')
        #path=Path('interp_mutualism_edge')
        competition(path=Path(path),
            mode=['dft','nested','bipartite'][0],
            rerun=0,measure=0,nsys=30,show=1,imit=False)
