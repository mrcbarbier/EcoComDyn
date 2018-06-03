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
    'structure':'bunin',
    'mean':3.,
    'std':1.,
    'symmetry':-1.,
    'dynamics':'nlv',
    },

}

dft_dyn={
    'nlv':{'type':'simplelv',
        'variables':[('n','com'),],
        },
    }

def test_net():
    kwargs={}
    prm=deepcopy(dft_prm)
    prm['community']['connectivity']=.2
    prm['community']['mode']='network'
    #prm['community']['structure']='barabasi'
    prm['n']['shape']=(100,)
    xs=np.linspace(-1,1,5)

    if 1:
        res=[]
        for x in xs:
            prm['community']['clustering']=x
            model=Model(parameters=prm,dynamics=dft_dyn)
            if 0:
                g=Graph(model.data['community'].matrix.T)
                g.plot(hold=1)
            mat=model.data['community'].matrix

            tri= np.sum( np.logical_and(mat!=0,np.dot(mat,mat)!=0 ) )
            deg=np.mean(np.sum(mat!=0,axis=1).astype('float'))
            res.append([ tri*1./deg,deg])

        res=np.array(res).T
        plt.subplot(121)
        plot(res[1],xs=xs,title='Degree',hold=1)
        plt.subplot(122)
        plot(res[0],xs=xs,title='Triangles')

    if 0:
        res=[]
        for x in xs:
            prm['community']['assortativity']=x
            model=Model(parameters=prm,dynamics=dft_dyn)
            if 0:
                g=Graph(model.data['community'].matrix.T)
                g.plot(hold=1)
            mat=model.data['community'].matrix
            degree=np.sum(mat!=0,axis=1).astype('float')

            links=np.where(mat!=0)
            corr= np.mean(degree[links[0]]* degree[links[1]] ) - np.mean(degree)**2
            corr /= np.mean(degree**2) - np.mean(degree)**2
            deg=np.mean(np.sum(mat!=0,axis=1).astype('float'))
            corr=np.corrcoef(degree[links[0]], degree[links[1]])[0,1]
            res.append([ corr,deg])
        res=np.array(res).T
        plt.subplot(121)
        plot(res[1],xs=xs,title='Degree',hold=1)
        plt.subplot(122)
        plot(res[0],xs=xs,title='Degree correlation')
def test_chi(path=Path('interp_chi'),rerun=0,measure=0,nsys=3,show=0):
    #Comparison between
    try:
        if not os.listdir(path):
            rerun=1
    except:
        rerun =1
        #If not already run before

    kwargs={}
    prm=deepcopy(dft_prm)
    kwargs['community_symmetry']=-1.
    kwargs['community_mean']=0.
    kwargs['n_shape']=(100,)

    sigmas=np.linspace(.3,1.5,50)
    #sigmas=np.linspace(.1,.8,60)
    sigmas=np.linspace(.0,1.,20)

    gammas=np.linspace(-1,1,3)
    axes=[
            ('community_std',sigmas),
            ('community_symmetry',gammas),
            ( 'sys',range(nsys) )
    ]


    if rerun:
        loop(path=path,axes=axes,check_key=0,print_msg=1,tmax=5000000,tsample=5000.,
            #keep='all',
            parameters=prm,dynamics=dft_dyn,
            **kwargs )

    if rerun or measure:
        rebuild_filelist(path)
        make_measures(path,use_measures=['usual','effective',measure_cavity,'matrix',
            'removal','testcavity'],death=10**-10)#,'trophic','assembly_corr'] )


    if not show:
        return

    axes=['n_shape','community_std']


    measures=open_measures(path)
    measures=measures[measures['avgN']!=0]


    plot_cavity(measures=measures, axes=axes[1:],
        split_by=['community_symmetry','n_shape'],)

    plot_data(measures=measures, axes=['community_std'],
        title='v by removal (dots) vs theory (line)',
        split_by=['community_symmetry','n_shape'],
        values=[('removal_v',{'style':'scatter'}),
            ('removal_v_std',{'style':'plot','linestyle':'None','marker':'x'}),
                ('v',{'style':'plot' } ),
             ],
        hold=1 )

    plot_data(measures=measures, axes=['community_std'],
        title='Difference made by one species: removal (dots) vs cavity (line)',
        split_by=['community_symmetry','n_shape'],
        values=[('removal_diff',{'style':'scatter'}),
                ('cavity_diff',{'style':'plot' } ),
             ],
        hold=1 )


    if 1:
        plot_data(measures=measures, axes=['community_std'],
            split_by=['community_symmetry','n_shape'],
            values=[('press',{'style':'plot','log':'y'}),
                    ('n_matrix_press',{'style':'scatter' } ),
                 ],
            hold=1)



    plot_data(measures=measures, axes=['community_std'],
        split_by=['community_symmetry','n_shape'],
        #title='Press and variability',
        values=[
            #('n_matrix_press',{'style':'scatter','log':'y','color':'b'}),
                ('n_matrix_var',{'style':'scatter'}),#,'color':'r' }),
                #('press',{'style':'plot','log':'y','color':'b'}),
                ('variability',{'style':'plot'}),#,'color':'r'}),
                #('n_pred_var',{'style':'plot','color':'g','linestyle':'None','marker':'^'}),
             ],
        hold=1 )


    if 0:
        plot_data(measures=measures, axes=['community_std'],
            split_by=['community_symmetry','n_shape'],
            values=[('n_matrix_eig',{'style':'scatter'}),
                    ('n_pool_symeig',{'style':'scatter','color':'r' } ),
                 ],
            hold=1 )

    plot_data(measures=measures, axes=['community_std'],
        split_by=['community_symmetry','n_shape'],
        values=[('phi',{'style':'plot'}), ('n_%alive',{'style':'scatter' } ) ] )


    plt.show()

def test_extinction(rerun=0,measure=0,nsys=10,show=0):
    #Comparison between competitive, trophic, disordered
    path=Path('test_extinction')
    try:
        if not os.listdir(path):
            rerun=1
    except:
        rerun =1
    if rerun:
        #If not already run before
        kwargs={}
        prm=deepcopy(dft_prm)



        loop(path=path,axes=axes,check_key=0,print_msg=1,tmax=500000,tsample=10000.,
            #keep='all',
            parameters=prm,dynamics=dft_dyn,
            **kwargs )

    if rerun or measure:
        rebuild_filelist(path)
        make_measures(path,use_measures=['usual','effective','abundance',
            'matrix',
            'removal',
            'testcavity',
            measure_cavity,
            ])

    if not show:
        return

    axes=['n_shape']

    measures=open_measures(path)
    measures=measures[measures['avgN']!=0]
    show_hist(measures=measures,axes=axes, values=['n_abundance'] )
    from cavity import cavity_distri
    xs=np.linspace(0,0.5,100)
    sigma,mu,gamma=[measures['n_couplings_'+z].mean() for z in ('sigma', 'mu','gamma')  ]
    dist=cavity_distri(S=2000,sigma=sigma,sigma_k=0,mu=mu,gamma=gamma)
    plot(xs,[dist(x) for x in xs],hold=1,linewidth=2 )
    plt.xlabel(r'Abundance $n_i$')


    plot_main(measures=measures,axes=axes)
    plot_cavity(measures=measures, axes=axes,)


def size_role(rerun=0,measure=0,nsys=3,show=0):
    #Comparison between
    path=Path('interp_size')
    try:
        if not os.listdir(path):
            rerun=1
    except:
        rerun =1
    if rerun:
        #If not already run before
        kwargs={}
        prm=deepcopy(dft_prm)

        #specs=[ ((16,),None,'nlv'),((5625,),'nlv',None) ][1:]
        specs=np.array((16,32,64,128,256,1024))[:4]
        mus=np.linspace(1,10,5)
        sigmas=np.linspace(.8,1.,3)
        gammas=np.linspace(-1,1,5)
        axes=[
                ( 'n_shape',[(n,) for n in specs] ),
                ('community_mean',mus),
                ('community_std',sigmas),
                ('community_symmetry',gammas),
                ( 'sys',range(nsys) )
        ]


        loop(path=path,axes=axes,check_key=0,print_msg=1,tmax=50000,tsample=1000.,
            #keep='all',
            parameters=prm,dynamics=dft_dyn,
            **kwargs )

    if rerun or measure:
        rebuild_filelist(path)
        make_measures(path,use_measures=['usual','effective',measure_cavity,
        'matrix',
        'testcavity',
        ])#,'trophic','assembly_corr'] )


    if not show:
        return

    axes=['n_shape','community_std']


    measures=open_measures(path)
    measures=measures[measures['avgN']!=0]
    plot_data(measures=measures, axes=['n_shape','community_std'],hold=1,
        #split_by=['community_symmetry'],
        values=[('phi',{'style':'wireframe'}), ('n_%alive',{'style':'scatter' } ) ] )
    #show_basic(measures=measures ,axes)
    #show_data(measures=measures,axes,  ['n_%alive','n_biomass','n_biomass_std','n_interactions_mu' ,'n_interactions_sigma'])


    plot_cavity(measures=measures,axes=['community_std'],
        split_by=['community_symmetry','n_shape'],)

    plt.show()
    plot_data(measures=measures, axes=['community_std'],hold=1,
        split_by=['community_symmetry','n_shape'],
        values=[('press',{'style':'plot'}), ('n_matrix_press',{'style':'scatter' } ) ] )


    plt.show()


def Rdistri_role(path=Path('interp_Rdistri'),rerun=0,measure=0,nsys=3,show=0):
    #Comparison between distributions for r_i

    try:
        if not os.listdir(path):
            rerun=1
    except:
        rerun =1

    kwargs={}
    prm=deepcopy(dft_prm)
    kwargs['n_shape']=(120,)

    sigmas=np.linspace(.0,1.,11)
    distris=['normal','exponential','diracs','uniform']
    #prm['growth']['values']
    axes=[
            ('growth_std',sigmas),
            ('growth_distribution',distris),
            ( 'sys',range(nsys) )
    ]


    if rerun:
        loop(path=path,axes=axes,check_key=0,print_msg=1,tmax=50000,tsample=1000.,
            #keep='all',
            parameters=prm,dynamics=dft_dyn,
            **kwargs )

    if rerun or measure:
        rebuild_filelist(path)
        make_measures(path,use_measures=['usual','effective','abundance',
            'matrix',
            measure_cavity,'removal','testcavity'])#,'trophic','assembly_corr'] )


    if not show:
        return

    axes=[
        'growth_std'
        #'n_capacity_std',
        ]
    split_by=['growth_distribution']


    measures=open_measures(path)
    measures=measures[measures['avgN']!=0]


    compare={'phi':'n_%alive','avgN':'n_biomass','stdN':'n_biomass_std','press':'n_matrix_press'}
    dico={'growth_distribution':r'$P(R_i)$',
                'growth_std':r'$\sigma_R$' ,
                }
    dico.update({j:cavity_dico.get(i,i) for i,j in compare.iteritems()})
    for i,j in compare.iteritems():
        plot_data(measures=measures, axes=axes,
            split_by=split_by,
            hold=1,
            dictionary=dico,
            values=[(i,{'style':'plot'}), (j,{'style':'scatter' } ) ] )


    plot_cavity(measures=measures, axes=axes,
        split_by=split_by,
    )

    plot_removal(measures=measures, axes=axes,
        split_by=split_by,log='')

    show_hist(measures=measures,axes=['growth_distribution'], values=['n_abundance'],
        log='y' ,bins=30)


    plt.show()


def test_connectivity(path=Path('test_connectivity'),rerun=0,measure=0,nsys=3,show=0):

    try:
        if not os.listdir(path):
            rerun=1
    except:
        rerun =1

    kwargs={}
    prm=deepcopy(dft_prm)
    prm['n']['shape']=(100,)
    prm['n']['death']=10.**-20
    prm['community']['symmetry']=-.5
    prm['community']['std']=0.4
    prm['community']['mean']=2.
    #prm['community']['structure']='random'

    sigmas=np.linspace(0,1.,51)[1:]
    axes=[
            ('community_mean',sigmas*10),
            #('community_std',sigmas*2),
            ('community_connectivity',[.1,.3,.5,.7,.9,1.][:1]),
            ( 'sys',range(nsys) )
    ]


    if rerun:
        loop(path=path,axes=axes,check_key=0,print_msg=1,tmax=500000,tsample=1000.,
            #keep='all',
            parameters=prm,dynamics=dft_dyn,
            **kwargs )

    if rerun or measure:
        rebuild_filelist(path)
        make_measures(path,use_measures=['usual','effective',#'abundance',
            'matrix',
            #measure_cavity,
            measure_funcresp_cavity,

            #'removal','testcavity'
            ])


    if not show:
        return

    axes=[ax[0] for ax in axes[:1]]
    split_by=['community_connectivity']
    log=''


    measures=open_measures(path)
    measures=measures[measures['avgN']!=0]

    if 1:
        test=measures.groupby('community_connectivity').min()
        idx=np.array(test.index)
        measures['n_couplings_max']=[np.round(-test.iloc[np.where(
            np.abs(idx-m)<0.01 )[0][0] ][ 'n_couplings_mean'],2)
             for m in measures['community_connectivity']]
        split_by=[ 'n_couplings_max']


    compare={'phi':'n_%alive','avgN':'n_biomass','stdN':'n_biomass_std','press':'n_matrix_press'}
    dico={'growth_distribution':r'$P(R_i)$',
                'growth_std':r'$\sigma_R$' ,
                 'n_couplings_max':'Coupling',
                 'community_mean':r'$\mu$'
                }
    dico.update({j:cavity_dico.get(i,i) for i,j in compare.iteritems()})
    for i,j in compare.iteritems():
        plot_data(measures=measures, axes=axes,
            split_by=split_by,
            hold=1,
            dictionary=dico,
            values=[(i,{'style':'plot'}), (j,{'style':'scatter' } ) ] )

    #plot_cavity(measures=measures, axes=axes,split_by=split_by,log=log)
    #plot_removal(measures=measures, axes=axes,split_by=split_by,log=log)

    plt.show()

def test_conn_simple(path=Path('test_conn_simple'),rerun=0,measure=0,nsys=3,show=0,
        **kwargs):
    S=200
    cs=np.linspace(0.1,1,15)[::-1]
    prm=deepcopy(dft_prm)
    if 0:
        prm['community']['mean']=-.05
        prm['community']['std']=0.05
        prm['community']['structure']='random'
    else:
        prm['community']['mean']=1.
        prm['community']['std']=0.4
        prm['community']['structure']='bunin'
        prm['community']['structure']='bunin'

    prm['n']['shape']=(S,)
    if 0:
        model=Model(parameters=prm,dynamics=dft_dyn)
        mat=model.data['community'].matrix.copy()
        simu=[]
        theo=[]
        xs=[]
        for conn in cs:
            mat2=mat.copy()
            for z in range(int( (1-conn)*S**2/2 ) ):
                i,j=np.random.randint(0,S,2)
                mat2[i,j]=0
                mat2[j,i]=0
                #print i,j
            #print mat2[0]
            model.data['community'].matrix=mat2
            model.evol(tmax=10000,print_msg=0)
            nalive=np.sum(model.results['n'][-1]>10**-10)*1. /S
            simu.append(nalive)

            trueconn=np.sum(mat2!=0)*1./(S*(S-1) )

            mu=-S*trueconn * np.mean( mat2[mat2!=0] )
            sigma=np.sqrt(S*trueconn * np.var( mat2[mat2!=0] ))
            measures={}
            measure_gen(model,measures,typ='effective')
            gamma=measures['n_couplings_gamma']
            print S,mu, sigma,gamma
            mu,sigma=measures['n_couplings_mu'],measures['n_couplings_sigma']
            theo.append(cavity_solve_props(S,mu,sigma,
                sigma_k=0,gamma=gamma)['phi'] )
            xs.append(measures['n_connectivity'])
        plot(xs,theo,hold=1)
        scatter(xs,simu)

    #rerun,measure=0,1
    axes=[
            ('community_connectivity',cs),
            ('sys',range(nsys))
    ]


    if rerun:
        loop(path=path,axes=axes,check_key=0,print_msg=0,tmax=50000,tsample=1000.,
            #keep='all',
            parameters=prm,dynamics=dft_dyn,
            **kwargs )

    if rerun or measure:
        rebuild_filelist(path)
        make_measures(path,use_measures=['usual','effective',#'abundance',
            'matrix',
            measure_cavity,
            ])

    #rerun,measure=0,1

    if not show:
        return

    axes=[ 'community_connectivity']
    split_by=None
    log=''


    measures=open_measures(path)
    measures=measures[measures['avgN']!=0]



    compare={'phi':'n_%alive','avgN':'n_biomass','stdN':'n_biomass_std','press':'n_matrix_press'}
    dico={'growth_distribution':r'$P(R_i)$',
                'growth_std':r'$\sigma_R$' ,
                }
    dico.update({j:cavity_dico.get(i,i) for i,j in compare.iteritems()})
    for i,j in compare.iteritems():
        plot_data(measures=measures, axes=axes,
            split_by=split_by,
            hold=1,
            dictionary=dico,
            values=[(i,{'style':'plot'}), (j,{'style':'scatter' } ) ] )

    plt.show()


def test_corrKA(path=Path('test_corrKA'),rerun=0,measure=0,show=1,nsys=1):


    try:
        if not os.listdir(path):
            rerun=1
    except:
        rerun =1
        #If not already run before

    kwargs={}
    from cavity_interp_run import dft_prm, dft_dyn
    prm=deepcopy(dft_prm)
    prm['n']['shape']=(100,)
    prm['community']['requires']=['growth']
    prm['community']['corr_KA']=0.
    prm['community']['symmetry']=-.5
    prm['growth']['std']=0.3
    prm['community']['mean']=3.


    rs=np.logspace(-4,-2,9)
    rs=np.concatenate((-rs[::-1],rs))
    sigmas=[ 0.3,0.5,0.9 ]
    axes=[
                ('community_corr_KA',rs),
                ('community_std',sigmas),
                ( 'sys',range(nsys) )
        ]
    Svar=('n_corr_KA','correlation',rs,rs )

    if rerun:
        loop(path=path,axes=axes,check_key=0,print_msg=1,tmax=500000,tsample=5000.,
            #keep='all',
            parameters=prm,dynamics=dft_dyn,
            **kwargs )

    if rerun or measure:
        rebuild_filelist(path)
        make_measures(path,use_measures=['usual',
            'effective',
            'matrix',
            measure_cavity,
            'testcavity'
            ],death=10**-10)#,'trophic','assembly_corr'] )
    axes=[ax[0] for ax in axes][:1]

    measures=open_measures(path)

    conv={'n_couplings_mu':'community_mean',
        'n_couplings_sigma':'community_std',
        'n_couplings_gamma':'community_symmetry',
        'n_capacity_mean':'growth_mean',
        'n_capacity_std':'growth_std',
        }

    theory={}
    dicos={}
    m=Model(parameters=prm)
    for sig in sigmas:

        theory[sig]=[]
        dicos[sig]=[]
        locprm=deepcopy(prm)
        for var in Svar[-1]:
            m.parameters['community']['corr_KA']=var
            m.parameters['community']['std']=sig
            dic={}
            dic.update(m.export_params())
            for i,j in conv.iteritems():
                dic[i]=dic[j]
            dic['n_corr_KA']=dic['community_corr_KA']
            measure_cavity(m,dic)
            #print 'yay', dic
            theory[sig].append(cavity_output(dic) )
            dicos[sig].append(dic)

    for trait in ('phi','avgN','stdN','press','variability'):
        plt.figure()
        plt.title(cavity_dico[trait])
        xlabel=Svar[1]
        for ir,sig in enumerate(sigmas):
            plt.xlim(xmin=min(Svar[-2]),xmax=max(Svar[-2]))
            plot(Svar[-2],[t[trait] for t in theory[sig]] ,hold=1,color=get_cmap(ir,len(sigmas)))

        plot_data(measures=measures, axes=axes,
            newfig=0,
            hold=1,
            split_by=['community_std'],
            dictionary={'r_std':r'$\sigma_\rho$'},
            values=[
                (trait,{'style':'plot','log':'','marker':'x','linestyle':'None'}),
                 (cavity_compare[trait],{'style':'scatter','log':'' ,'ylabel':'' ,'xlabel':xlabel} ) ]
            )


    plot_data(measures=measures, axes=axes,
        split_by=['community_std'],
        values=[
                ('n_corr_KA',{'style':'scatter' } ),
             ],
        hold=1)

    for i,j in conv.iteritems():
        plt.figure()
        #print i,j
        for ir,sig in enumerate(sigmas):
            plt.xlim(xmin=min(Svar[-2]),xmax=max(Svar[-2]))
            print sig
            plot(Svar[-2],[d[j] for d in dicos[sig]] ,hold=1,color=get_cmap(ir,len(sigmas)),ylabel=i)
            #print i,j,[d[i] for d in dicos[sig]],list(measures[j])[:5]

        plot_data(measures=measures, axes=axes,
            #split_by=['community_std'],
            newfig=0,
            values=[(j,{'style':'plot','marker':'x','linestyle':'None'}),
                    (i,{'style':'scatter' } ),
                 ],
            hold=1)

    plt.show()

def test_groupcalc(rerun=0,measure=1,nsys=20):
    #VERY BASIC TEST: DO I ACCOUTN FOR GAMMA PROPERLY WITH GROUPS?

    path=Path('test_groupcalc' )

    try:
        if not os.listdir(path) or (not measure and not 'files.csv' in os.listdir(path)) :
            rerun=1
    except:
        rerun =1

    prm={
            'n':{
        'type':'variable',
        'axes': ('com',),
        'death': 10.**-14,
        'shape':(80,),
        'mean':1.,
        'std':1,
        'sign':1,
        },

            'network':{
        'type':'matrix',
        'variables':[('n','com'),('n','com') ],
        'role':'edges',
        'save':False,
        'structure':'random',
        'connectivity':0.5,
        'mode':'network',
        },

            'growth':{
        'type':'matrix',
        'variables':['n'],
        'dynamics':'nlv',
        'role':'growth',
        'mean':1., 'std':0.4,
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

            'niche':{
        'type':'constant',
        'role':'niche',
        'axes': [('n','com')],
        'distribution':'diracs',
        'values':(0,1),
        'save':True,
        },


            'community':{
        'type':'matrix',
        'variables':[('n','com'),('n','com') ],
        'role':'interactions',
        'structure':'bunin',
        'edgefrom':'edges',
        'requires':['network'],
        'diagonal':0,
        'mean':.1,
        'stdrel':.5,
        'symmetry':-1.,
        'dynamics':'nlv',
        'diagonal':0,
        },

    }
    kwargs={}
    syst=range(nsys)

    axes=[('community_symmetry',np.linspace(-1,1,5)),
        ('sys',syst)]

    if rerun:
        #If not already run before

        loop(path=path,axes=axes,check_key=0,print_msg=1,tmax=500000.,tsample=10000,
            parameters=prm,dynamics=dft_dyn,
            **kwargs )

    if rerun or measure:
        rebuild_filelist(path)
        group_measures=make_measures(path,
                split_arrays={'contains':'n_groups','add':group_props_list+('n_abundance',) },
            use_measures=[
                measure_groups,
                measure_groups_cavity,

            ],
             )

        group_measures.to_csv(path+'group_measures.csv')
        make_measures(path,
            use_measures=[
                'usual',
                'effective',
                measure_cavity,
            ],
             )
    if not show:
        return

    measures=open_measures(path)
    group_measures=open_measures(path+'group_')

    axes=['community_symmetry']
    split_by=[]
    log=''
    dico={'community_symmetry':r'$\gamma$'}
    plot_main(measures=measures,axes=axes,split_by=split_by,log=log,dictionary=dico,subplots=1)
    plot_groups(measures=group_measures,axes=axes,split_by=split_by,log=log,dictionary=dico )
    plt.show()

if __name__=='__main__':

    if 1:
        #TEST THAT MY GROUP CALCULATION GETS CORRECT RESULTS WHEN DIVIDING IDENTICAL SPECIES
        test_groupcalc(rerun=0,measure=1,nsys=5)

    if 1:
        test_net()
    if 0:
        test_conn_simple(rerun=0,measure=0,show=1,nsys=10,path=Path('test_conn_simple2'))
    if 1:
        test_connectivity(rerun=0,measure=1,show=1,nsys=1,path=Path('test_connectivity'))

    #TEST OF ROLE OF CORRELATION B/W K AND A
    if 1:
        test_corrKA(rerun=0,measure=0,nsys=20)

    #STABILITY MEASURES
    if 1:
        test_chi(rerun=0,measure=1,show=1,nsys=1,path=Path('interp_chi_5'))

    #Role of distributons
    if 1:
        Rdistri_role(rerun=0,measure=0,show=1,nsys=10,path=Path('interp_Rdistri'))

    #DEPENDENCE IN S
    if 1:
        size_role(rerun=0,measure=0,show=1)#,nsys=10)