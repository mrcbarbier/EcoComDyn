# -*- coding: utf-8 -*-
## McArthur resource competition
from interp import *



dft_dyn={
    'nlv':{'type':'simplelv',
        'variables':[('n','com'),],
        },
    }

dft_prm={
        'n':{
    'type':'variable',
    'axes': ('com',),
    'death': 10.**-10,
    'shape':(30,),
    'mean':1.,
    'std':1.,
    'sign':1,
    },
        'r':{
    'type':'constant',
    'axes': ('com',),
    'role':'resource',
    'shape':(100,),
    'mean':1.,
    'std':0.4,
    'sign':1,
    #'scaling':{'std':-.5},
    },

        'consumption':{
    'type':'matrix',
    'variables':[('n','com'),('r','com') ],
    'mean':1., 'std':1./3.,
    'sign':1,
    'distribution':'exponential',
    'role':'consumption',
    'structure':'resource',
    'requires':['r'],
    },

        'mortality':{
    'type':'matrix',
    'variables':['n'],
    'role':'mortality',
    'mean':0.1, 'std':0,
    },


        'growth':{
    'type':'matrix',
    'variables':['n'],
    'dynamics':'nlv',
    'role':'growth',
    'structure':'resource',
    'requires':['r','consumption','mortality'],
    },
        'selfint':{
    'type':'matrix',
    'variables':['n'],
    'dynamics':'nlv',
    'role':'diagonal',
    'structure':'resource',
    'requires':['consumption'],
    },


        'community':{
    'type':'matrix',
    'variables':[('n','com'),('n','com') ],
    'role':'interactions',
    'dynamics':'nlv',
    'structure':'resource',
    'requires':['consumption'],
    },

    }

def measure_resource_spec(model,measure,**kwargs):
    r=model.data['r'].matrix
    xis=model.data['consumption'].matrix
    ms=model.data['mortality'].matrix
    c=np.mean(np.sum((xis>0),axis=0 )  )
    R=float(r.shape[0])
    S=float(xis.shape[0])
    xi=np.mean(xis[xis>0])

    measure['resource_availability']= np.sum( (r - np.dot(
        model.results['n'].matrix[-1],xis ))**2)
    Aij=model.data['community'].matrix
    Di=model.data['selfint'].matrix
    alpha=(Aij.T/Di)
    np.fill_diagonal(alpha,np.nan)
    alpha=nonan(alpha)
    from datatools import hist

    zs=np.sum((xis!=0),axis=1)
    zij=np.dot(xis,xis.T)[1,3]/xi**2
    #print zs[1],zs[3],np.dot(xis,xis.T)[1,3]/xi**2, zs[1]*zs[3]*xi**2/R, Aij[1,3], Aij[3,1]
    z=np.mean(zs)
    sigz=np.std(zs)
    #print '\n wait',z**2*np.mean(1./zs**2),(1-sigz**2/z**2) ,'\n\n'

    #print '\n',np.mean(r),(1-np.mean(ms))*(z* np.mean(xis[xis!=0])), np.mean(xis[xis!=0])**2 *z
    r=r/(1-np.mean(ms))/(z* xi)
    #print 'NOW',np.mean(r),np.std(r),'\n'
    #r=r/( np.mean(xis)*np.mean(r)-np.mean(ms))/z


    aij=(zij.T/zs).T
    #print 'mean', z/R, np.mean(aij),np.mean(alpha)
    #print 'std',sigz/R , np.std(aij),measure['n_couplings_sigma']/np.sqrt(S)
    #print 'std2',np.mean( np.sqrt(zs/R*(1-zs/R)/zs.reshape(1,S) ) ),np.sqrt((1-z/R)/R)  #,np.std(alpha)/np.sqrt(S)
    sigrho=np.std(r)

    measure['measured_R']=R
    measure['measured_sigrho']=sigrho
    measure['measured_z']=z
    measure['measured_sigz']=sigz

    #print S,R,z,sigz,sigrho
    mu,sigma,gamma,sigR,sigD,sigK,meanK,corrKA=convert_resources_spec(S,R,z,sigz,sigrho)
    #print '\n =>> meanK', meanK, measure['n_capacity_mean']
    #print '\n =>> stdK', meanK, sigrho*z*np.sqrt(np.mean(1./zs**2)), measure['n_capacity_std'],'\n'


    measure['predicted_mu']=mu
    measure['predicted_sigma']=sigma
    measure['predicted_gamma']=gamma
    measure['predicted_sigR']=sigR
    measure['predicted_sigD']=sigD
    measure['predicted_sigK']=sigK
    measure['predicted_meanK']=meanK
    measure['predicted_corrKA']=corrKA


    #if sigK>0.01:
        #hist(model.data['growth'].matrix/Di,normed=1,hold=1,bins=20)
        #rge=np.linspace(meanK-5*sigK,meanK+5*sigK,100)
        #plot(rge,np.exp(-(rge-meanK)**2/(2*sigK**2) )/np.sqrt(2*np.pi*sigK**2))

def measure_resource(model,measure,**kwargs):
    r=model.data['r'].matrix
    xis=model.data['consumption'].matrix
    ms=model.data['mortality'].matrix
    c=np.mean(np.sum((xis>0),axis=0 )  )
    R=float(r.shape[0])
    S=float(xis.shape[0])


    r=r/(1-np.mean(ms))/(R* np.mean(xis))
    measure['resource_availability']= np.sum( (r - np.dot(
        model.results['n'].matrix[-1],xis ))**2)

    measure['resource_%availability']=np.sqrt(measure['resource_availability']/np.sum(r**2))

    xip=xis[xis>0]
    xi=np.mean(xip )
    sigxi=np.std(xip)
    X=1./(1.+sigxi**2/xi**2)
    xi3=np.mean(( (xip-xi)/sigxi ) **3)
    xi4=np.mean(( (xip-xi)/sigxi )**4)
    rho=np.mean(r)
    sigrho=np.std(r)
    sigm=np.std(ms)
    m=np.mean(ms)
    measure['measured_X']=X
    measure['measured_xi']=xi
    measure['measured_sigma_xi']=sigxi
    measure['measured_sigma_rho']=sigrho

    measure['measured_xi3']=xi3
    measure['measured_xi4']=xi4
    distri=measure['consumption_distribution']
    pxi3,pxi4=moments_distri(distri)

    measure['predicted_xi3']=pxi3
    measure['predicted_xi4']=pxi4

    sig=np.sqrt(S * (1./R - R*xi**4) )
    mu,sigma,gamma,sigR,sigD,sigK,meanK,corrKA=convert_resources_gen(S,R,X,sigrho,
        xi3,xi4,m)

    measure['predicted_mu']=mu
    measure['predicted_bare_sigma']=sig
    measure['predicted_sigma']=sigma
    measure['predicted_gamma']=gamma
    measure['predicted_sigR']=sigR
    #print R* xi4+ (R**2-R)*np.mean(xip**2)**2, np.mean(model.data['selfint'].matrix**2 )
    measure['predicted_sigD']=sigD
    measure['predicted_sigK']=sigK
    measure['predicted_meanK']=meanK
    measure['predicted_corrKA']=corrKA
    #measure['predicted_sigDapprox']=np.sqrt(R* 2 * (xi**2 + sigxi**2)**2)

    #measure.update(cavity_solve_props(S=S,mu=mu,sigma=sigma,sigma_k=sigR+sigD,gamma=gamma))


def test_simple(profile=0):
    #TEST RESOURCE COMPETITION

    path=Path('tests/resource')
    tmax=200
    tsample=4.
    m=Model(parameters=dft_prm,dynamics=dft_dyn)

    if not profile:
        m.evol(tmax=tmax,print_msg=1,tsample=tsample)
    else:
        profiler.runcall(m.evol,tmax=tmax,tsample=tsample,print_msg=1)

        prolog('tmp.dat')

    m.save(path)
    m=Model.load(path,initialize=False)

    analyze_trajectory(m,hold=1,log='y')
    print np.mean(m.data['selfint'].matrix),np.mean(m.data['growth'].matrix),np.mean(m.data['community'].matrix)

    if 0:
        m.evol(tmax=tmax,print_msg=1,tsample=tsample)
        analyze_trajectory(m,hold=1,log='y')
    plt.show()


def imitate_model(inpath,outpath,use_measures=None,rerun=1,reproduce_corrKA=1,
        reproduce_connectivity=1,measure_func=None, **kwargs):
    #Imitate any arborescence of model realizations using closest Bunin-like model
    import pandas as pd
    COLS=set([])
    for content in extract_data(inpath,**kwargs):
        idx,filedata,model,props=content

        measure={}
        measure.update(props)
        measure_gen(model,measure,typ='effective')
        measure_func(model,measure)
        oldpath=filedata['path']
        newpath=Path(oldpath.replace(inpath,outpath) )
        newpath.mkdir()
        tmax,tsample=measure['tmax'],measure['tsample']
        from cavity_interp_run import dft_prm, dft_dyn
        prm=deepcopy(dft_prm)
        dyn=deepcopy(dft_dyn)

        sigma_K=measure['predicted_sigK']
        if reproduce_connectivity:
            conn=measure['n_connectivity']
        else:
            conn=1
        if reproduce_corrKA:
            #Correction added in matrix
            corrKA=measure['predicted_corrKA']
        else:
            #Correction integrated into sigma_K
            corrKA=0
            if sigma_K>10**-10:
                effS=model.parameters['n']['shape'][0]*measure['n_connectivity']
                #print sigma_K,  effS*measure['n_corr_KA']*2
                if not sigma_K**2>effS*measure['predicted_corrKA']*2:
                    print sigma_K**2,effS*measure['predicted_corrKA']*2, measure['n_connectivity']
                sigma_K=np.sqrt(sigma_K**2 - effS*measure['predicted_corrKA']*2)
                assert sigma_K>0

        prefix='predicted_'
        prm['n']['shape']=model.parameters['n']['shape']
        prm['growth']['std']=sigma_K
        prm['community'].update(
            {
                'mean':measure[prefix+'mu'],
                'std':measure[prefix+'sigma'],
                'symmetry':measure[prefix+'gamma'],
                'requires':['growth'],
                'corr_KA':corrKA,#measure['n_corr_KA'],
                'connectivity':conn,
                }
            )
        dic={}
        #dic.update(filedata)

        dic.update(model.export_params([p for p in model.parameters if not p in prm] ))
        dic.update({
            'tsample':tsample,
            'tmax':tmax,
            'path':newpath,
            })
        dic.update({m:measure[m] for m in measure if prefix in m})

        #print prm['growth']
        #print prm['selfint']
        newmodel=Model(parameters=prm,dynamics=dyn)


        dic.update(newmodel.export_params())
        COLS=COLS.union(dic.keys())
        if not rerun:
            continue

        success=newmodel.evol(tmax=tmax/50.,tsample=tsample,print_msg=1,**kwargs)
        if not success:
            print "ERROR!", oldpath
            continue
        newmodel.save(newpath,overwrite=1)

        #test={}
        #measure_gen(newmodel,test,typ='effective')
        #print [(measure['predicted_'+z], test['n_couplings_'+z]) for z in
            #('mu','sigma','gamma')]

        table=[dic]
        pd.DataFrame(table).to_csv(newpath+'files.csv')

    rebuild_filelist(Path(outpath))
    make_measures(Path(outpath),use_measures=use_measures,
        keep_cols=list(COLS))#,'trophic','assembly_corr'] )


def resource_spec (path=Path('interp_resource_spec'),mode=['R','Z'][0],rerun=0,
        measure=0,nsys=3,show=0,imit=0):
    #Resource specialist
    try:
        if not os.listdir(path):
            rerun=1
    except:
        rerun =1
        #If not already run before

    kwargs={}
    prm=deepcopy(dft_prm)
    kwargs['consumption_stdrel']=0.
    #prm['consumption']['save']=0
    prm['edges']={
        'type':'matrix',
        'variables':[('n','com'),('r','com') ],
        'mean':500, 'stdrel':.2,
        'distribution':'normal',
        'structure':'network',
        }

    prm['consumption']['requires']+=['edges']
    prm['consumption']['edgefrom']='edges'
    kwargs['n_shape']=(100,)

    rs=np.logspace(2,3,11)*5
    sigmas=[0,0.25,.5]#[:1]
    if mode=='R':
        axes=[
                ('r_shape',[(r,) for r in rs] ),
                ('r_std',sigmas),
                ( 'sys',range(nsys) )
        ]
        zs=rs
        Svar=('R',r'$A$',zs,zs )
    elif mode in ('Z','Z2'):

        kwargs['consumption_edgedist']='normal'
        if mode=='Z':
            prm['r']['shape']=(400,)
        else:
            prm['r']['shape']=(1600,)
        csigmas=np.linspace(0.,.4,11)[1:]*prm['r']['shape'][0]
        axes=[
                ('consumption_edgemean',csigmas),
                ('r_std',sigmas),
                ( 'sys',range(nsys) )
        ]
        #zmean=prm['consumption']['connectivity']*prm['r']['shape'][0]
        zs=np.linspace(csigmas[0],csigmas[-1],100)
        Svar=('z',r'$z$',zs, zs)

    elif mode=='sigZ':
        kwargs['consumption_edgedist']='normal'
        prm['r']['shape']=(400,)
        csigmas=np.linspace(0.,.5,11)[1:]
        axes=[
                ('consumption_edgestdrel',csigmas),
                ('r_std',sigmas),
                ( 'sys',range(nsys) )
        ]
        #zmean=prm['consumption']['connectivity']*prm['r']['shape'][0]
        zs=np.linspace(csigmas[0],csigmas[-1],100)
        Svar=('sigma_z',r'$sigz$',zs, zs)

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
            'abundance',
            measure_cavity,
            measure_resource_spec,
            'removal',
            'testcavity'],update_prev=0)

    if not show:
        return

    axes=[ax[0] for ax in axes][:1]
    split_by=['r_std']


    measures=open_measures(path)
    measures=measures[measures['avgN']!=0]

    measures_imit=None
    if imit:
        pathimit=Path(str(path).strip('\/')+'_imit')#model')
        try:
            measures_imit=open_measures(pathimit)
        except Exception as e:
            #raise e
            print e

    compare={'phi':'n_%alive','avgN':'n_biomass','stdN':'n_biomass_std','press':'n_matrix_press',
        'variability':'n_matrix_var'}
    theory={}
    for rstd in sorted(set(measures['r_std'])):
        dic={'S':np.mean([m[0] for m in measures['n_shape']]),
                'sigma_rho':rstd,
                'distribution':prm['consumption']['distribution'],
            'mortality':prm['mortality']['mean']}

        if mode =='R':
            dic['sigma_z']=(measures['measured_sigz']).mean()
            dic['z']=measures['consumption_edgemean'].mean()
        elif mode in ('Z','Z2','sigZ'):
            dic['R']=measures['measured_R'].mean()
            if mode =='sigZ':
                dic['z']=measures['consumption_edgemean'].mean()

            else:
                dic['sigma_z']=(measures['measured_sigz']).mean()

        theory[rstd]=[]
        for var in Svar[-1]:
            dic[Svar[0]]=var
            if Svar[0]=='z':
                #try:
                dic['sigma_z']=measures['consumption_edgestdrel'][0]*var
                #print "YAY"
                #except:
                    #pass
            #print dic
            theory[rstd].append(model_resource_specialists(**dic) )


    if mode=='R':
        log='x'
    else:
        log=''

    dico={'r_std':r'$\sigma_\rho$'}
    plot_main(measures=measures,axes=axes,measures_imit=measures_imit,
        theory=theory,Svar=Svar,split_by=split_by,log=log,dictionary=dico)
    plot_data(measures=measures,axes=axes,values=['resource_%availability'],
        hold=1,split_by=split_by)

    plot_data(measures=measures, axes=axes,
        hold=1,
        split_by=split_by,
        values=[
             ('cavity_meansum',{'style':'scatter','log':log,} ) ,
              ('cavity_corrsum',{'style':'plot','marker':'x'} ) ,
             ] )
    plot_cavity(measures=measures, axes=axes,split_by=split_by,log=log)
    plot_removal(measures=measures, axes=axes,split_by=split_by,log=log)

    #show_hist(measures=measures,axes=[axes[0]], values=['n_abundance'],log='y' ,bins=30)
    #plt.show()

    plot_modelprm(measures=measures, axes=axes,
        split_by=split_by,log=log)
    plt.show()


def resource_gen(path=Path('interp_resource_gen'),mode=['R','X'][0],rerun=0,measure=0,nsys=3,show=0,imit=0,hold=0):
    #Resource generalist
    try:
        if not os.listdir(path):
            rerun=1
    except:
        rerun =1
        #If not already run before

    kwargs={}
    prm=deepcopy(dft_prm)
    kwargs['consumption_std']=.3
    kwargs['r_std']=.4
    kwargs['n_shape']=(100,)
    if 'normal' in path:
        prm['consumption']['distribution']='normal'

    sigmas=[0,0.25,.5]
    rs=np.logspace(2,4,10)#[10:]
    if mode=='R':
        axes=[
                ('r_shape',[(r,) for r in rs] ),
                ('r_std',sigmas),
                ( 'sys',range(nsys) )
        ]
        cs=np.linspace(rs[0],rs[-1],100)
        Svar=('R',r'$R$',cs,cs )
    elif mode=='X':
        csigmas=np.linspace(0,1,15)
        axes=[
                ('consumption_std',csigmas),
                ('r_std',sigmas),
                ( 'sys',range(nsys) )
        ]
        cmean=prm['consumption']['mean']
        cs=np.linspace(csigmas[1],csigmas[-1],100)
        Svar=('X','Consumption variance',cs,1./(1+cs**2/cmean**2))
    elif mode=='RX':
        axes=[
                ('r_shape',[(r,) for r in rs] ),
                ('consumption_std',[0.2,0.4,0.6]),
                ( 'sys',range(nsys) )
        ]
        cs=np.linspace(rs[0],rs[-1],100)
        Svar=('R',r'$R$',cs,cs )
        prm['consumption']['save']=False
        prm['r']['save']=False


    if rerun:
        loop(path=path,axes=axes,check_key=0,print_msg=1,tmax=500000,tsample=10000.,converge=1,
            #keep='all',
            parameters=prm,dynamics=dft_dyn,
            **kwargs )

    def measure_tmp(model,measure,**kwargs):
        measure_groups(model,measure,groupby=np.ones(model.results['n'][-1].shape),)
        measure_groups_cavity(model,measure,prefix='group_',**kwargs)
        for z in measure:
            if 'group' in z and hasattr(measure[z],'__iter__'):
                measure[z]=measure[z][0]

    if rerun or measure:
        rebuild_filelist(path)
        make_measures(path,use_measures=['usual',
            #'abundance',
            'effective',
            'matrix',
            measure_cavity,
            #'removal',
            #'testcavity',
            ifelse(mode!='RX',measure_resource,None), #If mode == RX, I do not save the r matrix (too big)
            measure_tmp,
            ],update_prev=0)#,'trophic','assembly_corr'] )


    if not show:
        return

    axes=[ax[0] for ax in axes]

    split_by=axes[1:-1]
    axes=axes[:1]


    measures=open_measures(path)
    measures=measures[measures['avgN']!=0]


    measures_imit=None
    if imit:
        pathimit=Path(str(path).strip('\/')+'_imit')#model')
        try:
            measures_imit=open_measures(pathimit)
        except Exception as e:
            #raise e
            print e


    if 0:
        plot_data(measures=measures, axes=['r_shape'],
            split_by=split_by,
            values=[('n_matrix_eig',{'style':'scatter'}),
                    ('n_pool_symeig',{'style':'scatter','color':'r' } ),
                 ],
            hold=1 )

    theory={}
    for split in sorted(set(measures[split_by[0] ])):
        dic={'S':np.mean([m[0] for m in measures['n_shape']]),
            'distribution':list(measures['consumption_distribution'])[0],
            'mortality':measures['mortality_mean'].mean() }
        if mode =='R':
            if 'measured_X' in measures:
                dic['X']=measures['measured_X'].mean()
            else:
                dic['X']=1./(1.+measures['consumption_std'].mean()**2/measures['consumption_mean'].mean()**2)
        elif mode =='X':
            dic['R']=np.mean([m[0] for m in measures['r_shape']])

        if split_by[0]=='r_std':
            dic['sigma_rho']=split
        else:
            dic['sigma_rho']=measures['r_std'].mean()
            if mode=='RX':
                dic['X']=1./(1.+split**2/measures['consumption_mean'].mean()**2)


        theory[split]=[]
        for var in Svar[-1]:
            dic[Svar[0]]=var
            theory[split].append(model_resource_generalists(**dic) )


    dico={'r_std':r'$\sigma_\rho$','consumption_mean':r'$\xi$',
        'consumption_std':r'$\sigma_\xi$','r_shape':r'$R$'}
    plot_main(measures=measures,axes=axes,measures_imit=measures_imit,
        #theory=theory,Svar=Svar,
        split_by=split_by,log='x',dictionary=dico,traits=['totN','phi','simpsoninv','variability'],
        show_direct_calc=1,subplots=1)
    # plt.show()

    #plot_data(measures=measures,axes=axes,values=['resource_%availability'], split_by=split_by)
    #plt.figure()

    if Svar[0]=='X':
        plot(Svar[-2],Svar[-1],hold=1)
        plot_data(measures=measures, axes=axes,
                hold=1,
                newfig=0,
                split_by=split_by,
                values=[ ('measured_X',{'style':'scatter' } ) ] )

    #plot_cavity(measures,axes=axes,split_by=split_by,log='x')

    #plot_removal(measures=measures, axes=axes,
        #split_by=split_by,log='x')

    if 0:
        plot_modelprm(measures=measures, axes=axes,
            split_by=split_by,log='x')
    if 0:
        plot_data(measures=measures, axes=axes,
            hold=1,
            split_by=split_by,
            values=[('measured_xi3',{'style':'scatter','log':'x'}), ('predicted_xi3',{'style':'plot' } ),
                ('measured_xi4',{'style':'scatter','log':'x','color':'r'}), ('predicted_xi4',{'style':'plot' ,'color':'r'} ) ] )

    #show_hist(measures=measures,axes=['r_std'], values=split_by,log='y' ,bins=30)
    if not hold:
        plt.show()

    #TMP!!!!!!!!!!!!
    gpref='group_'
    plot_main(measures=measures,axes=axes,measures_imit=measures_imit,compare={gpref+i:j for i,j in GROUPS_COMPARE.items()},
        #theory=theory,Svar=Svar,
        split_by=split_by,log='x',dictionary=dico,traits=[gpref+i for i in ['totN','Phi']],
        show_direct_calc=1,subplots=1)

def intersection(rerun=1,measure=0,nsys=100,show=0):
    #Comparison between competitive trophic
    path=Path('interp_intersection')
        #(model_resource_generalists, {'vars':('R','X'), 'S':6, 'R':np.logspace(.8,.81,5), 'X':np.linspace(.9909,.991,5), 'sigma_rho':.25 }  ),
        #(model_trophic_satur,  {'vars':('m0','e'),'sigma_k':0.05,'S':2000, 'm0':np.linspace(12.1,12.2,6),'e':np.linspace(0.01,.03,11),'s0':2, } ),

    kwargs={}
    prm=deepcopy(dft_prm)
    kwargs['r_shape']=(6,)
    kwargs['n_shape']=(6,)
    kwargs['r_std']=0.25
    cmean=prm['consumption']['mean']
    prm['consumption']['save']=False
    prm['r']['save']=False
    X=0.991
    kwargs['consumption_std'] = np.sqrt(cmean**2  *( 1/X -1) )

    axes=[
                ( 'sys',range(nsys) )
        ]

    if rerun:
        loop(path=path,axes=axes,check_key=0,print_msg=1,converge=1,tmax=500000,tsample=15000.,
            #keep='all',
            parameters=prm,dynamics=dft_dyn,
            **kwargs )

    if rerun or measure:
        rebuild_filelist(path)
        make_measures(path,use_measures=['usual',
            'effective',
            'abundance',
            measure_cavity,],update_prev=0)

    measures=open_measures(path)
    kwargs={}
    show_hist(measures=measures,axes=['n_shape'], values=['n_abundance'],log='y' ,bins=30)
    from cavity import cavity_distri
    xs=np.linspace(0,0.5,100)
    sigma,mu,gamma=[measures['n_couplings_'+z].mean() for z in ('sigma', 'mu','gamma')  ]
    dist=cavity_distri(S=2000,sigma=sigma,sigma_k=0,mu=mu,gamma=gamma)
    plot(xs,[dist(x) for x in xs],hold=1,linewidth=2 )
    plt.xlabel(r'Abundance $n_i$')

    plt.show()


def resource_axis(path=Path('interp_resource_axis'),mode='niche',rerun=0,measure=0,
        nsys=3,show=0,imit=0,tmax=500000):
    #Resource generalist
    try:
        if not os.listdir(path):
            rerun=1
    except:
        rerun =1

    kwargs={}

    prm=deepcopy(dft_prm)
    prm.update({
        'niches':{
    'type':'constant',
    'role':'niche',
    'axes': [('n','com')],
    'mean':.5,
    'std':.5,
    'distribution':'uniform',
    'save':False,
    },

        'widths':{
    'type':'constant',
    'role':'niche_width',
    'axes': [('n','com')],
    'mean':.01,
    'std':.005,
    'distribution':'normal',
    'save':False,
    },
        })
    prm['consumption']['requires']+=['niches','widths']

    prm['consumption']['save']=False

    kwargs['n_shape']=(100,)
    kwargs['r_shape']=(1000,)
    kwargs['r_std']=.3
    width=kwargs['widths_mean']=1./kwargs['n_shape'][0]

    sigmas=np.logspace(-3.5,-2,60)
    ns=[(int(i),) for i in np.linspace(20,200,10) ]
    if mode=='niche':
        kwargs['widths_std']=0
        axes=[
                 ('widths_mean', sigmas),
                ('widths_stdrel',[0.0,0.3,0.5]),
                #('consumption_std',[0.0,0.1,0.3]),
                ( 'sys',range(nsys) )
        ]
    elif mode=='S':
        axes=[
                 ('n_shape', ns),#*width),
                ('widths_stdrel',[0.0,0.1,0.3]),
                ( 'sys',range(nsys) )
        ]
    elif mode=='randomize':
        axes=[
                 ('widths_mean', sigmas),
                 ('consumption_randomize', [0,1] ),#*width),
                ( 'sys',range(nsys) )
        ]


    if rerun:
        loop(path=path,axes=axes,check_key=0,print_msg=1,tmax=tmax,converge=0,tsample=5000.,
            #keep='all',
            parameters=prm,dynamics=dft_dyn,
            **kwargs )

    if rerun or measure:
        rebuild_filelist(path)
        make_measures(path,use_measures=['usual',
            'abundance',
            'effective',
            'matrix',
            measure_cavity,
            #'removal',
            #'testcavity',
            #measure_resource
            ],update_prev=0)

    if not show:
        return

    measures_imit=None
    if imit:
        pathimit=Path(path.strip()+'_imit')#model')
        try:
            measures_imit=open_measures(pathimit)
        except Exception as e:
            #raise e
            print e

    axes=[ax[0] for ax in axes]
    axes,split_by=axes[:1],axes[1:2]
    log='xy'


    measures=open_measures(path)
    measures=measures[measures['avgN']!=0]

    if 0:
        #measures=measures[measures['widths_stdrel']==0.0]
        split_by=None


    dico={'widths_mean':r'$w$',}

    plot_main(measures=measures,axes=axes,measures_imit=measures_imit,
        split_by=split_by,log='x',aggregator='median',dictionary=dico)

    plot_modelprm(measures=measures, axes=axes,
        split_by=split_by,log='x')

    plot_data(measures=measures, axes=axes,
            hold=1,
            split_by=split_by,
            values=[('n_connectivity',{}) ] )


    plt.show()



def resource_axis_prm(path=Path('interp_resource_axis_prm'),mode='niche',rerun=0,measure=0,
        nsys=3,show=0,imit=0,tmax=500000):
    #Resource generalist
    try:
        if not os.listdir(path):
            rerun=1
    except:
        rerun =1

    kwargs={}

    prm=deepcopy(dft_prm)
    prm.update({
        'niches':{
    'type':'constant',
    'role':'niche',
    'axes': [('n','com')],
    'mean':.5,
    'std':.5,
    'distribution':'uniform',
    'save':False,
    },

        'widths':{
    'type':'constant',
    'role':'niche_width',
    'axes': [('n','com')],
    'mean':.01,
    'std':.005,
    'distribution':'normal',
    'save':False,
    },
        })
    prm['consumption']['requires']+=['niches','widths']

    prm['consumption']['save']=False

    kwargs['n_shape']=(100,)
    kwargs['r_shape']=(1000,)
    width=kwargs['widths_mean']=1./kwargs['n_shape'][0]
    ns=[(int(i),) for i in np.linspace(20,200,10) ]
    sigmas=np.linspace(0,.5,10)
    axes=[
            ('n_shape', ns),
             ('widths_mean', sigmas[1:]*width*5),
            ('widths_stdrel',sigmas[1:] ),
            ('r_std',sigmas ),
            ( 'sys',range(nsys) )
    ]
    if rerun:
        table=PrmLooper().loop(axes=axes,check_key=0,print_msg=1,tmax=tmax,converge=0,tsample=5000.,
                #keep='all',
                parameters=prm,dynamics=dft_dyn,
                **kwargs )
        import pandas as pd
        pd.DataFrame(table).to_csv(path+'measures.csv')

    measures=open_measures(path)
    measures=measures[measures['n_shape']==(200,)]
    axes=[a[0] for a in axes if a[0]!='sys']
    axes=axes[1:3]
    split_by=['r_std']

    plot_modelprm(measures=measures, axes=axes,
        split_by=split_by,log='x',aggregator='median')
    plt.show()
#-------------------------------------------------------------------------
#test_simple(1)

#-------------------------------------------------------------------------

#intersection()

if __name__ == '__main__':

    path='interp_resource_axis'
    if 0:
        from cavity_interp_run import imitate
        imitate(path,path+'_imit',rerun=1,
            use_measures=['usual','effective',measure_cavity,'matrix'])


    if 0:
        resource_axis(path=Path(path),mode='niche',
            rerun=0,measure=0,nsys=5,show=1,imit=1,tmax=500000)

    if 0:
        path='interp_resource_axis_prm'
        resource_axis_prm(path=Path(path),mode='all',
            rerun=0,measure=0,nsys=10,show=1,imit=0,tmax=0)

    if 0:
        path='interp_resource_axis_rdm'
        resource_axis(path=Path(path),mode='randomize',
            rerun=0,measure=0,nsys=5,show=1,imit=1,tmax=500000)

    #-------------------------------------------------------------------------



    mode='RX'
    if mode=='R':
        path='interp_resource_gen_normal'
    elif mode=='X':
        path='interp_resource_gen_normal_X'
    elif mode=='RX':
        path='interp_resource_gen_normal_RX_small'
    if 0:
        from cavity_interp_run import imitate
        imitate_model(path,path+'_imit',rerun=1,measure_func=measure_resource,
            use_measures=['usual','effective',measure_cavity,'matrix'])

    if 1:
        resource_gen(path=Path(path),
            mode=mode,
            rerun=0,measure=0,nsys=70,show=1,imit=0)


    #-------------------------------------------------------------------------

    mode =['R','Z','Z2'][0]

    if mode =='R':
        path='interp_resource_spec_normal'
    if mode=='Z':
        path='interp_resource_spec_normal_Z'
    if mode=='Z2':
        path='interp_resource_spec_normal_Z2'

    if 0:
        from cavity_interp_run import imitate
        imitate_model(path,path+'_imit',rerun=1,measure_func=measure_resource_spec,
            reproduce_connectivity=0,
            use_measures=['usual','effective',measure_cavity,'matrix'])

    if 1:
        resource_spec(path=Path(path),
            mode=mode,
            rerun=0,measure=0,nsys=3,show=1,imit=0)


    #-------------------------------------------------------------------------
