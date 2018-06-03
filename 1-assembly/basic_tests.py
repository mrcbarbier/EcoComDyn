# -*- coding: utf-8 -*-


from datatools import *
dft_prm={
        'n':{
    'type':'variable',
    'axes': ('com',),
    'death': 10.**-10,
    'shape':(100,),
    'mean':1.,
    'std':1,
    'sign':1,
    },

        'growth':{
    'type':'matrix',
    'variables':['n'],
    'dynamics':'nlv',
    'role':'growth',
    'mean':1., 'std':0.2,
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

        'fluxes':{
    'type':'matrix',
    'variables':[('n','com'),('n','com') ],
    'role':'fluxes',
    'save':False,
    'mean':0.031, 'std':0,
    'sign':1,
    },

        'community':{
    'type':'matrix',
    'variables':[('n','com'),('n','com') ],
    'role':'interactions',
    # 'structure':'bunin',
    'mean':-2.3,
    'std':.005,
    'symmetry':-.5,
    'dynamics':'nlv',
    'diagonal':0,
    },

}

dft_dyn={
    'nlv':{'type':'simplelv',
        'variables':[('n','com'),],
        },
    }


def test_Rpackage():
    #NOTHING DEPENDS ON TSAMPLE IN DETERMINISTIC DYNAMICS
    prm=deepcopy(dft_prm)
    kwargs={}
    kwargs['community_std']=0.6
    m=Model(
        parameters=prm,dynamics=dft_dyn,
        **kwargs )
    tmax=50000
    m.evol(print_msg=1,tmax=tmax,tsample=tmax/100,)
    measure={}
    measure.update(m.export_params() )
    measure_gen(m,measure,typ=['usual','effective','matrix'])

    measure_cavity(m,measure)

    path=Path('/media/Storage/Recherche/Moulis/notes/genericpatterns/package/examples/01/')


    #### PRINTING OUT
    m.data['community'].matrix=-m.data['community'].matrix
    for i in sorted(measure):
        print i, measure[i]
    for i,j in {'growth':'r','selfint':'D','community':'A'}.iteritems():
        np.savetxt(path+j+'.txt',m.data[i].matrix)
    np.savetxt(path+'B.txt',m.results['n'].matrix[-1])


def test_smallS(nsys=20000):

    def compute(mode, pop):
        def offdiag(mat):
            diag=np.eye(mat.shape[0])
            return mat[diag ==0]
        if mode in ('zeta','K'):
            pop=[p[1] for p in pop]
            list1=np.concatenate(pop)
        else:
            pop=[p[0] for p in pop]
            list1=np.concatenate([offdiag(p) for p in pop]  )
        if mode=='gamma':
            list2 = np.concatenate([offdiag(p.T) for p in pop])
            return np.corrcoef(list1, list2)[0,1]
        if mode=='sigma':
            return np.std(list1)*np.sqrt(pop[0].shape[0])
        if mode=='mu':
            return -np.mean(list1)*(pop[0].shape[0])
        if mode=='zeta':
            return np.std(list1)
        if mode=='K':
            return np.mean(list1)

    Ss=[2,3,4,8,16]
    table=[]
    models={}
    populations={}
    for S in Ss:
        populations[S]=[]
        models[S]=[]
        loctable=[]
        for sys in range(max(1,nsys/S)):
            print S,sys
            prm=deepcopy(dft_prm)
            kwargs={}
            kwargs['n_shape']=(S,)
            kwargs['community_mean']=.3
            kwargs['community_std']=.4
            kwargs['community_symmetry']=.9
            m=Model(
                parameters=prm,dynamics=dft_dyn,
                **kwargs )
            #print m.data['community'].matrix
            tmax=5000
            m.evol(print_msg=0,tmax=tmax,tsample=tmax/3,converge=1)
            measure={}
            measure.update(m.export_params() )
            meas2=deepcopy(measure)
            models[S].append(m.copy())

            measure_gen(m,measure,typ=['usual','effective','matrix'])
            measure['S']=S
            measure_cavity(m,measure)
            loctable.append(measure)
            populations[S].append((m.find_data(role='interactions').matrix,m.find_data(role='growth').matrix) )
        table+=loctable
        for m,measure in zip(models[S],loctable):
            # code_debugger()
            #THEORETICAL PREDICTIONS WITH IMPOSED PARAMETER VALUES
            prefix = 'n_couplings_'
            meas2[prefix + 'mu'] = [ kwargs['community_mean'], compute('mu',populations[S]) ][-1]
            meas2[prefix+'sigma'] =[kwargs['community_std'] , compute('sigma',populations[S]) ][-1]
            meas2[prefix + 'gamma']=[kwargs['community_symmetry'], compute('gamma',populations[S]) ][-1]
            meas2['n_capacity_std']=[meas2['growth_std'] , compute('zeta',populations[S])][-1]
            meas2['n_capacity_mean']=compute('K',populations[S])
            measure_cavity(m,meas2)
            measure.update( {'THEO_'+i:j for i,j in meas2.iteritems() if not 'n_' in i}  )
            # code_debugger()
    table=pd.DataFrame(table)

    df=table.groupby('S').mean()
    #table=[]
    #for S in Ss:
        #measure=df.loc[S].to_dict()
        #measure['n_shape']=(S,)
        #measure_cavity(models.pop(0),measure)
        #table.append(measure)

    #df=pd.DataFrame(table)

    comp={'phi':'n_%alive','avgN':'n_biomass','stdN':'n_biomass_std','stdN':'n_biomass_std',
        'n_couplings_gamma':'n_couplings_gamma','n_couplings_sigma':'n_couplings_sigma',
          'n_couplings_mu':'n_couplings_mu'}


    sd=table[['S']+list(set(comp.keys()+comp.values()))].groupby('S').std()
    dico={'n_couplings_gamma':'gamma','n_couplings_sigma':'sigma'}
    for x in comp:
        plt.figure()
        for z in ('gamma','sigma','mu'):
            if z in x:
                res=[compute(z,populations[S]) for S in Ss]
                plot(Ss,res, color='r',hold=1)
        if not 'n_' in x:
            plot(Ss,df['THEO_'+x],hold=1,color='r')
        from datatools import errorbar
        errorbar(Ss,df[x],yerr=sd[x],hold=1 ,title=dico.get(x,x) )
        scatter(Ss,df[comp[x]],hold=1 )
    plt.show()


def test_tsample():
    #NOTHING DEPENDS ON TSAMPLE IN DETERMINISTIC DYNAMICS
    prm=deepcopy(dft_prm)
    kwargs={}
    m=Model(
        parameters=prm,dynamics=dft_dyn,
        **kwargs )
    tmax=5000
    m.evol(print_msg=1,tmax=tmax,tsample=tmax/100,)
    traj=m.results['n']
    plot(traj.index,traj.matrix,hold=1,log='y')
    m.evol(print_msg=1,tmax=tmax,tsample=tmax/10,reseed=0)
    traj=m.results['n']
    plot(traj.index,traj.matrix,hold=1,linestyle='None',marker='o')


    plt.show()


def test_noise():
    #NOISE EFFECT IS THOROUGHLY INDEPENDENT OF tsample EVER SINCE I SWITCHED TO dop853
    prm=deepcopy(dft_prm)
    prm['nnoise']={
        'type':'noise',
        'variables':[ 'n' ] ,
        'amplitude':1,
        'sign':1,
        'role':'noise',
        'dynamics':'noise',
        'rank':'full',
        #'direction':'random',
        }
    prm['n']['shape']=(3,)
    dyn=deepcopy(dft_dyn)
    dyn['noise']={'type':'noise',
        'variables':[('n','com'),],}
    kwargs={}
    m=Model(
        parameters=prm,dynamics=dyn,
        **kwargs )
    tmax=20
    nstep=10
    m.evol(print_msg=1,tmax=tmax,tsample=float(tmax)/nstep/10.)
    from datatools import plt
    traj=m.results['n']
    noise=m.data['nnoise']
    plt.subplot(121)
    plot(traj.index,traj.matrix,hold=1,log='y')
    plt.subplot(122)
    plot(noise.index,noise.matrix[:,:2],hold=1,log='y')
    m.evol(print_msg=1,tmax=tmax,tsample=float(tmax)/nstep,reseed=0)
    traj=m.results['n']
    noise=m.data['nnoise']
    plt.subplot(121)
    plot(traj.index,traj.matrix,hold=1,linestyle='None',marker='o')
    plt.subplot(122)
    plot(noise.index,noise.matrix[:,:2],hold=1,log='y',linestyle='None',marker='o')


    plt.show()

def test_noise_amplitude(**kwargs):
    #CONCLUSIONS: V ~ 1/2r works but requires very long tmax if r is not too large
    prm = {
        'n': {
            'type': 'variable',
            'axes': ('com',),
            'death': -10000,
            'shape': (1,),
            'mean': 0.,
            'std': 0,
        },
        'diag': {
            'type': 'matrix',
            'variables': ['n'],
            'role': 'diagonal',
            'mean': -2,
            'std': 0,
            'dynamics': 'lin',
        },
        # 'community': {
        #     'type': 'matrix',
        #     'variables': [('n', 'com'), ('n', 'com')],
            # 'role': 'interactions',
        #     'mean': 5.,
        #     'std': .5,
        #     'symmetry': .5,
        #     'dynamics': 'lin',
        # },
        'nnoise':{
            'type': 'noise',
            'variables': ['n'],
            'amplitude': 1,
            'role': 'noise',
            'dynamics': 'noise',
            'rank': 1,
            'direction':0,
        }

    }
    dyn={
        'lin':{'type':'linear',
            'variables':[('n','com'),],
            },
        'noise':{'type': 'noise',
              'variables': [('n', 'com'), ], }
    }
    m=Model(
        parameters=prm,dynamics=dyn,
        **kwargs )

    ts=np.linspace(0.02,3,50)
    dt=0.1
    tmax=30.
    from datatools import plot
    import scipy.integrate as scint
    for mode in (1,):
        #MODE = 0 : Hand integration
        #MODE = 1 : model.evol with pregenerated noise
        res=[]
        for t in ts:
            print t
            m.data['diag'][:]=-1./t
            # print tmax
            nstep=tmax/dt
            if mode:
                m.evol(print_msg=1,tmax=tmax,tsample=dt,death=-10000,dftol=-10000)
                traj=m.results['n'].matrix
            else:
                m.data['nnoise'].generate(0, 2 * tmax,dt)
                noise=np.array([m.data['nnoise'].get(t).matrix for t in np.linspace(0,tmax,nstep)])
                # xs=m.results['n'].index
                # plot(traj,xs=xs,hold=1)
                xs=np.linspace(0,tmax,nstep)[1:]
                traj=[0]
                t0=0
                def dxx(t,x):
                    # print m.data['nnoise'].get(t)[0]
                    return m.data['diag'][0]*x + np.random.normal(0,1)/np.sqrt(dt/100)
                for t in xs[1:]:
                    x=traj[-1]
                    for tt in np.linspace(t0,t,100):
                        x+=dxx(tt,x)* (t-t0)/100.
                    traj.append(x)
                    t0=t
                # plot(traj,xs=xs)
                # print '  RATIO   ', np.var(traj),np.var(noise)
            res.append( [np.var(traj)*2*np.abs(m.data['diag'][0]) ])#,np.var(noise)])
        res=np.array(res)
        plot(res,xs=ts,log='y',hold=1)
    plt.show()

def test_var():
    from scipy.linalg import solve_continuous_lyapunov as solve_lyapunov, norm, inv,eig
    prm=deepcopy(dft_prm)

    S=100
    gamma=-1.
    sigmas=np.linspace(0.01,.5,50)
    for sigma in sigmas:
        DOM=10000
        while DOM>0:
            #m=Model(prm)
            if gamma==0:
                asymm=0
            else:
                asymm=(1- np.sqrt(1-gamma**2) )/gamma

            aij= np.random.normal(0,1.,(S,S))/np.sqrt(S)
            aij= (aij+asymm * aij.T) / np.sqrt(1+ asymm**2 )
            np.fill_diagonal(aij,0)
            print np.std(aij)*np.sqrt(S),np.mean(aij*aij.T)*S

            J=- sigma*aij
            np.fill_diagonal(J,-1)

            DOM=np.max(np.real( eig(J)[0]))

        N=np.abs(np.random.normal(1,0.2,S ))
        CIJ=solve_lyapunov(J,-np.diag(N) )
        VarC=np.trace(CIJ)/S*2
        AA=sigma**2*( (aij**2).T *  (N/np.add.outer(N,N)  ) ) .T
        AG=sigma**2*( (aij*aij.T ) *   (N/np.add.outer(N,N)  ) )

        DEN=1 - np.sum(AG,axis=1)
        RESULT=np.dot( inv(np.diag(DEN) -AA  ), .5 * np.ones(S)  ) #/Slive

        VarP=1/(1- (1+gamma)/2 * sigma**2)
        #print VarP, sigma,RESULT
        #plot(np.logspace(-1,3),np.logspace(-1,3),hold=1)
        #scatter( np.sum(DEN*np.diag(CIJ)) ,np.sum(.5 + np.dot(AA,np.diag(CIJ) ) )  ,hold=1,log='xy')
        scatter(VarP,VarC,hold=1)
        scatter(VarP,np.mean(RESULT)*2,hold=1,color='r',alpha=.6)
    plt.show()



def test_chi_jeff(sigmas,nsys=100,Slive=30):
    def test(prm,Slive=30):
        measure={}
        var='n'
        shape=(Slive,)
        #print "\n\nA", prm['community']['std']
        Alive=BaseNetwork.make(shape=shape+shape,**prm['community']).matrix
        #rlive = BaseNetwork.make(shape=shape,**dft_prm['growth'])
        #print np.std(Alive)
        #print "D"
        Diag=np.eye(Slive)*BaseNetwork.make(shape=shape,**prm['selfint']).matrix
        Nlive=np.ones(Slive)

        print "SIGMA",prm['community']['std'], 'ASTD',np.std(Alive)/np.mean(np.abs(Alive))
        B=Alive-Diag#np.eye(Slive)*rlive/Klive
        D=np.diag(Nlive)
        J=np.dot(D,B)
        #print J
        from scipy.linalg import solve_continuous_lyapunov as solve_lyapunov, norm, inv
        from scipy.sparse.linalg import eigs

        def variance(J,D):
            #Jeff's variance
            return np.trace(solve_lyapunov(J,-np.dot(D,D) ))/float( D.shape[0] )
        def worst_var(A):
            I=np.eye(A.shape[0])
            return norm(inv(np.kron(A,I) + np.kron(I,A)) ,2)*2
        def worst_press(A):
            return norm(inv(A),2)

        try:
            measure['{}_matrix_eig'.format(var)]=eig=eigs(J,k=1,sigma=0)[0][0].real
        except:
            return None
        if eig>=0:
            return None
        measure['{}_matrix_press'.format(var)]=np.sqrt(norm(inv(B),ord='fro')**2/Slive)
        measure['{}_matrix_var'.format(var)]=variance(J,D**.5)*2
        measure['{}_worst_press'.format(var)]=worst_press(J)
        measure['{}_worst_var'.format(var)]=worst_var(J)
        #print n,v, variance(-np.eye(10), np.eye(10) )
        return measure



    sigmas=np.linspace(.3,1.5,50)
    prm=deepcopy(dft_prm)
    kwargs['community_symmetry']=.5
    plt.figure()
    ms=[]
    for sigma in sigmas:

        prm['community']['std']=sigma*5.
        measure={}
        for sys in range(nsys):
            m=None
            while m is None:
                m=test(prm,Slive)
            for k in m:
                measure.setdefault(k,0)
                measure[k]+=m[k]/float(nsys)
        ms.append(measure)
    scatter(sigmas,[m['n_matrix_press'] for m in ms],hold=1,log='y' )
    scatter(sigmas,[m['n_matrix_var'] for m in ms],c='r' ,hold=1)
    scatter(sigmas,[m['n_worst_press'] for m in ms],hold=1,log='y' ,color='g')
    scatter(sigmas,[m['n_worst_var'] for m in ms],c='k',hold=1 )
    scatter(sigmas,[-1./m['n_matrix_eig'] for m in ms],c='y' )

def test_cavity_funcresp():
    from datatools import *
    from cavity_conv import *
    dic ={'S':100.,'mu':10.,'sigma':.3,'sigma_k':1.,'gamma':.05,'Nc':50.,'avgK':1.}
    #S=dic['S']
    locals().update(dic)
    Nc=dic['Nc']

    print cavity_solve_props(**dic)
    N1,N2,v,phi,f= funcresp_cavity_solve(**dic)
    print N1,N2,v,phi, f#v*phi/2.,
    avgN=N1*phi
    q=N2*phi/avgN**2
    h= (dic.get('avgK',1)/avgN - mu)/sigma
    vv=v*sigma*phi
    print q,vv,h,phi
    print N1*phi,np.sqrt(N2*phi-(N1*phi)**2), N1*phi*S, N2/(S * phi * N1**2 )
    u=(1-mu*1./S)

    res=[]
    fs=np.linspace(0,1,20)
    for f in fs:
        utilde=u-gamma*v*(1-f)*phi*sigma**2
        effvar=(1-f)**2 *phi*sigma**2 * N2+sigma_k**2
        mean=(avgK-mu * (phi*N1*(1-f) + f * Nc/S  )  )/utilde
        var=effvar/utilde**2
        res.append( (mean,var) )
    res=zip(*res)
    #plot(res[0],res[1],xs=fs)
    #return


    res=[]
    comp=[]
    Ncs=np.logspace(0,2,30)
    for Nc in Ncs:
        dic['Nc']=Nc
        rs=[]
        for trials in range(3):
            r=funcresp_cavity_solve(**dic)
            rs.append(r)
        rs=np.array(rs)
        mean=(avgK-mu *   Nc/S  )  /u
        if (rs[:,-1]>0.8).any() and (rs[:,-1]<0.2).any():
            print rs.T[-2:]
        res.append(np.mean(rs,axis=0) )

        N1,N2,v,phi,f=res[-1]
        dic2={}
        dic2.update(dic)
        dic2['sigma']*=max(0.001,(1-f))
        dic2['gamma']/=np.clip((1-f),np.abs(dic2['gamma']) ,1.)
        dic2['mu']*= (1-f) + f * Nc/S /N1/phi
        c=cavity_solve_props(**dic2)
        comp.append( (c['avgN'],c['q'],c['v']/dic['sigma']/c['phi'],c['phi'])  )

    res=np.array(zip(*res))
    comp=np.array(zip(*comp))

    plot(res[0]*res[3],comp[0],xs=Ncs,hold=1)
    plt.figure()
    plot(res[1]/res[0]**2/res[3],comp[1],xs=Ncs,hold=1)
    plt.figure()
    plot(res[2],comp[2],xs=Ncs,hold=1)
    plt.figure()
    plot(res[3],comp[3],xs=Ncs,hold=1)
    plt.show()

    titles=['N1','N2','v','phi','f']
    for r in res:
        plt.figure()
        plt.title(titles.pop(0) )
        plot(Ncs,r,hold=1)
    plt.show()


def test_lognorm():
    from datatools import plt,plot,scatter
    prm=deepcopy(dft_prm)
    prm['n']['shape']=(100,)
    kwargs={}
    m=Model(
        parameters=prm,dynamics=dft_dyn,
        **kwargs )

    Nstar=np.random.lognormal(0,1,m.data['growth'].matrix.shape )
    m.data['growth'].matrix=np.dot(-m.data['community'].matrix+m.data['selfint'].matrix, Nstar )

    #print m.data['selfint'].matrix
    #print m.data['growth'].matrix
    #print m.data['community'].matrix
    tmax=2000
    m.evol(print_msg=1,tmax=tmax,tsample=tmax/100,converge=1)

    res=m.results['n'][-1]
    alive=res>10**-10
    print len(alive)
    #print res[alive]
    scatter(Nstar,res,log='xy')
    scatter(m.data['growth'].matrix[alive],res[alive],log='xy')

def test_convergence():
    prm=deepcopy(dft_prm)
    prm['n']['shape']=(300,)
    prm['n']['death']=10**-14
    prm['growth']['std']=0#10**-14
    kwargs={}
    m=Model(
        parameters=prm,dynamics=dft_dyn,
        **kwargs )

    # code_debugger()
    tmax=5000.
    m.evol(print_msg=1,tmax=tmax,tsample=tmax/30.,converge='force')
    trial=0
    while np.sum(m.results['n'][-1]>10**-10) <2:
        print trial
        m.evol(print_msg=0,tmax=tmax,tsample=tmax/30.,converge='force')
        trial+=1
    plot(m.results['n'].index,m.results['n'].matrix,log='xy',hold=1)
    code_debugger()

    m.evol(print_msg=1,tmax=tmax,reseed=0,tsample=tmax/30.,converge=0)

    plot(m.results['n'].index,m.results['n'].matrix,log='xy',hold=1,marker='o',linestyle='None')
    plt.show()
    code_debugger()
#test_lognorm()

#test_tsample()
# test_noise()
#test_noise_amplitude()
#test_var()


#test_cavity_funcresp()

# test_smallS()

test_convergence()