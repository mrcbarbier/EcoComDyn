# -*- coding: utf-8 -*-

from cavity_conv import *
from cavity import group_cavity_solve,group_variability,group_chi2,group_meanfield,group_press

group_props_list=('Phi','N1','N2','avgN','stdN','totN_joint','Phi_joint','stdN_joint','simpson_joint',
                  'V','chi2','press','press_rel','variability',
        'corr_AN','Nvar','totN','simpson','prod','totNmf')


GROUPS_COMPARE={'Phi':'n_groups_%alive','N1':'n_groups_biomass*','medianN':'n_groups_biomass_median',
        'stdN':'n_groups_biomass_std','press':'n_groups_matrix_press','variability':'n_groups_matrix_var',
        'totN':'n_groups_biomass_tot','simpson':'n_groups_simpson','simpsoninv':'n_groups_simpsoninv',
        'prod':'n_groups_productivity_ratio' }


def getmean(x):
    if len(x) > 0:
        try:
            return np.mean(x)
        except FloatingPointError:
            print x
            code_debugger()
    return 0


def getstd(x):
    if len(x) > 1:
        try:
            return np.std(x)
        except FloatingPointError:
            print x
            code_debugger()
    return 0

def aggregate_groups(measure):
    '''Transform predictions (phi,avgN,stdN) on groups into predictions on whole community'''
    #measure=retupleify(measure)
    S=measure['n_groups_S']
    Phi=measure['Phi']
    new=measure.copy()
    new['phi']=[np.sum(p*s)/np.sum(s) for p,s in zip(Phi,S)]

    #print measure['N1']
    for m in ('N1','N2') :
        new[m+'m']=[ np.sum(p*s*x)/np.sum(s) for p,s,x in zip(Phi,S,measure[m]) ]
    for m in ('simpson','press','variability'):#+tuple(z for z in new if 'n_groups_' in z):
        new[m.replace('groups_','')]=[ np.sum(s*x)/np.sum(s) for s,x in zip(S,measure[m]) ]


    new['totN']=[np.sum(x) for x in measure['totN']]
    new['avgN']=new['N1m']
    new['stdN']=[np.sqrt(n2-n1**2) for n2,n1 in zip(new['N2m'],new["N1m"])]

    for m in ('N1','N2') :
        del new[m+'m']

    return new

def group_interactions(Smean,Aij,mwhere,return_all=0,remove_diag=1):
    '''Compute mu,gamma... per group, given interaction matrix Aij and mwhere=[ [idxs of group 1], [...]]'''
    from itertools import product as iprod



    interactions=[[ [Aij[i,j] for i,j in iprod(m1,m2) if ifelse(remove_diag,i!=j,True) ] for m2 in mwhere] for m1 in mwhere]
    Amean=np.array([[getmean(i) for i in j] for j in interactions])
    Astd=np.array([[getstd(i) for i in j] for j in interactions])

    Arowstd=np.array([[ getstd([ getmean(Aij[i,m2])  for i in m1 ]) for m2 in mwhere] for m1 in mwhere])

    intsym=[[ [Aij[i,j]*Aij[j,i] for i,j in iprod(m1,m2) if i!=j or im1!=im2] for im2,m2 in enumerate(mwhere)] for im1,m1 in enumerate(mwhere)]
    Asym=np.array([[getmean(i) for i in j] for j in intsym])

    conn=np.array([[getmean(np.array(i)!=0) for i in j] for j in interactions])
    mu=-Smean*Amean
    sigma=np.sqrt(Smean) *Astd
    sigrow=Smean *Arowstd
    # gamma=(Asym-Amean*Amean.T )/(Astd*Astd.T +10**-15 )

    def getcorr(m1,m2,same=0):
        if len(m2)>1 and len(m1)>1 :
            if same:
                m=m1
                offdiag=(np.eye(len(m)) == 0)
                mat=Aij[np.ix_(m, m)]
                if np.std(mat[offdiag]):
                    return np.corrcoef(mat[offdiag], mat.T[offdiag])[0, 1]
            else:
                if np.std(Aij[np.ix_(m1, m2)]) and np.std(Aij[np.ix_(m2, m1)] ):
                    return np.corrcoef(Aij[np.ix_(m1, m2)].ravel(), Aij[np.ix_(m2, m1)].T.ravel())[0, 1]
        return 0

    gamma =np.array( [[ getcorr(m1,m2) for im2,m2 in enumerate(mwhere)]
            for im1,m1 in enumerate(mwhere)])
    #
    # def getcorr(m):
    #     if len(m):
    #         return np.corrcoef(Aij[np.ix_(m, m)][np.eye(len(m)) == 0],
    #                 Aij[np.ix_(m, m)].T[np.eye(len(m)) == 0])[0, 1]
    #     return 0
    gamma[np.eye(gamma.shape[0]) >0] = [ getcorr(m,m,same=1)  for m in mwhere ]

    try:
        assert not np.isnan(sigma).any()
    except:
        print 'SIGMA NAN'
        code_debugger()
    gamma[np.isnan(gamma)]=0
    try:
        assert np.max(np.abs(gamma))<=1
    except:
        print '|GAMMA| >1'
        code_debugger()
    if return_all:
        return conn,mu,sigma,sigrow,gamma,Amean,Astd,Asym
    return conn,mu,sigma,sigrow,gamma


def measure_tscore_fields():
    return ['trophic_scores','trophic_scores_max','trophic_depths','trophic_topos']

def measure_tscore(model,measure,**kwargs):
    if not 'trophic_scores' in measure:
        Nf=model.results['n'][-1].copy()
        r=model.data['growth'].matrix.copy()
        Aij=model.find_data(role='interactions').matrix.copy()

        if kwargs.get('remove_dead',0):
            alive=(Nf>-measure['n_death']*1.1 )
        else:
            alive=(Nf>-1)
        idxs=np.ix_(np.where(alive)[0],np.where(alive)[0])
        Alive=Aij[idxs]
        Nlive=Nf[alive]
        rlive=r[alive]
        tscore=np.zeros(r.shape[0])
        if not kwargs.get('remove_dead',0):
            try:
                Nlive+=0.001*np.min(Nlive[Nlive>10**-8])
            except:
                Nlive+=0.001

        tscore[alive]=trophic_score(Alive,
            pops=Nlive,
            influx= rlive*Nlive )
        tdepth=np.zeros(r.shape[0])

        #Top predator if has higher score than any of its own predator
        # toppred=[ not [j for j in np.where(Alive[i]<0)[0] if tscore[j]>tscore[i] ] for i in range(len(alive))]
        #OR!! influx is simply self-control:  species with much stronger self-control is effectively top
        toppred=model.find_data(role='diagonal').matrix.copy()[alive]*Nlive
        tdepth[alive]=trophic_score(-Alive,
            pops=Nlive,
            influx= np.array(toppred,dtype='float') ) #rlive*Nlive )
        measure['trophic_scores']=tscore
        measure['trophic_scores_max']=np.max(tscore)
        measure['trophic_depths']=tdepth
        # code_debugger()
    if not 'trophic_topos' in measure:
        r=model.data['growth'].matrix.copy()
        Aij=model.find_data(role='interactions').matrix.copy()
        tscore=trophic_score(Aij,
            pops=np.ones(r.shape) ,
            influx= r,ntries=0)
        measure['trophic_topos']=tscore



def group_coordinates(model,measure,groupby='random',cluster=0,refval=None,**kwargs):
    S=model.data['growth'].matrix.shape[0]
    nbins=kwargs.get('nbins',4)
    if not isinstance(groupby,basestring):
        assert groupby.shape[0]==S
        return groupby
    D= model.data['selfint'].matrix.copy().reshape(S, 1)
    if np.max(D)==0:
        D[:]=1
    A = (model.find_data(role='interactions').matrix.copy().T /D).T
    Ks =  model.data['growth'].matrix.copy().reshape(S, 1)/ D

    if '_' in groupby:
        groupby,mode=groupby.split('_')
    else:
        mode = ''

    if groupby=='random':
        rank=np.random.random(S)*nbins
    elif groupby == 'N':
        rank=model.results['n'][-1].copy()-Ks.squeeze()
    elif groupby in ('similarity','siminv','simlive','v'):
        X=None

        if groupby=='similarity':
            X= np.concatenate([A,A.T,Ks ],axis=1)
        elif groupby=='siminv':
            Ai= np.linalg.inv(np.eye(S )-A )
            X=np.concatenate([Ai,Ai.T,Ks ], axis=1 )

        else:
            N = model.results['n'][-1].copy()
            alive = (N > model.parameters['n'].get('death', 10 ** -8))
            Slive = sum(alive)
            alive = np.where(alive)[0]
            Alive = A[np.ix_(alive, alive)]
            Nlive = N[alive]
            Dlive = D[alive].squeeze()

            if groupby=='simlive':
                X=np.zeros((S,2*S))
                X[np.ix_(alive,alive )] = Alive
                X[np.ix_(alive,[a+S for a in alive] )] = Alive.T


            if groupby=='v':
                from scipy.linalg import inv
                B = (Alive - np.eye(Slive))
                diag = np.diag(Nlive * Dlive)
                J = np.dot(diag, B)
                invB = inv(J)
                Press = -invB  - np.eye(Slive)/Nlive

                if 0:
                    Press= Press*Nlive
                X=np.zeros((S,2*S))
                X[np.ix_(alive,alive )] = Press
                X[np.ix_(alive,[a+S for a in alive] )] = Press.T





        if cluster:
            print("Computing embedding")
            from sklearn import manifold
            from sklearn.cluster import AgglomerativeClustering
            X_red = manifold.SpectralEmbedding(n_components=2).fit_transform(X)
            clustering = AgglomerativeClustering(linkage='ward', n_clusters=nbins)
            rank=clustering.fit_predict(X_red)

            if not refval is None:
                #Order clusters by e.g. their average trophic position (since labels are arbitrarily ordered initially)
                ordered=np.unique(rank)
                ordered= np.argsort(np.argsort([ np.median(refval[rank==r]) for r in ordered ]))
                rank=np.array([ordered[r] for r in rank  ]  )
        else:
            rank=X
    elif groupby=='toptrophic':
        rank=measure['trophic_topos'].copy()
    elif groupby == 'trophic':
        rank=measure['trophic_scores'].copy()
    elif groupby == 'trophicd':
        rank=np.array((measure['trophic_scores'],-measure['trophic_depths'])).T
    elif 'niche' in groupby and model.find_data(role='niche'):
        rank=model.find_data(role='niche').matrix
        # print rank
    elif groupby in model.data:
        rank=model.data[groupby].matrix
    elif groupby in kwargs:
        rank=kwargs[groupby]
    else:
        rank=np.array(range(S))
    return rank

def group_agglomerate(rank,nbins):
    '''AGGLOMERATIVE CLUSTERING'''
    nbins=to_iterable(nbins)

    if len(rank.shape) < 2:
        rank = np.atleast_2d(rank).T

    grps = list(np.unique(rank, axis=0))
    import scipy.spatial.distance as sdist
    dist = sdist.squareform(sdist.pdist(np.array(grps)))

    pcdist = [np.linalg.norm(rank - g, axis=1) for g in grps]
    clusterlabel = np.argmin(pcdist, axis=0)

    pair = None
    clusterlabels = []
    clusterpos = []
    while len(grps) > 1:
        nbin = len(grps)
        if nbin < np.min(nbins) and len(clusterpos):
            break
        d2 = dist.copy()
        np.fill_diagonal(d2, 10 ** 15)  # Avoid diagonals
        pair = np.unravel_index(np.argmin(d2), dist.shape)

        i, j = int(min(pair)), int(max(pair))
        n = np.array([np.sum(clusterlabel == label) for label in range(nbin)])
        denominator = 1. * n[i] + n[j] + n
        ai = (n[i] + n) / denominator
        aj = (n[j] + n) / denominator
        b = (-n) / denominator
        clusterlabel[clusterlabel == j] = i
        clusterlabel[clusterlabel > j] -= 1
        # prevgp={x:[x] for x in range(j) }
        # prevgp.update({x-1:[x] for x in range(j+1,nbin)})
        # prevgp[i]=[i,j]
        grps[i] = (n[i] * grps[i] + n[j] * grps[j]) / (n[i] + n[j])
        grps.pop(j)
        assert (len(np.unique(clusterlabel)) == len(grps))

        newdist = ai * dist[i, :] + aj * dist[j, :] + b * dist[i, j]
        newdist[i] = 0
        dist[i, :] = newdist
        dist[:, i] = newdist
        # print 'AFFECT',i
        others = [k for k in range(nbin) if k != j]
        dist = dist[np.ix_(others, others)]

        if nbin in nbins:
            clusterlabels.append(clusterlabel.copy())
            clusterpos.append(list(grps))

    # code_debugger()
    return clusterlabels, clusterpos


def group_groupby(model,measure,grps=None,**kwargs):
    groupby=kwargs.get('groupby',None)

    members=None
    S=model.data['growth'].matrix.shape[0]
    nbins=kwargs.get('nbins',4)


    measure_tscore(model,measure)
    if 'ranks' in kwargs:
        rank=kwargs['ranks']
    else:
        if hasattr(groupby,'__iter__') and len(groupby)==S:
            rank=groupby
        else:
            rank=group_coordinates(model,measure,refval= measure['trophic_topos'].copy(),cluster=1,**kwargs )

    # print "YAY 2",grps, rank
    if grps is None:
        if len(np.array(rank).shape)==1:
            grps=sorted( set(rank) )
        else:
            grps=list(np.unique(rank,axis=0))
    else:
        kwargs.setdefault('remove_empty',0)

    # if 'trophic' in groupby:
    #     kwargs.setdefault('nbins',np.ceil(np.max(rank)-np.min(rank)))


    rshape=rank.shape
    if len(rshape)==1:
        rshape=(rshape[0],1)
    if ('nbins' in kwargs and len(grps)>nbins):
        #Obtain groups by discretization
        if 0:
            grps,S=put_in_bins(grps,nbins=nbins,binpos='left' )
        else:
            mins=np.atleast_1d(np.min(rank,axis=0))
            maxs=np.atleast_1d(np.max(rank,axis=0))
            def lmean(arr):
                return (arr[1:]+arr[:-1])/2.
            ranges=[ lmean(np.linspace(mins[i], maxs[i],nbins+1  ))
                 for i in range(rshape[1]) ]
            grps =list(itertools.product(*ranges))
    if members is None:
        #tmp=np.array(list(grps)+[10**16])-(10**-6)*np.mean(np.abs(grps))
        #members=[ np.logical_and((tmp[i+1]>rank),(rank>=tmp[i])) for i in range(len(grps)) ]

        dists= [np.linalg.norm(rank.reshape(rshape)-np.array(g),axis=1) for g in grps ]
        closest=np.argmin(dists,axis=0)
        members=[(closest==c) for c in range(len(grps)) ]
        #print grps
        #print [np.sum(m) for m in members ], np.sum(members)
        #print [ rank[m] for m in members]


    galive=[i for i,m in enumerate(members) if np.sum(m)>0 or not kwargs.get('remove_empty',1) ]

    def squeeze(g):
        if not hasattr(g,'__iter__'):
            return g
        if len(g)==1:
            return g[0]
        return tuple(g)


    if np.min(rank)<0:
        handled=1
        while not handled:
            handled=1
            code_debugger()

    return rank,[squeeze(grps[g]) for g in galive],[members[g] for g in galive]


def make_chain_zeros(lab,tol=0):
    '''If tol>0, keep links between levels if their distance is smaller than the minimal distance + tol'''
    if hasattr(lab[0], '__iter__'):
        bupos, tdpos = zip(*lab)
    else:
        bupos, tdpos = lab, lab
    bupos, tdpos = np.array(bupos), np.array(tdpos)
    budist = np.add.outer(bupos, -bupos)
    tddist = np.add.outer(-tdpos, tdpos)
    tmp = budist.copy()
    tmp[tmp <= 0] = np.max(tmp) + 1
    bumin = np.min(tmp, axis=1)
    tmp = tddist.copy()
    tmp[tmp <= 0] = np.max(tmp) + 1
    tdmin = np.min(tmp, axis=1)
    buzeros = np.logical_or(budist > bumin.reshape(bumin.shape + (1,))+tol, budist < 0)
    tdzeros = np.logical_or(tddist > tdmin.reshape(tdmin.shape + (1,))+tol, tddist < 0)
    # code_debugger()
    return np.logical_and(buzeros,tdzeros)


def measure_groups(model,measure,prefix='n_groups_',**kwargs):
    Nf=model.results['n'][-1].copy()
    # print "YAY YAY",kwargs.get('grps',None)
    alive=(Nf>10**-5 *np.mean(Nf) )#measure['n_death']*1.1 )
    if not np.sum(alive):
        return

    leveled,grps,members=group_groupby(model,measure,**kwargs)
    groupby=kwargs.get('groupby',None)
    if groupby is 'trophic':
        tscore=leveled
    else:
        tscore=None

    if kwargs.get('remove_dead',0) :
        galive=[i for i,m in enumerate(members)  if alive[m].any() ]
        members=[members[g] for g in galive]
        grps=[grps[g] for g in galive]
    S=np.array([np.sum(i) for i in members])

    Smean=np.mean(S)
    s=S/Smean

    measure[prefix.strip('_') ]=[str(x) for x in grps]
    measure[prefix+'members']=mwhere=[np.where(m)[0] for m in members]
    measure[prefix+'ranks']=leveled

    measure[prefix+'S']=S
    r=model.data['growth'].matrix.copy()
    Aij=model.find_data(role='interactions').matrix.copy()
    if 'selfint' in model.data:
        D=model.data['selfint'].matrix.copy()
    else:
        D=np.ones(r.shape )
    if np.max(D)==0:
        D[:]=1
    Aij=(Aij.T/D).T
    r/=D

    # print np.array([[np.sum([Aij[i, j] for i, j in iprod(m1, m2)]) / len(m1) for m2 in mwhere] for m1 in mwhere])

    #K=(model.data['growth']).matrix/D
    #print members

    conn,mu,sigma,sigrow,gamma,Amean,Astd,Asym=group_interactions(Smean,Aij,mwhere,return_all=1)

    trophic=kwargs.get('trophic_chain',0)
    if hasattr(trophic,'__iter__') or trophic:
        if  hasattr(trophic,'__iter__'):
            #EXPLICIT COORDINATES OF GROUPS TO CONSTRUCT CHAIN
            zeros=make_chain_zeros(trophic)

        else:
            zeros=np.triu(np.ones(mu.shape),k=2 ).astype('bool') +np.tril(np.ones(mu.shape),k=-2 ).astype('bool')
        mu[zeros]=0
        sigma[zeros]=0
        gamma[zeros]=0
        Amean[zeros]=0
        Astd[zeros]=0
        Asym[zeros]=0

    sigma[sigma==0]=10**-15
    assert not np.isnan(sigma).any()

    measure[prefix+'connectivity']=conn
    measure[prefix+'mu']=mu
    measure[prefix+'Amean']=Amean
    measure[prefix+'sigma']=sigma
    measure[prefix+'sigrow']=sigrow
    measure[prefix+'gamma']=gamma
    measure[prefix+'selfint']=D

    measure[prefix+'capacity_mean']=capa=np.array([np.mean(r[m]) if sum(m) else 0 for m in members])
    measure[prefix+'capacity_std']=np.array([np.std(r[m]) if sum(m)>1 else 0 for m in members])


    measure[prefix+'selfint_std']=np.array([np.std(D[m]) if sum(m)>1 else 0  for m in members])

    try:
        assert not np.isnan(capa).any()
    except:
        print 'CAPACITY NAN'
        code_debugger()


    measure[prefix+'#alive']=nalive=np.array([np.sum( (Nf[m]>measure['n_death']*1.1 ) ) for m in members])

    Splus=np.maximum(S,1)
    measure[prefix+'%alive']=Phi=nalive.astype('float')/Splus
    nalivetmp=nalive.copy()
    nalivetmp[nalive==0]=1


    #measure['n_abundance' ]=[ tuple(Nf[m].ravel()) for m in members]
    measure[prefix+'biomass_tot']=biom=np.array([np.sum( Nf[m] )   if sum(m) else 0  for m in members])
    measure[prefix+'biomass']=biom/Splus
    measure[prefix+'biomass_median']=np.array([np.median( Nf[m] )   if sum(m) else 0  for m in members])
    measure[prefix+'biomass*']=biom/nalivetmp

    measure[prefix+'biomass2_tot']=biom2=np.array([np.sum( Nf[m]**2 )  if sum(m) else 0  for m in members])
    measure[prefix+'biomass2']=biom2/Splus
    measure[prefix+'biomass*2']=biom2/nalivetmp
    #measure[prefix+'biomass_var']=
    var=biom2/Splus-biom**2/Splus**2
    measure[prefix+'biomass_std']=np.sqrt(np.clip(var,0,None))
    #measure[prefix+'biomass*_var']=
    varlive=biom2/nalivetmp-biom**2/nalivetmp**2
    measure[prefix+'biomass*_std']=np.sqrt(np.clip(varlive,0,None))
    measure[prefix+'simpson']=biom2/np.maximum(10**-15,biom)**2

    allmembers=np.concatenate(mwhere)
    measure[prefix + '%alive_joint'] = np.sum(Nf>measure['n_death']*1.1 )*1. / np.sum(S)
    measure[prefix+'biomass_tot_joint']=np.sum(Nf[allmembers])
    measure[prefix+'biomass_std_joint']=np.std(Nf[allmembers])#np.sqrt(np.clip(np.sum(biom2)-np.sum(biom)**2 ,0,None))/S
    measure[prefix+'simpson_joint']=np.sum(biom2)/np.sum(biom)**2
    #print np.array([np.sum( Nf[m]**2 ) for m in members]),np.array([np.sum( Nf[m] ) for m in members])**2

    measure[prefix+'productivity']=prod=np.dot(biom,Amean)
    # measure['n_relbiom']=np.sum(biom)/np.sum(r[r>0])
    # measure['n_relprod']=np.sum(prod)/np.sum(r[r>0])

    # print np.sum(biom) - np.sum(Nf[allmembers])

def measure_group_trophic(model,measure,prefix='n_groups_',**kwargs):
    Nf=model.results['n'][-1].copy()
    alive = (Nf> measure['n_death']*1.1 )
    if not np.sum(alive):
        return
    if not prefix+'members' in measure:
        measure_group_stability(model,measure,prefix=prefix,**kwargs)
    r=model.data['growth'].matrix.copy()
    Aij=model.find_data(role='interactions').matrix.copy()
    if 'selfint' in model.data:
        D=model.data['selfint'].matrix.copy()
    else:
        D=np.ones(r.shape )
    if np.max(D)==0:
        D[:]=1
    Aij=(Aij.T/D).T
    r/=D

    mwhere=deepcopy(measure[prefix+'members'] )
    walive=list(np.where(alive)[0])
    livemwhere=[ [walive.index(i) for i in m if i in walive] for m in mwhere  ]
    S=np.array([len(i) for i in mwhere])
    Smean=np.mean(S)


    ### NETWORK PROPERTIES
    tscore=measure['trophic_scores']
    xs=[np.logical_and(Aij[i]>0, alive)  for i in range(r.shape[0])]
    omnivory=np.array([ np.std(tscore[ x])  if x.any() else 0 for x in xs ])
    measure[prefix+'omnivory']=np.mean(omnivory)
    measure[prefix+'omnivory_max']=np.max(omnivory)
    measure[prefix+'omnivory_std']=np.std(omnivory)
    tmp=Aij.copy()
    tmp[tmp<0]=0
    preycomp=np.sum(np.dot(tmp,tmp.T)>np.mean(tmp[tmp>0])**2,axis=1)
    measure[prefix+'preycomp']=np.mean(preycomp)
    measure[prefix+'preycomp_std']=np.std(preycomp)
    biom=measure[prefix+'biomass_tot']
    def nonan(seq,dft=(1,) ):
        seq=np.array(seq)
        s=seq[np.isnan(seq)==False]
        if not s.shape:
            return np.array(dft)
        seq[np.isnan(seq)]=dft[0]
        seq[np.isinf(seq)]=dft[0]
        return seq
    measure[prefix+'biomass_ratio']=np.mean(np.log10( nonan(biom[1:]/biom[:-1] )))


    # lamb=np.log10(np.abs(Aij*Aij.T))
    # kappa=np.log10(np.abs(Aij/Aij.T))
    lamb=(np.abs(Aij*Aij.T))
    aij=Aij.copy()
    aij[aij==0]=10**-10
    kappa=1./(np.abs(aij/aij.T))

    def clean(a):
        a[np.isnan(a)] = 0
        a[np.isinf(a)] = 0
        return a

    clean(lamb)
    clean(kappa)

    dissip =lamb.copy()
    # clean(dissip)
    trophic=kwargs.get('trophic_chain')


    for name, mat in [ ('dissip',dissip), ('topfeed',lamb), ('inverted',kappa) ]:

        gpint=group_interactions(Smean,mat,mwhere,return_all=1,remove_diag=1)

        conn= measure[prefix+'connectivity']

        amean, astd, asym=[x/np.maximum(10**-10,measure[prefix+'connectivity']) for x in gpint[-3:] ]
        if name=='topfeed':
            # amean+=np.log10(np.multiply.outer(Smean,Smean)*conn*conn.T )
            amean*=(np.multiply.outer(Smean,Smean)*conn*conn.T )
        for a in (astd,amean):
            clean(a)
        amean=np.triu(amean)
        astd=np.triu(astd)
        if hasattr(trophic, '__iter__') or trophic:
            if hasattr(trophic, '__iter__'):
                # EXPLICIT COORDINATES OF GROUPS TO CONSTRUCT CHAIN
                zeros = make_chain_zeros(trophic)
            else:
                zeros=(np.triu(np.ones(amean.shape),k=2 ) +np.tril(np.ones(amean.shape),k=-2 ) + np.eye(amean.shape[0])  ).astype('bool')
        else:
            zeros= np.eye(amean.shape[0]).astype('bool')
        amean[zeros]=0
        astd[zeros]=0

        if np.sum(amean!=0):
            measure[prefix+name]= np.mean(amean[amean!=0])
            measure[prefix+name+'_std']=np.std(amean[amean!=0])
        else:
            measure[prefix+name]= 0
            measure[prefix+name+'_std']=0

        measure[prefix + name + '_intrastd'] = np.mean(astd[astd != 0])

    ### PATTERNS
    ## Utility variables
    Vars=[measure[prefix+'Varrel_{}'.format(iv)].copy() for iv in range(len(mwhere) )]
    Ntot=measure[prefix+'biomass_tot'].copy()
    dNtot=Ntot[1:]-Ntot[:-1]
    liveVars=[Var   for l,Var in zip(livemwhere,Vars) if len(l)>0  ]
    Vmean=measure[prefix+'Vmean'].copy()
    Vu,Vd = np.triu(Vmean,k=1),np.tril(Vmean,k=0)
    maxVars= [ np.argmax(np.diag(Var)) for Var in liveVars ]
    altsign=(-1) ** np.arange(0,len(Ntot))[::-1]
    zerothresh=-10**-9  #For things that are supposed to be zero but not numerically


    #TOP-DOWN SIGNS
    sprefix=prefix+'topdown_'
    tmps=np.zeros(4)
    if len(Ntot)>1:
        ln=int(np.floor(len(altsign)/2.))
        start= len(altsign)%2

        tmps=[
            # Largest CV in second-to-last level no matter which level is perturbed
            [m == len(maxVars) - 2 for i,m in enumerate(maxVars)  ],

            #Stronger topdown than bottomup press
            (np.sum(np.abs(Vu), axis=0) / (10 ** -15 + np.sum(Vu < 0, axis=0)) > zerothresh + np.sum(
                np.abs(Vd), axis=0) /
                    (10 ** -15 + np.sum(Vd != 0, axis=0)))[1:],
        ]

        if len(liveVars)>3:
            tmps+=[
            #2-step variance cascade
                ( liveVars[-2][start::2,start::2][ np.multiply.outer(altsign,altsign)[:ln,:ln ] <0] <=-zerothresh),
            ]
        else:
            tmps+=[0]

        tmps+=[
            # Alternance with largest on top,
             (dNtot * altsign[1:] >= zerothresh )[altsign[1:]>0]
        ]

    # ms=[np.mean(tmp)>0.5 for tmp in tmps]
    ms=[np.mean(tmp) for tmp in tmps]
    measure[sprefix + 'CV'],measure[sprefix+'V'],measure[sprefix + 'antiCV'],measure[sprefix+'biomass']=ms
    measure[sprefix+'total']=np.mean(ms)
    if np.isnan(ms).any():
        code_debugger()

    #BOTTOM-UP SIGN
    sprefix=prefix+'bottomup_'

    tmps=np.zeros(3)
    if len(Ntot)>1:
        tmps=[
        # Largest CV in perturbed level
            [m==i for i,m in  enumerate(maxVars)],
            #Stronger bottom up than topdown press
            ( np.sum(np.abs(Vu),axis=0) / (10**-15+np.sum(Vu<0,axis=0)) +zerothresh < np.sum(np.abs(Vd),axis=0)/
                                (10**-15+np.sum(Vd!=0,axis=0)) ),

            #Variance upflow (ALWAYS THE CASE!!)
            # ([np.min(V[i:,i: ]) for i,V in enumerate(Vars) ]>= zerothresh),

            # Pyramid
        ]
        if len(Ntot)>2:
            tmps+=[ (dNtot* np.sign(dNtot[0]) >=zerothresh)[altsign[1:]* np.sign(dNtot[0])<0]
            ]
        else:
            tmps+=[0]

    ms=[np.mean(tmp) for tmp in tmps]
    if np.isnan(ms).any():
        code_debugger()
    measure[sprefix + 'CV'],measure[sprefix+'V'],measure[sprefix+'biomass']=ms  #measure[sprefix + 'antiCV'],
    measure[sprefix+'total']=np.mean(ms)



def measure_group_stability(model,measure,prefix='n_groups_',**kwargs):
    Nf=model.results['n'][-1].copy()

    # alive=(Nf>10**-5 *np.mean(Nf) )
    alive = (Nf> measure['n_death']*1.1 )

    if not np.sum(alive):
        return

    if not prefix+'members' in measure:
        measure_groups(model,measure,prefix=prefix,**kwargs)
    r=model.data['growth'].matrix.copy()
    Aij=model.find_data(role='interactions').matrix.copy()
    if 'selfint' in model.data:
        D=model.data['selfint'].matrix.copy()
    else:
        D=np.ones(r.shape )
    if np.max(D)==0:
        D[:]=1
    Aij=(Aij.T/D).T
    r/=D

    mwhere=deepcopy(measure[prefix+'members'])
    S=np.array([len(i) for i in mwhere])

    Smean=np.mean(S)

    conn,mu,sigma,sigrow,gamma,Amean,Astd,Asym=group_interactions(Smean,Aij,mwhere,return_all=1)


    from scipy.linalg import inv,solve_continuous_lyapunov
    Alive=Aij[np.ix_(alive,alive)]
    Slive=np.sum(alive)
    Nlive=Nf[alive]
    Dlive=D[alive]
    B=(Alive-np.eye(Slive))
    DN=np.diag(Nlive)
    diag=np.diag(Nlive*Dlive)

    if  0 and not model.find_data(role='threshold'):
        J=np.dot(diag,B)
    else:
        #FUNCTIONAL RESPONSE
        # t=model.find_data(role='threshold')
        # q=model.find_data(role='exponent')
        # if not len(q):
        #     q=1
        ref=model.get_dx(0, model.get_labeled_result(idx=-1) )['n'].matrix
        refx=model.get_labeled_result(idx=-1)
        dN=np.diag(0.01*Nf)
        def adding(d,v):
            d2=deepcopy(d)
            d2['n']= d['n']+ v
            return d2

        J= np.array(  [(model.get_dx(0, adding( refx,dN[i]) )['n'].matrix - ref ) /( np.ones(len(Nf))*dN[i,i])
                       for i in range(len(Nf))  ] ).T
        J[np.isnan(J)]=0
        J=J[np.ix_(alive,alive) ]

    # Jrel=np.dot(np.dot(np.diag(1./Nlive),B),np.diag(Nlive))
    # Jrel = np.dot(J,np.diag(1./Nlive))
    Jrel=np.dot(np.dot(np.diag(1./Nlive**2/Dlive),J),np.diag(Nlive))
    # code_debugger()
    invJ=inv(J)
    Press=-invJ#*Nlive

    PressRel=-inv(Jrel)

    # code_debugger()
    measure[prefix+'Pressrel']=PressRel
    measure[prefix+'Press']=Press


    walive=list(np.where(alive)[0])
    livemwhere=[ [walive.index(i) for i in m if i in walive] for m in mwhere  ]
    Pgamma, Pmean, Pstd, Psym=group_interactions(Smean,Press,livemwhere,return_all=1,remove_diag=1)[-4:]
    measure[prefix+'Pressmean']=Pmean
    measure[prefix+'Pressstd']=Pstd
    measure[prefix+'Presssym']=Pgamma
    Vgamma, Vmean, Vstd, Vsym=group_interactions(Smean,PressRel,livemwhere,return_all=1,remove_diag=1)[-4:]
    measure[prefix+'Vmean']=Vmean
    measure[prefix+'Vstd']=Vstd
    measure[prefix+'Vsym']=Vgamma
    Vmean_self=Vmean.copy()
    # np.fill_diagonal(Vmean_self,0)
    Vmean_self+=np.diag([getmean(np.diag(PressRel)[m ])/len(m) if len(m) else 0 for m in livemwhere ])  #V adjusted so that diagonal reflects self-control
    measure[prefix+'Vmean_self']=Vmean_self

    def make_vec(i,S):
        vec=np.zeros(S)
        vec[i]=1
        return vec
    Var= np.array([solve_continuous_lyapunov(J, -diag * make_vec(m,Slive))  for m in livemwhere  ])
    for iv, V in enumerate(Var):
        measure[prefix+'Var_{}'.format(iv)]=np.array([[np.sum(V[np.ix_(m1,m2)]) for m2 in livemwhere] for m1 in livemwhere])

    Var_rel= np.array([solve_continuous_lyapunov(Jrel, -diag * make_vec(m,Slive))  for m in livemwhere  ])
    for iv, V in enumerate(Var_rel):
        measure[prefix+'Varrel_{}'.format(iv)]=np.array([[np.sum(V[np.ix_(m1,m2)]) for m2 in livemwhere] for m1 in livemwhere])



    conn,mu,sigma,sigrow,gamma,Amean,Astd,Asym=group_interactions(Smean,Aij,mwhere,return_all=1)

    def repl(A):
        slive=[len(i) for i in livemwhere]
        return np.concatenate( [ np.ones(si)*a for a,si in zip(A,slive)    ] )

    feedback=1-np.array([ np.dot(repl(Amean[i]),np.dot(Press, repl(Amean[:,i]))  )
                                 for i in range(len(livemwhere) )  ])
    # for i in range(len(livemwhere)):
        # measure[prefix + 'invasion_feedback_{}'.format(i)]=np.atleast_2d([feedback[i]])

def measure_cascade(model,measure,prefix='n_groups_',**kwargs):
    print 'MEASURING TROPHIC CASCADES'

    r=model.data['growth'].matrix.copy()
    Aij=model.find_data(role='interactions').matrix.copy()
    if 'selfint' in model.data:
        D=model.data['selfint'].matrix.copy()
    else:
        D=np.ones(r.shape )
    if np.max(D)==0:
        D[:]=1
    Aij=(Aij.T/D).T
    r/=D

    Nf=model.results['n'][-1].copy()

    alive=(Nf>10**-5 *np.mean(Nf) )#measure['n_death']*1.1 )

    #First compute jacobian
    Alive=Aij[np.ix_(alive,alive)]
    Slive=np.sum(alive)
    Nlive=Nf[alive]
    Dlive=D[alive]
    B=(Alive-np.eye(Slive))
    # B2=((Alive.T/Dlive).T-np.eye(Slive))
    from sicpy.linalg import inv

    DN=np.diag(Nlive)
    diag=np.diag(Nlive*Dlive)
    J=np.dot(diag,B)
    invJ=inv(J)
    Press=-invJ
    PressRel=-inv(B)

    # print J2

    tscore=measure['trophic_scores'].copy()

    #Then get matrix of variances
    from scipy.linalg import solve_continuous_lyapunov as solve_lyapunov,norm,inv
    Cij= solve_lyapunov(J,-DN )

    var='n'
    measure['{}_matrix_press'.format(var)]=norm(invB,ord='fro')**2/Slive
    measure['{}_matrix_var'.format(var)]=2*np.trace(Cij)/Slive

    print 'lyap',Cij

    #then get average of jacobian for species with Delta T ~ 1 and Delta T ~ 2
    Var=variance_interactions(J,lift=1)
    Var/=np.abs(np.mean(np.diag(Var) ))


    srcs=[('tscore',tscore)]

    niche=None
    if 'niche' in model.data:
        niche=model.data['niche'].matrix.copy()
        srcs+=[('niche',niche) ]

    print 'Var',Var
    print 'Press',Press
    print tscore
    print niche
    NF=Nf.copy()
    for tgt in range(len(tscore)):
        print '\n Hit ',tgt
        test=np.zeros(tscore.shape)
        test[tgt]=1
        print 'PRESS',np.dot(Press,np.dot(DN,test) )
        if 1:
            # THIS WORKS!
            m2=model.copy(initialize=False)
            # m2.data['growth'].matrix+=np.dot(DN,test)*0.05
            m2.parameters['immigration']={
                    'type': 'matrix',
                    'variables': ['n'],
                    'dynamics': 'nlv',
                    'role': 'influx',
                    'matrix':0.05*Nf[tgt]*test
            }
            m2.initialize(labels=['immigration'])
            m2.make_dynamics()
            m2.evol(converge=1,tmax=1000000,tsample=.1,print_msg=0,reseed=False,extend=True)
            print 'EMPIRICAL PRESS',(m2.results['n'][-1]-NF)/0.05
        print 'VAR\n',np.diag(solve_lyapunov(J,-DN*test ))
        if 0:
            plt.subplot(131)
            plt.imshow((solve_lyapunov(J, -DN * test)))
            plt.colorbar()
            plt.title('Theory')

            m3=model.copy(initialize=False)
            if 1:
                m3.parameters['nnoise'] ={
                'type':'noise',
                'variables':[ 'n' ] ,
                'amplitude':Nf[tgt]*0.05,
                'role':'noise',
                'dynamics':'noise',
                'rank':1,
                'direction':tgt,
                }

                m3.dynamics['noise']= {'type': 'noise',
                            'variables': [('n', 'com'), ], }
                m3.initialize(labels=['nnoise'])
                m3.make_dynamics()

            m3.evol(converge=0,tmax=20.,tsample=0.1,print_msg=1,extend=True,reseed=False)
            # print m3.results['n'].matrix
            from datatools import plot
            plt.subplot(133)
            plot( [np.var(m3.results['n'].matrix[:i,tgt]/0.05) for i in range(1,len(m3.results['n'].index) )],hold=1 )#m3.results['n'].matrix,xs=m3.results['n'].index ,hold=1)

            print 'EMPIRICAL VAR\n', np.diag(np.cov(m3.results['n'].matrix.T/0.05))

            plt.subplot(132)
            plt.imshow(np.cov(m3.results['n'].matrix.T/0.05))#/Nlive.reshape(len(Nlive),1)))
            plt.title('Empirical')
            plt.colorbar()
            # plt.subplot(133)
            # plt.imshow((solve_lyapunov(J, -DN * test))/(np.cov(m3.results['n'].matrix.T/Nlive.reshape(len(Nlive),1))) )
            # plt.title('Ratio')
            # plt.colorbar()
            plt.show()

    #print np.add.outer(tscore[alive],-tscore[alive])
    for dref,s,v,mode in iproduct( (1,2),srcs, [('press',Press),('var',Var )], ('' ,'loc')  ) :
        lab,src=s
        lval,val=v
        dist=np.add.outer(-src[alive],src[alive])

        std=1./5.
        if mode=='loc':
            #Include only the levels right below
            cyc=100000.
        else:
            #Go down the whole chain
            cyc=2.


        res=(val* np.exp(-(np.abs(dist-dref)%cyc)**2/(2*std**2) ) ) [dist>0]  #/np.sqrt(2*np.pi*std**2 )
        #res[dist<=0]=0
        if not res.shape:
            res=[0]
        #print lab,lval,dref,mode, res
        measure['n_cascade{}_{}_{}_{}'.format(mode,lab,lval, dref) ]=np.sum(res)
        measure['n_cascade{}_{}_{}_{}_std'.format(mode,lab,lval, dref) ]=np.std(res)*np.sqrt(len(res))

        if 0 and lval=='press' and lab=='tscore' and dref==1 and mode=='':
            from datatools import scatter3d

            #plt.subplot(121)
            #scatter(niche,tscore,hold=1)
            #plt.subplot(122)

            # kwargs['TMP'](model,measure,hold=1)
            plt.figure()

            val[dist<=0.1]=0
            print 'Nf=', Nf
            pos=np.add.outer(src[alive],np.zeros(len(src[alive])) )
            scatter3d(pos[val!=0],dist[val!=0],val[val!=0],hold=1)
            plt.xlabel('pos')
            plt.ylabel('dist')
            plt.show()

        if dref==1:
            resinv=val * (np.abs(dist)%cyc)
            measure['n_cascdist{}_{}_{}_pos'.format(mode,lab,lval) ]=np.sum(resinv[val>0])/np.sum(val[val>0])
            measure['n_cascdist{}_{}_{}_neg'.format(mode,lab,lval) ]=np.sum(resinv[val<0])/np.sum(val[val<0])

        if dref==2:
            if measure['n_cascade{}_{}_{}_1'.format(mode,lab,lval) ]<0:
                measure['n_cascade{}_{}_{}'.format(mode,lab,lval) ]=measure['n_cascade{}_{}_{}_2'.format(mode,lab,lval) ]/measure['n_cascade_{}_{}_1'.format(lab,lval) ]
            else:
                measure['n_cascade{}_{}_{}'.format(mode,lab,lval) ]=0

            #if np.abs(measure['n_cascade{}_{}_{}_1'.format(mode,lab,lval) ])>1.:

                #print '\n',Press
                #print Nf
                #print J
                #raise



def group_cavity_calc_props(prefix='',S=None,mu=None,sigma=None,mean_k=None,
            sigma_k=None,gamma=None,N1=None,N2=None,V=None,Phi=None,**kwargs):
    import pandas as pd
    props=[prefix+p for p in group_props_list]

    s=S/np.mean(S)
    base=list(np.zeros((len(props),len(S)) ) )

    if not N1.shape:
        N1=N1.reshape(base[1].shape)


    chi2=group_chi2(s,sigma**2,gamma,Phi,V)

    var= group_variability(s,sigma**2,gamma,Phi)

    Press=group_press(s,mu,Phi,V)
    Press_rel= np.multiply.outer(1./N1,N1) * Press

    N1joint=np.mean(Phi*N1*S)/np.mean(S)
    varN=Phi* (N2- Phi* N1**2)
    N2joint = np.mean(Phi*N2*S)/np.mean(S) #np.sum(Phi*N2*s)
    varNjoint=N2joint - N1joint**2

    #testvarNjoint= np.var( [ N1[i] + np.random.normal(0,np.sqrt(varN[i]), S[i])  for i in range(len(N1))  ]  )
    # code_debugger()
    Phijoint=np.sum(Phi*S)/np.sum(S)
    stdNjoint= np.sqrt(varNjoint)#np.sqrt(N2joint-N1joint**2)
    simpsonjoint=N2joint/N1joint**2 / np.sum(S)
    ones=np.ones(len(S))
    lst=[Phi,N1,N2,Phi*N1,np.sqrt(np.clip(varN,0,None) ),N1joint*np.sum(S),Phijoint,stdNjoint,simpsonjoint,V,chi2,Press,Press_rel,var]
    base[:len(lst)]=lst


    base[-5]=N2-N1**2
    base[-4]=S*Phi*N1
    base[-3]=N2/N1**2/S/Phi
    base[-3][base[-3]>1.]=1
    base[-3][N1<10**-5]=0
    base[-1]=0#group_meanfield(s,mu,mean_k,sigma_k=sigma_k,Ninit=)[0]*S

    #DO NOT SQUEEZE!
    return pd.Series(dict(zip( props ,
                base)))


def group_cavity_solve_props(S=None,mu=None,sigma=None,
            sigma_k=None,gamma=None,tau=None,**kwargs):
    #print 'RECEIVING',S,mu,sigma,sigma_k,gamma,kwargs

    N1,N2,V,Phi,success= group_cavity_solve(S=S,mu=mu,sigma=sigma,sigma_k=sigma_k,
        gamma=gamma,tau=tau,**kwargs)
    series= group_cavity_calc_props(S=S,mu=mu,sigma=sigma,sigma_k=sigma_k,gamma=gamma,
        N1=N1,N2=N2,V=V,Phi=Phi,tau=tau,**kwargs)
    series['success']=success
    return series


def group_cavity_output(prm,**kwargs):

    prefix=kwargs.get('dataprefix','n_groups_')
    S=prm[prefix+'S']
    mu=prm[prefix+'mu']
    sigma=prm[prefix+'sigma']
    if kwargs.get('use_sigrow',0):
        sigrow=prm[prefix+'sigrow']
    else:
        sigrow=None
    gamma=prm[prefix+'gamma']
    sigma_k=prm[prefix+'capacity_std']
    mean_k=prm[prefix+'capacity_mean']
    sigma_d=prm[prefix+'selfint_std']
    # N0 =None
    #
    # if len(S)>20:
    #     #For many groups, seed with empirical biomass
    #     N0=prm[prefix+'biomass*']

    #print S,mu,sigma,gamma,sigma_k,mean_k,sigma_d
    #raise
    return group_cavity_solve_props(S=S,mu=mu,sigma=sigma,sigma_k=sigma_k,gamma=gamma,
        corrKA=prm.get(prefix+'corr_KA',0),mean_k=mean_k,sigma_d=sigma_d,sigrow=sigrow,
                                   # N0=N0,
                                   **kwargs
        )

def measure_groups_cavity(model,measure,**kwargs):
    prefix=kwargs.get('prefix','groups_')
    kwargs['prefix']=prefix
    res= group_cavity_output(measure,**kwargs)

    measure.update(res)

    from datatools import nonan
    measure['cavity_compare']=np.mean([ np.abs(np.mean( nonan(
        measure[GROUPS_COMPARE[x]])/np.mean(measure[prefix+x]+0.0000001)-1) )
         for x in ('N1','stdN','Phi') ] )




def plot_groups(measures=None,axes=None,split_by=None,log='',dictionary=None,prefix='groups_',dataprefix='n_groups_',title='',
        traits=None,**kwargs):
    if dictionary is None:
        dictionary={}
    dictionary=dict(dictionary)
    compare={'Phi':'%alive','avgN':'biomass','totN':'biomass_tot',
        'totNmf':'biomass_tot',
        'stdN':'biomass_std','simpson':'simpson'}
    group_equiv={'Phi':'phi'}
    if not axes:
        axes=[]
    if not split_by:
        split_by=[]
    split_by=[z for z in split_by if z!='n_groups']
    if split_by:
        cur=split_by[0]

        for val in sorted(set(measures[cur])):
            plot_groups(measures=measures[measures[cur]==val],axes=axes,split_by=split_by[1:],
                log=log,dictionary=dictionary,prefix=prefix,title=title+'; {}:{}'.format(cur,val) ,
                **kwargs)
        return
    fig=plt.figure()
    fig.suptitle(title)
    split_by=ifelse(not 'n_groups' in axes, ['n_groups'],[])
    if traits is None:
        traits=sorted(compare)
    nbpanels=len(traits)
    panel=MutableInt(0)
    for  trait in traits:
        auto_subplot(plt,nbpanels,panel,rows=ifelse(len(traits)<4, 1,None) )
        ##if trait != 'Phi':
            #llog=log+'y'
        #else:
        llog=log
        plot_data(measures=measures,axes=axes,
            values=[(prefix+trait,{'style':'plot','linewidth':2}), (dataprefix+compare[trait],{'style':'scatter' } ) ] ,
                hold=1,newfig=0,
            split_by=split_by,log=llog,dictionary=dictionary,show_legends=((panel==1) and split_by),
            **kwargs
            )
        plt.title(cavity_dico.get(group_equiv.get(trait,trait),trait))
        plt.ylabel('')
