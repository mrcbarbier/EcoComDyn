# -*- coding: utf-8 -*-


from functions import *
from models import LabeledArray


def measure_gen(model,measure,typ='usual',**kwargs):
    '''Function called by most other measures'''
    model=model.copy(initialize=False)
    prm = model.parameters

    for var in model.results:
        #GENERAL SETUP
        comm=[i for i,j in prm.iteritems() if i in model.data and j.get('role',None)=='interactions'
            and var in j['variables'][0]]
        if not comm:
            continue
        comm=comm[0]

        growth=[i for i,j in prm.iteritems() if j.get('role',None)=='growth'
            and var in j['variables'][0]][0]
        A=model.get_labeled_data(comm)
        r=model.get_labeled_data(growth)
        axis=A.axes[0]

        try:
            capa=[i for i,j in prm.iteritems() if j.get('role',None)=='capacity'
                and var in j['variables'][0]][0]
            K=model.get_labeled_data(capa)
        except:
            capa=[i for i,j in prm.iteritems() if j.get('role',None)=='diagonal'
                and var in j['variables'][0]][0]
            if not capa in model.data:
                model.initialize()
            D=model.get_labeled_data(capa)
            K=r/D
            K.matrix=np.clip(K.matrix,0,100)
            measure['{}_capacity_std'.format(var) ]=np.std(K.matrix)
            measure['{}_capacity_mean'.format(var) ]=np.mean(K.matrix)
            measure['{}_selfint_std'.format(var) ]=np.std(D.matrix)

        try:
            traj=model.results[var]
            Nf=LabeledArray(traj.matrix[-1],axes=prm[var]['axes'])
        except:
            code_debugger()
        othax=tuple( i for i,a in enumerate(prm[var]['axes'])
             if not a == axis[-1] )
        death=kwargs.get('death',prm[var].get('death',0))
        if othax:
            alive=np.where(np.sum(Nf.matrix>death,axis=othax))[0]
        else:
            alive=np.where(Nf.matrix>death)[0]
        #print alive, Nf,np.sum(Nf.matrix>death,axis=othax),othax,axis,prm[var]['axes']

        Nlive=Nf.matrix.mean(othax)[alive]
        nlive=Nlive/np.mean(Nlive)

        def idxs(data):
            return np.ix_(*[alive if ax==axis else [0] for ax in data.axes ])
        S=A.matrix.shape[0]
        Alive=A[idxs(A)]
        rlive=r[idxs(r)]
        Klive=K[idxs(K)]
        Slive=len(alive)
        dx=model.get_dx(traj.index[-1],model.get_labeled_result(idx=-1))[var].matrix/Nf.matrix
        Slivestrict=S-np.sum(dx<-10**-8)


        #SPECIFIC MEASURES
        if 'abundance' in typ:
            measure['{}_abundance'.format(var) ]=tuple(Nf.matrix.ravel())

        if 'usual' in typ:
            '''Typical measurements: biomass, lyapunov function, number of species alive'''
            ax=prm[var]['axes'].index(axis[-1])
            measure['{}_#alive'.format(var)]=Slive
            measure['{}_%alive'.format(var)]=Slive/float(traj.shape[1+ax])
            measure['{}_%alive_strict'.format(var)]=Slivestrict/float(traj.shape[1+ax])
            #print traj.shape[1+ax], len(alive)
            measure['{}_biomass'.format(var)]=np.mean(traj[-1].sum( othax ))
            measure['{}_biomass2'.format(var)]=np.mean((traj[-1]**2).sum( othax ))
            measure['{}_biomass*'.format(var)]=measure['{}_biomass'.format(var)]/measure['{}_%alive'.format(var)]
            measure['{}_biomass_tot'.format(var)]=np.sum(traj[-1])
            measure['{}_biomass_std'.format(var)]=np.std(traj[-1].sum( othax ))
            measure['{}_simpson'.format(var)]=np.sum(traj[-1].sum( othax )**2)/measure['{}_biomass_tot'.format(var)]**2

            measure['{}_biomass_relstd'.format(var)]=measure['{}_biomass_std'.format(var)]/measure['{}_biomass'.format(var)]

            measure['{}_productivity'.format(var)]=np.sum(r.matrix*traj[-1])
            measure['{}_productivity_ratio'.format(var)]=measure['{}_productivity'.format(var)]/measure['{}_biomass_tot'.format(var)]
            measure['{}_diff'.format(var)]=np.median(traj[-1].std( othax ))
            #J= np.dot(np.diag(nlive),  )

        if 'degree' in typ:
            mat=A.matrix.copy()
            indeg=np.sum(mat!=0,axis=1)
            outdeg=np.sum(mat!=0,axis=0)
            measure['n_degree']=np.mean(indeg)
            measure['n_degree_std']=np.std(indeg)
            measure['n_alive_degree'.format(comm)]=np.mean(np.sum((A.matrix[np.ix_(alive,alive)]!=0),axis=1) )

        if 'effective' in typ:
            '''Measure effective value of interactions'''

            prefix='{}_'.format(var)
            try:
                rescaled_int=(A/D).matrix.copy()
            except:
                rescaled_int=(A*K/r).matrix.copy()

            matrix=A.matrix.copy()
            S=matrix.shape[0]
            if not np.isnan(rescaled_int).any():
                measure['n_connectivity']=np.sum(np.abs(rescaled_int)+np.abs(rescaled_int.T)>10**-15)*1./ (rescaled_int.size-S )
                #print measure['n_connectivity']
            else:
                measure['n_connectivity']=np.sum(np.abs(matrix)>10**-15)*1./ (matrix.size-S )
            np.fill_diagonal(matrix,np.nan)
            np.fill_diagonal(rescaled_int,np.nan)
            mats=[matrix,rescaled_int ]
            if 'finalprm' in typ:
                matlive=Alive.copy()
                np.fill_diagonal(matlive,np.nan)
                mats.append(matlive)

            for mat,pref in zip( mats,[prefix+'interactions_',prefix+'couplings_',prefix+'couplings_final_'] )[1:]:
                effS=np.sum((np.abs(mat)+np.abs(mat.T)) !=0)*1./mat.shape[0]#*measure['n_connectivity']
                ref=mat.copy()
                ref[np.isnan(ref)]=0
                mat=mat.copy()
                #mat[(np.abs(mat)+np.abs(mat).T) ==0]=np.nan #Ignore missing edges
                effS=mat.shape[0]
                mat2=mat.copy()
                #mat2[(mat) ==0]=np.nan #Ignore missing edges
                if not 'threshold' in model.data and 0:
                    #Remove edges that are too large
                    for i in range(mat.shape[0]):
                        if (mat[i]<-1).any():
                            mat[i,:]=np.nan
                            mat[:,i]=np.nan
                measure[pref+'mean'] = mu = np.mean(nonan(mat2))
                measure[pref+'std'] = std = np.std(nonan(mat2))
                measure[pref+'row_std']  =np.std(np.sum(ref,axis=1) ) # np.std([np.mean(nonan(mat[i]))
                       #  for i in range(S) if len(nonan(mat[i]))] )

                colstd=np.std(ref-np.mean(ref,axis=1),axis=1 )
                measure[pref+'col_std']  =np.std(colstd)*np.sqrt(effS)
                measure[pref+'col_sigma']  =np.mean(colstd)*np.sqrt(effS)
                measure[pref+'rowcol']  =np.mean(colstd*np.sum(ref,axis=1))
                tmp=np.abs(ref)+np.abs(ref.T)
                Nnei=np.dot((tmp!=0),Nf.matrix) / np.sum(tmp!=0,axis=1)
                measure[pref+'n_std']=np.std(Nnei)/np.mean(Nnei)
                measure[pref+'n_mean']=np.mean(Nnei)

                sym= (np.mean(nonan((mat-mu)*(mat-mu).T)))#-mu**2)
                #print sym, (np.mean(nonan(mat*mat.T))-mu**2)
                if 0:
                    measure[pref+'symmetry']=sym
                    measure[pref+'symstd']= np.std(nonan(mat*mat.T))
                mus=[-mu*effS,-np.mean(np.sum(ref,axis=1))]
                measure[pref+'mu']=mus[0]
                sigs=[std*np.sqrt(effS), ]
                measure[pref+'sigma']=sigs[0]
                if std>0:
                    gams=sym/std**2,#,np.mean(rm*rm.T)/np.mean(rm**2) ] #OPTION 2 IS WRONG: all the missing links contribute to the mean here
                    measure[pref+'gamma']=gams[0]
                else:
                    measure[pref+'gamma']=0

                #print ref[0],ref[:,0]
                #print mus,sigs,gams
                #print np.sum(ref,axis=1),  np.mean( nonan(mat)),np.mean(ref[ref!=0] )
                #print comm, measure[prefix+'symmetry'], std**2, mat[0,1],mat[1,0]
            measure['{}_growth_mean'.format(var)]=np.mean(r.matrix)
            measure['{}_growth_std'.format(var)]=np.std(r.matrix)
            measure['{}_growth_not'.format(var)]=np.sum((r.matrix<=0))
            adjusted=K.shape[0] *(K.shape[0]-1)
            measure['{}_corr_KA'.format(var)]=np.sum((-K*A/D).matrix)/adjusted -np.mean(
                K.matrix)*np.sum(-(A/D).matrix)/adjusted

        from datatools import hist, scatter
        if 'removal' in typ:
            remends=[]
            from trajectory import Trajectory
            difs=[]
            n0=[]
            vs=np.zeros((S,S))
            try:
                rescaled_int=(A/D).matrix.copy()
            except:
                rescaled_int=(A*K/r).matrix.copy()
            offdiag=zip(*[(i,j) for i in range(S) for j in range(S) if i!=j ] )
            a=-(rescaled_int- np.mean(rescaled_int[offdiag]) )/np.std(rescaled_int[offdiag])/np.sqrt(S)
            a[np.isnan(a)]=0
            for i in range(S):
                oth=range(S)
                oth.remove(i)
                if (not i in alive) or np.abs(Nf.matrix[i])<10**-12:
                    remends.append(Nf.matrix[oth]  )
                    continue
                #print "Removing",i
                mat=Nf.matrix.copy()
                mat[i]=0
                model.results[var]=Trajectory(init=mat)
                model.evol(tmax=50000,tsample=50000,reseed=0,keep='endpoint')

                nothf=Nf.matrix[oth]
                nothr=model.results[var][-1][oth]
                remends.append(nothr)
                difs.append(nothr-nothf)
                #tmp=a[oth,i]*measure['v']*Nf.matrix[i]/np.mean(Nf.matrix)
                #scatter(tmp,difs[-1]/np.mean(Nf.matrix))
                #vs[oth,i]=(nothf/np.mean(nothf) - nothr/np.mean(nothr) )*np.mean(Nf.matrix)/Nf.matrix[i]/a[oth,i]
                vs[oth,i]=(nothr - nothf )/Nf.matrix[i]/a[oth,i]
            vs[a==0]=0
            #print np.mean(np.abs(a.ravel() ) )
            #hist(vs.ravel(),log='y')

            measure['removal_diff']=np.mean([np.sum(np.abs(d)) for d in difs])/np.mean(Nf.matrix)
            if (vs!=0).any():
                measure['removal_v']=np.mean([np.mean(v[v!=0]) for v in list(vs)])
                measure['removal_v_std']=np.std([np.mean(v[v!=0]) for v in list(vs)])
            else:
                measure['removal_v']=measure['removal_v_std']=0
            #print measure['removal_v'],measure['removal_v_std']
            model.results[var]=Trajectory(init=Nf.matrix)

        if 'testcavity' in typ:
            n0=[]
            h=measure['h']
            v=measure['v']
            q=measure['q']
            phi=measure['phi']
            sigma=measure['n_couplings_sigma']
            mu=measure['n_couplings_mu']
            gamma=measure['n_couplings_gamma']
            sigma_k=measure['n_capacity_std']
            avgN=measure['avgN']
            meanN=np.mean(Nf.matrix)
            sigma_l=(sigma_k/sigma/avgN)
            l0s=[]
            try:
                rescaled_int=(A/D).matrix.copy()
            except:
                rescaled_int=(A*K/r).matrix.copy()
            #np.fill_diagonal(rescaled_int,np.nan)

            offdiag=zip(*[(i,j) for i in range(S) for j in range(S) if i!=j ] )
            a=-(rescaled_int+mu/S)/sigma

            sums=[]
            dN=[]
            difs=[]
            difnexts=[]
            n0rem=[]
            compv=[]
            for i in range(S):
                oth=range(S)
                oth.remove(i)
                diff=a[oth,i]*v*Nf.matrix[i]/meanN
                difs.append(diff)
                diffnext=np.dot(a[oth,:],a[:,i])*v**2 *Nf.matrix[i]/meanN
                difnexts.append(diffnext)

                nbefore= Nf.matrix[oth]/meanN + diff#+diffnext
                nbefore[Nf.matrix[oth]<10**-15]=0
                SUM=np.dot( a[i,oth],nbefore)
                l0=(K.matrix[i]-measure['n_capacity_mean'])/sigma/meanN
                l0s.append(l0)

                u=(1-mu/S)/sigma
                ut=(u-gamma*v)
                res=(h + l0 -SUM) /ut
                n0.append(res)
                sums.append(SUM)

                dd=(K.matrix[i]-Nf.matrix[i] + np.dot( rescaled_int[i,oth], Nf.matrix[oth]) )
                dN.append(dd)
                if 'removal' in typ:
                    nbefore2=remends[i]/meanN
                    #scatter(diff,nbefore2-Nf.matrix[oth]/meanN)
                    SUM2=np.dot( a[i,oth],nbefore2)
                    res2=(h + l0 -SUM2) /ut
                    n0rem.append(res2)
                    comp=nbefore-nbefore2
                    compv.append(comp)
            #scatter(n0,n0rem)
            n0rem=np.array(n0rem)
            n0=np.array(n0)
            dN=np.array(dN)
            measure['cavity_alive_truepos']=np.sum(n0[alive]>0 )*1./len(alive)
            measure['cavity_alive_falseneg']=np.sum(n0[alive]<0 )*1./len(alive)
            measure['cavity_alive']=np.sum(n0>0 )*1./n0.shape[0]
            measure['cavity_dN']=np.sum(dN>-10**-6 )*1./n0.shape[0]
            #print measure['cavity_dN']
            measure['cavity_diff']=np.mean([np.sum(np.abs(d)) for d in difs])
            measure['cavity_diffnext']=np.mean([np.sum(np.abs(d)) for d in difnexts])
            if 'removal' in typ:
                measure['cavity_v_vs_rem']=np.mean([np.sum(np.abs(d)) for d in compv])
                measure['cavity_alive_rem']=np.sum(n0rem>0 )*1./n0.shape[0]

            measure['cavity_meansum']=np.mean(SUM)
            measure['cavity_corrsum']=np.mean(np.array(SUM)*l0)-np.mean(SUM)*np.mean(l0)
            if 0 and rescaled_int.any():
                #PLOTS!
                from datatools import hist,plot,plt,scatter
                #if dN.any():
                    #print dN
                    #hist(dN)

                plt.subplot(224)
                hist(a[offdiag],hold=1,normed=1,log='y')
                hist(sums,hold=1,normed=1)
                plt.subplot(223)
                Nf.axes=[('n','com')]
                t=((r-D*Nf+A.dot(Nf,axes=[('n','com')])).matrix)/avgN
                hist(t[alive],hold=1,normed=0,bins=20,log='y')
                plt.subplot(221)
                xs=np.linspace(min(Nlive/avgN),max(Nlive/avgN),100)
                plt.ylim(ymin=0.0001,ymax=1)
                hist(Nlive/avgN,hold=1,normed=1,bins=20)
                var=(q-1)/phi
                plot(xs, np.exp(-(xs-1/phi)**2/2/var )/np.sqrt(var) ,log='y',hold=1,title='n+')

                plt.subplot(222)
                xs=np.linspace(min(n0),max(n0),100)
                plt.ylim(ymin=0.0001,ymax=1)
                hist(n0,hold=1,normed=1,bins=20,title='n0')
                #print 'SHOULD HAVE',S,mu,sigma,sigma_k,gamma
                print 'RECEIVING q v h phi',q,v,h,phi

                print 'phi',len([z for z in t if z>-10**-16])*1./S,len(alive)*1./S,phi,len([z for z in n0 if z>0])*1./S
                print '<n>',np.mean(Nf.matrix/avgN),1
                print '<n^2>',np.mean((Nlive/avgN)**2),q
                print '<n>+',np.mean(Nlive/avgN),1/phi
                print '<n^2>+',np.mean((Nlive/avgN)**2),q/phi
                print '<n0>',np.mean(n0),h /ut
                print 'var_n0',np.var(n0),(q+ sigma_l**2) /ut**2
                print '    (in which q',q/ut**2,'sigL', sigma_l**2/ ut**2,')'
                print 'sig_l',np.std(l0s),sigma_l,'sig_k',sigma_k, np.std(K.matrix)
                plot(xs, np.exp(-((ut* xs-h)**2/2/(q+ sigma_l**2) ) )/np.sqrt(q+ sigma_l**2) ,log='y')

        if 'finalprm' in typ:
            measure['{}_growth_final_mean'.format(var)]=np.mean(rlive)
            measure['{}_growth_final_std'.format(var)]=np.std(rlive)
            measure['{}_growth_final_not'.format(var)]=np.sum((rlive<=0))
            measure['{}_capacity_final_mean'.format(var)]=np.mean(Klive)
            measure['{}_capacity_final_std'.format(var)]=np.std(Klive)


        if 'matrix' in typ:
            def variance(J,D):
                #Jeff's variance
                #Jhat= lifted_matrix(J)
                from scipy.linalg import solve_continuous_lyapunov as solve_lyapunov, norm, inv
                #tmp=np.zeros(J.shape)
                #tmp[0,:]=Nlive
                #res=np.dot(inv(Jhat),tmp.ravel() )
                #print '1',res.reshape(J.shape)[:5,:5]

                #tmp=np.zeros(J.shape)
                #tmp[:,0]=Nlive
                #res=np.dot(inv(Jhat),tmp.ravel() )
                #print '2',res.reshape(J.shape)[:5,:5]


                res = solve_lyapunov(J,-np.dot(D,D) )
                #print 'X',(res)[:5,:5]
                return np.trace(res)/Slive

            Dlive=rlive/Klive
            Dlive[Klive==0]=1
            #DP=(r/K).matrix
            #DP[r.matrix==0]=1

            if 0:
                B=Alive-np.eye(Slive)*Dlive
                #BP=A.matrix-np.eye(S)*DP
            else:
                B=(Alive.T/Dlive).T-np.eye(Slive)

            D=np.diag(Nlive)
            J=np.dot(D,B)
            #2*variance(J,D**(1./2) )
            #raise

            from scipy.linalg import solve_continuous_lyapunov as solve_lyapunov, norm, inv
            from scipy.sparse.linalg import eigs
            if Slive>0:
                #measure['{}_empirical_press'.format(var)]=np.linalg.norm(np.linalg.inv(B),ord='fro')**2/Slive
                try:
                    measure['{}_matrix_press'.format(var)]=norm(inv(B),ord='fro')**2/Slive
                    measure['{}_matrix_var'.format(var)]=2*variance(J,D**(1./2) )
                except Exception as e:
                    print e
                    Slive=0
                #print CIJ[:3,:3]

            if Slive==0:
                #measure['{}_empirical_press'.format(var)]=0
                for suffix in ('press','var','eig','symeig'):
                    measure['{}_matrix_{}'.format(var,suffix)]=0
            elif 0:
                CIJ=solve_lyapunov(J,-D )
                from datatools import scatter,plot,hist
                prefix='n_couplings_'
                mu=measure[prefix+'mu']
                sigma=measure[prefix+'sigma']
                gamma=measure[prefix+'gamma']
                a=(-Alive-mu/S)/sigma
                ut=measure['utilde']
                n=Nlive
                np.fill_diagonal(a,0)
                AA=( (a**2).T *  (Nlive/np.add.outer(Nlive,Nlive)  ) ) .T
                AG=( (a*a.T ) *   (Nlive/np.add.outer(Nlive,Nlive)  ) )

                DEN=ut**2 - np.sum(AG,axis=1)
                RESULT=np.dot( inv(np.diag(DEN) -AA  ), ut/2 * np.ones(Slive)/sigma  ) #/Slive

                if 0:
                    #plot(np.diag(CIJ), np.diag(CIJ),hold=1)
                    #scatter( np.diag(CIJ) ,(ut/2/sigma + np.dot(AA,np.diag(CIJ) ) ) /DEN ,hold=1,log='xy')
                    CIJ2=CIJ.copy()
                    np.fill_diagonal(CIJ2,0)
                    test=np.dot(a,CIJ2)
                    np.fill_diagonal(test,0)
                    print 'a.(C-diagC),',test[:3,:3],np.mean(test)*Slive
                    print 'a.C',np.dot(a,CIJ)[:3,:3],np.mean(np.dot(a,CIJ))*Slive
                    print 'mean offdiag',np.mean(CIJ2), 'diag',np.mean(np.diag(CIJ))
                    print 'diag*u',np.mean(np.diag(CIJ))*ut
                #hist(CIJ2[CIJ2!=0],log=1)
                #plot(np.diag(CIJ), np.diag(CIJ),hold=1)
                #scatter( np.sum(DEN*np.diag(CIJ)) ,np.sum(ut/2/sigma + np.dot(AA,np.diag(CIJ) ) )  ,hold=1,log='xy')

                #DEN2=ut**2 - np.sum(AG,axis=1)
                #scatter( np.sum(DEN2*np.diag(CIJ)) ,np.sum(ut/2/sigma + np.dot(AA,np.diag(CIJ) ) )  ,color='r',hold=1,log='xy')


                measure['{}_pred_var'.format(var)]=  np.mean((ut/2/sigma + np.dot(AA,np.diag(CIJ) ) ) /DEN)*2 #np.mean(RESULT)*2
                print '{}_pred_var'.format(var), measure['{}_pred_var'.format(var)],measure['variability2']



        if 'matrix_eig' in typ:
                try:
                    measure['{}_matrix_eig'.format(var)]=eigs(J,k=1,sigma=0)[0][0].real
                    measure['{}_pool_symeig'.format(var)]=eigs( (BP+BP.T)/2.,k=1,sigma=0)[0][0].real
                except:
                    pass
                #test= np.sum(
                    #np.linalg.inv(B)**2)/Slive
                #assert abs(test / measure['{}_empirical_chi2'.format(var)] -1) <0.001



        if 'corr' in typ:
            prefix='{}_interactions_corr_'.format(var)
            def corr(form):
                return np.mean(nonan(np.einsum(form +'->ijk',mat,mat)) )
            measure[prefix+'ijkj']= (corr('ij,kj') - mu**2)/std**2
            measure[prefix+'ijik']=(corr('ij,ik') - mu**2)/std**2
            measure[prefix+'ijjk']=(corr('ij,jk') - mu**2)/std**2



        if 'trophic' in typ:
            '''Trophic score and fraction of basal species.
            NB: should I decide basal status from absence of outgoing trophic links,
            or from nonzero carrying capacity in the absence of trophic interactons?'''
            if not alive.shape[0]:
                Nlive=np.ones(1)
                tscore=np.zeros(1)
            else:
                tscore=trophic_score(Alive,
                    pops=Nlive,
                    influx= rlive*Nlive )


            prefix='{}_trophic_'.format(var)
            (measure[prefix+'max'],measure[prefix+'mean'],
                measure[prefix+'std'],measure[prefix+'20'],
                measure[prefix+'80'])=[
                    x(tscore)  for x in (
                    #lambda e: np.sum(e==0)/float(len(e)),
                    np.max,np.mean,np.std,
                    lambda e:np.percentile(e,20),
                    lambda e: np.percentile(e,80))
                     ]

            Asum=Alive.sum(axis=1)

            niches=[i for i,j in prm.iteritems() if j.get('role',None)=='niche']
            if niches:
                niche=model.data[niches[0]].matrix[alive]
                measure[prefix+'niche']=np.corrcoef(niche,tscore)[0,1]

            measure[prefix+'weighted_mean']=np.sum(tscore*Nlive)/np.sum(Nlive)
            measure[prefix+'corr_r']=np.mean(tscore*rlive) / (np.mean(tscore)*np.mean(rlive)) -1
            measure[prefix+'corr_N']=np.mean(tscore*Nlive)/ (np.mean(tscore)*np.mean(Nlive))-1
            measure[prefix+'corr_A']=np.mean(tscore*Asum)/ (np.mean(tscore)*np.mean(Asum))-1
            #print Asum, measure[prefix+'corr_A']

            if kwargs.get('trophic_score_pool',0):
                #DEACTIVATE FOR FASTER MEASUREMENTS
                pool_tscore=trophic_score(A.matrix,
                    influx= r.matrix )
                measure[prefix+'pool_max']=np.max(pool_tscore)
                measure[prefix+'pool_mean']=np.mean(pool_tscore)

            measure[prefix+'basal']=np.sum(tscore<1.1)/float(len(tscore))

        if 'assembly_corr' in typ:
            prefix='{}_corr_'.format(var)
            mask = np.ones(Alive.shape, dtype=bool)
            np.fill_diagonal(mask, 0)
            meanAlive=np.mean(Alive[mask])
            #meanK=np.mean(Klive)
            #meanNf=np.mean(Nf[alive] )
            #print meanA.shape, meanK.shape, Alive.shape,Klive.shape,(Alive*Klive).shape,np.mean((Alive*Klive)[mask]).shape
            if abs(meanAlive)<10**-5:
                meanAlive=1
                measure[prefix+'AK']=0
                measure[prefix+'AN']=0
            else:
                measure[prefix+'AK']=np.mean((Alive.T*Klive)[mask] )#/meanAlive/meanK -1
                measure[prefix+'AN']=np.mean( (Alive.T*nlive )[mask] )#/meanAlive/meanNf-1


def measure_usual(model,measure,**kwargs):
    return measure_gen(model,measure,typ='usual',**kwargs)

def measure_trophic(model,measure,**kwargs):
    return measure_gen(model,measure,typ='trophic',**kwargs)

def measure_assembly_corr(model,measure,**kwargs):
    return measure_gen(model,measure,typ='assembly_corr',**kwargs)

def measure_dead(model,measure,**kwargs):
    '''Clean up species that aren't extinct yet but are going to extinction and cannot invade.'''
    Nf=model.get_labeled_result('n',idx=-1).copy()
    S=(Nf.shape[0])
    death=kwargs.get('death',model.parameters['n'].get('death',10**-15))*1.1
    tf=model.results['n'].index[-1]
    # code_debugger()

    import scipy.optimize as sopt
    def get_dx(xx):
        x=10.**np.clip(xx,-20,15)
        N2 = Nf.copy()
        N2.matrix=x
        return model.get_dx(tf,{'n':N2})['n'].matrix/x

    r=sopt.root(get_dx, np.log10(Nf.matrix) )
    cleaned=np.where(np.logical_and(10**r.x<death,Nf.matrix>death) )[0]
    # model.results['n'].matrix[-1,cleaned]=death/2.
    measure['dead_removed']=tuple(cleaned)
    measure['#dead_removed']=len(cleaned)

    print 'Removed dead:', len(cleaned)#,cleaned, Nf[cleaned]
    model.results['n'].matrix[-1]=np.maximum(death/2,10**r.x)
    # code_debugger()
    return

    " OBSOLETE "

    print r.x
    code_debugger()

    dN=model.get_dx(tf,{'n':Nf})

    for i in range(S):
        if Nf[i]<death:
            continue
        if dN['n'][i]/Nf[i] < 0 and Nf[i] < max(death*10, np.mean(Nf.matrix)/np.sqrt(S)) :
            N2=Nf.copy()
            N2.matrix[i]=0
            dN2 = model.get_dx(tf, {'n':N2} )
            if dN2['n'][i]<=0:
                Nf.matrix[i]=0
    dNf = model.get_dx(tf, {'n':Nf})
    cleaned=[]
    for i in range(S):
        if Nf[i]==0 and dNf['n'][i] <=0 and model.results['n'][-1][i]>death:
            cleaned.append(i)
    model.results['n'][-1][cleaned] =death/1000.
    measure['dead_removed']=tuple(cleaned)
    measure['#dead_removed']=len(cleaned)

def measure_niche(model,measure,**kwargs):
    '''Measure assembled niche distribution'''
    prm = model.parameters
    model=model.copy(initialize=False)

    niches=[i for i,j in prm.iteritems() if j.get('role',None)=='niche']
    for niche in niches:
        axis=prm[niche]['axes'][0]
        if hasattr(axis,'__iter__'):
            var=axis[0]
        else:
            axis=(model.variables.keys()[0],axis)
            var=axis[0]

        #print 'meas_niche',niche, var
        comm=[i for i,j in prm.iteritems() if j.get('role',None)=='interactions'
            and var in j['variables'][0] and j.get('dynamics',None) ][0]
        A=model.get_labeled_data(comm)
        traj=model.results[var]

        Nf=LabeledArray(traj.matrix[-1],axes=prm[var]['axes'])
        othax=tuple( i for i,a in enumerate(prm[var]['axes'])
             if not a == axis[-1] )
        death=prm[var].get('death',0)
        alive=np.where(np.sum(Nf.matrix>death,axis=othax))[0]
        if alive.any():
            n=model.data[niche][alive]
        else:
            n=np.zeros(1)
        measure['{}_alive_mean'.format(niche)]=np.mean(n)
        measure['{}_alive_max'.format(niche)]=np.max(n)

############ OBSOLETE, MUST BE REDONE ##################

def measure_graph(model,measure,metrics=['all','fast'][1]):
    '''Graph measurements: degree, connectivtiy, betweenness...'''
    print 'WARNING: graphs not working for metacommunities'
    for traj in [ model.results['n'] ]:
        #plot(traj.index,traj.matrix,log='xy')
        alive=np.where(traj.matrix[-1]>0)[0]
        graph=nx.DiGraph(A.matrix[np.ix_(alive,alive)])
        allgraphs.append(graph)

    kmean=np.median([np.mean([g.degree(x) for x in g.nodes()] ) for g in allgraphs] )
    measure['kmean']=kmean
    measure['kmean/kmean0']=kmean/A.props['kmean']

    if metrics=='all':
        #SLOWER!
        betweenness=np.median([np.mean(nx.betweenness_centrality(g).values()) for g in allgraphs] )
        connectivity=np.median([nx.average_node_connectivity(g) for g in allgraphs] )
        measure['betweenness']=betweenness
        measure['connectivity']=connectivity


def measure_response(model,measure,**kwargs):
    allv=[]

    A=model.data['community']
    mu,S,sigma=[measure.get(x) for x in ('dissipation','N','strength')]
    u=(1-mu/S)/sigma
    for traj in [ model.results['n'] ]:
        alive=np.where(traj.matrix[-1]>0)[0]
        mat=(u-1)* np.eye( len(alive) ) + (A.matrix[np.ix_(alive,alive)]/S -mu)
        allv.append(np.mean(np.dot(np.linalg.inv(mat),np.ones(mat.shape[0]))))
    measure['response']=np.mean(allv)

def measure_varint(model,measure,**kwargs):
    '''Interaction between variances'''
    allt=[]
    A=model.data['community']
    for traj in [ model.results['n'] ]:
        alive=np.where(traj.matrix[-1]>0)[0]
        if 0:
            mat=A.matrix[np.ix_(alive,alive)]
            var=variance_interactions(mat,lift=True)
            allt.append(var)
        else:
            mat=variance_interactions(A.matrix,lift=True)
            var=mat[np.ix_(alive,alive)]
            allt.append((mat,var))

    res=np.array([[
            x(t[1])/x(t[0])-1 for x in (
                np.mean,np.var,
                lambda e:np.percentile(e,20),
                lambda e: np.percentile(e,80))
                 ]for t in allt])

    print 'varint', res

    if len(res.shape)>1:
        res=np.median(res,axis=0)

    (measure['varint_mean'],
        measure['varint_var'],measure['varint_20'],
        measure['varint_80'])=res
