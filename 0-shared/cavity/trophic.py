# -*- coding: utf-8 -*-
from datatools import *

from utility import *


def erfmom(mom,mean,var):
    thresh=0.0000001
    from numpy import sqrt,pi,exp
    from scipy.special import erf

    lowvar=(var<thresh)
    var=np.clip(var,thresh,None)
    present=(mean>0).astype('float')
    xx=mean/sqrt(2*var)

    mom0=.5 * (erf(xx)+1 )
    mom1=sqrt(var/2/pi) * exp(-xx**2)

    if mom==0:
        res= mom0
        res[lowvar]=present[lowvar]
    elif mom==1:
        res= mean* mom0 + mom1
        res[lowvar]=mean[lowvar]*present[lowvar]
    elif mom==2:
        res= (var+mean**2)*mom0 + mean*mom1
        res[lowvar]=mean[lowvar]**2*present[lowvar]
    return res

def cavity_w(k,delta):
    import scipy.integrate as sinteg
    res= sinteg.quad(lambda z:np.exp(-z**2/2)/np.sqrt(2*pi) *(z+delta)**k,-delta,np.inf)
    if res[1]>0.1:
        print 'WARNING: LARGE ERROR IN cavity_W'
    return res[0]

def trophic_cavity_solve_old(S=[100,100],mu=1,sigma=.1,mean_k=1., sigma_k=1,epsilon=.1,x0=None ,
        CALC_INTEGRAL=0,**kwargs):
    import scipy.optimize as sopt

    levels=len(S)
    S,mu,sigma,mean_k,sigma_k,epsilon=[np.array(x) for x in S,mu,sigma,mean_k,sigma_k,epsilon]
    if x0 is None:
        x0=np.ones(4*levels)
    avgN=1
    def eqs(vec,calc=CALC_INTEGRAL):
        meanN,varn,V,Phi=vec.reshape((4,levels))
        B=S*meanN/avgN
        loc={}
        loc.update( locals())
        up={z:np.concatenate((loc[z][1:],[0])) for z in ('varn','meanN','V','Phi','B') }
        down={z:np.concatenate(([0],loc[z][:-1])) for z in ('varn','meanN','V','Phi','B') }

        mean0= V*(mean_k/avgN + mu*(epsilon*up['B'] - down['B'] ) )
        var0= V**2*( sigma_k**2/avgN**2 + sigma**2  *(up['Phi']*up['varn']+
            epsilon**2 *down['Phi']*down['varn'] ) )

        eq4=1 - V*(1+sigma**2*epsilon *(up['Phi']*up['V']+down['Phi']*down['V'] ) )
        if calc:
            raise Exception("TODO: explicit calculation in trophic_cavity_solve")
        else:
            #Equivalent with directly computed moments of erf
            eq1=Phi-erfmom(0,mean0,var0)
            eq2=1-erfmom(1,mean0,var0)
            eq3=varn-erfmom(2,mean0,var0)
        res= np.concatenate( (eq1,eq2,eq3,eq4))

        zerocond=np.min(np.concatenate([varn.ravel(),Phi.ravel(),meanN.ravel()])) #Things that should be positive!
        if zerocond<0:
            res+=np.abs(zerocond )

        return res
    #root=sopt.newton_krylov
    #root=sopt.anderson
    root=sopt.root
    res= root(eqs,x0)

    result= res.x
    meanN,varn,V,Phi=result.reshape((4,levels))
    if not res.success or np.max(np.abs(eqs(result)))>10.**-4:
        print 'ERROR: {} {}'.format(res.message,list(res.fun.reshape((4,levels))) )

        print 'PARAMS: S {} mu {} sigma {} sigma_K {} epsilon {}'.format(S,mu,sigma,sigma_k,epsilon)
        meanN,varn,V,Phi=np.zeros((4,levels))
    return result.reshape((4,levels))



def trophic_cavity_solve(S=[100,100],mu=1,sigma=.0,mean_k=1., sigma_k=0,sigma_d=0,epsilon=.9,x0=None ,
        CALC_INTEGRAL=0,**kwargs):
    import scipy.optimize as sopt

    levels=len(S)
    if len(mu.shape)==2:
        mu=[0]+list(np.diag(mu,k=-1))
        sigma=[0]+list(np.diag(sigma,k=-1))
    S,mu,sigma,mean_k,sigma_k,epsilon=[np.array(x).astype('float')  if np.array(x).shape
        else np.array([float(x) for nb in range(levels)] )
         for x in S,mu,sigma,mean_k,sigma_k,epsilon]
    if x0 is None:
        x0=np.random.random((4,levels))
        #for i in range(levels):
            #x0[0,i]/=(mu[i]/epsilon[i] )**i
    x0=x0.ravel()

    s=S/np.mean(S)
    loc={}
    loc.update( locals())
    def eqs(vec,calc=CALC_INTEGRAL):
        N1,N2,V,Phi=vec.reshape((4,levels))
        loc.update( locals())
        locs=('N1','N2','V','Phi','s','mu','sigma')
        up={z:np.concatenate((loc[z][1:],[0])) for z in locs}
        down={z:np.concatenate(([0],loc[z][:-1])) for z in locs}

        #Correlations between N and D, to smallest nonzero order in sigma_d
        corrND=1- sigma_d**2
        corrND2=1 - 2* sigma_d**2

        #First moment
        meankeff=mean_k + epsilon*loc['mu']*down['s']*down['Phi']*down['N1'] -  up['mu']*up['s']*up['Phi']*up['N1']
        Veff=up['Phi']*up['s']*up['sigma']**2 *up['V']
        Veff+=down['Phi']*down['s']*loc['sigma']**2 *down['V']

        mean0= meankeff/( corrND - epsilon * Veff )

        #print N1,meankeff, mean_k[1],epsilon[1]*mu[1]*s[0]*N1[0]

        #Second moment
        varkeff=sigma_k**2 #rest doesn't vary!

        mvar0=varkeff+meankeff**2 + up['Phi']*up['s']*up['sigma']**2 *up['N2']
        mvar0+=epsilon**2 * down['Phi']*down['s']*loc['sigma']**2 * down['N2']


        mvar0/=( corrND2*(1+sigma_d**2) -2 * epsilon*Veff
                 + epsilon**2 *Veff**2 )
        var0=mvar0-mean0**2

        #print N1,mean0,var0,erfmom(1,mean0,var0)

        eq4=1 - V*(corrND - epsilon *Veff  )
        if calc:
            raise Exception("TODO: explicit calculation in trophic_cavity_solve")
        else:
            #Equivalent with directly computed moments of erf
            eq1=Phi-erfmom(0,mean0,var0)
            eq2=N1-erfmom(1,mean0,var0)
            eq3=N2-erfmom(2,mean0,var0)
        res= np.abs(np.concatenate( (eq1,eq2,eq3,eq4)) )**2

        zerocond=np.min(np.concatenate([N1.ravel(),Phi.ravel(),N2.ravel()])) #Things that should be positive!
        if zerocond<0:
            #print "ZEROCOND",res
            res+=np.abs(zerocond )

        return res
    #root=sopt.newton_krylov
    #root=sopt.anderson
    root=sopt.root
    res= root(eqs,x0)

    result= res.x
    N1,N2,V,Phi=result.reshape((4,levels))
    if not res.success or np.max(np.abs(eqs(result)))>10.**-4:
        trials=kwargs.pop('trials_left',20)
        if trials>0:
            return trophic_cavity_solve(S=S,mu=mu,sigma=sigma,mean_k=mean_k,
                sigma_k=sigma_k,sigma_d=sigma_d,epsilon=epsilon,x0=None ,
                CALC_INTEGRAL=CALC_INTEGRAL,trials_left =trials-1,**kwargs)
        print 'ERROR: {} {}'.format(res.message,list(res.fun.reshape((4,levels))) )
        print 'CURRENT',Phi,N1,N2,V
        print 'PARAMS: S {} mu {} sigma {} sigma_K {} epsilon {}'.format(S,mu,sigma,sigma_k,epsilon)
        N1,N2,V,Phi=np.zeros((4,levels))
    return N1,N2,V,Phi



def trophic_axis_cavity_calc(s,MU,SIGMA2,mean_k,sigma_k,sigma_d,epsilon,N1,N2,V,Phi,tau=1,debug=0):
    #Correlations between N and D, to smallest nonzero order in sigma_d
    corrND=1- sigma_d**2
    corrND2=1 - 2* sigma_d**2

    #First moment
    meankeff=mean_k + epsilon*np.dot(MU,s*Phi*N1)-  np.dot(s*Phi*N1,MU*tau)
    Veff=np.dot(SIGMA2,s*Phi*V)
    Veff+=np.dot(s*Phi*V,SIGMA2*tau**2)

    utilde= corrND - epsilon * Veff
    mean0= meankeff/utilde

    #print N1,meankeff, mean_k[1],epsilon[1]*mu[1]*s[0]*N1[0]

    #Second moment
    varkeff=sigma_k**2 #rest doesn't vary!

    mvar0=varkeff+meankeff**2 +  np.dot (N2*Phi*s,SIGMA2*tau**2) +epsilon**2 * np.dot(SIGMA2, N2*Phi*s)

    utilde2=( corrND2*(1+sigma_d**2) -2 * epsilon*Veff + epsilon**2 *Veff**2 )

    mvar0/=utilde**2 #
    var0=mvar0-mean0**2

    if debug:# and 0:
        print varkeff/var0

        #print Veff,utilde
        #print var0*utilde**2, np.dot (N2*Phi*s,SIGMA2*tau**2) +epsilon**2 * np.dot(SIGMA2, N2*Phi*s)

        #print utilde**2, utilde2 #Difference only if sigma_d !=0
        #print  varkeff, var0,#np.dot (N2*Phi*s,SIGMA2),epsilon**2 * np.dot(SIGMA2, N2*Phi*s)

    return mean0, var0, V*utilde

def trophic_axis_meanfield(s,MU,mean_k,epsilon,Phi=None,tau=1):
    if Phi is None:
        Phi=np.ones(s.shape[0])
    sPhi=s*Phi
    mat=np.diag(1/sPhi)- epsilon*MU + (MU*tau).T
    return np.dot(np.linalg.inv(mat),mean_k)/sPhi


def trophic_axis_variability(s,SIGMA2,epsilon,Phi=None):
    if Phi is None:
        Phi=np.ones(s.shape[0])
    sPhi=s*Phi
    mat=np.diag(1/sPhi)- epsilon*SIGMA2 + SIGMA2.T
    return np.dot(np.linalg.inv(mat),np.ones(s.shape[0]) )/sPhi


def trophic_axis_cavity_solve(S=[100,100],mu=1,sigma=.0,mean_k=1., sigma_k=0,
        sigma_d=0,epsilon=.9,tau=1,x0=None ,
        **kwargs):

    NTRIALS=50
    #tau = ratio of typical timescales (e.g. for metabolic)
    import scipy.optimize as sopt

    levels=len(S)
    S,mean_k,sigma_k,epsilon=[np.array(x).astype('float')  if np.array(x).shape
        else np.array([float(x) for nb in range(levels)] )
         for x in S,mean_k,sigma_k,epsilon]
    s=S/np.mean(S)

    #Making matrices
    if len(mu.shape)<2:
        MU=np.zeros( (levels,levels) )
        SIGMA2=np.zeros( (levels,levels) )
        structure=kwargs.get('structure','neighbor')
        rge=range(levels)
        MU[rge[1:],rge[:-1]]=mu[1:]
        SIGMA2[rge[1:],rge[:-1]]=sigma[1:]**2
    else:
        MU=mu
        SIGMA2=sigma**2
    #print MU
    #print SIGMA2,'\n'
    #Starting from MF soltion

    XMF=np.ones((4,levels))
    XMF[0]=np.clip(trophic_axis_meanfield(s,MU,mean_k,epsilon),0.00000000001,None)
    XMF[1]=XMF[0]**2
    if x0 is None:
        x0=np.random.random((4,levels))*XMF
        if kwargs.get('trials_left',NTRIALS)<NTRIALS*.75:
            x0+=np.random.random((4,levels))
        #for i in range(levels):
            #x0[0,i]/=(mu[i]/epsilon[i] )**i
    x0=x0.ravel()




    CALC_INTEGRAL=kwargs.get('CALC_INTEGRAL',0)
    def eqs(vec,calc=CALC_INTEGRAL):
        N1,N2,V,Phi=vec.reshape((4,levels))
        zerocond=np.min(np.concatenate([N1.ravel(),Phi.ravel(),N2.ravel()])) #Things that should be positive!
        N2=np.clip(N2,0,None)
        N1=np.clip(N1,0,None)
        Phi=np.clip(Phi,0,None)

        mean0,var0,Veq=trophic_axis_cavity_calc(s,MU,SIGMA2,mean_k,sigma_k,
            sigma_d,epsilon,N1,N2,V,Phi,tau)
        #print N1,mean0,var0,erfmom(1,mean0,var0)

        eq4=1 - Veq
        if calc:
            raise Exception("TODO: explicit calculation in trophic_cavity_solve")
        else:
            #Equivalent with directly computed moments of erf
            eq1=Phi-erfmom(0,mean0,var0)
            eq2=N1-erfmom(1,mean0,var0)
            eq3=N2-erfmom(2,mean0,var0)

        p=2
        if kwargs.get('trials_left',NTRIALS)<NTRIALS*.25:
            p=1
        res= np.abs(np.concatenate( (eq1,eq2,eq3,eq4)) )**p

        if zerocond<0:
            #print "ZEROCOND",res
            res+=np.abs(zerocond )

        return res
    #root=sopt.newton_krylov
    #root=sopt.anderson
    root=sopt.root
    res= root(eqs,x0)

    result= res.x
    N1,N2,V,Phi=result.reshape((4,levels))
    #mean0,var0,Veq=trophic_axis_cavity_calc(s,MU,SIGMA2,mean_k,sigma_k,sigma_d,epsilon,N1,N2,V,Phi,tau,debug=1)
    #print var0,N2, N2-N1**2
    if not res.success or np.max(np.abs(res.fun))>10.**-4:
        trials=kwargs.pop('trials_left',NTRIALS)
        if np.sum(res.fun**2)<np.sum(kwargs.get('best_fun',1000)**2):
            kwargs['best_fun']=res.fun
            kwargs['best_x']=res.x

        if trials>0:
            return trophic_axis_cavity_solve(S=S,mu=mu,sigma=sigma,mean_k=mean_k,
                sigma_k=sigma_k,sigma_d=sigma_d,epsilon=epsilon,x0=None ,
                trials_left =trials-1,**kwargs)


        N1,N2,V,Phi=kwargs.get('best_x',XMF).reshape((4,levels))
        print '\nERROR: {} {}'.format(res.message,list(kwargs.get('best_fun',res.fun).reshape((4,levels))) )
        print 'CURRENT',Phi,N1,N2,V
        print 'PARAMS: S {} mu {} sigma {} sigma_K {} epsilon {}'.format(S,mu,sigma,sigma_k,epsilon)
        mean0,var0,Veq=trophic_axis_cavity_calc(s,MU,SIGMA2,mean_k,sigma_k,
            sigma_d,epsilon,N1,N2,V,Phi,tau)
        print 'Mean,Var,Veq',mean0,var0,Veq
        print "MEANFIELD",XMF[0]
        #raise
        #N1,N2,V,Phi=XMF
    return N1,N2,V,Phi