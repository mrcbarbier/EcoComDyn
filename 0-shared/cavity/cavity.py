# -*- coding: utf-8 -*-
#from datatools import *

from utility import *
from numpy import sqrt,pi,exp,log,sin,cos
from scipy.special import erf,binom

#==================== BASICS ============================



def erfmom(mom,mean,var,lim=0):
    #lim = lower integration limit
    if lim!=0:
        if mom ==0:
            return erfmom(mom,mean-lim,var)
        moms=[erfmom(m,mean-lim,var)*lim**(mom-m)*binom(mom,m) for m in range(mom+1)]
        return np.sum(moms)
    if var<=0.001:
        var=0.001
    xx=mean/sqrt(2*var)

    mom0=.5 * (erf(xx)+1 )
    mom1=sqrt(var/2/pi) * exp(-xx**2)
    if mom==0:
        return mom0
    elif mom==1:
        return mean* mom0 + mom1
    elif mom==2:
        return (var+mean**2)*mom0 + mean*mom1

def cavity_w(k,delta):
    import scipy.integrate as sinteg
    res= sinteg.quad(lambda z:np.exp(-z**2/2)/np.sqrt(2*pi) *(z+delta)**k,-delta,np.inf)
    if res[1]>0.1:
        print 'WARNING: LARGE ERROR IN cavity_W'
    return res[0]

def cavity_solve(size=100,mu=1,sigma=1,sigma_k=1,gamma=1,x0=None,print_msg=1,
        CALC_INTEGRAL=0,**kwargs):
    import scipy.optimize as sopt
    NTRIALS=kwargs.get('NTRIALS',200)

    #print "SO",size,mu,sigma,sigma_k, gamma
    #if sigma==0:
        ##Mean field solution
        ##Ni= (Ki - mu avgN)/(1-mu/S)

        #sigma=0.01
        #x0=(1.,)
        #u=1-mu*1./size
        #def eq(x):
            #mean=(1-x)/u
            #var= sigma_k/u
            #eq1=phi-erfmom(0,mean,var)
            #eq2=1-erfmom(1,mean,var)
            #return eq1
        #root=sopt.root(eq,x0)
        #return q,v,h,phi

    u=(1-mu*1./size)/sigma

    avgK=kwargs.get('avgK',1)
    rowstd=kwargs.get('rowstd',0)
    #sigma=np.sqrt(max(0,sigma**2-rowstd**2/size) )
    #rowstd=np.sqrt(max(0,rowstd**2-sigma**2))
    #rowstd=0
    if x0 is None:
        x0=np.random.random(4)


    def eqs(vec,calc=CALC_INTEGRAL):
        q,v,h,phi=vec
        utilde=1./sigma-gamma*v
        avgN=1./(sigma*h + mu) * avgK
        sigma_lambda=np.sqrt(sigma_k**2 + rowstd**2 * avgN**2 )/sigma/avgN
        effvar=q+sigma_lambda**2
        correctionKA= size*(kwargs.get('corrKA',0))/sigma**2/avgN
        if correctionKA and not np.isnan(correctionKA) :
            effvar -=correctionKA #*2
            #print 'correction',correctionKA
        eq4=phi - v*utilde
        if calc:
            delta=h/np.sqrt(effvar)
            w0=cavity_w(0,delta)
            w1=cavity_w(1,delta)
            #w2=cavity_w(2,delta)
            w2=w0+delta*w1
            eq1=phi - w0
            eq2=utilde-np.sqrt(effvar)*w1
            eq3=utilde**2-effvar/q*w2
        else:
            #Equivalent with directly computed moments of erf
            mean=h/utilde
            var=effvar/utilde**2
            eq1=phi-erfmom(0,mean,var)
            eq2=1-erfmom(1,mean,var)
            eq3=q-erfmom(2,mean,var)
        res= np.array( (eq1,eq2,eq3,eq4))

        zerocond=min(q,phi,h+mu/sigma ) #Things that should be positive!
        if zerocond<0:
            res+=np.abs(zerocond )
        #print vec,res, w0,utilde
        return res
    #root=sopt.newton_krylov
    #root=sopt.anderson
    root=sopt.root
    res= root(eqs,x0)
    q,v,h,phi= res.x#np.abs(res.x)

    avgN=1./(sigma*h + mu)* avgK
    #print eqs((q,v,h,phi),0)
    if not res.success or np.max(np.abs(eqs((q,v,h,phi))))>10.**-4:
        trials=kwargs.pop('trials_left',NTRIALS)
        if trials>0:
            print 'trials',trials,x0
            return cavity_solve(size=size,mu=mu,sigma=sigma,
                sigma_k=sigma_k,gamma=gamma,x0=None ,
                trials_left =trials-1,**kwargs)

        if print_msg:
            print 'ERROR: {} {}'.format(res.message,res.fun)
            print 'PARAMS: S {} mu {} sigma {} sigma_K {} gamma {}'.format(size,mu,sigma,sigma_k,gamma)
            print 'VALS: q {} v {} h {} phi {} avgN {}'.format(q,v,h,phi,avgN)
        q=v=h=phi=0
    return q,v,h,phi



def cavity_AB(q,h,mu,sigma,sigma_k,gamma):
    avgN=1./(sigma*h + mu)
    varN=q* avgN**2

    B=(1./mu-avgN)/(varN + sigma_k**2/sigma**2)
    A=(1+gamma) (1./mu + gamma *avgN * B)
    A/=varN *(1+gamma) + sigma_k**2/sigma**2
    return avgN,varN,A,B

