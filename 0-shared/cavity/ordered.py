# -*- coding: utf-8 -*-
# -*- coding: utf-8 -*-
from cavity import *

import scipy.integrate as sinteg

import scipy.optimize as sopt

def obs_get_linear(x0=None,**kwargs):
    if x0 is None:
        x0=np.random.random(7)
        if mu<0:
            base=K#-mu *Nc/S
            x0[0]*=base*2
            x0[1]*=base**2*2
            x0[2]=x0[0]
            x0[3]=x0[1]
            x0[6]=.95

    refs={}

    def integr(mom,mean,var,exp=None,fnc=None):
        if fnc is None:
            fnc=lambda e:1
        if exp is None:
            exp=0
        res= sinteg.quad(lambda z:erfmom(mom,mean(z),var(z) ) * z**exp * fnc(z),0,1)
        if res[1]>0.1:
            print 'WARNING: LARGE ERROR IN ORDERED'
        return res[0]#*(1+exp)

    def utilde(x):
        pv=refs['pv']
        pvx=refs['pvx']
        return 1 -sigma**2 * ( (gamma + c_gamma *x) * pv - c_gamma *pvx  )
    def mean(x):
        pN1=refs['pN1']
        pN1x=refs['pN1x']
        return (K +c_K *x -(mu +c_mu*x) * pN1 + c_mu*pN1x  )  /utilde(x)
    def effvar(x):
        pN2=refs['pN2']
        pN2x=refs['pN2x']
        return(zeta**2 + c_zeta *x) + (sigma**2 + c_sigma*x)*pN2 - c_sigma*pN2x
    def var(x):
        return effvar(x)/utilde(x)**2

    def eqs(vec):
        pN1,pN2,pN1x,pN2x,pv,pvx,phi=vec
        refs.update({i:j for i,j in locals().iteritems() if i[0]=='p'} )
        eq1=phi-integr(0,mean,var)
        eq2=pN1-integr(1,mean,var)
        eq3=pN2-integr(2,mean,var)
        eq4=pv- integr(0,mean,var,fnc=lambda x: 1./utilde(x) )

        eq2x=pN1x-integr(1,mean,var,1)
        eq3x=pN2x-integr(2,mean,var,1)
        eq4x=pvx- integr(0,mean,var,1,fnc=lambda x: 1./utilde(x) )


        res= np.array( (eq1,eq2,eq3,eq4,eq2x,eq3x,eq4x))

        res=np.abs(res)
        zerocond=min(pN1,pN2,phi,pN1x,pN2x) #Things that should be positive!
        if zerocond<0:
            res+=np.abs(zerocond )
        if np.isnan(res).any():
            print res, vec, mu, sigma,gamma,mean,var,utilde,S
            raise
        return res

    root=sopt.root
    res= root(eqs,x0)
    pN1,pN2,pN1x,pN2x,pv,pvx,phi= res.x
    refs.update({i:j for i,j in locals().iteritems() if i[0]=='p'} )
    xmean=sinteg.quad(lambda z:erfmom(0,mean(z),var(z) ) * z,0,1)[0]/phi
    return res,pN1,pN2,pN1x,pN2x,pv,pvx,phi,xmean


def obs_get_exponential(x0=None,**kwargs):
    if x0 is None:
        x0=np.random.random(4)
        if mu<0:
            base=K#-mu *Nc/S
            x0[0]*=base*2
            x0[1]*=base**2*2
            x0[3]=.95

    refs={}

    def integr(mom,mean,var,exp=None,fnc=None):
        if fnc is None:
            fnc=lambda e:1
        if exp is None:
            exp=1
        res= sinteg.quad(lambda z:erfmom(mom,mean(z),var(z) ) * exp**z * fnc(z),0,1)
        if res[1]>0.1:
            print 'WARNING: LARGE ERROR IN ORDERED'
        return res[0]

    def utilde(x):
        pvx=refs['pvx']
        return 1 -sigma**2 * gamma  *pvx
    def mean(x):
        pN1x=refs['pN1x']
        return (K * c_K **x -mu  * c_mu**x * pN1x  )  /utilde(x)
    def var(x):
        pN2x=refs['pN2x']
        return (zeta**2 * c_zeta **x + sigma**2 *c_sigma**x *pN2x ) /utilde(x)**2
    def eqs(vec):
        if np.isnan(vec).any():
            return np.zeros(vec.shape)
        pN1x,pN2x,pvx,phi=vec
        refs['pN1x']=pN1x
        refs['pN2x']=pN2x
        refs['pvx']=pvx

        eq1=phi-integr(0,mean,var)
        eq2x=pN1x-integr(1,mean,var,1./c_mu)
        eq3x=pN2x-integr(2,mean,var,1./c_sigma)
        eq4x=pvx- integr(0,mean,var,c_gamma,fnc=lambda x: 1./utilde(x) )


        res= np.array( (eq1,eq2x,eq3x,eq4x))

        res=np.abs(res)
        zerocond=min(phi,pN1x,pN2x) #Things that should be positive!
        if zerocond<0:
            res+=np.abs(zerocond )
        return res

    root=sopt.root
    res= root(eqs,x0)
    if np.isnan(res.x).any():
        code_debugger()
        # print res, vec, mu, sigma,gamma,mean,var,utilde,S
        raise
    pN1x,pN2x,pvx,phi= res.x
    refs['pN1x']=pN1x
    refs['pN2x']=pN2x
    refs['pvx']=pvx

    pN1=integr(1,mean,var)
    pN2=integr(2,mean,var)
    pv= integr(0,mean,var,fnc=lambda x: 1./utilde(x) )
    xmean=sinteg.quad(lambda z:erfmom(0,mean(z),var(z) ) * z,0,1)[0]/phi

    return res,pN1,pN2,pN1x,pN2x,pv,pvx,phi,xmean

def ordered_cavity_solve(x0=None,**kwargs):
    NTRIALS=200

    mode=kwargs.get('mode','linear')
    st=globals()
    S=float(kwargs['S'])
    for z in ('mu','sigma','gamma','K','zeta'):
        st[z]=float(kwargs[z])
        st['c_'+z]=float(kwargs['c_'+z])


    if mode =='linear':
        res,pN1,pN2,pN1x,pN2x,pv,pvx,phi,xmean=obs_get_linear(**kwargs)

    elif mode=='exponential':
        res,pN1,pN2,pN1x,pN2x,pv,pvx,phi,xmean=obs_get_exponential(**kwargs)


    if ( np.max(np.abs(res.fun))>10.**-4 ):
        trials=kwargs.pop('trials_left',NTRIALS)
        #print 'REDOING,',np.max(np.abs(res.fun))
        if trials>0:
            if np.sum(res.fun**2)<np.sum(kwargs.get('best_fun',1000)**2):
                kwargs['best_fun']=res.fun
                kwargs['best_x']=pN1,pN2,pN1x,pN2x,pv,pvx,phi,xmean
            return ordered_cavity_solve(x0=None ,
                trials_left =trials-1,**kwargs)

        pN1,pN2,pN1x,pN2x,pv,pvx,phi,xmean=kwargs.get('best_x',np.zeros(7))
        print '\nERROR: {} {}'.format(res.message,kwargs.get('best_fun',res.fun))

        print 'PARAMS: S {} mu {} sigma {} zeta {} gamma {}'.format(S,mu,sigma,zeta,gamma)
        print 'VALS: N1 {} N2 {} v {} phi {} '.format(pN1/phi,pN2/phi,pv/phi,phi)

    N1,N2,N1x,N2x,v,vx=[z/phi for z in (pN1,pN2,pN1x,pN2x,pv,pvx)]
    #print 'TRIALS',NTRIALS-kwargs.get('trials_left',NTRIALS)
    return N1,N2,N1x,N2x,v,vx,phi,xmean
