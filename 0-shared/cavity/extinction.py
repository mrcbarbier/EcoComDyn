# -*- coding: utf-8 -*-
from cavity import *

def x_integral(mom,xmean,xvar,Nmean,Nvar,phi,DIRECT_CALC=False):
    import scipy.integrate as sinteg

    def xfunc(mom,N):

        return erfmom(mom,xmean,xvar,lim=-N)

    if DIRECT_CALC:
        res= sinteg.quad(lambda N:np.exp(-(N-Nmean)**2/2/Nvar)/np.sqrt(2*pi*Nvar)*
            sinteg.quad(lambda x:np.exp(-(x-xmean)**2/2/xvar)/np.sqrt(2*pi*xvar) *x**mom,
            -N,np.inf)[0],0,np.inf)/phi
    else:
        res= sinteg.quad(lambda N:np.exp(-(N-Nmean)**2/2/Nvar)/np.sqrt(2*pi*Nvar)*
            xfunc(mom,N),0,np.inf)

    if res[1]>0.1:
        print 'WARNING: LARGE ERROR IN CAVITY EXTINCTION x_integral'
    return res[0]/phi

def cavity_extinction_simple(N0=1,S=100,mu=1,sigma=1,gamma=1,v=1,phi=1,avgN=1,varN=1,**kwargs):

    u=(1-mu*1./S)

    utilde=u-phi*sigma**2*gamma*v

    xmean=(N0*mu/S )/utilde
    xmean /= 1+mu* phi/utilde

    xvar=( mu**2 *((N0/S)**2 + phi**2 * xmean**2
        - 2 * N0 * phi * xmean/S ) )/utilde**2 - xmean**2
    xvar/= 1-phi*sigma**2 /utilde**2

    return xmean,xvar

def cavity_extinction(N0=1,S=100,mu=1,sigma=1,gamma=1,v=1,phi=1,avgN=1,varN=1,x0=None,**kwargs):
    import scipy.optimize as sopt

    if phi==0 or varN==0:
        return 0,0,0

    base=cavity_extinction_simple(N0=N0,S=S,mu=mu,sigma=sigma,gamma=gamma,v=v,
                phi=phi,avgN=avgN,varN=varN,**kwargs)
    return (max(0,1-x_integral(0,base[0],base[1],avgN,varN,phi)),)+base

    u=(1-mu*1./S)

    if x0 is None:
        x0=np.array(base)*(1+.3*(np.random.random(2)-.5) )

    utilde=u-phi*sigma**2*gamma*v

    #print '\nyay',avgN,varN,phi,erfmom(2,0,2,lim=0),erfmom(2,0,2,lim=-100)#, x_integral(1, 1.,0.01, 1.,1.,1.)


    def eqs(vec):
        xmean,xvar=vec

        effmean=(N0*mu/S - mu* phi * xmean)/utilde
        effvar=(phi*sigma**2*xvar + mu**2 *((N0/S)**2 + phi**2 * xmean**2
            - 2 * N0 * phi * xmean/S ) )/utilde**2 - effmean**2
        effvar=max(0,effvar)

        eq1=xmean- x_integral(1,effmean,effvar,avgN,varN,phi)
        eq2=xvar - x_integral(2,effmean,effvar,avgN,varN,phi)
        res= np.array( (eq1,eq2))

        zerocond=min(0,xvar ) #Things that should be positive!
        #res=res**2
        if zerocond<0:
            res+=np.abs(zerocond )
            #print 'damn',xvar
        #print vec,res, w0,utilde
        return res
    #root=sopt.newton_krylov
    #root=sopt.anderson
    root=sopt.root
    res= root(eqs,x0)
    xmean,xvar= res.x#np.abs(res.x)

    effmean=(N0*mu/S - mu* phi * xmean)/utilde
    effvar=(phi*sigma**2*xvar + mu**2 *((N0/S)**2 + phi**2 * xmean**2
            - 2 * N0 * phi * xmean/S ) )/utilde**2 - effmean**2
    #print '!',mu,N0,S,phi,xmean,utilde,effmean,effvar
    frac=max(0,1-x_integral(0,effmean,effvar,avgN,varN,phi))
    if not res.success or np.max(np.abs(eqs((xmean,xvar))))>10.**-4:
        trials=kwargs.pop('trials_left',20)
        if trials>0:
            return cavity_extinction(N0=N0,S=S,mu=mu,sigma=sigma,gamma=gamma,v=v,
                phi=phi,avgN=avgN,varN=varN,trials_left=trials-1,**kwargs)
        print 'EXTINCT ERROR: {} {}'.format(res.message,res.fun)

        if np.max(np.abs(eqs((xmean,xvar))))>10.**-4:
            frac,xmean,xvar=0,0,0
    else:
        print 'EXTINCT SUCCESS!',frac,xmean,xvar
    return frac,xmean,xvar
