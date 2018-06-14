# -*- coding: utf-8 -*-
from cavity import *
import scipy.optimize as sopt
from itertools import product as iproduct

def erfmom(mom,mean,var):
    thresh=0.0000001
    from numpy import sqrt,pi,exp
    from scipy.special import erf

    lowvar=(var<thresh)
    var=np.clip(var,thresh,None)
    present=(mean>0).astype('float')
    xx=np.clip(mean/sqrt(2*var),-np.sqrt(40),np.sqrt(40))

    mom0=.5 * (erf(xx)+1 )
    mom1=sqrt(var/2/pi) * exp(-xx**2)

    res=0
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

def group_cavity_calc(S,MU,SIGMA2,GAMMA,mean_k,sigma_k,sigma_d,N1,N2,V,Phi,tau=1,debug=0,SIGROW2=None):
    #Correlations between N and D, to smallest nonzero order in sigma_d
    s=S/np.mean(S)


    corrND=1- sigma_d**2
    corrND2=1 - 2* sigma_d**2

    #First moment
    meankeff=mean_k -np.dot(MU,s*Phi*N1)
    Veff=np.dot(np.sqrt(SIGMA2*SIGMA2.T)*GAMMA,s*Phi*V)

    utilde= corrND -   Veff
    mean0= meankeff/utilde

    #print N1,meankeff, mean_k[1],epsilon[1]*mu[1]*s[0]*N1[0]

    #Second moment
    varkeff=sigma_k**2  #rest doesn't vary!
    if SIGROW2 is 0:
        SIGCOL2=SIGMA2
    else:
        SIGCOL2=SIGMA2-SIGROW2/np.mean(S)
        varkeff+= np.dot(SIGROW2, N1**2 * Phi**2 * s**2)


    mvar0=varkeff+meankeff**2 +   np.dot(SIGCOL2, N2*Phi*s)

    #utilde2=( corrND2*(1+sigma_d**2) -2  *Veff + Veff**2

    mvar0/=utilde**2
    var0=mvar0-mean0**2

    if debug:
        print 'groups.py DEBUG',varkeff/var0

        #print Veff,utilde
        #print var0*utilde**2, np.dot (N2*Phi*s,SIGMA2*tau**2) +epsilon**2 * np.dot(SIGMA2, N2*Phi*s)

        #print utilde**2, utilde2 #Difference only if sigma_d !=0
        #print  varkeff, var0,#np.dot (N2*Phi*s,SIGMA2),epsilon**2 * np.dot(SIGMA2, N2*Phi*s)

    return mean0, var0, V*utilde

def group_meanfield(s,MU,mean_k,Phi=None,sigma_k=None,Ninit =None):
    if Phi is None:
        Phi=np.ones(mean_k.shape)
        if 0:
            if not Ninit is None:
                Phi = erfmom(0, Ninit, sigma_k ** 2)
            elif sigma_k is None:
                Phi=np.ones(s.shape[0])
            else:
                Phi=erfmom(0,mean_k,sigma_k**2)
        Phi=np.clip(Phi,0.01,None)
    #mat=np.diag(1/sPhi)+ MU
    root=sopt.root
    phi=Phi
    import scipy.linalg as la
    if Ninit is None:
        if not sigma_k is None:
            Ninit=group_meanfield(s,MU,mean_k)[0]
        else:
            Ninit=np.dot(la.inv(np.eye(s.shape[0])+(s*Phi*MU.T).T ),mean_k)
            if (Ninit<0).any():
                Ninit=mean_k - np.dot(MU,s*Phi* np.clip(mean_k,0.1,None) )
            #print Ninit
            Ninit=np.clip(Ninit,0.1,None)


    #return np.dot(np.linalg.inv(mat),mean_k)/sPhi
    traceback=[]
    def eqs(x):
        traceback.append(x)
        if x.shape==mean_k.shape:
            N1=x
            Phi=phi
            roof = np.max(mean_k)*10**6
        else:
            N1,Phi=x.reshape(2,x.shape[0]/2 )
            roof=np.max( np.sqrt(mean_k**2+sigma_k**2) )*10**6

        eqN1=(N1 -np.clip(mean_k - np.dot(MU,np.clip(s*Phi*N1,0,None) ) ,0,roof))#/np.max(N1)

        if sigma_k is None:
            res= eqN1
        else:
            eqPhi= Phi- erfmom(0,N1,sigma_k**2)
            res=np.concatenate([eqN1,eqPhi] )
        return res
    if sigma_k is None:
        res=root(eqs,Ninit  )
        N1=res.x
    else:
        res=root(eqs, np.concatenate([Ninit,Phi] ),options={'maxfev':10000,
            #'factor':0.1
            }
            )
        N1,Phi=res.x.reshape(2,res.x.shape[0]/2 )
    #print '\n =>',s, MU, mean_k, res.x, res.fun
    #print res.x,eqs(res.x),mean_k,sPhi,mean_k + sPhi*s np.dot(MU,res.x )
    #raise
    if not res.success and np.max(np.abs(res.fun) )>10**-1:
        if 0 and np.max(np.abs(la.eigvals(-s*Phi*MU.T)))>1:
            return np.zeros((2,s.shape[0]) )
        print 'Meanfield failure',np.max(np.abs(res.fun))

        handled=1
        if len(s)>20 and False:
            #FOR SOME REASON, THIS BASICALLY NEVER WORKS
            print "Bootstrap using Lotka Volterra"
            N1=N0=np.clip(mean_k.copy(),.1,None)
            import scipy.integrate as scint
            def LOCDX(Xs):
                dx= np.clip(Xs,0,None)*(mean_k - Xs - np.dot(MU,s*Xs))/np.max(Xs)
                dx[Xs<10**-10]= np.clip(dx[Xs<10**-10],0,None)
                return dx


            integrator=scint.ode(LOCDX).set_integrator('dop853',nsteps=100000)
            integrator.set_initial_value(N1,0)
            N1=integrator.integrate(1000)
            Phi=(N1>10**-10* np.mean(N1) )
            if np.max(N1)> np.sum(mean_k[mean_k>0] )*100:
                N1[:]=0
                Phi[:]=0
            if sigma_k is None:
                bs=np.max(np.abs(eqs( N1 )))
            else:
                bs=np.max(np.abs(eqs(  np.concatenate([N1,Phi]) )))
            print "     bootstrap: ERROR",bs
            if bs>np.max(N1):
                handled=1

        while not handled:
            code_debugger()
    return np.clip(N1,10**-15,None),np.clip(Phi,10**-15,None)


def group_variability(s,SIGMA2,GAMMA,Phi=None):
    if Phi is None:
        Phi=np.ones(s.shape[0])
    sPhi=s*Phi
    mat=np.diag(1/sPhi) + SIGMA2
    return np.dot(np.linalg.inv(mat),np.ones(s.shape[0]) )/sPhi

def group_press(s,mu,Phi,V):
    from scipy.linalg import inv
    alive=np.where(Phi>0.001)[0]
    press=np.zeros(( V.shape[0],V.shape[0]) )
    if not len(alive):
        return press
    press[np.ix_(alive,alive)]= inv(np.dot(np.diag(( s * V*Phi)[alive] ),
                            np.eye(len(alive)) + np.dot(mu[np.ix_(alive,alive)], np.diag((s * V*Phi)[alive] ) )))
    return press

def group_chi2(s, SIGMA2,GAMMA,Phi,V):
    diag= ( 1- np.dot(np.sqrt(SIGMA2*SIGMA2.T)*GAMMA,s*Phi*V) )**2
    mat=SIGMA2*SIGMA2.T*s*Phi
    if sum(diag.shape)>2:
        assert mat[1,2]==SIGMA2[1,2]*SIGMA2[2,1]*s[2]*Phi[2]
    import numpy.linalg as la
    return la.inv(np.diag(diag)-mat)


def root_helper(eqs,x0,bounds=None,tol=0.01):
    import numpy.linalg as la,scipy.optimize as sopt
    def fun(x):
        return la.norm(eqs(x))
    res=sopt.minimize(fun, x0,bounds=bounds,tol=tol )
    if res.success or fun(res.x)<fun(x0):
        return res.x
    return x0


def group_cavity_solve(S=(100,100),mu=1,sigma=.0,gamma=0,mean_k=1., sigma_k=0,
        sigma_d=0,tau=1,x0=None ,NTRIALS=3,EPS=10**-15,logmode=0, ftol=0.003,
        SIGMA_THRESH=0.01, #Threshold for recursing
        varmode=1,sigrow=None,
        **kwargs):
    USE_RECURS=kwargs.get('USE_RECURS',1)
    USE_ROOT_HELPER=kwargs.get('USE_ROOT_HELPER',1)

    #tau = ratio of typical timescales (e.g. for metabolic)


    levels=len(S)
    S,mean_k,sigma_k=[np.array(x).astype('float')  if np.array(x).shape
        else np.array([float(x) for nb in range(levels)] )
         for x in S,mean_k,sigma_k]
    s=S/np.mean(S)

    #Making matrices
    if len(mu.shape)<2:
        MU=(np.ones( (levels,levels) )*mu).T
    else:
        MU=mu
    if len(sigma.shape)<2:
        SIGMA2=(np.ones( (levels,levels) )*sigma**2).T
    else:
        SIGMA2=sigma**2
    if len(gamma.shape)<2:
        GAMMA=(np.ones( (levels,levels) )*gamma).T
    else:
        GAMMA=gamma#sigma**2

    if not sigrow is None:
        SIGROW2= sigrow**2
    else:
        SIGROW2=0

    if 0 and (MU<0).any():
        test=-(s*MU.T)
        import numpy.linalg as la
        if np.max(np.abs(la.eigvals(test)))>1:
            return np.zeros((4,levels) )+EPS



    CALC_INTEGRAL=kwargs.get('CALC_INTEGRAL',0)
    traceback=[]
    def eqs(vec,calc=CALC_INTEGRAL,logmode=logmode):
        if logmode:
            logvec=vec
            vec=np.exp(np.clip(vec,None,30)  )
        N1,N2,V,Phi=vec.reshape((4,levels))
        if varmode:
            N2=N2*N1**2
        if logmode:
            # Phi=np.clip(Phi,0,1)
            pass
        else:
            zerocond = np.min(np.concatenate([N1.ravel(), Phi.ravel(), N2.ravel()]))  # Things that should be positive!

        if 0:
            traceback.append(vec)
        #N2=np.clip(N2,0,None)
        #N1=np.clip(N1,0,None)
        #Phi=np.clip(Phi,0,None)

        mean0,var0,Veq=group_cavity_calc(S,MU,SIGMA2,GAMMA,mean_k,sigma_k,
            sigma_d,N1,N2,V,Phi,tau,SIGROW2=SIGROW2)
        #print N1,mean0,var0,erfmom(1,mean0,var0)

        eq4=1 - Veq
        if calc:
            raise Exception("TODO: explicit calculation in group_cavity_solve not done yet")
        else:
            #Equivalent with directly computed moments of erf

            if logmode:
                flor=EPS
            else:
                flor=0
            roof=np.max( np.sqrt( mean_k**2 + sigma_k**2 ))*10**7
            theophi=np.clip(erfmom(0,mean0,var0),flor,1)
            theoN1=np.clip(erfmom(1,mean0,var0),flor,roof)
            theoN2=np.clip(erfmom(2,mean0,var0),flor,roof**2 )
            #ref[theophi<10**-5]=0


            def refscale(x):
                if 0:
                    return 1
                else:
                    return np.mean(x)

            if  not logmode:
                #Putting all the equations on the same scale, while regularizing zero denominators
                # tmp=refscale(Phi)
                # eq1=1 -(Phi+tmp)/(theophi+tmp)
                eq1 = Phi-theophi
                phi=np.clip(theophi,0.001,None)
                tmp=refscale(theoN1)
                # eq2=1-(np.abs(Phi)*N1  +tmp)/(theoN1+tmp)
                eq2 = (phi*N1-theoN1)/tmp
                tmp=refscale(theoN2)
                # eq3=1-(np.abs(Phi)*N2 +tmp)/(theoN2+tmp)
                eq3 = (phi*N2-theoN2)/tmp
            else:
                if 0:
                    eq1=(Phi-theophi)/max(theophi)
                    eq2= (Phi*N1- theoN1 )/max(theoN1)
                    eq3=( Phi*N2-theoN2 )/max(theoN2)
                else:
                    logN1,logN2,logV,logPhi=logvec.reshape((4,levels))
                    tmp = refscale(theophi)
                    eq1=1 -(Phi+tmp)/(theophi+tmp)
                    tmp = refscale(np.log(theoN1))
                    eq2=1-(logN1  +tmp)/(np.log(theoN1)-np.log(theophi)+tmp)
                    tmp = refscale(np.log(theoN2))
                    eq3=1-(logN2 +tmp)/(np.log(theoN2)-np.log(theophi)+tmp)


        res= np.concatenate( (eq2,eq3,eq4,eq1))
        # if logmode:
            # res*=max(1,1+np.max(vec)/20.)

        return res
    #root=sopt.newton_krylov
    #root=sopt.anderson


    tries=0
    class Tmp():
        success=False
        fun =1000

    res=Tmp()
    trials_left=kwargs.pop('trials_left',NTRIALS)
    recurs_depth=kwargs.pop('recurs_depth',0)
    X0=None
    XMF=None
    while not res.success and np.max(np.abs(res.fun))>ftol and trials_left -tries >0 :
        newx0=0
        if (recurs_depth>5 or  tries >1) and np.max(sigma)>SIGMA_THRESH and USE_RECURS:
            kw={}
            kw.update(kwargs)
            #print 'SIGMA', sigma
            #print kwargs.get('trials_left')
            factor=max(1.03, np.exp( np.log(np.max(sigma)/0.01)/10. ) )
            print '====> RECURSE FROM', np.max(sigma**2), 'TRIES LEFT',trials_left - tries, 'DEPTH',recurs_depth
            x0=np.concatenate(group_cavity_solve(S=S,mu=mu,sigma=sigma/factor,mean_k=mean_k,gamma=gamma,
                sigma_k=sigma_k,sigma_d=sigma_d,x0=None ,logmode=False,
                NTRIALS=2,recurs_depth=recurs_depth+1,**kwargs)[:4] )
            tries=trials_left
            newx0=1

        elif x0 is None:
            #Starting from Mean Field solution

            XMF=np.ones((4,levels))
            Phi=erfmom(0,mean_k,sigma_k**2)
            Phi[Phi==0]=1
            if not kwargs.get('N0',None) is None:
                print 'INITIALIZING WITH EMPIRICAL Ns'
            XMF[0],XMF[3]=group_meanfield(s,MU,mean_k,sigma_k=sigma_k,Ninit=kwargs.get('N0',None))
            if ( (XMF[0],XMF[3])<=10**-14+np.zeros((2,levels))).all():
                print "ZEROS FROM MEANFIELD"
                return np.zeros((4,levels))+EPS
            XMF[1]=XMF[0]**2
            XMF[2]=1
            #XMF[3]=np.clip(Phi,0.01,1)
            mean0,var0,Veq=group_cavity_calc(S,MU,SIGMA2,GAMMA,mean_k,sigma_k,
                sigma_d,XMF[0],XMF[1],XMF[2],Phi,tau,SIGROW2=SIGROW2)
            XMF[3]=np.clip(np.minimum(XMF[3],erfmom(0,mean0,var0)),0.01,1)

            x0=XMF
            x0=x0.ravel()
            newx0=1


        if varmode and newx0:
            x0[levels:2*levels]=np.clip(x0[levels:2*levels]/x0[:levels]**2,0,100)
            x0[levels:2 * levels][x0[:levels]==0]=1

        if logmode:
            X0=np.log(x0+EPS)
        else:
            X0=x0
        funinit=np.abs(eqs(X0))
        if np.max(np.abs(funinit ))< ftol/10:
            res.x=x0
            res.success=True
            res.fun = funinit
        if XMF is None:
            print "    x0:",np.sum(funinit),"maxerr",np.max(funinit)


        root=sopt.root
        if logmode:
            bounds=None
        else:
            bounds=[(0,None) for i in range(3*levels)]+ [ (0,1) for i in range(levels) ]

        methods = ['hybr']
        if np.max(sigma)<SIGMA_THRESH:
            #Last chance, no recursion, so better try methods
            if logmode:
                methods=['Krylov','df-sane','hybr']
            else:
                methods=['hybr','Krylov','df-sane'] #df-sane never seems to work better than hybr; Krylov rarely (but very slow)

        if ( tries>0 or recurs_depth>0 ) and USE_ROOT_HELPER:
            print '    Using root helper'
            X1=root_helper(eqs,X0,bounds=bounds,tol=ftol)
            import numpy.linalg as la
            nm=la.norm(X1-X0)
            mx=np.argmax(np.abs(X1-X0))
            print '    Root helper pushed',nm, ['{}_{}'.format(i,j) for i,j in iproduct(['N1','N2','V','Phi'],[str(z) for z in range(levels) ] ) ][mx], X0[mx],'->',X1[mx]
            X0=X1
        while not res.success and methods:
            method=methods.pop(0)
            try:
                res= root(eqs,X0,method=method,tol=ftol )
            except Exception as e:
                print e
            if not res.success:
                print '    Root finding failure:', method
            # if not methods and not logmode:
            #     methods=['Krylov']
            #     logmode=1
        tries+=1


    print MU, SIGMA2, GAMMA, SIGROW2
    if not res.success:
        code_debugger()

    result= np.clip(res.x,EPS,None)
    if logmode:
        result=np.exp(np.clip( result,np.log(EPS),42  ))
    #print eqs(XMF.ravel()),res.fun
    success= 1
    error=np.max(np.abs(res.fun))

    if res.success and error>ftol:
        res.success=False
        print '     WARNING: groups.py claims success but error={}>ftol, REJECTED!'.format(error)

    if not res.success and error>ftol:
        success=0
        DEBUG=kwargs.get('DEBUG',False)
        if DEBUG:
            print '\nERROR: {}'.format(res.message)#,list(kwargs.get('best_fun',res.fun).reshape((4,levels))) )
            x0error=eqs(X0)
            x0better= (np.abs(res.fun)>np.abs(x0error) )
            labels=['{}_{}'.format(i,j) for i,j in iproduct(['N1','N2','V','Phi'],[str(z) for z in range(levels) ] ) ]
            ERRORS= ['{}: {} ({}) -- prev {} ({})'.format(*i) for i in zip(labels, result,res.fun,x0,x0error) if np.abs(i[2])>ftol]
            nberr=15
            print ERRORS[:nberr],
            if len(ERRORS)>nberr:
                print 'and {} others'.format(len(ERRORS)-nberr)
            else:
                print ''
            print 'RECURSION DEPTH',recurs_depth, 'NTRIALS',trials_left, "SIGMA2",np.max(sigma**2)
            result[x0better]=x0[x0better]
            error=min(error, np.max(np.abs(x0error)) )

        handled= not DEBUG
        while not handled:
            N1,N2,V,Phi=result.reshape((4,levels))
            if varmode:
                N2=N2*N1**2
            handled=1
            # code_debugger()
    else:
        if tries>1:
            print '    Root finding: finally, success'

        #mean0,var0,Veq=group_cavity_calc(s,MU,SIGMA2,GAMMA,mean_k,sigma_k,
            #sigma_d,XMF[0],XMF[1],XMF[2],XMF[3],tau)
        #print mean0
        #print np.sqrt(var0)
        #print np.clip(erfmom(0,mean0,var0),0,1)
        #print Veq
        #Veff=np.dot(np.sqrt(SIGMA2*SIGMA2.T),s*Phi*V)
        #print 1-Veff
    #else:
        #print res.fun


    N1,N2,V,Phi=result.reshape((4,levels))
    N1[Phi<EPS]=EPS
    N2[Phi<EPS]=EPS
    if varmode:
        N2=N2*N1**2
    return N1,N2,V,Phi,success/error