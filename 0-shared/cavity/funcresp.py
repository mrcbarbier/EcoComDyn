# -*- coding: utf-8 -*-
from cavity import *


def funcresp_cavity_solve(S=100,mu=1,sigma=1,sigma_k=1,gamma=1,Nc=10,x0=None,
        CALC_INTEGRAL=0,**kwargs):
    import scipy.optimize as sopt

    import scipy.integrate as sinteg
    NTRIALS=2
    #Old method

    f=1
    S,mu,sigma,sigma_k,gamma,Nc=(float(z) for z in (S,mu,sigma,sigma_k,gamma,Nc))

    avgK=kwargs.get('avgK',1)

    rowstd=kwargs.get('rowstd',0)
    colstd=kwargs.get('colstd',0)
    #rowstd=np.sqrt(max(0,rowstd**2-sigma**2))

    def eqs(vec,calc=CALC_INTEGRAL):
        N1,N2,v,phi=vec
        Ac=mu/S*Nc/phi

        #sigmu=1+(kwargs.get('Kcol',0)-kwargs.get('rowcol',0)*phi*N1 )/(avgK-mu*phi*N1)

        def calcmom(mom,mean,var,momk=None):
                #return 0


            if momk is None:
                momk=mom
            PSVG=phi*sigma**2*v*gamma
            #PSVG=np.clip(PSVG,None,.99)
            def Kmin(z):
                return z/(1+np.abs(z/Ac ))
            #print sinteg.quad(lambda z:np.exp(-(z-mean)**2/2/var)
                #/np.sqrt(2*pi*var),-np.inf,np.inf),mean,var,N1,N2,v,phi
            #std=np.sqrt(var)
            if Nc>5*S*N1 or N1>100*avgK:
                if momk != mom:
                    res=[phi/(1-PSVG) ,0 ]
                else:
                    u=1-PSVG
                    res=[ erfmom(mom,(avgK-mean)/u, (sigma_k**2 + var)/u**2 ),0 ]
            else:
                if var <0.01:
                    res=(1./ (1- PSVG/(1+np.abs(mean/Ac) )**2)**mom *
                         erfmom(momk,avgK-Kmin(mean),sigma_k**2) )
                else:
                    res= sinteg.quad(lambda z:np.exp(-(z-mean)**2/2/var)
                        /np.sqrt(2*pi*var) / (1- PSVG/(1+np.abs(z/Ac) )**2)**mom *
                         erfmom(momk,avgK-Kmin(z),sigma_k**2),
                        -np.inf,np.inf)
                #mean-4*std,mean+4*std )
            #print mean,std,res
            if res[1]>0.1:
                print 'WARNING: LARGE ERROR IN FUNCRESP'

            return res[0]

        mean=mu*N1*phi
        var=np.clip(sigma**2*N2*phi,0.001,None)
        eq1=phi-calcmom(0,mean,var)
        eq2=phi*N1-calcmom(1,mean,var)
        eq3=phi*N2-calcmom(2,mean,var)

        eq4=v-calcmom(1,mean,var,momk=0)/phi

        #eqf=f-calcf()
        res= np.array( (eq1,eq2,eq3,eq4))

        res=np.abs(res)
        zerocond=min(N1,N2,phi) #Things that should be positive!
        if zerocond<0:
            res+=np.abs(zerocond )
        if np.isnan(res).any():
            print res, vec, mu, sigma,gamma,mean,var,S
            raise
        #print vec,res
        return res
    #root=sopt.newton_krylov
    #root=sopt.anderson


    class Tmp():
        success=False
    res=Tmp()
    tries=0
    def meanfield():
        v = avgK/(1+mu*(1-1/S) )
        fi=erfmom(0,v,sigma_k**2)
        return np.array( (v/fi,(v**2+sigma_k**2)/fi**2,1,fi) )

    while not res.success :# and np.sum(kwargs.get('best_fun',1000))>10**-:
        if tries>NTRIALS:
            #q,vv,h,fi= cavity_solve( S=S,mu=mu,sigma=sigma,print_msg=0,
                #sigma_k=sigma_k,gamma=gamma,x0=None,trials_left=1,**kwargs)
            #avgN=1./(sigma*h + mu)* avgK
            res.x= np.zeros(4)+0.0001#kwargs.get('best_x',np.zeros(4))
            res.x[0]=0
            #if avgN <=0:
                #res.x=np.zeros(4)+0.0001
                #break
            #print '\nFR-ERROR: {} {}'.format(res.message,kwargs.get('best_fun',res.fun))
            #print 'PARAMS: S {} mu {} sigma {} sigma_K {} gamma {} Nc {}'.format(S,mu,sigma,sigma_k,gamma,Nc)
            #print 'VALS: N1 {} N2 {} v {} phi {} f {}'.format(N1,N2,v,phi,f)
            #print 'x0',x0#, meanfield()
            #print 'other', avgN/fi, avgN**2*(1 +q )/fi**2, vv,fi
            break

        if tries >0:
            #print 'failing at life',S,mu,sigma,sigma_k,gamma,Nc
            x0=funcresp_cavity_solve(S=S,mu=mu,sigma=sigma/1.1,
                sigma_k=sigma_k/1.01,gamma=gamma,x0=None ,Nc=Nc,**kwargs)[:-1]
            if x0[0]==0:
                res.x=x0
                break
        else:
            x0=meanfield()

        if sigma<=0.01:
            res.x=x0
            break
        root=sopt.root
        res= root(eqs,x0)
        tries+=1
        if np.sum(res.fun**2)<np.sum(kwargs.get('best_fun',1000)**2):
            kwargs['best_fun']=res.fun
            kwargs['best_x']=res.x

    N1,N2,v,phi= res.x
    #print res.x
    Ac=mu/S*Nc/phi
    mean=phi*mu*N1
    var=np.clip(phi*sigma**2*N2,0.001,None)
    f=sinteg.quad(lambda z:np.exp(-(z-np.abs(mean) )**2/2/var)/np.sqrt(2*pi*var),np.abs(Ac),np.inf)[0]

    return N1,N2,v,phi,f


def funcresp_cavity_solve_old(S=100,mu=1,sigma=1,sigma_k=1,gamma=1,Nc=10,x0=None,
        CALC_INTEGRAL=0,**kwargs):
    import scipy.optimize as sopt
    NTRIALS=200


    S,mu,sigma,sigma_k,gamma,Nc=(float(z) for z in (S,mu,sigma,sigma_k,gamma,Nc))
    avgK=kwargs.get('avgK',1)
    if x0 is None:
        x0=np.random.random(5)
        if mu<0:
            base=avgK-mu *Nc/S
            x0[0]*=base*2
            x0[1]*=base**2*2
            x0[3]=.95
    #print x0
    rowstd=kwargs.get('rowstd',0)
    colstd=kwargs.get('colstd',0)
    #rowstd=np.sqrt(max(0,rowstd**2-sigma**2))

    def eqs(vec,calc=CALC_INTEGRAL):
        N1,N2,v,phi,f=vec

        sigmu=1+(kwargs.get('Kcol',0)-kwargs.get('rowcol',0)*phi*N1 )/(avgK-mu*phi*N1)

        utilde=1 -(v*phi* (gamma* sigma**2 )  )*(1-f) *sigmu  #


        effvar=(1-f)**2 *phi*( (sigma **2  + colstd**2) * N2 +rowstd**2 *phi*N1**2) +sigma_k**2
        correctionKA= S*(kwargs.get('corrKA',0))
        if correctionKA and not np.isnan(correctionKA) :
            effvar -=correctionKA #*2
            #print 'correction',correctionKA
        eq4=1 - v*utilde

        mean=(avgK-mu * min(phi*N1, (phi*N1*(1-f) + f * Nc/S  ))  )/utilde
        #print mean, N1

        var=effvar/utilde**2
        eq1=phi-erfmom(0,mean,var)
        eq2=phi*N1-erfmom(1,mean,var)
        eq3=phi*N2-erfmom(2,mean,var)

        def calcf():
            if Nc>10000*N1:
                return 0
            mean=np.abs(mu)*phi*N1
            thresh=np.abs(mu)*Nc/S
            var=sigma**2 * phi * N2
            var=np.clip(var,0.000001,None)
            if 0:
                return erfmom(0,mean-thresh,var )
            else:
                import scipy.integrate as sinteg

                norm=erfmom(0,mean,var )
                if norm < 10**-10:
                    return 0
                res= sinteg.quad(lambda z:np.exp(-(z-mean)**2/2/var)
                    /np.sqrt(2*pi*var)/(1+thresh/z ),
                    0,np.inf)/norm
                if res[1]>0.1:
                    print 'WARNING: LARGE ERROR IN FUNCRESP'
                return res[0]

        eqf=f-calcf()
        res= np.array( (eq1,eq2,eq3,eq4,eqf))

        res=np.abs(res)
        zerocond=min(N1,N2,phi,f) #Things that should be positive!
        if zerocond<0:
            res+=np.abs(zerocond )
        if np.isnan(res).any():
            print res, vec, mu, sigma,gamma,mean,var,utilde,S
            raise
        return res
    #root=sopt.newton_krylov
    #root=sopt.anderson
    root=sopt.root
    res= root(eqs,x0)
    N1,N2,v,phi,f= res.x
    #raise
    #print f*mu*Nc/S,N1,f,eqs(res.x),res.fun
    #raise
    if not res.success or np.max(np.abs(res.fun))>10.**-4:
        trials=kwargs.pop('trials_left',NTRIALS)
        if trials>0:
            if np.sum(res.fun**2)<np.sum(kwargs.get('best_fun',1000)**2):
                kwargs['best_fun']=res.fun
                kwargs['best_x']=res.x
            return funcresp_cavity_solve(S=S,mu=mu,sigma=sigma,
                sigma_k=sigma_k,gamma=gamma,x0=None ,Nc=Nc,
                trials_left =trials-1,**kwargs)

        N1,N2,v,phi,f=kwargs.get('best_x',np.zeros(5))
        print '\nERROR: {} {}'.format(res.message,kwargs.get('best_fun',res.fun))

        print 'PARAMS: S {} mu {} sigma {} sigma_K {} gamma {} Nc {}'.format(S,mu,sigma,sigma_k,gamma,Nc)
        print 'VALS: N1 {} N2 {} v {} phi {} f {}'.format(N1,N2,v,phi,f)
    return N1,N2,v,phi,f
