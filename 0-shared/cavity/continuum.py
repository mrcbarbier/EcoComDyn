# -*- coding: utf-8 -*-

from cavity import *
import scipy.integrate as sinteg
import scipy.interpolate as sinterp


def continuum_quad(xy,*args):
    if len(args)==3:
        x=xy
        a,b,c=args
        return a+b*x+c*x**2
    elif len(args)==7:
        a,bx,cx,by,cy,bxy,cxy=args
        x=xy[0]
        y=xy[1]
        return a+(bx*x+cx*x**2)  + (bxy*x+cxy*(x-y)**2) + (by*y+cy*y**2)



def continuum_func(opt,sgn=1,mode='exp',thresh=10**-5,dim=2):

    func=continuum_quad
    if mode=='exp':
        mod=np.exp
    else:
        mod=lambda x:x

    sgn=np.array(sgn)
    if not sgn.shape:
        sgn=sgn.reshape((1,))
    if len(sgn)<len(opt):
        sgn=np.concatenate([(1,),sgn])

    if dim==2:
        return np.vectorize(lambda xs,ys:np.sum([sgn[i]*mod( np.maximum(np.log(thresh),
            func(np.array((xs,ys)),*o))) for i,o in enumerate(opt) ] ))
    else:
        return np.vectorize(lambda xs:np.sum([sgn[i]*mod( np.maximum(np.log(thresh),
            func(xs,*o))) for i,o in enumerate(opt) ] ))
    raise Exception("SHOOT")


def continuum_mf(**kwargs):
    ranks=kwargs['ranks']
    avgK=kwargs.get('avgK')(ranks)
    rank2=np.tile(ranks,(len(ranks),1) )
    mu=kwargs.get('mu')( rank2.T,rank2  )
    import scipy.linalg as la
    n1= np.clip( np.dot(la.inv(1+mu),avgK),0.1,None)

    return np.concatenate( (n1,n1**2, np.ones(len(ranks)), np.ones(len(ranks)) )  )

def continuum_cavity_solve(**kwargs):
    import scipy.optimize as sopt

    fnc=continuum_func
    for n in ('S','avgK','varK' , 'mu','sigma','gamma' ):
        if not n in kwargs:
            kwargs[n] = fnc(kwargs[n+'prm'],**dict(kwargs.get(n+'opt', () )) )
    if not 'ranks' in kwargs:
        ranks=kwargs['ranks']=np.linspace(kwargs.get('ranks_min',0),
            kwargs.get('ranks_max',1),kwargs.get('resolution',20))
    else:
        ranks=kwargs['ranks']
    resolution=len(ranks)

    def gaussian(x,mean,var):
        return np.exp(-(x-mean)**2/(2*var) )/np.sqrt(2*np.pi*var )

    def calcprm(vec,**kwargs):
        print 'tmp'
        n1,n2,v,phi=[ sinterp.interp1d(ranks,el,fill_value='extrapolate') for el in vec ]

        mean0= [kwargs['avgK'](r) + sinteg.quad(lambda y: kwargs['mu'](r,y)*n1(y)*phi(y),
            -np.inf,np.inf)[0]  for r in ranks]
        var0= [kwargs['varK'](r) + sinteg.quad(lambda y: kwargs['sigma'](r,y)*n2(y)*phi(y),
            -np.inf,np.inf)[0]  for r in ranks]
        v0=[1./(1 + sinteg.quad(lambda y: np.sqrt(kwargs['sigma'](r,y)*kwargs['sigma'](r,y))
            *phi(y)*kwargs['gamma'](r,y)*v(y),
            -np.inf,np.inf)[0] )  for r in ranks]
        print 'tmp2'
        return  mean0,var0,v0

    def calcmom(mom,mean0,var0):
        mean0,var0=[sinterp.interp1d(ranks,z,fill_value='extrapolate' ) for z in (mean0,var0) ]
        res= [sinteg.quad(lambda x: x**mom *gaussian(x,mean0(r),var0(r) ),0,np.inf)[0]
            for r in ranks]
        return res

    def model(x):
        N1,N2,V,Phi=x.reshape((4,resolution))
        mean0,var0,v=calcprm([N1,N2,V,Phi],**kwargs)

        eq1=Phi-calcmom(0,mean0,var0)
        eq2=Phi*N1-calcmom(1,mean0,var0)
        eq3=Phi*N2-calcmom(2,mean0,var0)
        eq4=V-v

        res= np.concatenate( (eq1,eq2,eq3,eq4))
        print res
        return res

    tries=0
    class Tmp():
        success=False
    res=Tmp()
    while not res.success:
        if tries >0:
            kw={}
            kw.update(kwargs)
            kw['sigma']=lambda x,y: kwargs['sigma'](x,y)/1.01
            x0=continuum_cavity_solve(**kw)
        else:
            #x0=continuum_mf(**kwargs)
            x0=np.random.random( 4*resolution )
            print x0
            #raise

        root=sopt.root
        res= root(model,x0)
        tries+=1

    result= res.x
    N1,N2,V,Phi=result.reshape((4,resolution))
    return N1,N2,V,Phi
