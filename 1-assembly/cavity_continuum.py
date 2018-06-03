# -*- coding: utf-8 -*-
from cavity_conv import *
from group_library import group_groupby,group_cavity_solve

import scipy.interpolate as inter
import scipy.integrate as sinteg
from scipy.optimize import curve_fit

ORDERED_COMPARE={'phi':'n_%alive','N1':'n_biomass*','avgN':'n_biomass',
        'stdN':'n_biomass_std','press':'n_matrix_press','variability':'n_matrix_var',
        'totN':'n_biomass_tot','simpson':'n_simpson','simpsoninv':'n_simpsoninv',
        'prod':'n_productivity_ratio' ,'relstdN':'n_biomass_relstd' }

def continuum_interactions(measure,rank,A,K):



    func=lambda x,a,b,c:continuum_quad(x,a,b,c)
    hist,dens=np.histogram(rank,bins=50,normed=1)
    bins=np.array((dens[1:]+dens[:-1]))/2.
    Stot=len(rank)

    Sprm=np.array([curve_fit(func,bins[hist>0],np.log(hist[hist>0]*Stot) )[0]])
    measure['Sprm']=Sprm
    measure['Sopt']=(('dim',1), )

    S=continuum_func(Sprm,**dict(measure['Sopt']))


    Kpos,Kneg=(K>0),(K<=0)

    avgKprm=np.array([curve_fit(func,np.array(rank),
         np.log(np.abs( ks )) )[0] for ks in (np.clip(K,0.00001,None),np.clip(K,None,-0.00001) ,)])

    measure['avgKprm']=avgKprm
    measure['avgKopt']=(('sgn',(1,-1) ),('dim',1) )

    avgK=continuum_func(avgKprm,**dict(measure['avgKopt']))
    varKprm=np.array([curve_fit(func,np.array(rank[Ax]),
         np.log( (K[Ax] -avgK(rank[Ax]) )**2 ) )[0] for Ax in (Kpos,Kneg) if sum(Ax) ])


    measure['varKprm']=varKprm
    measure['varKopt']=(('dim',1), )


    if 0:
        Rank=sorted(rank)
        print sinteg.quad(S,Rank[0],Rank[-1]),Stot
        from datatools import plot,scatter,plt
        plot(Rank,S(Rank),hold=1)
        scatter(bins,hist*Stot,xlabel='rank',ylabel="S",title='S'  )
        plot(Rank,avgK(Rank),hold=1)
        scatter(rank,K)

    func=lambda x,a,bx,cx,by,cy,bxy,cxy:continuum_quad(x,a,bx,cx,by,cy,bxy,cxy)


    measure['ranks_min']=np.min(rank)
    measure['ranks_max']=np.max(rank)

    ranks=np.tile(rank,(len(rank),1) )
    shape=ranks.shape
    xs=ranks.ravel()
    ys=ranks.T.ravel()
    A[np.abs(A)<10**-5]=0
    Apos=(A.ravel()>0)
    Aneg=(A.ravel()<0)


    muprm=np.array([curve_fit(func,np.array((xs[Ax],ys[Ax])),
         np.log(np.abs(A).ravel()[Ax]) )[0] for Ax in (Apos,Aneg) if np.sum(Ax)>7 ])

    #print mupopt,mumopt
    mu=continuum_func(muprm,sgn=ifelse(np.sum(Aneg)>7,-1,1) )

    measure['muprm']=muprm
    measure['muopt']=(('sgn',-1), )


    if 1:
        from datatools import scatter3d,plt
        Ranks=np.tile(sorted(rank),(len(rank),1) )
        XS=Ranks
        YS=Ranks.T

    if 0:
        scatter3d(xs[::],ys[::],np.log(np.abs(A).ravel()[::]),hold=1,alpha=0.3,c='r')
        plt.gca().plot_wireframe(XS,YS,
            np.log(np.abs(mu(XS.ravel(),YS.ravel())) ).reshape(shape))
        plt.show()

    sigs= np.clip( (A.ravel()-mu(xs,ys))**2, 10**-5,None)
    sigprm=np.array([curve_fit(func,np.array((xs[Ax],ys[Ax])),
         np.log(  sigs[Ax] )   )[0]
          for Ax in ((A.ravel()!=0),)  ]) #(Apos,Aneg)])


    measure['sigmaprm']=sigprm
    sigma=continuum_func(sigprm,sgn=1)

    if 0:
        scatter3d(xs[:],ys[:],np.log( sigs[:] ),hold=1,alpha=0.3,c='r')
        plt.gca().plot_wireframe(XS,YS,
            np.log(sigma(XS.ravel(),YS.ravel()) ).reshape(shape))
        plt.show()

    gam=(A.ravel()-mu(xs,ys))*(A.T.ravel()-mu(ys,xs)) /np.sqrt(sigma(xs,ys),sigma(ys,xs))
    gammprm=np.array([ curve_fit( func,np.array((xs[Ax],ys[Ax])),np.log(np.abs(gam[Ax]) ),
        method='trf',loss='soft_l1' ) [0]
        for Ax in (Apos,Aneg)])


    measure['gammaprm']=gammprm
    gamma=continuum_func(gammprm,sgn=[-1,-1],mode='exp')

    if 0:
        gam= np.clip(gam, -1,1)
        #scatter3d(xs[::],ys[::], gam,hold=1,alpha=0.3,c=gam)
        scatter3d(0,0,0,hold=1)
        plt.gca().plot_wireframe(XS,YS,
            gamma(XS.ravel(),YS.ravel() ).reshape(shape))
        plt.show()


def measure_continuum(model,measure,FIT_COEFFICIENTS=1,**kwargs):

    rank=group_groupby(model,measure,**kwargs)[0]
    measure['continuum_ranks']=tuple(rank)


    Aij=(model.find_data(role='interactions')[0]#/model.data['selfint']
        ).matrix
    r=model.find_data(role='growth')[0].matrix

    if 'selfint' in model.data:
        D=model.data['selfint'].matrix
    else:
        D=np.ones(r.shape )

    rescrank=(rank-np.min(rank))/(np.max(rank)-np.min(rank))

    alive=np.where(model.results['n'].matrix[-1]>10**-5)[0]

    measure['n_mean_rank']=np.mean( rescrank[alive] )

    continuum_interactions(measure,rank,-(Aij.T/D).T,r/D)

def continuum_calc_props(prefix='',S=None,N1=None,N2=None,V=None,Phi=None,conn=1,**kwargs):
    #Without Bunin scaling

    props=('v','phi','chi2','press','variability','varcavity','utilde',
        'avgN','avgN2','relstdN','stdN','corr_AN','totN','simpson','simpsoninv','prod','extinction')


    rmin,rmax= kwargs['ranks_min'],kwargs['ranks_max']

    Stot=sinteg.quad(S, rmin,rmax)[0]
    #
    v=sinteg.quad(lambda x:V(x)*Phi(x), rmin,rmax)[0]
    phi=sinteg.quad(lambda x:Phi(x), rmin,rmax)[0]
    chi2=0
    press=0
    variability=0
    C=0
    utilde=0
    avgN=sinteg.quad(lambda x:N1(x)*Phi(x), rmin,rmax)[0]
    avgN2=sinteg.quad(lambda x:N2(x)*Phi(x), rmin,rmax)[0]
    stdN=sinteg.quad(lambda x:(N2(x)-N1(x)**2)*Phi(x), rmin,rmax)[0]
    relstdN=sinteg.quad(lambda x:(N2(x)/N1(x)**2-1)*Phi(x), rmin,rmax)[0]
    meanAN=0
    simpson=avgN2/(avgN)**2/Stot/phi
    prod=0
    extinct=0


    props=[prefix+p for p in props]
    return pd.Series(dict(zip( props ,
        (v,phi,chi2,press,variability,C,utilde,avgN,avgN2,relstdN,stdN,-meanAN,
             avgN*Stot,simpson,1./simpson,prod,Stot*extinct) )))

    #if avgN<=0 or phi<=0:
        #base=np.zeros(len(props))
        ##set press and variability to 1 by default
        #ref=7
        #for p in props:
            #if p in default:
                #base[props.index(p)]=default[p]

        #props=[prefix+p for p in props]
        #return pd.Series(dict(zip( props ,
            #base)))

    ##sigma_l=(sigma_k/sigma/avgN)
    ##meanAN= mu - (1+ gamma)/(h + mu/sigma) * (1+ h / ( phi * (q + sigma_l**2) ) )
    ##meanAN/=S * phi
    #meanAN=0
    #u=(1.-mu/Seff)
    #ut=(u-phi* gamma*v*sigma**2)
    #press=1./( ut**2 - phi*sigma**2  )
    #chi2=phi*sigma**2*press
    #utilde= ut

    #variability=varJeff
    #extinct=0
    #simpson=N2/S/phi/N1**2
    ##extinct,deltaN,deltaN2=cavity_extinction(N0=avgN,S=S,mu=mu,sigma=sigma,gamma=gamma,v=v,phi=phi,avgN=avgN,varN=stdN**2)



def continuum_cavity_output(prm,**kwargs):
    dic={}
    dic.update(kwargs)
    dic.update(prm)
    Stot=prm['n_shape'][0]

    if kwargs.get('groups',1):
        resolution=dic.get('resolution',10)
        rank=dic['ranks']=np.linspace(dic.get('ranks_min',0),
            dic.get('ranks_max',1),resolution)

        fnc=continuum_func
        functions={}
        conv={'avgK':'mean_k','varK':'sigma_k'}
        fnconv={'varK':np.sqrt,'sigma':np.sqrt}
        for n in ('S','avgK','varK' , 'mu','sigma','gamma' ):
            functions[n]=fnc(dic[n+'prm'],**dict(dic.get(n+'opt', () )) )
            if n in ('mu','sigma','gamma'):
                ranks=np.tile(rank,(len(rank),1) )
                xs=( ranks.T, ranks )
            else:
                xs=[rank]
            dic[conv.get(n,n)] = fnconv.get(n,lambda x:x)(functions[n](*xs ))

        dic['S']/= np.sum(dic['S'])/Stot
        Smean=np.mean(dic['S'])
        dic['mu']*=Smean
        dic['sigma']*=np.sqrt(Smean)


        #return group_cavity_solve_props(**dic)
    #else:

        #print dic
        N1,N2,V,Phi,success=group_cavity_solve(**dic)
        #from datatools import plot,scatter
        #plot(rank,N1*Phi)

        plot(N1,N2,Phi,xs=rank,legend=['N1','N2','Phi'],hold=1)
        N1,N2,V,Phi= [inter.interp1d(rank,x) for x in (N1,N2,V,Phi) ]#continuum_cavity_solve(**dic)

        for n in functions:
            dic[n]=functions[n]
        res=continuum_calc_props(N1=N1,N2=N2,V=V,Phi=Phi,**dic)
    else:
        raise

    return res


def measure_continuum_cavity(model,measure,**kwargs):
    #from datatools import plot,scatter
    #scatter(model.data['niches'].matrix, model.results['n'].matrix[-1],hold=1)
    res= continuum_cavity_output(measure,**kwargs)
    measure.update(res)

    scatter(measure['continuum_ranks'],model.results['n'][-1],hold=1 )
    scatter(measure['continuum_ranks'],(model.results['n'][-1]>10**-10),hold=1,c='g' )
    plt.show()
    #prefix=kwargs.get('prefix','')
    #measure['cavity_compare']=np.mean([ np.abs(measure[ORDERED_COMPARE[x]]/measure[prefix+x]-1)
         #for x in ('avgN','simpson','phi') ] )
#