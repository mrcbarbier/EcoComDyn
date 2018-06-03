# -*- coding: utf-8 -*-

### MODELS FOR interp.py

from datatools import ADD_FOLDER_TO_PATH
ADD_FOLDER_TO_PATH('/../')
from imports import *
from cavity_conv import *
from group_library import *
from cavity_ordered import *
from cavity_continuum import *
from numpy import sqrt
from SBM import SBM, group_imit


CAVITY_COMPARE={'phi':'n_%alive','avgN':'n_biomass',
        'stdN':'n_biomass_std','press':'n_matrix_press','variability':'n_matrix_var',
        'totN':'n_biomass_tot','simpson':'n_simpson','simpsoninv':'n_simpsoninv',
        'prod':'n_productivity_ratio' ,'relstdN':'n_biomass_relstd'}

def measure_cavity_compare(model,measure,**kwargs):
    #orders=['nestedness','partition','edgeorder','order']

    #results=[]
    #for order in orders:
        #label='community_{}'.format(order)
        #if label in measure:
            #results.append(measure[label])
    #measure['community_order']=np.mean(results)

    #Comparison to cavity
    # print [(measure[CAVITY_COMPARE[x]],measure[x])
    #      for x in ('avgN','simpson','phi') ]

    measure['cavity_compare']=np.mean([ np.abs( measure[x]- measure[CAVITY_COMPARE[x]]) / ( measure[x]+ measure[CAVITY_COMPARE[x]])
         for x in ('avgN','simpson','phi') ] )


class PrmLooper(Looper):

    def loop_core(self,use_measures=['effective'],**kwargs):
        for i,j in self.attrs.iteritems():
            kwargs.setdefault(i,j)


        table=[]
        dic={}

        model=self.Model(**kwargs)
        model.evol(tmax=0,converge=0)

        dic.update(model.export_params())
        measure_gen(model,dic,typ=use_measures)
        if dic['n_couplings_sigma']>100:
            #for i,j in model.data.iteritems():
                #print i
                #print j.matrix
            raise
        table.append(dic)
        return table


def cavity_props_from_matrix(A,r,K):
    '''Compute analytically predicted properties from empirical matrices'''
    S=A.shape[0]
    mat= np.einsum('i,ij->ij',1/(r/K),A)
    np.fill_diagonal(mat,np.nan)
    mu=np.mean(nonan(mat)) * S
    sigma=np.std(nonan(mat)) * sqrt(S)
    gamma=S*(np.mean(nonan(mat*mat.T) ) -np.mean(nonan(mat))**2 )/sigma**2
    sigma_k=np.std(K)

    return cavity_solve_props(S=S,mu=mu,sigma=sigma,
            sigma_k=sigma_k,gamma=gamma)

def cavity_props_from_matrix_props(S=100,mean=1,std=1,Tmean=1,sigma_k=1):
    '''Compute analytically predicted properties from empirical matrix props'''
    mu=mean* S
    sigma=std * sqrt(S)
    gamma=S*(Tmean -mean**2 )/sigma**2
    print mu,sigma,gamma, sigma_k
    #if sigma_k==0:
        #sigma_k=0.0001

    return cavity_solve_props(S=S,mu=mu,sigma=sigma,
            sigma_k=sigma_k,gamma=gamma)



def model_sample(S=100,nmat=100,sampler=None,**kwargs):
    props=[]
    for m in range(nmat):
        A,r,K=sampler(**kwargs)
        mat= np.einsum('i,ij->ij',1/(r/K),A)
        np.fill_diagonal(mat,np.nan)
        mean=np.mean(nonan(mat))
        std=np.std(nonan(mat))
        Tmean=np.mean(nonan(mat*mat.T) )
        sigma_k=np.std(K/np.mean(K))
        props.append({'S':S,'mean':mean,'std':std,'Tmean':Tmean,'sigma_k':sigma_k})
        #props.append(cavity_props_from_matrix(A,r,K))
    props=pd.DataFrame(props)
    props= props.mean()
    attrs={}
    attrs.update(kwargs)
    attrs.update(props)
    return cavity_props_from_matrix_props(**props)



def model_macarthur_discrete(S=100,nmat=100,
        #R=5., #Ratio
        R=1000,
        mean_R=10,std_R=4,
        mean_xi=1,std_xi=0.3,
        mean_m=.5,std_m=.1,
        **kwargs):

    def sampler(**kw):
        B=np.random.normal(mean_R,std_R,R)
        xi=np.random.normal(mean_xi,std_xi, (S,R))/R
        try:
            m=np.abs(np.random.normal(mean_m,std_m,S))
        except:
            m=np.ones(S)*mean_m
        A=np.dot(xi,xi.T)
        r= np.dot(xi,B)-m
        K= r / np.diag(A)
        return A,r,K
    return model_sample(S=S,sampler=sampler,nmat=nmat,**kwargs)


def model_macarthur_axis(S=100,nmat=100,
        delta=1.,#niche separation
        sigma=.52,#niche width
        KR=5.,#resource capacity
        m=1., #mortality
        mean_xi=.2, #consumption
        std_xi=0.0,
        **kwargs):
    S=int(S)

    def sampler(**kw):
        i=np.array(range(S)).astype('float')+np.random.normal(0,0.05,int(S))
        try:
            xi=np.random.normal(mean_xi,std_xi, S)/sigma
        except:
            xi=np.ones(S)*mean_xi/sigma
        r=KR*xi-m
        K=2*sqrt(np.pi)*sigma*(1-m/xi/KR )/xi
        #print r.shape, K.shape
        A=np.einsum('i,ij->ij',(r/K),
                np.exp( - np.add.outer(i,-i)**2 *delta**2 / 4/ sigma**2 ) )
        return A,r,K
    return model_sample(S=S,sampler=sampler,nmat=nmat,**kwargs)


def model_indep(S=100,mu=0.01,sigma=0.005,sigma_k=.3,gamma=1,**kwargs):
    '''Symmetrical interactions ~ S^0 '''
    mu=mu*S
    sigma=sigma*sqrt(S)

    return cavity_solve_props(S=S,mu=mu,sigma=sigma,
            sigma_k=sigma_k,gamma=gamma,**kwargs)

def model_satur(S=100,mu=2,sigma=1.5,sigma_k=.3,gamma=1,**kwargs):
    '''Symmetrical interactions ~ S^-1'''
    mu=mu
    sigma=sigma/sqrt(S)
    return cavity_solve_props(S=S,mu=mu,sigma=sigma,
            sigma_k=sigma_k,gamma=gamma,**kwargs)

def model_trophic_indep(S=100,m=0.01,s=.005,e=.3,sigma_k=.3,**kwargs):
    m=m*S
    s=s*S
    mu,sigma,gamma=conv_trophic(S,m,s,e)
    return cavity_solve_props(S=S,mu=mu,sigma=sigma,
            sigma_k=sigma_k,gamma=gamma,**kwargs)

def model_trophic_satur(S=100,m0=2.,s0=1.,e=.3,sigma_k=.3,**kwargs):
    m=m0
    s=s0
    mu,sigma,gamma=conv_trophic(S,m,s,e)

    res= cavity_solve_props(S=S,mu=mu,sigma=sigma,
            sigma_k=sigma_k,gamma=gamma,**kwargs)
    res['mu']=mu
    res['sigma']=sigma
    res['gamma']=gamma
    res['zeta']=sigma_k

    return res

def convert_resources_spec(S,R,z,sigz,sigrho):
    S,R=float(S),float(R)
    mu=S*z/R
    Z=sigz/z
    #print "YAY", (1-z/R), sigrho**2*(1.-z/R)
    z=min(z,R*0.999)
    sigma = np.sqrt(S)* np.sqrt( sigz**2/R**2 + (R-z)**2/R**3 )
    gamma= 1./(1+ R* sigz**2/(R-z)**2 )
    sigR=np.sqrt( Z**2 + sigrho**2*(1.-z/R)) #*z/z**2/(1.-m)**2 )
    sigD=Z
    varK= sigrho**2#/z*(1-z/R) #-( 1-2*z+ (sigz**2+z**2)/R + (1-z/R)**2 )  /(z**2 * (1-m)**2 * R)
    sigK= np.sqrt(varK)#np.sqrt(sigrho**2 *R/z**2 )#
    #if sigK>2:
        #print "YAY",z,sigK,sigrho
        #raise
    meanK=1
    corrKA=0
    return mu,sigma,gamma,sigR,sigD,sigK,meanK,corrKA

def convert_resources_gen(S,R,X,sigrho,xi3,xi4,m):
    #sigR2=np.sqrt(sigrho**2+ R*rho**2 *sigxi**2 +sigm**2)
    #print '===<',sigR, sigR2
    #xi=np.sqrt(X/float(R))
    Y=1-X
    sigrho*=np.sqrt(Y/X)
    sigR=np.sqrt(sigrho**2  + (1-X)/(X*R))  # *(1-m)**2+ (1-X)/X/R )
    #mu=S*R*xi**2 = S*X
    mu=S*(X - 2*X*Y/R-Y*np.sqrt(X*Y)*xi3/R   )
    gamma=(1 + X**2 * (-5 + 4 *X) +
         2* (-X *np.sqrt(-(-1 + X) *X) + np.sqrt(-(-1 + X)* X**5 )) *xi3)/(1 +
         X**2 * (-6 - 5 *(-2 + X)*X) + (-2 *X *np.sqrt((1 - X) * X) +
            6 *np.sqrt(-(-1 + X)* X**5) - 4 * np.sqrt((1 - X) * X**7)  ) * xi3 + (-1 + X)**2 * X**2 * xi4)
    sigma=np.sqrt(S/R*(1 + X**2 *(-6 - 5* (-2 + X) *X) +
        (-2 *X* np.sqrt(-(-1 + X) *X) +
        6 *np.sqrt(-(-1 + X)* X**5) -
         4 *np.sqrt(-(-1 + X)* X**7 ))* xi3
         + (-1 + X)**2 *X**2 *xi4) )
    #print sig2,xi3,xi4,gamma

    #sigD2=np.sqrt(R* (np.mean(xip**4)- (xi**2 + sigxi**2)**2))
    sigD=np.sqrt((-1 + X)* (1 - 5 *X -  4 *np.sqrt((1 - X)* X)*xi3 + (-1 + X)* xi4))/np.sqrt(R)

    #sigK=np.sqrt ( sigrho**2 * (1-m)**2 -((-1 + X)* (1 +
       #5 *(-1 + X)* X - (-1 + m)**2 * X *(-1 + 9 *X)* sigrho**2 -
       #2 *(np.sqrt( (1 - X)* X) - 2 *np.sqrt((1 - X)* X**3) +
      #4 *(-1 + m)**2* np.sqrt(-(-1 + X)* X**3) *sigrho**2)* xi3 + (-1 + X)* X *(-1 + 2 *(-1 + m)**2 *sigrho**2) *xi4))/(X* R) )


    sigK=np.sqrt(sigrho**2 + (sigrho**2*X*Y - 9*sigrho**2*X**2*Y + (1 + 5*(-1 + X)*X)*Y +
        (1 - 2*sigrho**2)*X*xi4*Y**2 - 2*(1 + (-2 + 4*sigrho**2)*X)*xi3*Y*(X*Y)**0.5)/(R*X) )

    meanK=1 + (-2  + 2* X - (1 - X)**(1.5)/X**.5 *xi3)/(R)

    #corrKA=((-1 + X)* (-1 - 5 *(-1 + X)* X +
           #2* (np.sqrt(-(-1 + X) *X) - 2* np.sqrt((1 - X)*X**3 )) *xi3 + (-1 + X)* X*xi4))/R /2.
    #TODO: Understand factor 1/2 that I had to add to corrKA for perfect matching

    corrKA=(Y*(1 - 5*X*Y + X*xi4*Y + 2*(-1 + 2*X)*xi3*(X*Y)**0.5))/R
    return mu,sigma,gamma,sigR,sigD,sigK,meanK,corrKA

def moments_distri(distri):
    if distri=='exponential':
        from math import factorial

        momxi,momxi2,momxi3,momxi4=[factorial(n) #*(xi/sigxi)**n
                                    for n in range(1,5)]
        #momxi2=1
        #print '\n          ',xi, sigxi, xi3,xi4
    elif distri=='normal':
        momxi,momxi2,momxi3,momxi4=0.,1.,0.,3.
    xi3=momxi3 - 3*momxi*momxi2 +3* momxi**2 * momxi -momxi**3
    xi4=momxi4 - 4*momxi*momxi3 + 6*momxi**2*momxi2 +4* momxi**3 * momxi -momxi**4
    return xi3,xi4

def model_resource_generalists(S=100,R=200.,X=.5,sigma_rho=0.0,**kwargs):
    ''''''
    R=float(R)
    mu=S*X
    #sigma=(1-X**2)*S/float(R)
    #print "YAY",S,R,X,sigma_rho

    m=kwargs.get('mortality',0)

    xi=np.sqrt(X/R)
    sigxi=np.sqrt((1.-X)/R)
    #mu=S*R*xi**2

    #moments of (xi - <xi>)/sigxi
    distri=kwargs.get('distribution','exponential')
    xi3,xi4=moments_distri(distri)
    rho=1/(1.-m)/(R *xi)


    mu,sigma,gamma,sigR,sigD,sigK,meanK,corrKA=convert_resources_gen(S,R,X,sigma_rho*rho,xi3,xi4,m)
    #print sigK
    res= cavity_solve_props(S=S,mu=mu,sigma=sigma,
            sigma_k=sigK,gamma=gamma,avgK=meanK,corrKA=corrKA,**kwargs)


    res['mu']=mu
    res['sigma']=sigma
    res['gamma']=gamma
    res['zeta']=sigK
    return res

def model_resource_specialists(S=100,R=200.,z=50.,sigma_z=10.,sigma_rho=0.0,**kwargs):
    ''''''
    m=kwargs.get('mortality',0)

    z=float(z)
    R=float(R)
    rho=1/np.sqrt(z)/(1-m)

    #print S,R,z,sigma_z,sigma_rho*rho

    mu,sigma,gamma,sigR,sigD,sigK,meanK,corrKA=convert_resources_spec(S,R,z,sigma_z,sigma_rho*rho)

    res= cavity_solve_props(S=S,mu=mu,sigma=sigma,
            sigma_k=sigK,gamma=gamma,avgK=meanK,corrKA=corrKA,**kwargs)
    #return cavity_solve_props(S=S,mu=mu,sigma=sigma,
            #sigma_k=sigma_R,gamma=1,**kwargs)

    res['mu']=mu
    res['sigma']=sigma
    res['gamma']=gamma
    res['zeta']=sigK
    return res

