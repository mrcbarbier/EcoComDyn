# -*- coding: utf-8 -*-


from imports import *
from cavity import *
from numpy import sqrt

cavity_compare={'phi':'n_%alive','avgN':'n_biomass','stdN':'n_biomass_std','press':'n_matrix_press',
        'variability':'n_matrix_var','totN':'n_biomass_tot'}


def parametrization(param='biomass',prm=None):
    if param=='biomass':
        prm['interactions']['formula']='community'
    elif param=='rate':
        prm['interactions']['formula']='community*growth'
    elif param=='bunin':
        prm['interactions']['formula']='community*growth/capacity'
    elif param=='macarthur':
        #TODO!
        raise Exception("WIP")


def cavity_output(prm,**kwargs):
    """Properly formats output of cavity_solve"""
    S=prm['n_shape'][0]
    try:
        prefix='n_couplings_'
        mu=prm[prefix+'mu']
        sigma=prm[prefix+'sigma']
        #sigma =prm[prefix+'sigma']
        #print sigma,prm[prefix+'sigma']
        gamma=prm[prefix+'gamma']
    except:
        print 'COULD NOT FIND MEASURED COUPLINGS'
        mu=prm['community_mean']
        sigma=prm['community_std']
        gamma=prm['community_symmetry']
    try:
        sigma_k=prm['n_capacity_std']
    except:
        sigma_k=prm['growth_std']

    #print prm['n_corr_KA']
    #print '!!!!',prm.get('n_capacity_mean',1)
    return cavity_solve_props(S=S,mu=mu,sigma=sigma,sigma_k=sigma_k,gamma=gamma,
        corrKA=prm.get('n_corr_KA',0),avgK=prm.get('n_capacity_mean',1),
        conn=prm.get('n_connectivity',1),
        **kwargs
        )

def cavity_rowstd_output(prm,**kwargs):
    """Properly formats output of cavity_solve, separating col and row std"""
    S=prm['n_shape'][0]
    prefix='n_couplings_'
    mu=prm[prefix+'mu']
    sigma=prm['n_couplings_col_sigma']
    gamma=prm[prefix+'gamma']
    sigma_k=prm['n_capacity_std']

    return cavity_solve_props(S=S,mu=mu,sigma=sigma,sigma_k=sigma_k,gamma=gamma,
        corrKA=prm.get('n_corr_KA',0),avgK=prm.get('n_capacity_mean',1),
        rowstd=prm.get('n_couplings_row_std',0),
        conn=prm.get('n_connectivity',1),
        **kwargs
        )

def measure_cavity(model,measure,**kwargs):
    res=cavity_output(measure)
    #print res
    measure.update(res )

def measure_rowstd_cavity(model,measure,**kwargs):
    res=cavity_rowstd_output(measure)
    #print res
    measure.update(res )


def measure_funcresp_cavity(model,measure,**kwargs):
    def conv_funcresp(*arg,**kw):
        kw['Nc']=measure.get('threshold_mean',10**10)/2
        N1,N2,V,phi,f=funcresp_cavity_solve(*arg,**kw)
        #print N1
        measure['saturation_fraction']=f
        measure['saturation_flux']=kw['mu']*min(1,((1-f) + f * kw['Nc']/kw['S'] /N1/phi))
        avgN=N1*phi
        mu,sigma=kw['mu'],kw['sigma']
        v=V*sigma*phi
        q=N2*phi/avgN**2
        h= (kw.get('avgK',1)/avgN - mu)/sigma
    #avgN=1/(sigma * h + mu)*kwargs.get('avgK',1)
        return q,v,h,phi
    kwargs['solver']=conv_funcresp
    res=cavity_rowstd_output(measure,**kwargs)
    measure.update(res )

def conv_trophic(S=100,m=1,s=1,e=1):
    '''Given interactions of the form A_ij= (m+s N(0,1) )/S , A_ji=-e A_ij
    compute Bunin parameters mu, sigma, gamma'''
    S,m,s,e=[max(0.0001,float(x)) for x in (S,m,s,e)]
    e=np.clip(e,0,1)
    mu= -m*(e-1)/2
    sigma= (1+e)/sqrt(2 * S ) * sqrt(m**2/2 + s**2 )
    gamma = (-mu**2 - e*(m**2 + s**2) ) / sigma**2 /S
    return mu,sigma,gamma


def invert_conv_trophic(S=400,mu=1,sigma=0.05,gamma=-1,x0=None):
    from scipy.optimize import root
    if x0 is None:
        x0=np.ones(4)
        x0[0]=S

    def eqs(x):
        S,m,s,e=x

        mu2,sigma2,gamma2=conv_trophic(S=S,m=m,s=s,e=e)
        print x,mu2,sigma2,gamma2,mu,sigma,gamma
        cost=max(0,e-1)**2 + max(0,-e)**2 + max(0,-S)**2 + max(0,-m)**2+max(0,-s)**2
        return (mu2-mu,sigma2-sigma,gamma2-gamma,cost)#( S+mu2 +sigma2+gamma2-np.abs(S +mu2 +sigma2 +gamma2) )**2 )
    res= root(eqs,x0)
    if res.success:
        return res.x
    else:
        print res.message, res.fun
        return False

def cavity_trophic(prm):
    """Properly formats output of cavity_solve, first calculating mu and sigma from trophic"""
    from numpy import sqrt
    #print type(prm['n_shape'])

    S=prm.get('community_degree',prm['n_shape'][0])
    m=prm['community_mean']*S
    s=prm['community_std']*S
    e=prm['community_efficiency']
    mu,sigma,gamma=conv_trophic(S,m,s,e)

    try:
        print S,'|',m ,s,e,'|', mu, sigma, gamma, '|', S*prm['n_interactions_mean'],np.sqrt(S)* prm['n_interactions_std']
        mu=-S*prm['n_interactions_mean']
        sigma=np.sqrt(S)* prm['n_interactions_std']
    except:
        pass

    sigma_k=prm['capacity_std']
    return cavity_solve_props(S=S,mu=mu,sigma=sigma,sigma_k=sigma_k,gamma=gamma)

def cavity_cellular(prm):
    from numpy import sqrt
    #print type(prm['n_shape'])

    S=prm['n_shape'][0]
    sm=prm['interaction_self_mean']
    mu=-prm['interaction_other_mean']*S/sm
    sigma=prm['interaction_other_std']*sqrt(S)/sm
    gamma=prm['community_symmetry']
    sigma_k=prm['capacity_std']/prm['capacity_mean']

    ntries=10
    handled=0
    while ntries and not handled:
        res=cavity_solve_props(S=S,mu=mu,sigma=sigma,sigma_k=sigma_k,gamma=gamma,x0=np.random.random(4))
        if res['phi']>0:
            handled=1
        ntries-=1
    return res




if __name__ =='__main__':
    resolution=128
    for gamma in (1,0,-1):
        mat=cavity_show_loop(title='',sigma=np.linspace(.2,1.5,resolution),
            mu=np.linspace(-2,12,resolution),sigma_k=0.3,gamma=gamma,S=1000,rows=1 )

    plt.show()