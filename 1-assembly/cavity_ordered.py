
from cavity_conv import *
from group_library import group_interactions,group_groupby,group_coordinates

ORDERED_COMPARE={'phi':'n_%alive','N1':'n_biomass*','avgN':'n_biomass',
        'stdN':'n_biomass_std','press':'n_matrix_press','variability':'n_matrix_var',
        'totN':'n_biomass_tot','simpson':'n_simpson','simpsoninv':'n_simpsoninv',
        'prod':'n_productivity_ratio' ,'relstdN':'n_biomass_relstd' }



def measure_ordered(model,measure,FIT_COEFFICIENTS=1,maxgroup=1000,**kwargs):

    rank,grps,members=group_groupby(model,measure,**kwargs)

    Aij=(model.find_data(role='interactions')#/model.data['selfint']
        ).matrix
    r=model.data['growth'].matrix

    mode=kwargs.get('mode',model.parameters['community'].get('mode','linear'))

    if 'mode' in model.parameters['growth'] and model.parameters['growth']['mode']!=mode:
        print 'WARNING: DIFFERENT MODE FOR GROWTH AND COMMUNITY, NOT ALLOWED YET'

    if 'selfint' in model.data:
        D=model.data['selfint'].matrix
    else:
        D=np.ones(r.shape )

    from scipy.optimize import curve_fit
    def f(x,a,b):
        return a+b*x

    S=np.array([np.sum(i) for i in members])
    mwhere=[np.where(m)[0] for m in members]
    conn,mu,sigma,sigrow,gamma=group_interactions(np.sum(S),Aij/D,mwhere)


    grank=np.array(grps)
    if np.max(rank)!=0:
        grank/=np.max(rank)
    Ks=[(r/D)[m] for m in mwhere]
    K=[np.mean(k) for k in Ks]
    zeta=[np.std(k) for k in Ks]


    alive=np.where(model.results['n'].matrix[-1]>10**-5)[0]

    rescrank=(rank-np.min(rank))/(np.max(rank)-np.min(rank))
    measure['n_mean_rank']=np.mean( rescrank[alive] )

    rankdiff = np.add.outer(grank, -grank)
    xs = grank
    if  FIT_COEFFICIENTS:
        #fit coefficients from matrices

        for k in ('K','zeta','mu','sigma','gamma'):
            #print locals()[k]
            if k=='mu':
                xs=rankdiff.ravel()
            ys=np.array(locals()[k]).ravel()
            if mode=='linear':
                popt,pcov=curve_fit(f, xs, ys)
            if mode=='exponential':
                if len(ys)<2:
                    popt=[1,1]
                elif 0 in ys:
                    popt=[0,1]
                else:
                    try:
                        popt,pcov=curve_fit(f, xs, np.log(np.abs(ys) ) )
                    except Exception as e:
                        code_debugger()
                    popt=np.exp(popt)
                    popt[0]=(2*(np.mean(ys)>0)-1) * popt[0]

            # print  popt
            measure['n_order_cons{}'.format(k) ]=popt[0]
            measure['n_order_slope{}'.format(k) ]=popt[1]

            # code_debugger()
            #print np.mean(locals()[k]),popt[0]
            #print popt[0]+popt[1]*rankdiff

    else:
        #Get coefficient from theoretical values
        print "DEBUG MODE in measure_ordered: get coefficients from theoretical values"
        S=r.shape[0]
        measure['n_order_consgamma']=measure['community_symmetry']
        measure['n_order_slopegamma']=0
        measure['n_order_consmu']=-measure['community_mean']*S
        measure['n_order_conssigma']=measure['community_std']*np.sqrt(S)
        measure['n_order_consK']=measure['capacity_mean']
        measure['n_order_slopeK']=measure['capacity_meanslope']
        measure['n_order_conszeta']=measure['capacity_std']
        measure['n_order_slopezeta']=measure['capacity_stdslope']
        if mode=='exponential':
            measure['n_order_slopemu']=measure['community_meanslope']
            measure['n_order_slopesigma']=measure['community_meanslope']
            for z in measure.keys():
                if 'n_order_slope' in z:
                    measure[z]=np.exp(measure[z])
                    if  'K' in z or 'zeta' in z:
                        measure[z]=measure[z]**2
            measure['n_order_consK']=measure['capacity_mean']/np.exp(measure['capacity_meanslope'])

        elif mode=='linear':

            measure['n_order_slopemu']=-measure['community_meanslope']*S
            measure['n_order_slopesigma']=measure['community_meanslope']*np.sqrt(S)
            measure['n_order_consK']=measure['capacity_mean']-measure['capacity_meanslope']
            measure['n_order_slopeK']=measure['capacity_meanslope']*2
            measure['n_order_conszeta']=measure['capacity_std']
            measure['n_order_slopezeta']=measure['capacity_stdslope']*2



def ordered_cavity_output(prm,**kwargs):

    kwargs['S']=prm['n_shape'][0]

    for k in ('K','zeta','mu','sigma','gamma'):
        kwargs[k ]=prm['n_order_cons{}'.format(k)]
        kwargs['c_{}'.format(k) ]=prm['n_order_slope{}'.format(k)]

    N1,N2,N1x,N2x,V,Vx,Phi,xmean= ordered_cavity_solve(**kwargs)

    res= basic_calc_props(sigma_k=kwargs['zeta'],N1=N1,v=V,N2=N2,phi=Phi,conn=1,**kwargs)
    res['mean_rank']=xmean
    return res


def measure_ordered_cavity(model,measure,**kwargs):
    kwargs['mode']=kwargs.get('mode',model.parameters['community'].get('mode','linear'))
    res= ordered_cavity_output(measure,**kwargs)
    measure.update(res)
    prefix=kwargs.get('prefix','')
    # measure['cavity_compare']=np.mean([ np.abs(measure[ORDERED_COMPARE[x]]/measure[prefix+x]-1)
    #      for x in ('avgN','simpson','phi') ] )
