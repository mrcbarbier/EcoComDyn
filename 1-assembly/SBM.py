import graph_tool.all as gt
from imports import *
from group_library import *



COMPARE = {
    'Phi': 'n_groups_%alive',
    # 'avgN':'n_groups_biomass',
    'totN': 'n_groups_biomass_tot',
    # 'totNmf':'n_groups_biomass_tot',
     'stdN':'n_groups_biomass_std',
    # 'simpson': 'n_groups_simpson',
}
COMPARE_MATS={
    # 'press_rel': 'n_groups_Vmean_self'
}


def SBM_state(Aij):
    g = gt.Graph()
    lst = zip(*np.where(Aij != 0))
    for z in np.where(np.sum(np.abs(Aij)+np.abs(Aij.T),axis=0)==0 )[0]:
        lst+=[(z,z)]
    g.add_edge_list(np.array(lst))
    tmp = g.new_edge_property('double')
    tmp.get_array()[:] = [Aij[i, j] for i, j in lst]
    g.edge_properties["weight"] = tmp
    state = gt.minimize_nested_blockmodel_dl(g, state_args=dict(recs=[g.ep.weight], rec_types=["real-normal"]))
    return [np.array(s) for s in state.get_bs()]


def SBM(Aij,lenmax=100,return_all=0):
    state=SBM_state(Aij)
    gps=[]
    for s in state:
        if not gps:
            gps.append(s)
        else:
            gps.append(np.array([s[x] for x in gps[-1]]  ) )
    gps=[ g for g in gps if len(set(g))<=lenmax ]
    if return_all:
        return gps
    return gps[0]


imitprm = {
    'niche':{
        'type': 'matrix',
        'variables': ['n'],
        'role': 'niche',
    },

    'n': {
        'type': 'variable',
        'axes': ('com',),
        'death': 10. ** -14,
        'mean': 1.,
        'std': 1,
        'sign': 1,
    },

    'growth': {
        'type': 'matrix',
        'variables': ['n'],
        'dynamics': 'nlv',
        'role': 'growth',
        'structure': 'mixture',
        'requires': ['niche'],
    },
    'selfint': {
        'type': 'matrix',
        'variables': ['n'],
        'dynamics': 'nlv',
        'role': 'diagonal',
        'mean': 1, 'std': 0.,
    },

    'community': {
        'type': 'matrix',
        'variables': [('n', 'com'), ('n', 'com')],
        'role': 'interactions',
        'structure': 'mixture',
        'dynamics': 'nlv',
        'diagonal': 0,
        'requires': ['niche'],
    },

}


def errorcalc(theo,obs,gselect=None,mode=0,remove_diag=0):
    if len(theo.shape)>1:
        if not gselect is None and theo.shape[0] ==len(gselect):
            theo = theo[np.ix_(gselect,gselect) ]
            obs = obs[np.ix_(gselect,gselect)]
        if remove_diag:
            # print theo
            # code_debugger()
            theo=theo[np.eye(theo.shape[0])!=1 ]
            obs = obs[np.eye(obs.shape[0]) != 1]
    elif not gselect is None:
        theo = theo[gselect]
        obs = obs[gselect]


    if np.sum(theo)<10**-10 and np.sum(obs)<10**-10:
        return 0
    if mode==0:
        #LIN ERROR
        return np.sum(np.abs(theo - obs)) /np.sum(np.abs(theo + obs))
    elif mode==1:
        #LOG ERROR
        return np.mean(np.abs(np.log10(theo) - np.log10(obs)) )  # /np.sum(meas[j])
    elif mode ==2:
        #WEIGHTED LOG ERROR
        weights=obs/np.max(obs)
        w = weights**.25 * np.exp(-0.1/weights**.5  ) #Seriously discounts values below 1% of max
        w[weights<10**10]=0
        if not w.any():
            w[:]=1
        weights = w/np.sum(w)
        good=np.logical_and(obs>0,theo>0 )
        return  np.sum(np.abs(np.log(theo[good] ) - np.log(obs[good] ))*weights[good] )


def group_imit(model,meas,*args,**kwargs):
    # =============== SIMULATION TESTS ====================

    compare=kwargs.get('compare',COMPARE)
    compare_mats=kwargs.get('compare_mats',COMPARE_MATS)
    gamma = meas['n_groups_gamma']
    sigma = meas['n_groups_sigma']
    S = model.parameters['n']['shape'][0]
    death = model.parameters['n']['death']

    x0 = np.concatenate([meas['n_groups_biomass*'], meas['n_groups_biomass*2'], meas['groups_V'], meas['n_groups_%alive']])

    Smean = np.mean(meas['n_groups_S'])
    clusterlabel=meas['n_groups_ranks']
    gamma[sigma < 10 ** -4] = 0

    prm = deepcopy(imitprm)
    prm['niche'].update({
        'matrix': clusterlabel, })

    prm['growth'].update({
        'rmean': meas['n_groups_capacity_mean'],
        'rstd': meas['n_groups_capacity_std'], })

    prm['n'].update({
        'shape': (S,), })

    prm['community'].update({
        'Amean': -meas['n_groups_mu'] / Smean,
        'Astd': sigma / np.sqrt(Smean),
        'Asym': gamma, })

    tests = range(1)
    locmeasures = []
    locmodels = []
    for test in tests:
        locmeas = {}
        locmeasures.append(locmeas)
        mtest = Model(parameters=prm, dynamics=model.dynamics)
        mtest.evol(tmax=100000, tsample=10000, converge=1, print_msg=1)
        locmodels.append(mtest)

        locmeas.update(mtest.export_params())
        # for field in measure_tscore_fields():
            # if field in prevmeas:
            #     locmeas[field] = meas[field]
        measure_groups(mtest, locmeas, groupby=clusterlabel)
        measure_groups_cavity(mtest, locmeas, x0=x0, prefix='')
        measure_group_stability(mtest, locmeas)

    locmeasures = pd.DataFrame(locmeasures)

    gselect = (meas['n_groups_biomass_tot'] > death * 10)
    for i in compare:
        j = compare[i]
        meas['simu_' + i] = tmp = np.mean(locmeasures[j], axis=0)
        meas['simu_error_{}'.format(i)] = errorcalc(tmp, meas[j], gselect)
        meas['simu_errorjoint_{}'.format(i)] = errorcalc(np.mean(locmeasures[j + '_joint'], axis=0),
                                                         meas[j + '_joint'])
        meas['simutheo_' + i] = np.mean(locmeasures[i], axis=0)

    for i in compare_mats:
        j = compare_mats[i]
        for m in (meas, locmeasures):
            if not j + '_joint' in m:
                m[j + '_joint'] = [np.array([np.sum(mj)]) for mj in m[j].values]
        meas['simu_' + i] = tmp = np.mean(locmeasures[j], axis=0)
        meas['simu_error_{}'.format(i)] = errorcalc(tmp, meas[j], gselect)
        meas['simu_errorjoint_{}'.format(i)] = errorcalc(
            np.mean(locmeasures[j + '_joint'], axis=0), meas[j + '_joint'])
