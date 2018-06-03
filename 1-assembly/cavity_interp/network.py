# -*- coding: utf-8 -*-

from run import *
from funcresp import measure_funcresp



def ordering(path=Path('ordering'),mode='dft',struct='all',
        rerun=0,measure=0,nsys=3,sysmin=0,show=0,funcresp=False,**kw):

    if struct == 'all':
        struct=['triangular','nested','bipartite','scalefree','assortativity','clustering','connectivity']
    if hasattr(struct,'__iter__'):
        for st in struct:
            ordering(path=path,mode=mode,struct=st,
            rerun=rerun,measure=0,nsys=nsys,sysmin=sysmin,show=0,funcresp=funcresp)
        ordering(path=path,mode=mode,struct=None,
            rerun=0,measure=measure,nsys=nsys,sysmin=sysmin,show=show,funcresp=funcresp)
        return


    if mode != 'dft':
        path=Path(path.strip()+'_{}'.format(mode))
    try:
        if not os.listdir(str(path)):
            rerun=1
    except OSError as e:
        print e
        rerun =1
    kwargs={}
    prm=deepcopy(dft_prm)
    prm['growth']['std']=0.25
    prm['n']['shape']=(120,)
    prm['community']['stdrel']=.3
    prm['community']['diagonal']=0
    prm['community']['structure']='random'
    prm['edges']={
        'type':'matrix',
        'variables':[('n','com'),('n','com') ],
        'mode':'network',
        'connectivity':.2,
        'role':'edges',
        'save':False,
        'structure':'random',
        }

    if funcresp:
        prm['threshold']={
            'type':'matrix',
            'variables':['n'],
            'dynamics':'nlv',
            'role':'threshold',
            'mean':5., 'std':0.,
            'sign':1,
        }
        dyn={
            'nlv':{'type':'simplefr',
                'variables':[('n','com'),],
                },
            }
    else:
        dyn=deepcopy(dft_dyn)


    if mode=='competition':
        prm['community']['mean']=-0.04
        prm['community']['symmetry']=.9
        prm['growth']['sign']=1

    elif mode=='predation':

        prm['final']=prm.pop('community')

        prm['community']={
            'type':'matrix',
            'variables':[('n','com'),('n','com') ],
            'role':'fluxes',
            'mean':0.05,'stdrel':.3,
            'save':False,
        }

        prm['final']['requires']=['community']
        prm['final']['efficiency']=0.2
        prm['final']['structure']='trophic'
        prm['final']['ordered']=0
        prm['final'].pop('symmetry',None)
        prm['growth']['sign']=0

    elif mode=='mutualism':
        if not funcresp:
            prm['community']['rowmax']=0.5
        else:
            prm['community']['mean']=0.05

        prm['community']['symmetry']=0.3
        prm['growth']['sign']=0
        #prm['growth']['std']=0.25



    prm['community']['requires']=['edges']
    prm['community']['edgefrom']='edges'

    order= np.linspace(0,1,11)
    axes=[]
    if struct =='nested':
        prm['edges']['structure']='nested'
        #prm['community']['connectivity']=.6
        axes=[
                ( 'edges_nestedness',order ),
            ]
    elif struct =='triangular':
        prm['edges']['structure']='triangular'
        #prm['community']['connectivity']=.6
        axes=[
                ( 'edges_triangular',order[6:] ),
            ]

    elif struct=='bipartite':
        prm['edges']['structure']='multipartite'
        axes=[
                ( 'edges_partition',order ),
            ]
    elif struct=='scalefree':
        prm['edges']['structure']='barabasi'
        axes=[
                 ('edges_order', order ),
            ]
        prm['edges']['connectivity']=prm['edges']['connectivity']/3.
        prm['community']['mean']=prm['community']['mean']*3.
    elif struct=='assortativity':
        axes=[
                 ('edges_assortativity', np.linspace(-1,1,2*len(order)) ),
            ]
    elif struct=='clustering':
        axes=[
                 ('edges_clustering',  np.linspace(-1,1,2*len(order)) ),
            ]
    elif struct=='connectivity':
        axes=[
                 ('edges_connectivity', np.linspace(0.05,1,len(order)) ),
            ]


    elif not struct is None:

        raise Exception('Unknown structure: {}'.format(struct) )

    axes+=[
                ( 'sys',range(sysmin,nsys+sysmin) ) ]

    if not struct is None:
        axes[0]= (axes[0][0],'community_order'), zip(axes[0][1],np.abs(axes[0][1]))
        kwargs['community_label']=struct


    if rerun:
        loop(path=path,axes=axes,check_key=0,print_msg=1,tmax=500000,ftol=0.0002,tsample=10000.,
            parameters=prm,dynamics=dyn,converge=1,
            **kwargs )


    if rerun or measure:
        # rebuild_filelist(path)
        make_measures(path,use_measures=['usual',
            'effective',
            ifelse(funcresp, measure_funcresp, 'usual'),
            ifelse(funcresp, measure_funcresp_cavity, measure_cavity),
            #'matrix',
            #'degree',
            #measure_ordered,
            #measure_ordered_cavity,
                                         # measure_groups,measure_groups_cavity,
            measure_cavity_compare,
            #'abundance',
            ],

            # read_filter=ifelse( struct is None,None,('community_label','==',struct)),
                death=10**-10)

    if not show:
        return

    axes=['community_order']
    split_by=['community_label']

    measures=open_measures(path)


    dico={'community_edgemean':r'$s$','n_shape':r'$S$','community_nestedness':'nestedness',
        'community_partition':'partition','edges_structure':'network',
        'community_structure':'structure','cavity_compare':'Distance to reference',
        'community_order':'Ordering','multipartite':'bipartite','random':'scale-free'}

    plot_modelprm(measures=measures,axes=axes,split_by=split_by,dictionary=dico)

    #plt.show()

    plot_main(measures=measures,axes=axes,split_by=split_by,dictionary=dico,aggregator='mean',log='')

    plot_data(measures=measures, axes=axes,
        split_by=split_by,
        hold=1,
        values=[('n_couplings_n_std',{'style':'plot'}),('n_couplings_n_mean',{'style':'plot','color':'r'}) ]
        )


    plot_data(measures=measures, axes=axes,
        split_by=split_by,
        hold=1,
        values=['n_connectivity']
        )

    if 0:
        plot_data(measures=measures, axes=axes,
            split_by=split_by,
            hold=1,
            values=['n_degree_std']
            )
        plot_data(measures=measures, axes=axes,split_by=split_by,hold=1,
            values=[('n_biomass_tot',{'style':'scatter'}),('totN',{'style':'plot'})],
            dictionary=dico,title='{}: Total biomass'.format(mode.capitalize()),)

    plot_data(measures=measures, axes=axes,
        split_by=split_by,
        hold=1,
        values=['cavity_compare'],
        dictionary=dico,
        title=mode.capitalize(),
        )
    plt.ylim(ymin=0,ymax=.5)

    if 0:
        #ABUNDANCE HIST
        if struct!='all':
            measures=measures[measures['community_structure']==struct]
        show_hist(measures=measures[measures['community_order']==1],axes=axes, values=['n_abundance'],dictionary=dico ,log='y')

        sigma,mu,gamma=[measures['n_couplings_'+z].mean() for z in ('sigma', 'mu','gamma')  ]
        zeta=measures['n_capacity_std'].mean()
        from cavity import cavity_distri
        print sigma, mu ,gamma
        dist=cavity_distri(S=measures['n_shape'][0][0],sigma=sigma,sigma_k=zeta,mu=mu,gamma=gamma)
        xs=np.linspace(0,np.max( [np.max(z) for z in measures['n_abundance']]),100)
        plot(xs,[dist(x) for x in xs],hold=1,linewidth=2 ,log='y')
        plt.show()
    plt.show()

def show_ordering(path='ordering'):
    path=Path(path)
    modes=['competition','predation','mutualism',]

    axes=['community_order']
    split_by=['community_label']


    dico={'community_edgemean':r'$s$','n_shape':r'$S$','nested':r'nested','znested':r'nested',
        'community_structure':r'structure','cavity_compare':r'Prediction error (%)',
        'community_order':r'Control parameter','zbipartite':r'bipartite','bipartite':r'bipartite','scalefree':r'scale-free',
        'connectivity':r'connectance','triangular':r'cascade','assortativity':r'assortative','clustering':r'clustered',
        'community_label':''}
    panel=MutableInt(0)
    for mode in modes:
        lpath=Path('{}_{}'.format(path,mode))
        auto_subplot(plt,len(modes),panel,rows=1)

        measures=open_measures(lpath)#,fname='ordered')
        measures['cavity_compare']=100./(0.5+ 1./ measures['cavity_compare'])
        # measures['cavity_compare']=100* measures['cavity_compare']
        # print measures['cavity_compare']
        print np.min(measures['cavity_compare'])

        measures['community_label']=[v if v != 'bipartite' else 'zbipartite' for v in measures['community_label'].values]
        measures['community_label']=[v if v != 'nested' else 'znested' for v in measures['community_label'].values]
        plot_data(measures=measures, axes=axes,
            split_by=split_by,
            hold=1,
            newfig=0,
            values=['cavity_compare'],
            # log='y',
            dictionary=dico,
            aggregator='median',
            title= r'({}) {} '.format('abc'[panel.val-1],mode.capitalize()),
            )
        from mpl_toolkits.axes_grid1.inset_locator import inset_axes
        ax=plt.gca()
        # plt.axhline(5)
        if panel.val>1:
            plt.ylabel('')
            ax.set_yticklabels(['' for item in ax.get_xticklabels()])
        if panel.val!=2:
            plt.xlabel('')
        if panel.val<len(modes):
            ax.legend_ = None
        plt.ylim(ymin=0,ymax=50)
        inset_axes = inset_axes(ax,
                                width="40%",  # width = 30% of parent_bbox
                                height=1.,  # height : 1 inch
                                loc=2)
        # plt.fill_between([-1,2],0.01,5,color='b',alpha=0.1)
        # plt.fill_between([-1,2],5,100,color='r',alpha=0.1)
        plot_data(measures=measures, axes=axes,
            split_by=split_by,
            hold=1,
            newfig=0,
            values=['cavity_compare'],
            log='y',show_legends=False,show_labels=False,
            dictionary=dico,
            aggregator='median',
            )
        plt.ylim(ymin=0.01,ymax=100)
        plt.xlim(xmin=0.,xmax=1)
        if panel.val<len(modes):
            plt.yticks([0.01,1,100])
        else:
            plt.yticks([])
        plt.xticks([])
        ax=plt.gca()
        ax.yaxis.tick_right()
        # plt.axhline(5,color='k')

    plt.tight_layout()
    plt.show()

if __name__=='__main__':
    RERUN='rerun' in sys.argv

    path=Path('ordering')

    ordering(path=Path(path),
        funcresp=False,
        mode=['competition','mutualism','predation'][1],
        struct='all',#['triangular','nested','bipartite','scalefree','connectance','assortativity','clustering'],
        rerun=RERUN,measure=1,nsys=100,sysmin=0,show=0)


    show_ordering(path='ordering')
