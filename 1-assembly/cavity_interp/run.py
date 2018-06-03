# -*- coding: utf-8 -*-

from interp import *
dft_prm={
        'n':{
    'type':'variable',
    'axes': ('com',),
    'death': 10.**-14,
    'shape':(80,),
    'mean':1.,
    'std':1,
    'sign':1,
    },

        'growth':{
    'type':'matrix',
    'variables':['n'],
    'dynamics':'nlv',
    'role':'growth',
    'mean':1., 'std':0.,
    'sign':1,
    },
        'selfint':{
    'type':'matrix',
    'variables':['n'],
    'dynamics':'nlv',
    'role':'diagonal',
    'mean':1, 'std':0.,
    'save':1,
    'sign':1,
    },


        'community':{
    'type':'matrix',
    'variables':[('n','com'),('n','com') ],
    'role':'interactions',
    'structure':'bunin',
    'mean':3.,
    'std':1.,
    #'symmetry':-1.,
    'dynamics':'nlv',
    'diagonal':0,
    },

}

dft_dyn={
    'nlv':{'type':'simplelv',
        'variables':[('n','com'),],
        },
    }



def imitate(inpath,outpath,use_measures=None,rerun=1,reproduce_corrKA=1,
        reproduce_connectivity=1, **kwargs):
    #Imitate any arborescence of model realizations using closest Bunin-like model
    import pandas as pd
    COLS=set([])
    for content in extract_data(inpath,**kwargs):
        idx,filedata,model,props=content

        measure={}
        measure.update(props)
        measure_gen(model,measure,typ='effective')
        #print [(v,measure[v]) for v in measure if 'bunin' in v]
        oldpath=filedata['path']
        newpath=Path(oldpath.replace(inpath,outpath) )
        newpath.mkdir()
        tmax,tsample=measure['tmax'],measure['tsample']
        prm=deepcopy(dft_prm)
        dyn=deepcopy(dft_dyn)
        prefix='n_couplings_cavity_'

        if reproduce_connectivity:
            conn=measure['n_connectivity']
        else:
            conn=1
        try:
            sigma_K=measure['n_capacity_std']
            meanK=measure['n_capacity_mean']
        except:
            sigma_K=prm['growth']['std']=measure['growth_std']
            meanK=measure['growth_mean']
        if reproduce_corrKA:
            #Correction added in matrix
            corrKA=measure['n_corr_KA']
        else:
            #Correction integrated into sigma_K
            corrKA=0
            if sigma_K>10**-10:
                effS=model.parameters['n']['shape'][0]*measure['n_connectivity']
                #print sigma_K,  effS*measure['n_corr_KA']*2
                if not sigma_K**2>effS*measure['n_corr_KA']*2:
                    print sigma_K**2,effS*measure['n_corr_KA']*2, measure['n_connectivity']
                sigma_K=np.sqrt(sigma_K**2 - effS*measure['n_corr_KA']*2)
                assert sigma_K>0

        prm['n']['shape']=model.parameters['n']['shape']
        prm['growth']['std']=sigma_K
        prm['growth']['mean']=meanK
        prm['community'].update(
            {
                'mean':measure[prefix+'mu'],
                'std':measure[prefix+'sigma'],
                'symmetry':measure[prefix+'gamma'],
                'requires':['growth'],
                'corr_KA':corrKA,#measure['n_corr_KA'],
                'connectivity':conn,
                }
            )
        #print measure
        dic={}
        dic.update(filedata)
        dic.update(model.export_params([p for p in model.parameters if not p in prm] ))
        dic.update({
            'tsample':tsample,
            'tmax':tmax,
            'path':newpath,
            })
        COLS=COLS.union(dic.keys())
        if not rerun:
            break

        newmodel=Model(parameters=prm,dynamics=dyn)


        if 0:
            A=newmodel.get_labeled_data('community')
            K=newmodel.get_labeled_data('growth')
            S=A.matrix.shape[0]

            adjusted=S*(S-1)

            #print '!!!', (K*A).matrix[1,5], K.matrix[1]*A.matrix[1,5], K.matrix[5]*A.matrix[1,5]
            print '!!!'
            print np.mean(A.matrix),np.std(A.matrix), np.std(K.matrix),np.std(newmodel.data['selfint'].matrix)
            oldA,oldK=(model.data['community'].matrix.T/model.data['selfint'].matrix).T,model.data['growth'].matrix/model.data['selfint'].matrix
            print np.mean(oldA), np.std(oldA),np.std(oldK)
            print measure['n_couplings_mu']/(model.parameters['n']['shape'][0]*measure['n_connectivity'])

            print '!!!',newmodel.parameters['community']['corr_KA'],np.sum((K*A).matrix)/adjusted -np.mean(
                    K.matrix)*np.sum(A.matrix)/adjusted



        success=newmodel.evol(tmax=tmax/50.,tsample=tsample,print_msg=1,**kwargs)
        if not success:
            print "ERROR!", oldpath
            continue
        newmodel.save(newpath,overwrite=1)


        dic.update(newmodel.export_params())

        table=[dic]
        pd.DataFrame(table).to_csv(newpath+'files.csv')

    rebuild_filelist(Path(outpath))
    make_measures(Path(outpath),use_measures=use_measures,death=10**-10,
        keep_cols=list(COLS))#,'trophic','assembly_corr'] )




def gamma_role(rerun=0,measure=0,nsys=10,show=0):
    #Comparison between competitive, trophic, disordered
    path=Path('interp_gamma_cavity')
    try:
        if not os.listdir(path):
            rerun=1
    except:
        rerun =1
    if rerun:
        #If not already run before
        kwargs={}
        prm=deepcopy(dft_prm)

        prm['fluxes']={
            'type':'matrix',
            'variables':[('n','com'),('n','com') ],
            'role':'fluxes',
            'save':False,
            'mean':50./1600., 'std':0,
            'sign':1,
            }
        prm['commtroph']={
            'type':'matrix',
            'variables':[('n','com'),('n','com') ],
            'role':'interactions',
            'structure':'trophic',
            'requires':['fluxes'],
            #'efficiency':0.724,
            'efficiency':0.6,
            'ordered':0,
            'dynamics':None,
            }

        prm['community'].update(
            {
                #'distribution':'exponential'
                }
            )

        #specs=[ ((16,),None,'nlv'),((5625,),'nlv',None) ][1:]
        specs=[ ((300,),None,'nlv'),((400,),None,'nlv'),((600,),None,'nlv'),((1600,),'nlv',None) ]#[:-1]
        axes=[
                ( ('n_shape','commtroph_dynamics','community_dynamics'),specs ),
                ( 'sys',range(10) )
        ]


        loop(path=path,axes=axes,check_key=0,print_msg=1,tmax=500000,tsample=10000.,
            #keep='all',
            parameters=prm,dynamics=dft_dyn,
            **kwargs )

    if rerun or measure:
        rebuild_filelist(path)
        make_measures(path,use_measures=['usual','effective','abundance',
            'matrix',
            measure_cavity,'testcavity'
            ])#,'trophic','assembly_corr'] )

    if not show:
        return

    axes=['n_shape']

    measures=open_measures(path)
    measures=measures[measures['avgN']!=0]
    show_hist(measures=measures,axes=axes, values=['n_abundance'] )
    from cavity import cavity_distri
    xs=np.linspace(0,0.5,100)
    sigma,mu,gamma=[measures['n_couplings_cavity_'+z].mean() for z in ('sigma', 'mu','gamma')  ]
    dist=cavity_distri(S=2000,sigma=sigma,sigma_k=0,mu=mu,gamma=gamma)
    plot(xs,[dist(x) for x in xs],hold=1,linewidth=2 )
    plt.xlabel(r'Abundance $n_i$')


    plot_main(measures=measures,axes=axes)
    plot_cavity(measures=measures, axes=axes,)

    if 0:
        plot_data(measures=measures, axes=axes,
            #split_by=['community_symmetry'],
            values=[('press',{'style':'plot','linewidth':2}), ('n_matrix_press',{'style':'scatter' } ) ] ,
            hold=1)
        plt.ylabel(r'Response to press$')
        plt.xlabel(r'Pool size $S$')
        #show_data(measures=measures,axes=axes, values= ['n_%alive','phi','avgN','stdN','n_biomass','n_biomass_std','n_interactions_mu' ,'n_interactions_sigma'])


    if 0:
        compare={'phi':'n_%alive','avgN':'n_biomass','stdN':'n_biomass_std','press':'n_matrix_press',}
        for trait in ('phi','avgN','stdN','press'):
            plot_data(measures=measures, axes=axes,
                hold=1,newfig=1,
                title=trait,
                values=[
                    (trait,{'style':'plot','log':'x'}),
                     (compare[trait],{'style':'scatter' } ) ]
                )

    plt.show()


def model_funcresp(S=100,mu=-.1,sigma=.3, gamma=-0.05,zeta=.3,Nc=10,**kwargs):
    dic= {'sigma':sigma,'mu':mu,'gamma':gamma,'sigma_k':zeta,'Nc':10**10 }
    rs=[]
    #for trials in range(3):
    rs=funcresp_cavity_solve(**dic)
        #rs.append(r)
    #rs=np.array(rs)
    #rs=np.mean(rs,axis=0)

    N1,N2,v,phi,f=rs
    #print f,mu
    dic2={}
    dic2.update(dic)
    #print dic
    import scipy.optimize as sopt
    #print N1,N2,v,phi
    print zeta,mu,sigma,gamma
    def inv(x):
        dici= {'sigma':np.clip(x[2]/1.5,0.001,None),'mu':np.clip(x[1],-1,None),'gamma':np.clip(x[3],-1,1),
            'sigma_k':np.clip(x[0],0.001,None),'Nc':10**10}
        res=funcresp_cavity_solve(**dici)
        if res[0]==0:
            return np.array(rs[:-1])
        #print '  TEST -',x
        #print '  NOW --',res
        return np.array(rs[:-1])-res[:-1]
    opt=sopt.root(inv, (zeta,mu,sigma,gamma) )
    x=opt.x
    print 'FINALLY',x
    #raise
    dic2= {'sigma':np.clip(x[2],0,None),'mu':x[1],'gamma':np.clip(x[3],-1,1),'sigma_k':np.clip(x[0],0,None) }
    #dic2['sigma']*=max(0.001,(1-f))
    #dic2['gamma']/=np.clip((1-f),np.abs(dic2['gamma']) ,1.)
    #dic2['mu']*= (1-f) + f * Nc/S /N1/phi
    #dic2['zeta']=dic2['sigma_k']
    #print dic2
    return dic2

def model_nei(S=100,c=1,mu=10,sigma=.5, gamma=-0.05,zeta=.3,Nc=10,**kwargs):
    k=S*c
    return {'sigma':sigma*np.sqrt(k)/S,'mu':mu*k/S,'gamma':gamma,'sigma_k':zeta }

def compare_prms_models(models,mode='trace'):
    import pandas as pd

    #spaces=[('mu','sigma'),('sigma','gamma'),('sigma','zeta') ]
    spaces=[('mu','sigma'),('gamma','zeta') ]
    nbpanels=len(spaces)
    panel=MutableInt(0)
    if mode=='hull':
        plt.figure()

    for space in spaces:
        xrge,yrge=np.array([0,0]),np.array([0,0])
        auto_subplot(plt,nbpanels,panel,rows=1)
        X=space[0]
        Y=space[1]
        plt.xlabel(r'$\{}$'.format(X ) )
        plt.ylabel(r'$\{}$'.format(Y) )

        for im,m in enumerate(models):
            color=get_cmap(im,len(models))
            model,mdict=m
            plotopts=mdict.get('options',{})
            var1,var2=mdict['vars']
            ms=[]
            for v2 in mdict[var2]:
                table=[]
                for v1 in mdict[var1]:
            #for var in itertools.product(*[mdict[l] for l in vlabels] ) :
                    dico={}
                    dico.update(mdict)
                    dico[var1]=v1
                    dico[var2]=v2
                    #for lab,val in zip(vlabels,var):
                        #dico[lab]=val
                    table.append(dict(model(**dico)))
                    #print table[-1]
                measure=pd.DataFrame(table)
                x=measure[X]
                y=measure[Y]
                ms.append((x,y))
                if mode=='trace':
                    plot(x,y,hold=1,color=color,**plotopts)
                elif mode=='scatter':
                    scatter(x,y,hold=1,c=color,alpha=0.1,edgecolors='none',**plotopts)


            if mode=='fill':
                x,y=zip(*ms)
                plt.fill_between(x[0], np.min(y,axis=0), np.max(y,axis=0),facecolor=color)
            elif mode=='hull':
                from scipy.spatial import ConvexHull
                from matplotlib.patches import Polygon
                points=np.concatenate(ms,axis=1).T
                if Y=='zeta' and model==model_trophic_satur:
                    points[:,1]=np.random.randint(0,2,points.shape[0])

                #points=np.array((np.array(x).ravel(),np.array(y).ravel())).T
                hull = ConvexHull(points)
                #print points[hull.simplices]
                #scatter(points[:,0],points[:,1],hold=1,c=color,alpha=1.,edgecolors='none',**plotopts)

                cent = np.mean(points, 0)
                pts = []
                for pt in points[hull.simplices]:
                    pts.append(pt[0].tolist())
                    pts.append(pt[1].tolist())

                pts.sort(key=lambda p: np.arctan2(p[1] - cent[1],
                                                p[0] - cent[0]))
                pts = pts[0::2]  # Deleting duplicates
                pts.insert(len(pts), pts[0])
                k = 1.1
                poly = Polygon(k*(np.array(pts)- cent) + cent ,alpha=0.1)
                               #facecolor=color, alpha=1)
                poly.set_capstyle('round')
                plt.gca().add_patch(poly)

            xrge=[f(points[:,0]) for f in (np.min,np.max)]
            yrge=[f(points[:,1]) for f in (np.min,np.max)]
            xspan=xrge[1]-xrge[0]
            yspan=yrge[1]-yrge[0]
            plt.xlim(xmin=xrge[0]-xspan/2.,xmax=xrge[1]+xspan/2.)
            plt.ylim(ymin=yrge[0]-yspan/2.,ymax=yrge[1]+yspan/2.)
            plt.tight_layout()
    plt.show()



def compare_prms(paths=[],options={}):
    axes={'sigma':'n_couplings_sigma','mu':'n_couplings_mu',
        'zeta':'n_capacity_std','gamma':'n_couplings_gamma'}
    measures=[]
    for path in paths:
        measure=open_measures(Path(path))
        if path in options:
            print options[path]
        measures.append(measure)

    spaces=[('mu','sigma')]#,('sigma','gamma'),('sigma','zeta') ]
    nbpanels=len(spaces)
    panel=MutableInt(0)
    for space in spaces:
        auto_subplot(plt,nbpanels,panel,rows=1)
        X=space[0]
        Y=space[1]
        plt.xlabel(r'$\{}$'.format(X ) )
        plt.ylabel(r'$\{}$'.format(Y) )
        for measure in measures:
            plotopts={}
            for z in (X,Y):
                plotopts.update(options.get(z,{}) )
            x=measure[axes[X]]
            y=measure[axes[Y]]
            scatter(x,y,hold=1,**plotopts)
    plt.show()



def intersection(rerun=0,nsys=10,measure=0,show=0):
    #Comparison between competitive, trophic
    path=Path('interp_intersection')
    if 0:
        compare_prms_models(
            [
                (model_trophic_satur,  {'vars':('m0','e'),'m0':np.linspace(20,20.1,10),'e':np.linspace(0,0.001,10)[1:-1],
                        'S':100,'s0':0.,'sigma_k':.3 } ),
            ],
            mode='scatter'
            )
    try:
        if not os.listdir(path):
            rerun=1
    except:
        rerun =1
    if rerun:
        #If not already run before
        kwargs={}
        prm=deepcopy(dft_prm)

                #(model_trophic_satur,  {'vars':('m0','e'),'m0':np.linspace(20,20.1,10),'e':np.linspace(0,0.001,10)[1:-1],
                        #'S':100,'s0':1.,'sigma_k':.3 } ),
        prm['fluxes']={
            'type':'matrix',
            'structure':'random',
            #'distribution':'exponential',
            'variables':[('n','com'),('n','com') ],
            'role':'fluxes',
            'save':False,
            'mean':20./100., 'std':1./100,
            'sign':1,
            }
        prm['commtroph']={
            'type':'matrix',
            'variables':[('n','com'),('n','com') ],
            'role':'interactions',
            'structure':'trophic',
            'requires':['fluxes'],
            'efficiency':0.,
            'ordered':0,
            'dynamics':None,
            'diagonal':0,
            }

        prm['community'].update(
            {
            'mean':10.,
            'std':1.,
            'symmetry':-1.,
            'diagonal':0,
            })

        prm['growth']['std']=0
        prm['n']['shape']=(100,)

        specs=[ (None,'nlv'),('nlv',None) ]
        axes=[
                ( ('commtroph_dynamics','community_dynamics'),specs ),
                ( 'sys',range(nsys) )
        ]


        loop(path=path,axes=axes,check_key=0,print_msg=1,tmax=500000,tsample=10000.,converge=1,
            parameters=prm,dynamics=dft_dyn,
            **kwargs )

    if rerun or measure:
        rebuild_filelist(path)
        make_measures(path,use_measures=['usual','effective','abundance',
            #'matrix',
            measure_cavity,#'testcavity'
            ])#,'trophic','assembly_corr'] )

    if not show:
        return

    axes=['commtroph_dynamics']

    dico={'commtroph_dynamics':'community',0:'competition',1:'predation' , 'n_abundance':'Species abundance'}
    measures=open_measures(path)
    measures['commtroph_dynamics'][measures['commtroph_dynamics'].isnull()]=0
    measures['commtroph_dynamics'][measures['commtroph_dynamics']=='nlv']=1

    show_hist(measures=measures,axes=axes, values=['n_abundance'],dictionary=dico ,log='y')
    from cavity import cavity_distri
    xs=np.linspace(0,np.max( [np.max(z) for z in measures['n_abundance']]),100)
    sigma,mu,gamma=[measures['n_couplings_cavity_'+z].mean() for z in ('sigma', 'mu','gamma')  ]
    #sigma,mu,gamma=1,10,-1
    dist=cavity_distri(S=2000,sigma=sigma,sigma_k=measures['n_capacity_std'].mean(),mu=mu,gamma=gamma)
    plot(xs,[dist(x) for x in xs],hold=1,linewidth=2 ,log='y')
    #plt.xlabel(r'Abundance $n_i$')


    #plot_main(measures=measures,axes=axes)
    #plot_cavity(measures=measures, axes=axes,)



    plt.show()

if __name__=='__main__':
    if 0:
        intersection(rerun=0,nsys=100,measure=0,show=1)
    if 0:
        compare_prms( paths=[ 'interp_mutualism_nested', 'interp_resource_gen_normal_RX',
        'trophic_order_stupid'
        ])#,options={ 'mu':{'log':'x'} }  )

    if 1:
        compare_prms_models(
            [
                (model_resource_generalists, {'vars':('R','X'), 'S':20, 'R':np.logspace(1,4,15), 'X':np.linspace(.3,.9,15), 'sigma_rho':.25 }  ),
                (model_trophic_satur,  {'vars':('m0','e'),'m0':np.linspace(0.1,25,10),'e':np.linspace(0,1,10)[1:-1] } ),
                #(model_resource_specialists, {'vars':('sigma_z','S'), 'sigma_z':np.linspace(1,100,100), 'S':np.linspace(1,100,100), 'R':300,'z':200,}  ),
                #(model_funcresp,{'vars':('mu','Nc'),'mu':np.linspace(0,1,40)[1:],'Nc':np.logspace(0.5,2.5,6) }  )
                #(model_funcresp,{'vars':('sigma','Nc'),'sigma':np.linspace(0,1,5)[1:],'mu':-.5,'Nc':np.logspace(0,2,3) }  )
                # (model_nei,{'vars':('S','c'),'c':np.logspace(-2,0,10),'S':np.logspace(0,2,10) }  )
                ],

            mode='hull',
            #mode='fill',
            )

    if 0:
        compare_prms_models(
            [

                (model_resource_generalists, {'vars':('R','X'), 'S':6, 'R':np.logspace(.8,.81,5), 'X':np.linspace(.9909,.991,5), 'sigma_rho':.25 }  ),
                (model_trophic_satur,  {'vars':('m0','e'),'sigma_k':0.05,'S':2000, 'm0':np.linspace(12.1,12.2,6),'e':np.linspace(0.01,.03,11),'s0':2, } ),
                ],

            mode='scatter',
            #mode='fill',
            )


    #ROLE OF GAMMA (DIFFERENT INTERACTION NATURES
    if 0:
        print conv_trophic(400.,3,0,0.33333)

        #Possibility 1: mu=60, sigma=2, S=100
        print invert_conv_trophic(mu=60,sigma=2,gamma=-1,x0=(600.,3,0.01,0.33333) )
        print conv_trophic(12368.,282,0.,0.575)

        #Possibility 2: mu=24, sigma=2, S=16
        print invert_conv_trophic(mu=24,sigma=2,gamma=-1,x0=(600.,3,0.01,0.33333) )
        print conv_trophic(5625.,174,0.,0.724)

        print invert_conv_trophic(mu=10.,sigma=1.,gamma=-1,x0=(800.,200,0.01,0.6) )
        print conv_trophic(1600.,50,0.,0.6)

        success=0
        sigma=0
        while not success:
            success=invert_conv_trophic(mu=3,sigma=sigma,gamma=-1)
            sigma+=0.1
        print success
    if 0:
        gamma_role(rerun=0,measure=0,show=1,nsys=10)