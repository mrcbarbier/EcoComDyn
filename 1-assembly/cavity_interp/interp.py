# -*- coding: utf-8 -*-



from base import *

'''Plots for the cavity method interpretation paper'''



def S_trajectory(model=model_indep,var=('S',r'$S$'),**kwargs):
    '''Plot trajectory of bunin state variables as a function of variable (dft S) '''
    res=[]
    variable,name=var
    Smin,Smax=kwargs.get('Smin',10),kwargs.get('Smax',1000)
    if kwargs.get('log',True):
        Ss=np.logspace(np.log10(Smin),np.log10(Smax),100)
        log='x'
    else:
        Ss=np.linspace(Smin,Smax,100)
        log=''
    for S in Ss:
        print S
        kwargs[variable]=S
        r=model(**kwargs)
        kwargs['x0']=tuple(r[z] for z in ('q','v','h','phi'))
        if kwargs['x0'][-1]<=0 or np.isnan(kwargs['x0']).any():
            del kwargs['x0']
        r[variable]=S
        res.append(r)
    traj=pd.DataFrame(res)

    dico=kwargs.get('dico',cavity_dico)
    nbpanels=len([j for j in traj if j in dico])
    panel=MutableInt(0)
    for j in dico:
        #if j in cavity_dico:
            auto_subplot(plt,nbpanels,panel,rows=kwargs.get('rows',None))
            #plt.figure()
            plot(traj[variable],traj[j],title=dico[j],hold=1,log=log,xlabel=name)
            plt.xlim(xmin=min(traj[variable]),xmax=max(traj[variable]))
    if not kwargs.get('hold',0):
        plt.show()

def prm_plot(model=model_indep,axes=[],style='heatmap',log='',**kwargs):
    res=[]
    xlabel,xs=axes[0]
    ylabel,ys=axes[1]
    for x in xs:
        for y in ys:
            r={}
            kwargs[xlabel]=x
            kwargs[ylabel]=y
            r.update(model(**kwargs))
            #print r,kwargs
            r[xlabel]=x
            r[ylabel]=y
            res.append(r)
    traj=pd.DataFrame(res)

    nbpanels=len([j for j in traj if j in cavity_dico])
    panel=MutableInt(0)
    for j in sorted(traj):
        if j in cavity_dico:
            if style=='heatmap':
                ax=auto_subplot(plt,nbpanels,panel)
            else:
                ax=auto_subplot(plt,nbpanels,panel,projection='3d')
            #plt.figure()
            #plot([],[],title=cavity_dico[j],
                #hold=1,xlabel=xlabel,ylabel=ylabel)
            #print j, traj[j]
            plt.title(cavity_dico[j])
            plt.xlabel(xlabel)
            plt.ylabel(ylabel)
            if 'x' in log:
                plt.xscale('log')
            if 'y' in log:
                plt.yscale('log')
            X=traj[xlabel]
            Y=traj[ylabel]
            Z=traj[j]
            shape=len(xs),len(ys)
            if style=='heatmap':
                #plt.title('{}{}'.format(title,values[0]))
                plt.pcolor(X.reshape(shape),Y.reshape(shape),Z.reshape(shape),
                    vmin=max(min(Z),0),vmax=min(max(Z),10))
                plt.colorbar()
                plt.xlim(xmin=min(X),xmax=max(X))
                plt.ylim(ymin=min(Y),ymax=max(Y))
            else:
                ax.scatter(X,Y,Z,c=Z)
    plt.show()


def plot_main(measures=None,measures_imit=None,axes=None,split_by=None,log='',traits=None,
        theory=None,Svar=None,show_direct_calc=True,dictionary={},subplots=False,prefix='',
        compare=CAVITY_COMPARE,**kwargs ):


    if subplots:
        plt.figure()
        plt.suptitle(kwargs.pop('title',''))
    if traits is None:
        traits=('phi','avgN','totN','stdN','simpsoninv','relstdN')#,'press','variability')#,'prod')
    #measures[prefix+'simpsoninv']=1/measures[prefix+'simpson']
    measures['n_simpsoninv']=1/measures['n_simpson']
    panel=MutableInt(0)
    for trait in traits:
        if not subplots:
            plt.figure()
        else:
            ax,nbrows,panels_per_row=auto_subplot(plt,len(traits),panel,rows=kwargs.get('rows',None),return_all=1)

        plt.title(cavity_dico[trait])

        if theory:
            #comparison to theory
            xlabel=Svar[1]
            rstds=sorted(set(measures[split_by[0]]))
            for ir,rstd in enumerate(rstds):
                plt.xlim(xmin=min(Svar[-2]),xmax=max(Svar[-2]))
                plot(Svar[-2],[t[trait] for t in theory[rstd]] ,hold=1,color=get_cmap(ir,len(rstds)))
        else:
            xlabel=axes[0]


        if not measures_imit is None:
            plot_data(measures=measures_imit, axes=axes,
                newfig=0,
                split_by=split_by,
                hold=1,
                dictionary=dictionary,
                values=[
                    (compare[trait],{'style':'plot','marker':'^','linestyle':'None'}),
                     ],**kwargs
                )

        values=[ (compare[trait],{'style':'scatter','log':log ,
            'ylabel':'' ,'xlabel':dictionary.get(xlabel,xlabel)} ) ]

        if show_direct_calc or not theory:
            #Show trait as directly computed from measurements
            values+=[ (prefix+trait,{'style':'plot'}),]
            if theory:
                values[-1][-1].update({'marker':'x','linestyle':'None'})


        plot_data(measures=measures, axes=axes,
            newfig=0,
            split_by=split_by,
            hold=1,
            dictionary=dictionary,
            show_legends=(trait==traits[1] or not subplots),
            values=values,**kwargs
            )
        plt.ylabel('')
        if subplots:
            if panel.val<= panels_per_row*(nbrows-1):
                plt.xlabel('')
                plt.xticks([])



def plot_cavity(measures=None,axes=None,split_by=None,log=''):

    vals=[
             ('n_%alive',{'style':'scatter','log':log} ) ,
              ('cavity_alive',{'style':'plot',} ) ,
              #('cavity_alive_next',{'style':'plot','linestyle':'dotted','lw':4, } ) ,
              #('n_%alive_strict',{'style':'plot','linestyle':'--'} ) ,
              ('cavity_dN',{'style':'plot','linestyle':'--','marker':'+'} ) ,
              #('cavity_alive_truepos',{'style':'plot','linestyle':'None', 'marker':'x'} ) ,
              ('phi',{'style':'plot','linestyle':'None', 'marker':'x'} ) ,
             ]
    legs=['cavity_alive','cavity_dN','phi']
    if 'cavity_alive_rem' in measures:
        vals.append(
            ('cavity_alive_rem',{'style':'plot','linestyle':'None', 'marker':'^'})
              )
        legs.append('cavity_alive_rem')

    plot_data(measures=measures, axes=axes,
        hold=1,
        split_by=split_by,
        title='Cavity test (dots=%alive)',
        values=vals )
    plt.legend(legs)



def plot_removal(measures=None,axes=None,split_by=None,log=''):

    plot_data(measures=measures, axes=axes,
        title='v by removal (dots) vs theory (line)',
        split_by=split_by,
        values=[('removal_v',{'style':'scatter','log':log}),
            ('removal_v_std',{'style':'plot','linestyle':'None','marker':'x'}),
                ('v',{'style':'plot' } ),
             ],
        hold=1 )

    plot_data(measures=measures, axes=axes,
        title='Difference made by one species: removal (dots) vs cavity (line)',
        split_by=split_by,
        values=[('removal_diff',{'style':'scatter','log':log}),
                ('extinction',{'style':'plot' } ),
                #('cavity_diff',{'style':'plot' } ),
                #('cavity_diffnext',{'style':'plot','linestyle':'dotted' } ),
             ],
        hold=1 )


    plot_data(measures=measures, axes=axes,
        title='Previous abundance: |nprev_cavity - nprev_removal|',
        split_by=split_by,
        values=[
                ('cavity_v_vs_rem',{'style':'plot' ,'log':log} ),
             ],
        hold=1 )

def plot_modelprm(measures=None,axes=None,split_by=None,log='',measures_imit=None,
    **kwargs):
    conv={'n_couplings_mu':'predicted_mu',
        'n_couplings_sigma':'predicted_sigma',
        'n_couplings_gamma':'predicted_gamma',
        'n_capacity_mean':'predicted_meanK',
        'n_capacity_std':'predicted_sigK',
        'n_growth_std':'predicted_sigR',
        'n_corr_KA':'predicted_corrKA',
        }

            #('n_selfint_std',{'style':'scatter','log':'x'}),
            #('n_capacity_std',{'style':'scatter'}),
             #('predicted_sigD',{'style':'plot' } ) ,

    for i,j in conv.iteritems():
        values=[
                 (i,{'style':'scatter','log':log} ) ,
                (j,{'style':'plot'  } ) ,
                 ]
        if not j in measures:
            values=values[:1]
        plot_data(measures=measures, axes=axes,
            title=j.split('_')[1],
            hold=1,split_by=split_by,
            values=values,**kwargs )
        if not measures_imit is None:
            plot_data(measures=measures_imit, axes=axes,
            hold=1,newfig=0,split_by=split_by,
            values=[
                (i,{'style':'plot','marker':'^','linestyle':'None'}),
                 ] ,**kwargs)


if __name__=='__main__':
    #S_trajectory(model=model_indep)
    #S_trajectory(model=model_satur)
    if 0 :
        S_trajectory(model=model_resource_generalists,var=('R',r'$A$'),Smin=100,Smax=15000,hold=1)
        S_trajectory(model=model_resource_generalists,var=('R',r'$A$'),sigma_rho=0.2,Smin=100,Smax=15000,hold=1)
        S_trajectory(model=model_resource_generalists,var=('R',r'$A$'),sigma_rho=0.5,Smin=100,Smax=15000)
    if 0 :
        S_trajectory(model=model_resource_specialists,var=('sigma_z',r'$A$'),R=300,z=200,Smin=1.,Smax=100,hold=0)
    #S_trajectory(model=model_trophic_indep)
    #S_trajectory(model=model_trophic_satur)
    #S_trajectory(model=model_macarthur_discrete,Smax=100)
    #S_trajectory(model=model_macarthur_axis,Smax=100)
    if 0:
        prm_plot(model=model_resource_generalists,axes=[
            ('sigma_rho',np.linspace(0.0,0.5,64)),('R',np.linspace(100,400,64))   ],
                log='',X=.5)
    if 1:
        prm_plot(model=model_resource_specialists,axes=[
            ('sigma_rho',np.linspace(0.0,0.5,32)),('R',np.linspace(100,400,32))   ],
                log='',sigma_z=0.3)

    if 0:
        prm_plot(model=model_trophic_indep,axes=[ ('m',np.linspace(0.01,0.1,32)),
            ('e',np.linspace(0.1,0.9,32)  ) ])

    if 0:
        prm_plot(model=model_macarthur_axis,axes=[ ('sigma',np.linspace(0.4,0.7,32)),
            ('m',np.linspace(0,1,32)  ) ])