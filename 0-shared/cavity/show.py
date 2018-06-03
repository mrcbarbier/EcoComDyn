# -*- coding: utf-8 -*-
from cavity import *
from extinction import *
from collections import OrderedDict

#===========
cavity_dico=OrderedDict([
    ('phi',r'Fraction of survivors'),
    ('avgN',r'$\langle N \rangle$'),
    ('stdN',r'$\sigma_N $'),
    ('relstdN',r'$\sigma_n $'),
    ('totN',r'Total biomass'),
    ('totNmf',r'Total biomass (mean-field)'),
    ('simpson',r'Simpson index'),
    ('simpsoninv',r'Simpson diversity'),
    ('v','Species independence'),
    #('corr_AN',r'$\langle N_i A_{ij} \rangle$'),
    ('chi2','Phase parameter'),
    ('press','Response to press'),
    ('variability','Variability'),
    ('prod','Productivity ratio')
    ])



def basic_calc_props(prefix='',S=None,mu=None,sigma=None,
            sigma_k=None,gamma=None,N1=None,N2=None,v=None,phi=None,avgK=1,conn=1,**kwargs):
    #Without Bunin scaling


    props=('v','phi','chi2','press','variability','varcavity','utilde',
        'avgN','avgN2','relstdN','stdN','corr_AN','totN','simpson','simpsoninv','prod','extinction')

    zeta=np.clip(sigma_k,0.0001,None)
    default={
        'variability':1,
        'prod':1,
        'press':1,
        }
    #q=N2/N1**2/Phi
    #h=(1./(N1*Phi) - mu)/sigma
    Seff=S*conn
    sol=N1,N2,v,phi
    import pandas as pd
    avgN=N1*phi
    avgN2=N2*phi
    stdN=np.sqrt(avgN2-avgN**2)
    relstdN=stdN/avgN

    if avgN<=0 or phi<=0:
        base=np.zeros(len(props))
        #set press and variability to 1 by default
        ref=7
        for p in props:
            if p in default:
                base[props.index(p)]=default[p]

        props=[prefix+p for p in props]
        return pd.Series(dict(zip( props ,
            base)))

    #sigma_l=(sigma_k/sigma/avgN)
    #meanAN= mu - (1+ gamma)/(h + mu/sigma) * (1+ h / ( phi * (q + sigma_l**2) ) )
    #meanAN/=S * phi
    meanAN=0
    u=(1.-mu/Seff)
    ut=(u-phi* gamma*v*sigma**2)
    press=1./( ut**2 - phi*sigma**2  )
    chi2=phi*sigma**2*press
    utilde= ut

    w=(1 -np.sqrt(1-gamma*sigma**2*phi))/(gamma*sigma**2*phi)
    C=1/(1-mu/Seff -(gamma+1)*sigma**2 *phi/(2-sigma**2 *gamma *w *phi) - 1./S *
         (2 * phi * mu**2)/(2-sigma**2*phi*gamma*w - 2*mu*phi ) )

    U=1-mu/Seff
    IMF=(U+phi*mu)/(2+ (mu*phi-1)/U )
    ISIG=sigma**2 * phi*(gamma+1.)/(2*U)
    varJeff= 1/(IMF - ISIG)

    variability=varJeff
    extinct=0
    simpson=N2/S/phi/N1**2
    #extinct,deltaN,deltaN2=cavity_extinction(N0=avgN,S=S,mu=mu,sigma=sigma,gamma=gamma,v=v,phi=phi,avgN=avgN,varN=stdN**2)

    # <N><K> + <N>^2 sig (u n^2 - h n)/(1+q/sigl^2)
    prod= (ut*avgN2 + mu*avgN**2 +sigma**2/zeta**2 *avgN2 * avgK *avgN
        )/ (1+sigma**2/zeta**2*avgN2)

    props=[prefix+p for p in props]
    return pd.Series(dict(zip( props ,
        (v,phi,chi2,press,variability,C,utilde,avgN,avgN2,relstdN,stdN,-meanAN, avgN*S,simpson,1./simpson,prod,S*extinct) )))

def cavity_calc_props(prefix='',S=None,mu=None,sigma=None,
            sigma_k=None,gamma=None,q=None,v=None,h=None,phi=None,conn=1,**kwargs):

    props=('q','v','h','phi','chi2','press','variability','varcavity','utilde',
        'avgN','avgN2','relstdN','stdN','corr_AN','totN','simpson','simpsoninv','prod','extinction')


    default={
        'variability':1,
        'prod':1,
        'press':1,
        }

    Seff=S*conn
    sol=q,v,h,phi
    import pandas as pd
    avgK=kwargs.get('avgK',1)
    avgN=1./(sigma * h + mu)
    if avgK>0:
        avgN*=avgK
    else:
        print 'DANGER: avgK<=0!'
    relstdN=sqrt(max(0,q-1.))
    stdN=relstdN*avgN

    avgN2=q*avgN**2

    if avgN<=0 or phi<=0:
        base=np.zeros(len(props))
        #set press and variability to 1 by default
        ref=7
        for p in props:
            if p in default:
                base[props.index(p)]=default[p]

        props=[prefix+p for p in props]
        return pd.Series(dict(zip( props ,
            base)))

    sigma_l=(sigma_k/sigma/avgN)
    meanAN= mu - (1+ gamma)/(h + mu/sigma) * (1+ h / ( phi * (q + sigma_l**2) ) )
    meanAN/=S * phi
    u=(1-mu/Seff)/sigma
    ut=(u-gamma*v)
    chi2=phi/( ut**2 - phi )
    press=chi2/phi/sigma**2
    press=1/( ut**2*sigma**2 - phi*sigma**2 )
    utilde= ut
    #variability = ut*sigma /(ut**2*sigma**2 - phi*sigma**2 *(1+gamma)/2 )

    w=(1 -np.sqrt(1-gamma*sigma**2*phi))/(gamma*sigma**2*phi)
    #Cold=1/(1-mu/S -(gamma+1)*sigma**2 *phi/(2-sigma**2 *gamma *w *phi) )
    C=1/(1-mu/Seff -(gamma+1)*sigma**2 *phi/(2-sigma**2 *gamma *w *phi) - 1./S *
         (2 * phi * mu**2)/(2-sigma**2*phi*gamma*w - 2*mu*phi ) )

    U=1-mu/Seff
    IMF=(U+phi*mu)/(2+ (mu*phi-1)/U )
    ISIG=sigma**2 * phi*(gamma+1.)/(2*U)
    varJeff= 1/(IMF - ISIG)

    variability=varJeff
    extinct=0
    #extinct,deltaN,deltaN2=cavity_extinction(N0=avgN,S=S,mu=mu,sigma=sigma,gamma=gamma,v=v,phi=phi,avgN=avgN,varN=stdN**2)

    # <N><K> + <N>^2 sig (u n^2 - h n)/(1+q/sigl^2)
    prod= avgN*kwargs.get('avgK',1) + avgN**2 * sigma * (ut *q -h )/(1+ q/sigma_l**2)

    props=[prefix+p for p in props]
    return pd.Series(dict(zip( props ,
        sol+(chi2,press,variability,C,utilde,avgN,avgN2,relstdN,stdN,-meanAN, avgN*S, q/S,S/q,prod/avgN,S*extinct) )))

def cavity_solve_props(S=None,mu=None,sigma=None,
            sigma_k=None,gamma=None,conn=1,**kwargs):
    #print 'RECEIVING',S,mu,sigma,sigma_k,gamma,kwargs

    solver=kwargs.get('solver',cavity_solve)
    q,v,h,phi=solver(S=S*conn,mu=mu,sigma=sigma,sigma_k=sigma_k,gamma=gamma,**kwargs)
    return cavity_calc_props(S=S,mu=mu,sigma=sigma,sigma_k=sigma_k,gamma=gamma,q=q,v=v,h=h,phi=phi,conn=conn,**kwargs)

#=================== ANALYTICAL PLOTS ===============================


def cavity_loop(cols=None,**attr):
    axes={}
    log=attr.pop('log',['avgN','totN'])
    for i in attr.keys():
        if i in ('cols',):
            continue
        if np.array(attr[i]).shape:
            axes[i]=attr[i]
            del attr[i]
    labx,laby=axes.keys()
    xs,ys=axes[labx],axes[laby]
    if cols is None:
        cols=cavity_dico.keys()

    prev=(1.,1.,1.,1.)
    mat = np.zeros( (len(xs),len(ys), len(cols) ) )
    for ix,x in enumerate(xs):
        for iy,y in enumerate(ys):
            attr[labx]=x
            attr[laby]=y
            #print attr
            if not 'S' in attr:
                attr['S']=attr['size']
            props=cavity_solve_props(x0=prev,**attr )
            for p in props:
                if p in log:
                    props[p]=np.log10(props[p])
            phi,avgN=props['phi'], props['avgN']
            if phi==0 or avgN<=0:
                props['avgN']=100
                props['totN']=100
                ##ERROR: Divergence in the code
                props['chi2']=-4
                props['prod']=1
            else:
                prev=[props[g] for g in ('q','v','h','phi')]

            mat[ix,iy]=[props[col] for col in cols]
            #if avgN<0:
                #mat[ix,iy]=0
    for i, col in enumerate(cols):
        if col in ('avgN','totN'):
            y=mat[:,:,i]
            y[y>=100]=np.max(y[y!=100])
        #if col==None: # I DON"T REMEMBER WHAT THIS WAS SUPPOSED TO BE
            #z[z==-100]=min(0,np.min(z[z!=-100]))
    return mat

def cavity_show_loop_old(**attr):
    if attr.get('newfig',1):
        plt.figure()

    axes={}
    title=[]
    for i in sorted(attr.keys()):
        if i in ('cols',):
            continue
        if np.array(attr[i]).shape:
            axes[i]=attr[i]
        else:
            title.append('{}: {}'.format(i,attr[i]))
    title=attr.get('title','; '.join(title))
    if title:
        plt.suptitle(title)
    labx,laby=axes.keys()
    xs,ys=axes[labx],axes[laby]


    cols =attr.pop('cols',None)
    if cols is None:
        cols=cavity_dico.keys()
    attr['cols']=cols
    xs,ys=np.array([xs for y in ys]).T,np.array([ys for x in xs])
    mat=cavity_loop(**attr)

    nbpanels=mat.shape[-1]
    panel=MutableInt(0)
    for i in range(nbpanels):

        auto_subplot(plt,nbpanels,panel,rows=2)

        #plt.subplot(221+i)
        plt.xlabel(r'$\{}$'.format(labx))
        plt.ylabel(r'$\{}$'.format(laby))
        plt.title(cavity_dico[cols[i]]   )
        plt.ylim(ymin=np.min(ys),ymax=np.max(ys))
        plt.pcolor(xs,ys,mat[:,:,i])
        plt.colorbar()
    #plt.tight_layout(pad=0)
    plt.subplots_adjust(left=0.03,  right=0.97)



def cavity_show_loop(cols=['phi',#'avgN',#'stdN','press',
            'variability','totN','simpson'],maxN=2,minchi2=-8,**attr):
    if attr.get('newfig',1):
        plt.figure()

    log=attr.pop('log',['avgN','totN'])
    axes={}
    title=[]
    for i in sorted(attr.keys()):
        if i in ('cols',):
            continue
        if np.array(attr[i]).shape:
            axes[i]=attr[i]
        else:
            title.append('{}: {}'.format(i,attr[i]))
    title=attr.get('title','; '.join(title))
    if title:
        plt.suptitle(title)
    labx,laby=axes.keys()
    xs,ys=axes[labx],axes[laby]


    mat={c:np.zeros( (len(xs),len(ys)) ) for c in cols  }
    x0=None
    for ix,x in enumerate(xs):
        for iy,y in enumerate(ys):
            prm={}
            prm.update(attr)
            prm[labx]=x
            prm[laby]=y
            if not x0 is None:
                prm['x0']=x0
            if not 'S' in prm:
                prm['S']=float(prm['size'])
            res=cavity_solve_props(**prm)


            for p in res.keys():
                if p in log and res[p]>0:
                    res[p]=np.log10(res[p])


            if res['avgN']>maxN:
                res['avgN']=maxN
                res['totN']=maxN#*prm['S']
            if res['chi2']==0:
                res['chi2']=minchi2
            res['chi2']=np.clip(res['chi2'],minchi2,-minchi2)
            res['variability']=np.clip(res['variability'],0,3 )
            for c in mat:
                mat[c][ix,iy]=res[c]
            x0=[res[z] for z in ('q','v','h','phi')]


    xs,ys=np.array([xs for y in ys]).T,np.array([ys for x in xs])


    nbpanels=len(cols)
    panel=MutableInt(0)
    for i in range(nbpanels):

        auto_subplot(plt,nbpanels,panel,rows=attr.get('rows',None) )
        plt.xlabel(r'$\{}$'.format(labx))
        plt.ylabel(r'$\{}$'.format(laby))
        plt.title(cavity_dico[cols[i]])
        plt.ylim(ymin=np.min(ys),ymax=np.max(ys))
        res=mat[cols[i]]
        if np.min(res)>=0 and np.abs( np.max(res)/np.min(res[res>0]) )>100:
            res[res==0]=np.min(res[res!=0])
            res=np.log10(res)
            print cols[i]
        vmin=np.min(res)
        vmax=np.max(res)
        if cols[i]=='chi2':
            vmax=-minchi2
            vmin=minchi2
        plt.pcolor(xs,ys,res,vmin=vmin,vmax=vmax,cmap='plasma')
        plt.colorbar()
        #plt.figure()
    #plt.tight_layout(pad=0)
    plt.subplots_adjust(left=0.03,  right=0.97)

def cavity_distri(S=100,gamma=.9, sigma_k=2., mu=2.,sigma=.1):
    #Predicted abundance distribution
    q,v,h,phi=cavity_solve(S,mu,sigma,sigma_k,gamma )

    u=(1-mu/float(S))/sigma
    avgN=1./(sigma*h + mu)
    sigma_lambda=sigma_k/sigma/avgN
    import scipy.stats as SS

    #print q,v,h,phi
    utilde=u-gamma*v
    mean=h/utilde
    var=(q+sigma_lambda**2)/utilde**2
    gauss=SS.norm(loc=mean,scale=np.sqrt(var))
    dist= lambda x:gauss.pdf(x/avgN)/avgN

    return dist




def phi_vs_sigma(S=100,gamma=.9, sigma_k=2., mu=2.,sigma=.1,hold=0,do_plot=1,do_scale=0,
        sizes=np.logspace(.5,3.,80), sigmas=np.linspace(0.2,2.4,50)):
    from datatools import plt,plot
    plt.figure()

    #sigmas=sigma*np.sqrt(sizes)
    surv=[]
    resp=[]
    q,v,h,phi=(1.,1.,1.,1.)

    #print q,v,h,phi
    #for S in sizes:
    for sigma in sigmas:
        if do_scale:
            q,v,h,phi=cavity_solve(S,mu*S,sigma*np.sqrt(S),sigma_k,gamma,x0=(q,v,h,phi) )
        else:
            q,v,h,phi=cavity_solve(S,mu,sigma,sigma_k,gamma,x0=(q,v,h,phi) )

        #print S,phi
        surv.append(S*phi)
        #surv.append(phi)
        resp.append(v)
        #u=(1-mu/float(S))/sigma
        #resp.append(phi/( (u-gamma*v )**2 - phi) )
    if do_plot:
        plt.xlabel(r'$S$')
        plt.ylabel(r'$S^*$')
        plot(sizes,surv,hold=hold)
        #plot(sigmas,surv,hold=hold)
    #return sizes,surv,resp
    return sigmas,surv,resp

