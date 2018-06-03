# -*- coding: utf-8 -*-
import numpy as np, scipy, networkx as nx, scipy.linalg as la
from math import *
import bisect as bisct,random
import matplotlib.pyplot as plt
import cProfile,pstats,sys,time,itertools,re
import inspect
from copy import deepcopy

###======================== DATA ANALYSIS TOOLSET ======================###


profiler=cProfile.Profile()


def code_debugger(skip=0):
    import code
    import inspect
    stack=inspect.stack()
    dic = {}
    dic.update(stack[1+skip][0].f_globals)
    dic.update(stack[1+skip][0].f_locals)
    code.interact(local=dic)


def tupleify(measure):
    def totuple(a):
        if not hasattr(a,'__iter__'):
            return a
        try:
            return tuple(totuple(i) for i in a)
        except TypeError:
            return a
    if not hasattr(measure,'iteritems'):
        if hasattr(measure,'__iter__'):
            return measure.__class__([tupleify(m) for m in measure])
        return measure
    for i,j in measure.iteritems():
        if isinstance(j, np.ndarray):
            measure[i]='array({})'.format(totuple(j))
        else:
            measure[i]=totuple(j)
    return measure

def retupleify(df,recursive=1,verbose=0):
    '''Tuples are recovered as strings by read_csv, convert them to tuples'''
    import ast, pandas as pd
    pandas_mode= isinstance(df,pd.DataFrame    )
    for z in df:
        #print z, df[z][0],isinstance(df[z][0],basestring) and '(' in df[z][0]
        if verbose:
            print "Retupleifying",z
        if pandas_mode:
            val=[x for x in df[z] if isinstance(x,basestring) and '(' in x]
            if val:
                val=val[0]
        else:
            val=df[z]
        if isinstance(val,basestring) and '(' in val:
            def safeval(v):
                dic={'nan':np.nan,'inf':np.inf}
                if not isinstance(v,basestring):
                    return v
                if 'array' in v or 'list' in v:
                    v=v.replace('array','list')
                    ls= eval(v,{},dic)
                    ls= [safeval(l) for l in ls ]
                    return np.array(ls)
                return eval(v,{},dic)

            try:
                if pandas_mode:
                    df[z]=df[z].map(safeval)
                else:
                    df[z]=ast.safeval(df[z])
                #print "NOW",type(df[z][0])
            except Exception as e:
                print 'retupleify:',e
        if recursive:
            if hasattr(val,'keys') or isinstance(val,pd.DataFrame):
                df[z]=retupleify(df[z])
    return df


def kth_diag_indices(a, k):
    rows, cols = np.diag_indices_from(a)
    if k < 0:
        return rows[-k:], cols[:k]
    elif k > 0:
        return rows[:-k], cols[k:]
    else:
        return rows, cols


class MutableInt(object):
    def __init__(self,val):
        self.val = val
    def __add__(self,i):
        return self.val+i
    def __eq__(self,i):
        return self.val==i


def auto_subplot(plt,nbpanels,panel,projection=None,rows=None,return_all=0):
    i=panel.val
    if rows is None:
        panels_per_row=np.ceil(np.sqrt(nbpanels) )
    else:
        panels_per_row=np.ceil(nbpanels/rows).astype('int')
    nbrows=np.ceil(nbpanels*1./panels_per_row)
    ax=plt.subplot(100*nbrows+panels_per_row*10+i +1 ,projection=projection)
    panel.val+=1
    if return_all:
        return ax,nbrows,panels_per_row
    return ax



def gaussian(x,mean,var):
    return np.exp(-(x-mean)**2/2/var)/np.sqrt(2*np.pi*var)

def lognormal(x,mean,var):
    return np.exp(-(np.log(x)-mean)**2/2/var)/np.sqrt(2*np.pi*var)/x

def ifelse(cond,res,els):
    if cond:
        return res
    return els



# Public Domain, i.e. feel free to copy/paste
# Considered a hack in Python 2
def debug_caller_name(skip=2):
    """Get a name of a caller in the format module.class.method
    `skip` specifies how many levels of stack to skip while getting caller
    name. skip=1 means "who calls me", skip=2 "who calls my caller" etc.
    An empty string is returned if skipped levels exceed stack height
    """
    stack = inspect.stack()
    start = 0 + skip
    if len(stack) < start + 1:
        return ''
    parentframe = stack[start][0]
    name = []
    module = inspect.getmodule(parentframe)
    # `modname` can be None when frame is executed directly in console
    if module:
        name.append(module.__name__)
        # detect classname
    if 'self' in parentframe.f_locals:
        # I don't know any way to detect call from the object method
        # XXX: there seems to be no way to detect static method call - it will
        # be just a function call
        name.append(parentframe.f_locals['self'].__class__.__name__)
    codename = parentframe.f_code.co_name
    if codename != '<module>': # top level usually
        name.append( codename ) # function or a method
    del parentframe
    return ".".join(name)


def prolog(fname,col=2):
    stats=[[i.code ,i.totaltime,i.inlinetime,i.callcount,i.reccallcount] for i in profiler.getstats()]
    stats=sorted(stats,key=lambda e:e[col],reverse=1)
    with open(fname,'w') as prolog:
        for i in stats:
            if not i:
                continue
            st=' '.join([str(z) for z in  i])
            prolog.write(st+'\n')


def ADD_FOLDER_TO_PATH(folder,relative=1,depth=1):
    import os, sys, inspect
    if relative:
        cmd_folder = os.path.realpath(os.path.abspath(os.path.split(inspect.getfile(
            inspect.currentframe(depth) ))[0])) +folder
    else:
        cmd_folder = os.path.realpath(os.path.abspath(folder))
    if cmd_folder not in sys.path:
        sys.path.insert(0, cmd_folder)



class Path(str):
    '''Strings that represent filesystem paths.
    Overloads __add__:
     - when paths are added, gives a path
     - when a string is added, gives a string'''
    def __add__(self,x):
        import os
        if isinstance(x,Path):
            return Path(os.path.normpath(os.path.join(str(self),x)))
        return os.path.normpath(os.path.join(str(self),x))

    def norm(self):
        import os
        return Path(os.path.normpath(str(self)))

    def osnorm(self):
        """Deal with different separators between OSes."""
        import os
        if os.sep=='/' and "\\" in str(self):
            return Path(os.path.normpath(str(self).replace('\\','/' )))
        elif os.sep=='\\' and "/" in str(self):
            return Path(os.path.normpath(str(self).replace('/','\\' )))
        else:
            return self.norm()

    def prev(self):
        import os
        lst=self.split()
        path=os.path.join(lst[:-1])
        return path.osnorm()

    def split(self):
        """"""
        import os
        lst=[]
        cur=os.path.split(self.norm())
        while cur[-1]!='':
            lst.insert(0,cur[-1])
            cur=os.path.split(cur[0])
        return lst

    def mkdir(self,rmdir=False):
        """Make directories in path that don't exist. If rmdir, first clean up."""
        import os
        if rmdir:
            os.rmdir(str(self))
        cur=Path('./')
        for intdir in self.split():
            cur+=Path(intdir)
            if not os.path.isdir(cur):
                os.mkdir(cur)

    def copy(self):
        return Path(self)

    def strip(self):
        '''Return string without final / or \\ to suffix/modify it.'''
        return str(self).strip('\/')

def wplot(typ,*args,**kwargs):
    USE_MYHIST=kwargs.get('USE_MYHIST',1)
    if kwargs.pop('newfig',0):
        plt.figure()
    save=kwargs.pop('save',0)
    hold=kwargs.pop('hold',False)
    title=kwargs.pop('title',False)
    legs=kwargs.pop('legends',False)
    invert=kwargs.pop('invert',0)
    projection=kwargs.pop('projection','2d')
    if not legs:
        legs=kwargs.pop('legend',False)
    if not legs:
        legs=kwargs.pop('leg',False)
    lkw=kwargs.pop('legopt',{})
    if 'xlabel' in kwargs:
        plt.xlabel(kwargs.pop('xlabel'))
    if 'ylabel' in kwargs:
        plt.ylabel(kwargs.pop('ylabel'))
    if 'labels' in kwargs:
        labs=kwargs.pop('xlabel')
        plt.xlabel(labs[0])
        plt.ylabel(labs[1])

    ax=plt.gca()

    if 'log' in kwargs:
        lg=kwargs.pop('log')
        if not isinstance(lg,str):
            if typ!= 'hist':
                kwargs['log']='xy'
            else:
                kwargs['log']=lg
        else:
            if 'x' in lg:
                ax.set_xscale('log')
            if typ!= 'hist' or USE_MYHIST:
                if 'y' in lg :
                    ax.set_yscale('log')
                if 'z' in lg:
                    ax.set_zscale('log')
            if typ=='hist':
                if 'y' in lg:
                    kwargs['log']=1
                if 'x' in lg:
                    if min(args[0])<=0:
                        args=([ar for ar in args[0] if ar>0],)+args[1:]
                    if not 'bins' in kwargs or isinstance(kwargs['bins'],int):
                        kwargs['bins']=np.logspace(log10(min(args[0])),
                            log10(max(args[0])),kwargs.get('bins',100),10)
    elif typ=='hist' and min(args[0])>0 and max(args[0])/min(args[0])>10**3:
        kwargs['log']=1
    if 'xs' in kwargs :
        xs=kwargs.pop('xs')
        args=list(itertools.chain.from_iterable([[xs[:len(a)],a] for a in args]))
    if not args:
        return
    handle =None
    force_zlim=kwargs.pop('force_zlim',False)

    while handle is None:
        try:
            if typ=='plot':
                handle=plt.plot(*args,**kwargs)
            if typ =='hist':
                if not 'bins' in kwargs:
                    kwargs['bins']=100
                if kwargs.pop('cleverbins',False):
                    if isinstance(kwargs['bins'],int):
                        kwargs['bins']=min(kwargs['bins'],len(args[0])/20)
                kwargs['hold']=1
                kwargs['density']=kwargs.pop('density',kwargs.pop('normed',1))
                if USE_MYHIST:
                    def myhist(a,**kws):
                        kw={}
                        plt.plot([],[])
                        kw.update(kws)
                        h=None
                        while h is None:
                            try:
                                h,b=np.histogram(a,**kw)
                            except TypeError as e:
                                # print e
                                del kw[str(e).strip().split(' ')[-1].replace("'",'').replace('"','') ]
                    #(b[1:]+b[:-1])/2,
                        handle=None

                        kw={}
                        kw.update(kws)
                        histtype=kw.get('histtype',None)

                        if kwargs.get('cumulative',0):
                            h=np.cumsum(h)
                            h/=h[-1]
                        while handle is None:
                            try:
                                if histtype=='step' and 0:
                                    x=b[:-1]
                                    y=h[:]
                                    if kw.get('log',0):
                                        x,y=x[x>0],y[x>0]
                                        x,y=x[y>0],y[y>0]
                                    handle=plt.plot( x, y, **kw )
                                else:
                                    handle=plt.bar(  b[:-1], h,b[1:]-b[:-1], **kw )
                            except AttributeError as e:
                                # print e
                                del kw[str(e).strip().split(' ')[-1].replace("'",'').replace('"','') ]
                        # plt.ylim(ymin=np.min(h[h>0]) )
                        # plt.autoscale(enable=1,axis='y')
                        return handle
                    handle=[myhist(a,**kwargs) for a in args ]
                else:
                    handle=[plt.hist(a,**kwargs)[2][0] for a in args]

            if typ=='scatter':
                nbaxes=2
                if projection=='3d':
                    nbaxes=3
                if len(args)>nbaxes:
                    handle=[]
                    for idx in range(len(args)/nbaxes):
                        handle.append(ax.scatter(*args[idx:idx+nbaxes],**kwargs))
                else:
                    handle=ax.scatter(*args,**kwargs)
                if force_zlim:
                    Z=args[2]
                    ax.set_zlim(bottom=min(Z), top=max(Z))

            if typ=='error':
                handle=plt.errorbar(*args,**kwargs)
            if typ=='errorfill':
                from mpltools import special as pltspecial
                handle=pltspecial.errorfill(*args,**kwargs)
            if typ=='contour':
                handle=plt.contour(*args,**kwargs)
            if typ=='contourf':
                handle=plt.contourf(*args,**kwargs)
            if typ=='bar':
                handle=plt.bar(*args,**kwargs)
        except AttributeError as e:
            #print e
            del kwargs[str(e).strip().split(' ')[-1] ]

        #except e:
            #print "Cannot plot:", sys.exc_info()[0],e
            #return

    if kwargs.get('colorbar',0):
        plt.colorbar(handle)
    if title:
        plt.title(title)
    if legs:
        plt.legend(handle,legs,**lkw)
    if invert:
        if 'x' in invert:
            plt.xlim(plt.xlim()[::-1])
        elif 'y' in invert:
            plt.ylim(plt.ylim()[::-1])
    if save:
        plt.savefig(save)
        plt.clf()
        plt.close()
        return
    elif not hold:
        plt.show()
    return handle

def errorbar(*args,**kwargs):
    return wplot('error',*args,**kwargs)
def errorfill(*args,**kwargs):
    return wplot('errorfill',*args,**kwargs)
def plot(*args,**kwargs):
    return wplot('plot',*args,**kwargs)
def hist(*args,**kwargs):
    return wplot('hist',*args,**kwargs)
def scatter(*args,**kwargs):
    return wplot('scatter',*args,**kwargs)
def scatter3d(*args,**kwargs):
    ax=plt.gcf().add_subplot(111,projection='3d')
    return wplot('scatter',*args,projection='3d',**kwargs)
def contour(*args,**kwargs):
    return wplot('contour',*args,**kwargs)
def contourf(*args,**kwargs):
    return wplot('contourf',*args,**kwargs)
def bar(*args,**kwargs):
    return wplot('bar',*args,**kwargs)

def cumulplot(*args,**kwargs):
    x=sorted(args[0])
    cx = np.linspace(0, 1, len(x))
    return wplot(kwargs.pop('plottype','plot'),x,cx,*args[1:],**kwargs)

def make_double_loghist(xs,ys,**opts):
    '''Plot histogram where xs is in logscale both for xs<0 and xs>0 in two adjacent panels'''
    maxxs=np.max(np.log10(np.abs(xs)))
    minxs=np.log10(np.max([np.min(xs[xs>0]),np.min(-xs[xs<0]) ]))

    ysp,ysm=ys[xs>0],ys[xs<0]
    minys,maxys=np.max([np.min(ysp[ysp>0]), np.min(ysm[ysm>0])]),np.max(ys)#*1.2
    fig=plt.subplot(121)
    plt.subplots_adjust(wspace=0, hspace=0)
    fig.set_xlim(xmin=-maxxs,xmax=-minxs)
    fig.set_ylim(ymax=maxys,ymin=minys)
    plot(-np.log10(-xs[xs<0]),ysm,hold=1,log='y',**opts)

    fig.axes.set_xticklabels([r'$-10^{{{}}}$'.format(-float(x)) for x in fig.axes.get_xticks()])

    fig=plt.subplot(122)
    fig.set_xlim(xmax=maxxs,xmin=minxs)
    fig.set_ylim(ymax=maxys,ymin=minys)
    fig.axes.get_yaxis().set_visible(False)
    plot(np.log10(xs[xs>0]),ysp,hold=1,log='y',**opts)
    fig.axes.set_xticklabels([r'$10^{{{}}}$'.format(float(x)) for x in fig.axes.get_xticks()])


def mhist(*data,**kwargs):
    if len(data)>1:
        return [mhist(x,**kwargs) for x in data]
    else:
        data=data[0]
    kwargs.setdefault('log',0)
    if kwargs['log']:
        if min(data)<=0:
            data=[ar for ar in data if ar>0]
        if not 'bins' in kwargs or isinstance(kwargs['bins'],int):
            rg=kwargs.get('range', ( max(10**-50,min(data)) ,max(10**-50,max(data))) )
            kwargs['bins']=np.logspace(log10(rg[0] ),
                log10(rg[1]),kwargs.get('bins',100),10)
    bhist,bbins=np.histogram(data,bins=kwargs.pop('bins',200),density=kwargs.pop('normed',1),
        range=kwargs.get('range',None))
    wbins=(bbins[1:]-bbins[:-1])
    bbins=(bbins[:-1]+bbins[1:])/2
    if kwargs.get('cumulative',0):
        if kwargs['cumulative']>0:
            bhist=np.cumsum(bhist)
        else:
            bhist=np.sum(bhist)-np.cumsum(bhist)
    if kwargs.get('typ',False):
        wplot(kwargs['typ'],bbins,bhist,**kwargs)
    bhist[bhist<10**-50]=0
    if kwargs.get('width',0):
        return bbins,wbins,bhist
    return bbins,bhist


def powlaw(expo,xmin,xmax,shape=None):
    #distribution proportional to x**expo
    xx=expo+1
    if not xx:
        return exp(log(xmin) + log(xmax/xmin)*np.random.random(shape))
    return (xmin**xx + (xmax**xx-xmin**xx)*np.random.random(shape))**(1./xx)


def rint(nb):
    return int(round(nb))


def make_density(*data,**kwargs):
    import scipy.stats as stats

    from numpy import array
    if len(data)==1:
        data=data[0]
    mat=array(data,dtype="float")
    if mat.shape[1]==2 :
        #data is couples (x,y)
        mat=mat.T
    normed=kwargs.get('normed',False)
    nbin=kwargs.get('bins',100)
    logy=kwargs.get('logy',False)
    remove0=1-kwargs.get('include_zeroes',1)
    bintype='linear'
    xmin,xmax=np.min(mat[0]),np.max(mat[0])
    bins=None
    if kwargs.get('logx',False):
        numat=mat[0][mat[0]>10**(-100)]
        xmin,xmax=min(nozero(nonan(numat))),max(nonan(numat))
        bintype='log'
        bins=np.logspace(log10(xmin),log10(xmax),nbin,False)
        binw=[bins[i+1]-bins[i] for i in xrange(nbin-1)]
        binw.append(xmax-bins[-1])
    xspan=xmax-xmin
    if bins is None:
        #linear spacing case
        binw=xspan/nbin
        bins=np.linspace(xmin,xmax,nbin,False)
    binnage=[[] for i in xrange(nbin)]
    for x,y in mat.T :
        if x<xmin or x>xmax:
            continue
        if bintype=='linear':
            xbin=int(floor(float(x-xmin)/binw))
        else :
            xbin=bisct.bisect_left(bins,x)
        if xbin ==nbin: #maxvalue
            xbin=nbin-1
        if not remove0 or abs(y)>10**(-40):
            binnage[xbin].append(y)
    res=array([stats.describe(i)[2:]+(min(i),max(i)) if i else stats.describe([0])[2:]+(0,0) for i in binnage])
    sspercen=scipy.stats.scoreatpercentile
    if kwargs.get('relative',1):
        quantile=array([array([sspercen(i,50),sspercen(i,50)-sspercen(i,5),sspercen(i,95)-sspercen(i,50)]) if i else array([0,0,0]) for i in binnage])
        res2=array([-res[:,-2]+res[:,0],res[:,-1]-res[:,0]])
    else:
        quantile=array([array([sspercen(i,50),sspercen(i,5),sspercen(i,95)]) if i else array([0,0,0]) for i in binnage])
        res2=array([res[:,-2],res[:,-1]])
    quantile=quantile.T
    if normed :
        if bintype=='linear':
            res[:,0]/=sum(res[:,0])*binw
        else :
            sm=np.dot(res[:,0],np.array(binw))
            res[:,0]/=binw/sum(binw)#sum(sm)
    return bins,res[:,0],res[:,1],res2,quantile[0],quantile[1:],array([len(i) for i in binnage])

def plot_density(*data,**kwargs):
    bins,res,std=make_density(*data,**kwargs)[:3]
    kwargs.pop('bins',None)
    kwargs.pop('relative',None)
    kwargs.pop('include_zeroes',None)
    kwargs.pop('normed',None)
    kwargs.pop('logx',None)
    kwargs.pop('logy',None)
    errorbar(bins,res,yerr=std,**kwargs)

def nonan(mat):
    return mat[np.isnan(mat)==False]
def nozero(mat):
    return mat[mat!=0]
def noinf(mat):
    return mat[np.logical_and(mat!=np.inf,mat!=-np.inf)]

