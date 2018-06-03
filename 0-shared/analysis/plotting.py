# -*- coding: utf-8 -*-
import os, sys, inspect

from imports import *
import numpy as np
from datatools import plot,hist, scatter
import matplotlib.pyplot as plt
from itertools import product as iproduct



def data_to_matrix(measures=None,path=None,axes=None,values=None ,
        split_by=None,select=None,
        mode=None,**kwargs):
    """Create index+matrix from measures in path. Take indices from 'axes',
    and matrix content from 'values'.
    If split_by is a tuple of column names, instead returns two dictionaries where
    keys are every tuple of possible values of each column to split by, associated
    with respective indices and matrices."""
    values=[val if not isinstance(val,tuple) else val[0] for val in values ]

    aggregator=kwargs.pop('aggregator','mean')
    if isinstance(aggregator,basestring):
        aggregator=lambda e,a=aggregator: getattr(e,a)()


    def binner(idx,measures):
        bins=[]

        for i in idx :
            m=measures
            if hasattr(i, '__iter__'):
                for col,val in zip(axes,i):
                    m=m[m[col]==val]
            else:
                m=m[m[axes[0]]==i]
            bins.append({k:m[k] for k in values } )
        return bins

    if measures is None:
        measures=open_measures(path)
    if not select is None:
        for i,j in select.iteritems():
            measures=measures[ measures[i] == j]

    if hasattr(axes[0],'__iter__'):
        print "AXES",axes
        raise

    measures=measures[ [ax for axlist in (axes,values,split_by)
         if hasattr(axlist,'__iter__')  for ax in axlist ]]


    #Remove tuples for aggregation; retain old measure
    oldmeasures=measures[:]
    for col in measures:
        if isinstance(measures.iloc[0][col],tuple):
            #try:
                #measures[col]=map(lambda x:np.product(x),measures[col])
            #except:
                measures[col]=map(lambda x:x[0],measures[col])


    vals=[v if not isinstance(v,tuple) else v[0] for v in values]

    if split_by is None:

        grp= aggregator(measures.groupby(axes))[vals]#.as_matrix()
        idx=np.array(grp.index.tolist())
        if mode =='bin':
            return idx,binner(idx,oldmeasures)
        else:
            return idx, grp.as_matrix()
    else:
        result={}
        index={}
        keys=[tuple(zip(split_by,val)) for val in iproduct(*[np.unique(measures[j]) for j in split_by]) ]
        #print 'data_to_matrix', keys,[val for val in iproduct(*[measures[j].unique() for j in split_by])]
        for k in keys:
            meas=measures.copy()
            for col,val in k:
                try:
                    meas=meas[meas[col]==val]
                except:
                    meas = None
                    break
            if not meas is None:
                grp=aggregator(meas.groupby(axes))[vals]
                idx=np.array(grp.index.tolist())
                index[k]=idx
                if mode =='bin':
                    result[k]=binner(idx,oldmeasures)
                else:
                    result[k]=grp.as_matrix()
        return index, result


def show_hist(axes=None,values=None,**kwargs):
    idx,binned=data_to_matrix(axes=axes,values=values,mode='bin',**kwargs)
    if kwargs.get('newfig',1):
        plt.figure()
    from datatools import mhist
    nbpanels=len(values)
    panel=MutableInt(0)
    dico=kwargs.get('dictionary',{})
    def get_dico(val):
        return dico.get(val,val)
    results={}
    for val in values:
        if nbpanels>1:
            auto_subplot(plt,nbpanels,panel)
        plt.xlabel(get_dico(val))
        plt.ylabel('Frequency')
        if kwargs.get('split_by',0):
            raise Exception("Not done yet")
        else:
            coords=[]
            legends=[]
            for i,binn in zip(idx,binned):
                print 'showhist',i, len(binn), len(binn[val])
                lst=[ll for l in binn[val] for ll in l]
                #for lst in binn[val]:
                    #print len(lst)
                xs,ys=mhist(lst,bins=kwargs.get('bins',30))
                coords.append( (xs,ys) )
                legends.append('{}'.format(get_dico(i)) )
            for c in coords:
                plot(c[0],c[1],hold=1,linestyle='None',marker='o',**{k:kwargs[k] for k in kwargs if
                    k in ('log',) })
            plt.legend(legends)
            results[val] = coords
    return results

def make_figure(idx,mat,title='',newfig=1,axes=None,values=None,runavg=None,
        style='heatmap',log='',zrange=None,logcutoff=10**-10,**kwargs):
    show_labels=kwargs.get('show_labels',1)
    if len(values)>1 or isinstance(values[0],tuple):
        #print mat.shape,idx.shape
        for  iv,val in enumerate(values):
            if isinstance(val,tuple):
                val,prm=val
            else:
                prm={}
            kw={}
            kw.update(kwargs)
            kw['style']=style
            kw['log']=log
            kw['zrange']=zrange
            kw['title']=title
            kw.update(prm)
            #print val, kw
            make_figure(idx,mat[[slice(None) for x in mat.shape[:-1] ]+[iv]  ],
                newfig=newfig,axes=axes,values=[val],**kw)
            newfig=False
        return

    dico=kwargs.pop('dictionary',{})
    def get_dico(val):
        return dico.get(val,val)

    if newfig:
        fig=plt.figure()
    else:
        fig=plt.gcf()
    color=kwargs.pop('color',kwargs.get('c','b'))
    if len(axes)==1:
        if title:
            plt.title(title)
        X=idx
        Y=mat
        if runavg:
            #Running average
            from datatools import runavg as average
            X=average(X,runavg)
            Y=average(Y,runavg)

        if 'y' in log:
            X=X[Y.squeeze()>logcutoff]
            Y=Y[Y>logcutoff]
        if show_labels:
            plt.xlabel(get_dico(axes[0]))
            plt.ylabel(get_dico(values[0]))
        if style =='plot':
            mk=kwargs.pop('marker',None)
            plot(X,Y,hold=1,color=color,log=log,**kwargs)
        elif style=='scatter':
            #print color, kwargs
            if 'marker' in kwargs and kwargs['marker'].get_fillstyle()=='none':
                    kwargs['c']='none'
                    kwargs['edgecolor']=color
            else:
                kwargs['c']=color
            scatter(X,Y,hold=1,log=log, **kwargs)

    elif len(axes)==2:
        if title:
            title='{}: '.format(get_dico(title))
        X,Y=idx.T
        Z=mat
        if runavg:
            print 'WARNING: running average not ready for 3d.'
            from datatools import runavg as average
            X=average(X,runavg)
            Y=average(Y,runavg)
            Z=average(Z,runavg)

        xnb=len(set(X.ravel()))
        ynb=len(set(Y.ravel()))
        shape=(xnb,ynb)
        if xnb*ynb!= X.shape[0]:
            if style=='wireframe':
                print 'Changing style to scatter',values, shape,X.shape
                color='k'
            style='scatter'
        if style =='heatmap':
            #X,Y,Z=xs,ys,H
            if 'x' in log:
                plt.xscale('log')
            if 'y' in log:
                plt.yscale('log')
            if show_labels:
                plt.xlabel(get_dico(axes[0]))
                plt.ylabel(get_dico(axes[1]))
            plt.title('{}{}'.format(title,values[0]))
            plt.pcolor(X.reshape(shape),Y.reshape(shape),Z.reshape(shape))
            plt.colorbar()
            #scatter(X,Y,c=Z,s=300,log='x')
        elif style =='wireframe':
            if newfig:
                ax=fig.add_subplot(111,projection='3d')
            else:
                ax=plt.gca()
            if show_labels:
                plt.xlabel(get_dico(axes[0]))
                plt.ylabel(get_dico(axes[1]))
            plt.title('{}{}'.format(title,values[0]))
            ax.set_zlim(bottom=min(Z), top=max(Z))
            if zrange:
                ax.set_zlim(bottom=zrange[0], top=zrange[1])
            #X,Y,Z=xs,ys,H
            try:
                ax.plot_wireframe(X.reshape(shape),Y.reshape(shape),Z.reshape(shape),
                    )
            except Exception as e:
                print 'COULD NOT PLOT',values,title, X.shape, Y.shape,shape
                print e
        else:
            if newfig:
                ax=fig.add_subplot(111,projection='3d')
            else:
                ax=plt.gca()
            if show_labels:
                plt.xlabel(get_dico(axes[0]))
                plt.ylabel(get_dico(axes[1]))
            plt.title('{}{}'.format(title,values[0]))
            ax.set_zlim(bottom=np.min(Z), top=np.max(Z))
            if zrange:
                ax.set_zlim(bottom=zrange[0], top=zrange[1])
            if log:
                ax.set_zscale('log')
            #X,Y,Z=xs,ys,H
            if color is None:
                color=Z
            ax.scatter(X,Y,Z,c=color)

def plot_function(path=None,measures=None,axes=None,function=None,values=None,split_by=None,
        select=None ,hold=0,**kwargs):
    '''Passed function should return pandas.Series objects'''

    import pandas as pd
    if measures is None:
        measures=open_measures(path)

    newcols=measures.apply(function,axis=1)
    if values is None:
        values=list(newcols)[:1]
    elif not hasattr(values,'__iter__'):
        values=[values]
    #print newcol
    measures=measures.join(newcols)
    #print measures
    idx,mat=data_to_matrix(measures=measures,axes=axes,values=values,
        split_by=split_by,select=select,**kwargs)
    kwargs.setdefault('style','wireframe')
    if split_by:
        for key in sorted(idx):
            make_figure(idx[key],mat[key],title=key,axes=axes,values=values,**kwargs)
    else:
        make_figure(idx,mat,axes=axes,values=values,**kwargs)

    if not hold:
        plt.show()




def plot_data(path=None,measures=None,axes=None,values=None,split_by=None,select=None ,
        hold=0,**kwargs):
    """Create matrix using 'data_to_matrix', then plots in either 2d and 3d."""
    if None in (axes,values):
        raise Exception( 'plot_data is missing arguments.')

    if len(axes)>2:
        print('ERROR: plot_data cannot plot more than 2 axes. Taking only first 2 axes.')
        axes=axes[:2]



    if measures is None:
        measures=open_measures(path)
    functions=set([v[1].pop('function') for v in values if isinstance(v,tuple) and 'function' in v[1]  ])
    for function in functions:
        newcols=measures.apply(function,axis=1)
        if values is None:
            values=list(newcols)[:1]
        elif not hasattr(values,'__iter__'):
            values=[values]
        #print newcol
        measures=measures.join(newcols)

    values =to_iterable(values)
    if not hasattr(values,'__iter__'):
        values=[values]

    kwargs.setdefault('style','scatter')
    idx,mat=data_to_matrix(measures=measures,axes=axes,values=[v[0] if isinstance(v,tuple)
        else v for v in values],split_by=split_by,select=select,**kwargs)

    dico=kwargs.get('dictionary',{})
    def get_dico(val):
        return dico.get(val,val)
    title=kwargs.pop('title','')
    if not split_by is None:
        newfig=kwargs.pop('newfig',1)
        for ik,key in enumerate(sorted(idx)):
            kw={}
            kw.update(kwargs)
            kw['newfig']=newfig
            kw.setdefault('color',get_cmap(ik,len(idx)))
            from matplotlib.markers import MarkerStyle
            from itertools import product as iprod
            markers= sorted(list( iprod(['o',  'X', '^', 's', 'p', 'v','P'][:(len(idx)+1)/2], ('full','none')  )), key=lambda e:e[1],reverse=1)*50
            kw.setdefault('marker',MarkerStyle(marker =markers[ik][0],fillstyle=markers[ik][1]  ) )
            make_figure(idx[key],mat[key],title=title,axes=axes,values=values,**kw)
            newfig=False
        if kwargs.get('show_legends',True):
            plt.legend([', '.join([r'{}={}'.format(get_dico(j[0]),get_dico(j[1]))
                    if get_dico(j[0]) else r'{}'.format(get_dico(j[1])) for j in i ])
                    for i in sorted(idx)])
    else:
        #print split_by, mat.keys()
        make_figure(idx,mat,axes=axes,values=values,title=title,**kwargs)

    if not hold:
        plt.show()


def show_data(path=None,measures=None,axes=None,values=None,subplots=1,**kwargs):
    if kwargs.get('newfig',1):
        plt.figure()
    if subplots:
        panel=MutableInt(0)
    for val in values:
        prm={}
        prm.update(kwargs)
        prm['newfig']=0
        if isinstance(val,tuple):
            prm.update(val[1])
            val=val[0]

        if subplots:
            auto_subplot(plt,len(values),panel)
            if val != values[0]:
                plt.legend([])
        if 'function' in prm:
            plot_function(path=path,measures=measures,axes=axes,values=[val],
                hold=1,**prm)

        else:
            plot_data(path=path,measures=measures,axes=axes,values=[val],
                hold=1,**prm)

def data_xy(path,xvar='delay',yvar='harvest',relative='xy',**kwargs):
    '''Return one of the model.data variables as a function of another,
    aggregating over all the files in "path". '''
    results={}
    keyorder=None
    groupby=kwargs.get('groupby',[])
    aggregate_over=kwargs.get('aggregate_over',[])
    for idx,f,model,measure in extract_data(path,**kwargs):
        if keyorder is None:
            keyorder=tuple(x for x in sorted(measure) if x!='replica')
        key=tuple( (k,measure[k]) for k in keyorder if ( (groupby and k in groupby)
             or not k in aggregate_over) )
        results.setdefault(key,[])
        xs=model.data[xvar][1:]
        ys=model.data[yvar][1:]
        if 'x' in relative:
            xs/=np.mean(xs)
        if 'y' in relative:
            ys/=np.mean(ys)
        results[key].append( (xs.ravel(),ys.ravel())  )
    #Simplify the keys
    diffcols=set([])
    keyref=results.keys()[0]
    for key in results:
        for i,tup in enumerate(key):
            k,m = tup
            if m!= keyref[i][1]:
                diffcols.add(i)
    final={}
    for key in results:

        simkey=tuple(key[c] for c in diffcols)
        print 'data_xy',simkey
        final.setdefault(simkey,[[],[]])
        for xs, ys in results[key]:
            final[simkey][0]+=list(xs)
            final[simkey][1]+=list(ys)
    return final


def data_xy_plot(path,xvar=None,yvar=None,relative=(),**kwargs):
    """Calls data_xy, then scatterplots."""
    bins=kwargs.get('bins',40)
    for key,tup in data_xy(path,xvar=xvar,yvar=yvar,
                relative=relative,**kwargs).iteritems():
        plt.figure()
        if bins:
            from datatools import hist_binner, errorbar
            binlist,res=hist_binner(tup[0],tup[1],bins=bins)
            errorbar(binlist, [np.mean(res.get(b,0)) for b in binlist],
                yerr=[np.std(res.get(b,0)) for b in binlist] ,
                hold=1,title=str(key),xlabel=xvar,ylabel=yvar)
        else:
            scatter(tup[0],tup[1],hold=1,title=str(key),xlabel=xvar,ylabel=yvar)
    plt.show()



def analyze_trajectory(model,hold=0,log='y'):
    if isinstance(model,basestring):
        from models import BaseModel
        m=BaseModel.load(model,initialize=False)
    else:
        m=model

    from datatools import plot,plt
    for var in m.results:
        plt.figure()
        res=m.results[var].matrix
        plot(np.array([res[t].ravel() for t in range(res.shape[0])]),
            hold=1,xs=m.results[var].index,
            log=log,title=var)
    if not hold:
        plt.show()





def split_plot(measures,col,axis,splitby,**kws):
    '''Internal use, for 'show'. '''
    print 'COL',col
    plt.figure()
    kwargs={}
    kwargs.update(kws)

    plot_fnc=kwargs.pop('plot_fnc',plot)
    aggregator=kwargs.pop('aggregator','median')
    if isinstance(aggregator,basestring):
        aggregator=lambda e,a=aggregator: getattr(e,a)()
    kwargs['hold']=True

    if not splitby:
        lst= aggregator(measures.groupby(axis))
        plot_fnc(lst.index,lst[col],**kwargs)
        return lst
    lst=[]
    legends=[]
    ip=0
    grouped=measures.groupby(splitby)
    for name,group in grouped:
        grp= aggregator(group.groupby(axis))
        lst.append( grp[col] )
        legends.append('{}:{}'.format(splitby,name))
        kwargs['hold']=1
        kwargs['alpha']=.5
        kwargs.setdefault('linewidth',3)

        kwargs['color']=get_cmap(ip,len(grouped))
        plot_fnc(grp.index,grp[col],**kwargs)
        ip+=1
    plt.legend(legends)
    return lst


def show(path,plots=(),axes='',plotdim=2,labels=None,use_measures=(),hold=0,
        recompute=None,keep_cols=[],plot_fnc=None,aggregator='median',
        aggregate_over=None,read_filter=None):
    '''plotdim = 2 or 3 (3 requires axes to have at least two items)
    axes: can be a single label or a list. If list, then grouping axes beyond
        plot dimension are used to separate data into different curves/surfaces
    keep_cols: specify which colums from files.csv should be kept into measures.csv
        (when unusual properties must be kept track of)
    plot_fnc: specifiy which function (plot,scatter) to use in split_plot
    aggregator: function to call on data once it has been separated using grouby? by default: take median
    '''
    import pandas as pd
    path = Path(path)

    if not recompute:
        try:
            measures=open_measures(path)
            if recompute is None and modification_date(path+'measures.csv')<modification_date(path+'files.csv'):
                print 'Outdated measures.csv, recomputing...'
                recompute=True
        except IOError as e:
            if recompute is False:
                print 'Cannot find measures.csv, recomputing...'
            recompute=True

        if not read_filter is None:
            key,operator,val =read_filter
            measures=measures[eval('measures[key] {} val'.format(operator))]
            #print measures

    if recompute:
        measures=make_measures(path,use_measures=use_measures,keep_cols=list(keep_cols),
            read_filter=read_filter,aggregate_over=aggregate_over)


    if aggregate_over:
        dico={}
        grouped=measures[m].groupby('group')
        mgroup=grouped.mean()
        sgroup=grouped.std()
        for m in measures:
            try:
                dico[m+'-mean']=mgroup[m]
                dico[m+'-std']=sgroup[m]

            except:
                pass
        measures=pd.DataFrame(dico)

    if labels is None:
        labels=axes[:]
    if not isinstance(axes,basestring):
        axis,splitby=axes[:plotdim-1],axes[plotdim-1:]
    else:
        if plotdim>2:
            raise Exception("Error: Cannot do 3d plot with single axis")
        splitby=None
        axis=axes

    if plotdim==2:
        #2d plot
        xlabel=labels[0]
        for top in plots:
            try:
                kw=plots[top]
            except TypeError as e:
                kw={}
            if not 'title' in kw:
                kw['title']=top
            if plot_fnc:
                kw['plot_fnc']=plot_fnc
            if aggregator:
                kw['aggregator']=aggregator
            split_plot(measures,top,axis,splitby,xlabel=xlabel,**kw)
    else:
        #3d plot
        print '3D PLOT NOT IMPLEMENTED YET'
        xlabel,ylabel=labels[:2]
        plt.figure()


    if not hold:
        plt.show()

