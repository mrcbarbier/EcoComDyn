# -*- coding: utf-8 -*-

from imports import *
import os
import datetime
import matplotlib.pyplot as plt
import matplotlib.cm as cmx
import matplotlib.colors as colors
import numpy as np, pandas as pd

from copy import deepcopy

def code_debugger(skip=0):
    import code
    import inspect
    stack=inspect.stack()
    dic = {}
    dic.update(stack[1+skip][0].f_globals)
    dic.update(stack[1+skip][0].f_locals)
    code.interact(local=dic)

class MirrorDict(dict):
    def __missing__(self, key):
        return key


def dumps(fil,obj):
    fil.write(str(obj))

def loads(fil):
    txt=''
    for l in fil:
        txt+=l.strip()
    return eval(txt,{},{'array':np.array,'nan':np.nan})



def save_table(table,path=None,index=None,footer=None):
    '''Save data in table form (labelled columns sharing an index)'''
    if not table:
        return False
    df=None
    for name,data in table.iteritems():
        if data is None:
            print "save_table: No data for",name
            continue
        if hasattr(data,'index'):
            #Trajectory
            index=data.index
            matrix=data.matrix
        else:
            #Normal matrix
            matrix=data
            if index is None:
                index=np.array(range(matrix.shape[0]))

        dic=[]
        for t in range(index.shape[0]):
            if matrix[t].shape:
                dic.append(
                {(name,i):j for i,j in mat_to_dict(matrix[t]).iteritems() })
            else:
                dic.append({(name,(0,)):matrix[t] })

        if df is None:
            df=pd.DataFrame(dic,index=index)
        else:
            df=df.join(pd.DataFrame(dic,index))
    if df is None:
        return False
    df.to_csv(path,sep='\t')
    if not footer is None:
        fil=open(path,'a')
        fil.write('#{}'.format(footer))
        fil.close()
    return True

def load_table(path=None,index_col=0):
    '''Load data in table form (labelled columns sharing an index)'''
    from StringIO import StringIO
    fil=open(path,'r')
    comments,uncommented='',''
    for l in fil:
        if '#' in l:
            l,sep,comment=l.partition('#')
            comments+=comment
            if l:
                uncommented+=l+'\n'
        elif l.strip():
            uncommented+=l


    df=pd.read_csv(StringIO(uncommented),sep='\t',header=0,index_col=index_col)
    index=np.array(df.index)
    matrix={}
    for idx,r in df.iterrows():
        dr={}
        for i,j in r.iteritems():
            name,col=eval(i)
            dr.setdefault(name,{})
            dr[name][col]=j
        for name in dr:
            matrix.setdefault(name,[])
            matrix[name].append(dict_to_mat(dr[name]) )
    for name,mat in matrix.iteritems():
        mat=np.array(mat)
        mat[np.isnan(mat)]=0
        matrix[name]=mat

    if comments:
        return index,matrix,comments
    return index,matrix

def to_iterable(iterable,klass=list):
    if not hasattr(iterable,'__iter__'):
        return klass( (iterable,) )
    return klass(iterable)

def put_in_bins(xs,ys=None,nbins=10,binpos='default',bintype='percentile',**kwargs ):
    #use xs as coordinates and ys as things to put in bins
    xs=np.array(xs)
    if ys is None:
        ys=range(xs.shape[0])
    ys=np.array(ys)
    bins=[]
    xbins=None
    if bintype=='percentile':
        xbins=np.unique(np.array([np.percentile(xs,p) for p in np.linspace(0,100,nbins+1)]))
        nbins=len(xbins)-1
    elif bintype=='lin':
        xmin,xmax=np.min(xs),np.max(xs)
        xbins= np.linspace(xmin,xmax,nbins+1)
    elif bintype=='log':
        xmin,xmax=np.min(xs[xs>0]),np.max(xs)
        if xmax<0:
            raise Exception('put_in_bins: xs<0, log impossible')
        xbins=np.logspace(np.log10(xmin),np.log10(xmax),nbins+1)

    xbins[-1]+= 0.001*np.abs(xbins[-1])
    nonemptybins=[]
    for n in range(nbins):
        b =ys[np.logical_and(xbins[n]<=xs,xs<xbins[n+1])]
        if len(b)==0:
            continue
        bins.append(b)
        nonemptybins.append(n)
    nonemptybins.append(-1)
    xbins=xbins[nonemptybins]
    if binpos =='default':
        #Center of bin except edges for first and last bin
        xs=(xbins[1:]+xbins[:-1])/2.
        xs[0]=xbins[0]
        xs[-1]=xbins[-1]
    elif binpos =='center':
        xs=(xbins[1:]+xbins[:-1])/2.
    elif binpos =='left':
        xs=xbins[:-1]
    elif binpos =='right':
        xs=xbins[1:]
    return xs,bins

def binpercentile(xs,ys,nbins=10,percs=(20,50,80),binpos='default',bintype='percentile',**kwargs ):
    #Returns two things: bin coordinates (center of bin, by default), and array of shape ( len(percs), len(data) ) with percentiles percs as rows

    data=np.array(xs)
    xs,oldbins=put_in_bins(xs,ys,nbins=nbins,percs=percs,binpos=binpos,bintype=bintype,**kwargs )
    bins=[]
    for b in oldbins:
        bins.append([ np.percentile(b,p) for p in percs ]+ [float(len(b))] )
    bins= np.array(bins).T
    if len(bins)==0:
        print 'ERROR!',xs,xbins
        raise
    return xs, bins

def binfill(xs,ys,alpha_density=0,**kwargs):
    kw={k:kwargs.pop(k) for k in ('nbins','percs','binpos','bintype') if k in kwargs}
    xbins, bins=binpercentile(xs,ys,**kw)
    print xbins.shape,bins.shape
    plt.plot(xbins,bins[1],**kwargs)
    if alpha_density:
        alpha=kwargs.pop('alpha',.6)* ( bins[-1]/np.max(bins[-1]) )**.5
        for i in range(len(xbins)-1):
            print alpha[i]
            kwargs['linewidth']=0
            plt.fill_between(xbins[i:i+2], bins[0,i:i+2], bins[2,i:i+2],alpha=alpha[i], **kwargs)
    else:
        plt.fill_between(xbins, bins[0], bins[2], **kwargs)

    if not kwargs.get('hold',0):
        plt.show()


def kth_diag_indices(a, k):
    rows, cols = np.diag_indices_from(a)
    if k < 0:
        return rows[-k:], cols[:k]
    elif k > 0:
        return rows[:-k], cols[k:]
    else:
        return rows, cols


def tupleify(measure):
    def totuple(a):
        if not hasattr(a,'__iter__'):
            return a
        try:
            return tuple(totuple(i) for i in a)
        except TypeError:
            return a
    if not hasattr(measure,'keys'):
        if hasattr(measure,'__iter__'):
            return tuple([tupleify(m) for m in measure])
        return measure
    for i,j in measure.items():
        if isinstance(j, np.ndarray):
            measure[i]='array({})'.format(totuple(j))
        else:
            measure[i]=totuple(j)
    return measure

#def tupleify(df,recursive=1):
    #'''Arrays are transformed into strings'''
    #pandas_mode= isinstance(df,pd.DataFrame    )
    #if pandas_mode:
        #df=df.copy()
    #for z in df:
        ##print z, df[z][0],isinstance(df[z][0],basestring) and '(' in df[z][0]
        #if pandas_mode:
            #val=df[z][0]
        #else:
            #val=df[z]
        #if isinstance(val,np.ndarray):
            #try:
                #if pandas_mode:
                    #df[z]=df[z].map(lambda x: to_iterable(x,klass=tuple) )
                #else:
                    #df[z]=ast.safeval(df[z])
                ##print "NOW",type(df[z][0])
            #except Exception as e:
                #print 'tupleify:',e
        #if recursive:
            #if hasattr(val,'keys') or isinstance(val,pd.DataFrame):
                #df[z]=tupleify(df[z])
    #return df


def retupleify(df,recursive=1,verbose=0):
    '''Tuples are recovered as strings by read_csv, convert them to tuples'''
    import ast
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
                    v=v.replace('array','')
                    ls= eval(v,{},dic)
                    if hasattr(ls,'__iter__'):
                        ls= [safeval(l) for l in ls ]
                    else:
                        ls=[ls]
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

def read_csv(*args,**kwargs):
    # try:
    df=pd.read_csv(*args,**kwargs)
    # except Exception as e:
        # df=pd.read_csv(*args)
        # print 'error in read_csv',e
        # code_debugger()
    df =df[ [k for k in df if not 'Unnamed' in k] ]
    retupleify(df)
    return df

def open_measures(path,fname='measures',**kwargs):
    if not '.' in fname:
        fname=fname+'.csv'


    return read_csv(Path(path)+fname,**kwargs )

def get_cmap(i,N):
    import matplotlib.cm
    return matplotlib.cm.get_cmap('jet')( (i*1./N ) )

def old_get_cmap(N):
    '''Returns a function that maps each index in 0, 1, ... N-1 to a distinct
    RGB color.'''
    color_norm  = colors.Normalize(vmin=0, vmax=N-1)
    scalar_map = cmx.ScalarMappable(norm=color_norm, cmap='hsv')
    def map_index_to_rgb_color(index):
        return scalar_map.to_rgba(index)
    return map_index_to_rgb_color

def modification_date(filename):
    t = os.path.getmtime(filename)
    return datetime.datetime.fromtimestamp(t)

def mat_to_dict(mat):
    if not mat.shape:
        return {(0,):mat }
    idx=np.indices(mat.shape)
    idx=np.swapaxes(np.array([idx[x].ravel() for x in xrange(idx.shape[0]) ] ),0,-1)
    return {tuple(i):j for i,j in zip(idx,mat.ravel()) }


def mat_to_cols(mat,idx=None,drop_zeros=0):
    index=np.indices(mat.shape)
    if not idx is None:
        mat=mat[np.ix_(idx,idx)]
        index=np.array([index[x][np.ix_(idx,idx)] for x in xrange(index.shape[0])])
    index=np.swapaxes(np.array([index[x].ravel() for x in xrange(index.shape[0]) ] ),0,-1)

    mat=mat.ravel()
    if drop_zeros:
        xs=(mat!=0)
        index=index[xs]
        mat=mat[xs]
    return [tuple(i) for i in index ],mat

def dict_to_mat(dic,shape=None,keytype='int'):
    #Works if dictionary has tuples of numbers as keys
    keys = [to_iterable(k,klass=tuple) for k in dic.keys()]
    if keytype=='any':
        minlen=min([len(k) for k in keys])
        skeys=[sorted(set([k[j] for k in keys]))   for j in range(minlen)]
        keydic={i:tuple([skeys[j].index(i[j] ) for j in range(minlen)]) for i in dic.keys() }
        keys=[np.argsort(s) for s in skeys ]
    else:
        keydic={i:i for i in dic.keys()}

    if shape is None:
        shape=tuple(np.array(np.max(keys,axis=0)+1))
        elsh=np.max([np.array(d).shape  for d in dic.values( )],axis=0)

    mat=np.zeros(shape)
    for i,j in dic.iteritems():
        mat[keydic[i]]=j

    if mat.shape==(1,):
        mat=mat.squeeze()
    return mat

def anydict_to_mat(matrix,shape=None):
    array=np.array
    #Works for any keys
    if hasattr(matrix,"__iter__") and not hasattr(matrix,"keys"):
        return array(matrix)
    if hasattr(matrix.keys()[0],'__iter__'):
        keys=sorted(matrix)
        if hasattr(matrix.values()[0],'__iter__'):
            elsh=array(matrix.values()[0]).shape
        else:
            elsh=()
        sh=(sorted(matrix)[-1][0]+1,sorted(matrix,key=lambda e:e[1])[-1][1]+1)+elsh
        sh=tuple(array(sh))

        if not shape is None:
            sh=shape
        mat=np.zeros(sh)
        for i,j in sorted(matrix):
            mat[i,j]=matrix[(i,j)]
        return mat
    return np.array([[matrix[i][j] for j in sorted(matrix[i])] for i in sorted(matrix)])

def nozero(mat):
    return mat[mat!=0]
def nonan(mat):
    return mat[np.isnan(mat)==False]


def auto_subplot(plt,nbpanels,panel,projection=None,rows=None,return_all=0):
    i=panel.val
    if rows is None:
        panels_per_row=np.ceil(np.sqrt(nbpanels) )
    else:
        panels_per_row=np.ceil(nbpanels/rows).astype('int')
    nbrows=np.ceil(nbpanels*1./panels_per_row)
    ax=plt.subplot(nbrows,panels_per_row,i +1 ,projection=projection)
    panel.val+=1
    if return_all:
        return ax,nbrows,panels_per_row
    return ax

class MutableInt(object):
    def __init__(self,val):
        self.val = val
    def __add__(self,i):
        return self.val+i
    def __eq__(self,i):
        return self.val==i


import networkx as nx

class Graph(nx.MultiDiGraph):
    pos=None
    _label=None
    def __str__(self):
        if self._label:
            return self._label
        return nx.MultiDiGraph.__repr__(self)
    def levels_layout(self,ranks=None,xpos=None):
        pos={}
        refpos=nx.spring_layout(self)
        if not isinstance(xpos,dict):
            xpos={i:j for i,j in zip(self.nodes(),xpos ) }
        maxx,maxy=np.max(refpos.values(),axis=0)
        for i,d in self.nodes(data=True):
            if not xpos is None and i in xpos:
                x=xpos[i]
            else:
                x= np.random.random()*maxx
            if ranks is None:
                y=d['height']*maxy
            else:
                y=ranks[i]
            pos[i]=x,y
        self.pos=pos

    def plot(self,newfig=True,hold=False,labels=None,edge_labels=None,nscale=1,minsize=0.001,**kwargs):
        '''Use matplotlib to plot the network'''
        import matplotlib.pyplot as plt
        if self.pos is None:
            pos=self.pos=nx.spring_layout(self)
        else:
            pos=self.pos
        if newfig:
            plt.figure()
        node_size=np.ones(len(self.node) ) *kwargs.get('node_size',1)
        node_size[node_size<minsize]=0
        node_size/=np.max(node_size)/2.
        node_size[node_size>0]=np.clip(node_size[node_size>0],0.1,None)*nscale


        node_size[np.random.random(node_size.shape )<kwargs.get('subsample',0)  ]=0
        nodidx={n:i for i,n in enumerate(self.nodes()) }

        edge_width=kwargs.get('edge_width',1)
        if edge_width:
            sign=kwargs.get('sign',0)
            if not sign or sign<0:
                nx.draw_networkx_edges(self,pos=pos,edgelist=[(i,j) for i,j,d in
                    self.edges(data=True) if d.get('weight',0)<0 and node_size[nodidx[i]]*node_size[nodidx[j]]>0],
                    edge_color='r',width=edge_width )
            if not sign or sign>0:
                nx.draw_networkx_edges(self,pos=pos,edgelist=[(i,j) for i,j,d in
                    self.edges(data=True) if d.get('weight',0)>=0 and node_size[nodidx[i]]*node_size[nodidx[j]]>0.],
                    edge_color='b',width=edge_width )
        nx.draw_networkx_nodes(self,pos=pos ,node_size=node_size*200.,linewidths=0,
                node_color= ifelse('node_color' in kwargs, np.array(kwargs.get('node_color',[])),'blue') )
        if labels is None:
            labels={n:str(n) for n in self.nodes() }
        elif not isinstance(labels,dict):
            labels={n:labels[i] for i,n in enumerate(self.nodes()) }
        nx.draw_networkx_labels(self,pos=pos,labels={n:l for n,l in labels.items() if node_size[nodidx[n]]>0} )
        if not edge_labels is None:
            if edge_labels=='weights':
                #print [self.edge[j][i][0] for i,j,d in
                    #self.edges(data=True)]
                edge_labels={(i,j):'{:.2f},{:.2f}'.format(self.edge[i][j][0]['weight'],self.edge[j][i][0]['weight']) for i,j,d in
                    self.edges(data=True) if # self.edge[i][j][0]['weight']>-self.edge[j][i][0]['weight'] and
                      node_size[nodidx[i]]*node_size[nodidx[j]]>0}
            nx.draw_networkx_edge_labels(self,pos=pos,edge_labels= edge_labels )
        if 'title' in kwargs:
            plt.title(kwargs['title'])
        if 'xlabel' in kwargs:
            plt.xlabel(kwargs['xlabel'])
        if 'ylabel' in kwargs:
            plt.ylabel(kwargs['ylabel'])
        if not hold:
            plt.show()



if 0:
    from collections import OrderedDict, Mapping, Container
    from pprint import pprint
    from sys import getsizeof


    def deep_getsizeof(o, ids):
        """Find the memory footprint of a Python object

        This is a recursive function that rills down a Python object graph
        like a dictionary holding nested ditionaries with lists of lists
        and tuples and sets.

        The sys.getsizeof function does a shallow size of only. It counts each
        object inside a container as pointer only regardless of how big it
        really is.

        :param o: the object
        :param ids: should be set()
        :return:
        """
        d = deep_getsizeof
        if id(o) in ids:
            return 0

        r = getsizeof(o)
        ids.add(id(o))
        if not isinstance(r, int):
            print o
            raise

        if isinstance(o, str) or isinstance(o, unicode) or isinstance(o, np.ndarray):
            return r

        if isinstance(o, Mapping):
            return r + sum(d(k, ids) + d(v, ids) for k, v in o.iteritems())

        if isinstance(o, Container):
            try:
                return r + sum(d(x, ids) for x in o)
            except Exception as e:
                # print e
                pass

        return r
