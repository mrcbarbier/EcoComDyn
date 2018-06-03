# -*- coding: utf-8 -*-
#Generates variables from a given distribution
from utility import *
import scipy.stats
from itertools import product as iproduct





class Distribution(object):
    dft={
        'mean':0,
        'var':1,
        'std':1,
        'resolution':100,
        'distribution':'normal',
        'freq':None, #
        'vanish':None, #

        }

    def __init__(self,**attrs):
        for i,j in self.dft.iteritems():
            attrs.setdefault(i,j)
        self.attrs=deepcopy(attrs)
        mean=attrs['mean']
        std=attrs.get('std',np.sqrt(attrs['var']))
        distribution=attrs['distribution']
        resolution=attrs['resolution']

        self.dist=dist=None
        if distribution =='normal':
            dist=scipy.stats.norm(loc=mean,scale=std)
        if distribution=='uniform':
            dist=scipy.stats.uniform(loc=mean-std,scale=2*std)
        if distribution=='gamma':
            dist=scipy.stats.gamma(a=skew,loc=weight,scale=width)


        if dist:
            self.dist=np.vectorize(dist.ppf)(np.linspace(0.001,.999,resolution) )
            self.ppf=np.vectorize(self.get_ppf)
        else:
            self.ppf = self.sample

    def sample(self,val,t=0):
        """Convert uniformly distributed val to variable distributed according
        to own distribution, when not defined with ppf."""

        if not self.attr.get('freq',None) is None:
            val=np.cos((self.attrs['freq']*t*+val)*2*np.pi )

        if self.attrs['distri']=='dirac':
            res= self.attrs['mean']
        elif self.attrs['distri']=='diracs':
            #print 'diracs'
            if val>(self.attr.get('skew',0)+1.)/2.:
                res= self.attrs['mean']*(1+self.attrs['std'])
            else:
                res=  self.attrs['mean']*(1-self.attrs['std'])
            #print self.attrs['weight']*(1+self.attrs['width']/2.),

        if  not self.attr.get('vanish',None) is None:
            res=res*np.exp(-self.attrs['freq']*t)
        return res


    def get_ppf(self,val):
        return self.dist[int(val*self.dist.shape[0]) ]

    def get(self,val,**kwargs):
        if self.dist is None:
            res= val
        else:
            res=self.ppf(val,**kwargs)
        sign=self.attrs.get('sign',None)
        if not sign is None:
            res=np.abs(res)
            if sign<0:
                res=-res
        return res



class BaseSampler(object):
    def __init__(self,*args,**kwargs):
        self.attrs=kwargs

    def sample(self,**attrs):
        for i in self.attrs:
            attrs.setdefault(i,self.attrs[i])
        if 'use' in attrs:
            return attrs['use']
        distribution=attrs.get('distribution','normal')
        shape=tuple(int(x) for x in to_iterable(attrs.get('shape',100),tuple) )
        for att in attrs.get('scaling',[]):
            attrs[att]*=shape[0]**attrs['scaling'][att]
            #print att,attrs[att]

        mean=attrs.get('mean',0.)
        if 'stdrel' in attrs:
            std=attrs['stdrel']*np.abs(mean)
        else:
            std=attrs.get('std',np.sqrt(attrs.get('var',1.)))
        if distribution== 'normal':
            if np.min(np.atleast_1d(std))<=0 or np.isnan(std).any():
                mat= np.ones(shape)*mean
            else:
                #if hasattr(std,'__iter__'):
                #mat=scipy.stats.multivariate_normal.rvs(size=shape)
                mat= np.random.normal(mean,std,shape)
        if distribution== 'exponential':
            #Mean simply gives a translation
            if std==0:
                mat= np.ones(shape)*mean
            else:
                #if hasattr(std,'__iter__'):
                #mat=scipy.stats.multivariate_normal.rvs(size=shape)
                mat= std*np.random.exponential(1,shape)-std+mean
        if distribution=='uniform' :
            if 'min' in attrs or 'max' in attrs:
                mat=np.random.uniform(attrs.get('min',0),attrs.get('max',1),shape)
            else:
                mat=np.random.uniform(mean-std*np.sqrt(3),mean+std*np.sqrt(3),shape)
        if distribution=='randint' :
            if 'min' in attrs or 'max' in attrs:
                mn,mx=round(attrs.get('min',0)),round(attrs.get('max',1)+1)
            else:
                mn,mx=round(mean-std*np.sqrt(3)),round(mean+std*np.sqrt(3)+1)
            mn,mx=int(mn),int(mx)
            weights=attrs.get('weights',None)
            if 'weightfactor' in attrs:
                weights=attrs['weightfactor']**np.arange(mn,mx)

            if weights is None:
                mat=np.random.randint(mn,mx,shape)
            else:
                weights=np.array(weights) / np.sum(weights)
                mat=np.random.choice(range(mn,mx),p=weights ,size=shape)
            absent=list(range(mn,mx))
            while set(absent)!=set(mat.ravel()):
                if len(absent) > len(mat.ravel())+1:
                    break
                a=np.random.choice(absent)
                if not a in mat:
                    mat[np.unravel_index( np.random.randint(0,mat.size),mat.shape) ]=a

        if distribution=='normals':
            mat=np.zeros(shape)
            if 'values' in attrs:
                mean=attrs['values']
            mean,std=[np.array(z,ndmin=1) for z in (mean,std)]
            if mean.shape!=std.shape:
                mean,std=zip(*list(iproduct(mean,std)))
            i=0
            pos=np.random.randint(0,len(mean),shape[0])
            for m,s in zip(mean,std):
                mat[pos==i]= np.random.normal(m,s,shape)[pos==i]
                i+=1
        if distribution=='diracs':
            if not 'values' in attrs:
                tot=(mean**2+std**2)
                attrs['values']=[0,tot/mean]
                attrs['weights']=[std**2/tot,mean**2/tot ]
            mat=np.random.choice(attrs['values'],size=shape,p=attrs.get('weights',None) )
            if len(attrs['values'])==shape[0] and len(shape)==1:
                tmp=np.array(attrs['values'])
                np.random.shuffle(tmp)
                mat=tmp
        if distribution=='beta':
            A, sigA = np.abs(mean), std
            a = -(A * ((-1 + A) * A + sigA ** 2)) / sigA ** 2
            b = -1 + A + ((-1 + A) ** 2 * A) / sigA ** 2
            mat = np.sign(mean)*np.random.beta(a, b, shape)
        if distribution=='powlaw':
            exp=-attrs.get('exp',2)
            from datatools import powlaw
            mat= mean*powlaw(exp,attrs.get('min',0.01),attrs.get('max',1), shape )

        diag=attrs.get('diagonal',None)
        if not diag is None:
            if hasattr(diag,'keys'):
                mat[np.eye(mat.shape[0]).astype('bool')] = self.sample(shape=(mat.shape[0],),**diag)
            else:
                if hasattr(diag,'__iter__'):
                    mat[np.eye(mat.shape[0]).astype('bool')]=np.diag(diag)
                else:
                    np.fill_diagonal(mat,diag)
        sign=attrs.get('sign',None)
        if not sign is None:
            if sign>0:
                mat=np.abs(mat)
            if sign<0:
                mat=-np.abs(mat)
        return mat
