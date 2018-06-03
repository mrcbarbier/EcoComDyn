# -*- coding: utf-8 -*-
from utility import *
from sampler import Distribution
from labeledarray import LabeledArray

class Trajectory(object):

    def __init__(self,init=None,matrix=None,t=0,index=None,label='traj',**kwargs):
        if matrix is None:
            self.matrix=np.array([init])
        else:
            self.matrix=matrix.copy()
        if index is None:
            self.index=np.array([t]) #TIME
        else:
            self.index=index.copy()
        self.label=label
        self.prm=kwargs

    def append(self,array,index=None):
        self.matrix=np.concatenate((self.matrix,array.reshape((1,)+array.shape) ))
        if index is None:
            index=self.index[-1]+1
        self.index= np.concatenate((self.index,[index]))
        return self

    def add(self,other,index=None):
        if index is None:
            try:
                dt=self.index[-1]-self.index[-2]
            except:
                dt=1
            index=np.linspace(self.index[-1]+dt,
                self.index[-1]+other.shape[0]*dt ,other.shape[0])
        self.matrix= np.concatenate((self.matrix,other.matrix),axis=0)
        self.index= np.concatenate((self.index,index))


    def __getitem__(self,*args,**kwargs):
        return self.matrix.__getitem__(*args,**kwargs)

    def __setitem__(self,*args,**kwargs):
        val= self.matrix.__setitem__(*args,**kwargs)
        return val

    @property
    def shape(self):
        return self.matrix.shape

    def reshape(self,shape):
        return self.matrix.reshape( (len(self.index),)+tuple(shape) )

    def set(self,array,index=None):
        if not index is None:
            self.index=index.copy()
        self.matrix=array.copy()
        return self

    def view(self):
        return self.index, self.matrix

    def save(self,path,fname=None):
        if fname is None:
            fname=self.label
        return save_table({self.label:self},path=path+fname ,footer=str(self.prm) )

    @classmethod
    def load(klass,path,fname=None):
        index,tab,footer=load_table(path+fname)
        from ast import literal_eval
        kwargs= literal_eval(footer)
        return klass(matrix=tab.get(fname, tab.values()[0] ),index=index,**kwargs )

    def copy(self):
        t=self.__class__(**deepcopy(self.prm))
        if not self.matrix is None:
            t.set(self.matrix.copy(), self.index.copy())
        return t



    def get(self,t,interpolate=True):
        '''Find value whose index is left of t. If 'interpolate', interpolate linearly between increments.'''
        from bisect import bisect_left
        xbin=max(bisect_left(self.index,t)-1,0)
        if interpolate:
            return self.matrix[xbin] + ( t-self.index[xbin] ) * (self.matrix[xbin+1]-self.matrix[xbin]) / (
                self.index[xbin+1]-self.index[xbin])
        else:
            return self.matrix[xbin]



class Noise(Trajectory):
    '''Pregenerated stochastic time series for use in dynamics.'''

    dft={
        'rank':1,
        'direction':None,
        }

    def __init__(self,label='noise',index=None,matrix=None,**kwargs):
        self.label=label
        self.distribution=Distribution(**kwargs)
        for i,j in self.dft.iteritems():
            kwargs.setdefault(i,j)
        self.prm=kwargs
        self.index=index
        self.matrix=matrix
        self.dt=None

    def copy(self):
        t=Trajectory.copy(self)
        t.dt=self.dt.copy()
        return t

    def generate(self,t1,t2,dt,shape=None ):
        '''Generate uniform noise between t1 and t2 with timestep dt'''
        self.index=np.linspace(t1,t2,(t2-t1)/dt )
        if shape is None:
            shape=self.get_rank()
        shape=to_iterable(shape,tuple)
        self.matrix=np.random.random( self.index.shape + shape )
        self.dt=Trajectory(matrix=dt*np.ones(self.index.shape),index=self.index)
        return self

    def extend(self,t,dt=None):
        '''Generate more noise if t is outside of range.'''
        if self.index[0]<= t <=self.index[-1]:
            return
        if dt is None:
            dt=self.index[-1]-self.index[-2]
        prev=len(self.index)
        handled= self.add(self.__class__().generate(self.index[-1]+dt,t+10*dt,dt,self.matrix.shape[1:]) )
        self.dt.add(dt*np.ones(len(self.index[prev:]) ),index=self.index[prev:] )
        return handled

    def get_rank(self):
        rank=self.prm['rank']
        if rank is 'full':
            return self.prm['shape']
        return to_iterable(rank,tuple)

    def get(self,t):
        self.extend(t)
        res= self.distribution.get(Trajectory.get(self,t))/self.dt.get(t)**0.5 #NB: I have tested that 1/sqrt(dt), it is correct here
        res *=self.prm.get('amplitude',1)
        direction = self.prm['direction']
        shape=self.prm['shape']
        rank=self.get_rank()
        if np.array(rank).squeeze()==1:
            rank=1
        if not direction is None:
            if direction == 'random':
                if not rank is 1:
                    direction=np.zeros(shape+rank )
                else:
                    direction=np.zeros(shape)
                direction[tuple(np.random.randint(0,s) for s in shape) ]=1
            elif isinstance(direction,int):
                d=direction
                if not rank is 1 :
                    direction=np.zeros(shape+rank )
                else:
                    direction=np.zeros(shape)
                direction[d]=1

            if not rank is 1:
                res=np.dot(direction,res)
            else:
                res=direction*res
        return res

class LabeledNoise(Noise):
    axes=None
    def __init__(self,axes=None,*args,**kwargs):
        self.axes=axes
        return Noise.__init__(self,*args,**kwargs)

    def get(self,t):
        return LabeledArray(Noise.get(self,t),axes=self.axes )

    @classmethod
    def from_noise(self,noise,axes=None):
        n= LabeledNoise(label=noise.label,index=noise.index,
            matrix=noise.matrix,axes=axes,**n.prm)
        n.distribution=noise.distribution
        return n
