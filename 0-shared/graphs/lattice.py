# -*- coding: utf-8 -*-

from utility import *
import itertools


def neighbors(arr,x,y,n=3):
    #From unutbu, stackoverflow
    ''' Given a 2D-array, returns an nxn array whose "center" element is arr[x,y]'''
    arr=np.roll(np.roll(arr,shift=-x+1,axis=0),shift=-y+1,axis=1)
    return arr[:n,:n]


class NeighborGraph(object):
    '''Generator for species'''
    _name='species'

    def __init__(self,neighbors=None,**kwargs):
        if neighbors is None:
            neighbors={}
        self.neighbors=neighbors
        self.props=kwargs
        self.name=kwargs.get('name',self._name)

    def __getitem__(self,*args,**kwargs):
        return self.neighbors.__getitem__(*args,**kwargs)

    def __repr__(self):
        return "NeighborGraph(neighbors={},**{})".format(self.neighbors,self.props)

    def save(self,path,fname=None):
        if fname is None:
            fname=self.name
        dumps(open(path+fname,'w'),self.props)

    @classmethod
    def load(klass,path,fname=None):
        if fname is None:
            fname=self.name
        return klass.make(**loads(open(path+fname,'r')))

    def copy(self):
        return self.__class__(deepcopy(self.neighbors),**deepcopy(self.props))

    @classmethod
    def make(klass,**attrs):
        """Generate neighbor list"""
        #struct=attrs['structure']
        shape=attrs['shape']
        bc=attrs['boundary']
        neigh={}
        k=0
        for coords in itertools.product(*[range(i) for i in shape] ) :
            nei=[]
            coords=list(coords)
            for d in range(len(shape)):
                new=coords[:]
                #print 'yes',coords,d
                if coords[d]>0:
                    new[d]=new[d]-1
                elif bc=='periodic':
                    new[d]=shape[d]-1
                else:
                    new =None
                if new:
                    nei.append( tuple(new) )

                new=coords[:]
                if coords[d]<shape[d]-1:
                    new[d]=new[d]+1
                elif bc=='periodic':
                    new[d]=0
                else:
                    new =None
                if new:
                    nei.append( tuple(new) )

            neigh[tuple(coords)]=nei
            #print coords,nei
            k+=len(nei)
        attrs['k']=k/float(np.prod(shape))
        obj=klass(neighbors=neigh,**attrs)
        return obj
