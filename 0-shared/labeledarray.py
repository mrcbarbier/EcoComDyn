# -*- coding: utf-8 -*-
import numpy as np

from utility import save_table,load_table
import types

class LabeledArray(object):
    '''Array where axes are given labels, to allow automatic reshaping
    and alignment of axes and Einstein summation.'''

    def __init__(self,matrix=None,axes=None,**kwargs):
        if hasattr(matrix,'shape') and not matrix.shape:
            matrix=matrix.reshape((1,))
        self.matrix=matrix
        self.axes=axes

    def __getitem__(self,*args,**kwargs):
        return self.matrix.__getitem__(*args,**kwargs)

    def __setitem__(self,*args,**kwargs):
        val= self.matrix.__setitem__(*args,**kwargs)
        return val

    def __repr__(self):
        return 'LabeledArray(matrix={},axes={})'.format(self.matrix,self.axes)

    def __str__(self):
        return '{}:{}'.format(self.axes,self.matrix)

    @property
    def shape(self):
        return self.matrix.shape

    def reshape(self,*args):
        return self.matrix.reshape(*args)

    def ravel(self):
        return self.matrix.ravel()



    def save(self,path,fname=''):
        return save_table({'values':self},path=path+fname ,footer=str(self.axes) )

    @classmethod
    def load(klass,path,fname=''):
        index,tab,footer=load_table(path+fname)
        from ast import literal_eval
        axes= literal_eval(footer)
        #print tab['values'].shape,axes
        return klass(matrix=tab['values'] ,axes=axes )

    def copy(self):
        return self.__class__(matrix=self.matrix.copy(),axes=tuple(self.axes) )

    @property
    def T(self):
        return self.__class__(matrix=self.matrix.T,axes=tuple(self.axes) )

    def dot(self,array, axes=None,**kwargs):
        '''Tensor dot product, summing over products for axes in 'axes'
        (or by default each axis that occurs in both self and array)
        and multiplying together axes with same name that were not summed over.'''
        saxes=self.axes
        oaxes=array.axes
        if axes is None:
            axes=[]
            for ax in saxes[::-1]:
                 if ax in oaxes and not ax in axes:
                     axes.append(ax)

        letters='abcdefghijklmnopqrstuvwxyz'
        slet,letters=letters[:len(saxes)],letters[len(saxes):]
        olet=letters[:len(oaxes)]

        rlet=''
        raxes=[]
        tmp=axes[:]
        for i, ax in enumerate(saxes[::-1]):
            sindex=len(saxes)-1-i
            if ax in oaxes:
                oindex=oaxes.index(ax)
                if not olet[oindex] in slet:
                    olet=olet.replace(olet[oindex],slet[sindex])
            if ax in tmp:
                tmp.remove(ax)
            else:
                rlet=slet[sindex]+rlet
                raxes.insert(0,saxes[sindex])

        for i,let in enumerate(olet):
            if not let in slet:
                rlet=rlet+let
                raxes.append(oaxes[i])

        form='{},{}->{}'.format(slet,olet,rlet)
        #print form
        result=np.einsum(form,self.matrix,array.matrix )
        return LabeledArray(result,axes=raxes )

    def abs(self):
        return LabeledArray(np.abs(self.matrix),self.axes)

    def sum(self,axes=None):
        iaxes=tuple(self.axes.index(ax) for ax in axes)
        #print iaxes
        return LabeledArray(np.sum(self.matrix,axis=iaxes ),
            [ax for i,ax in enumerate(self.axes) if not i in axes])

    def operation(self,val,method):
        '''Applies numpy array operation, aligning axes if val is another
        LabeledArray.'''
        if not hasattr(val,'axes'):
            return LabeledArray(getattr(self.matrix,method)(val),self.axes)
        if val.axes==self.axes:
            return LabeledArray(getattr(self.matrix,method)(val.matrix),self.axes)

        smatrix=self.matrix.copy()
        omatrix=val.matrix.copy()
        #print smatrix,omatrix,self.axes,val.axes

        def unicize(axes):
            new=[]
            for ax in axes:
                if ax in new:
                    if hasattr(ax,'__iter__'):
                        ax=tuple(ax)
                    else:
                        ax=(ax,)
                    ax=ax+('$2',)
                i=3
                while ax in new:
                    ax=ax[:-1]+('${}'.format(i),)
                    i+=1
                new.append(ax)
            return new
        def deunicize(axes):
            return [ax[:-1] if isinstance(ax,tuple) and '$' in ax[-1]
                else ax for ax in axes]

        saxes,oaxes=unicize(self.axes),unicize(val.axes)

        letters='abcdefghijklmnopqrstuvwxyz'
        slet,letters=letters[:len(saxes)],letters[len(saxes):]
        olet=letters[:len(oaxes)]

        oadd=''
        for sindex, ax in enumerate(saxes):
            if ax in oaxes:
                oindex=oaxes.index(ax)
                olet=olet.replace(olet[oindex],slet[sindex])
            else:
                oadd+=slet[sindex]

        sadd=''
        addaxes=[]
        for oindex,ax in enumerate(oaxes):
            if not ax in saxes:
                sadd+=olet[oindex]
                addaxes.append(ax)

        raxes=deunicize(saxes+addaxes)


        smatrix=smatrix.reshape(smatrix.shape
            +tuple(1 for l in sadd) )
        #print form, olet,oadd,slet,sadd

        if olet+oadd == slet+sadd:
            omatrix=omatrix.reshape(omatrix.shape
                +tuple(1 for l in oadd) )
        else:
            form='{}->{}'.format(olet+oadd,slet+sadd)
            omatrix=np.einsum(form,omatrix.reshape(omatrix.shape
                +tuple(1 for l in oadd) ) )

        return LabeledArray(getattr(smatrix,method)(omatrix),raxes)


    def unary(self,method):
        return LabeledArray(getattr(self.matrix,method)(),self.axes)

def la_make_operation(method):
    def operation(self,val):
        return self.operation(val,method)
    return operation

def la_make_assignation(method):
    def assignation(self,val):
        self.matrix= self.operation(val,method).matrix
    return assignation

def la_make_unary(method):
    def unary(self):
        return self.unary(method)
    return unary

for meth in ('add','sub','mul','div','truediv','floordiv','pow'):
    for prefix in ('','r'):
        method='__{}{}__'.format(prefix,meth)
        setattr(LabeledArray,method,la_make_operation(method))
    method='__i{}__'.format(prefix,meth)
    setattr(LabeledArray,method,la_make_assignation('__{}__'.format(meth)))

for meth in ('pos','neg'):
    meth='__{}__'.format(meth)
    setattr(LabeledArray,meth,la_make_unary(meth))

if __name__ =='__main__':
    #Tests
    A=LabeledArray(np.random.random((4,4,5)),('x','x','z') )
    v=LabeledArray(np.random.random((5,4,3)),('z','x','y'))
    #print v.matrix,(v+1).matrix,(2/v).matrix
    g=A.dot(v,axes=['x'])


    h=np.tensordot(A.matrix,v.matrix,axes=[[1],[1]])
    print 'TEST_LA', g[1,3],h[1,3,3]

    i=A.dot(v,axes=[])
    print i.shape, i.axes


    A=LabeledArray(np.random.random((2,3)),('x','y') )
    v=LabeledArray(np.random.random((1,2)),('y','x'))

    print 'A',A.matrix
    print 'v',v.matrix
    print 'A+v',(A+v).matrix
    #print A[0,:],v[0,:],(A+v)[0,0]