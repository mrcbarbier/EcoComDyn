# -*- coding: utf-8 -*-
import numpy as np
from labeledarray import LabeledArray

class Module(object):

    def __init__(self,parent=None,*args,**kwargs):
        if parent:
            self.parent=parent
    def __str__(self):
        return "Module"

    def E(self,*args,**kwargs):
        return 0

    def copy(self):
        copy= self.__class__(self.parent)
        return copy

    def get_data(self,label,**kwargs):
        return self.parent.data[kwargs[label]],self.parent.parameters[kwargs[label]]


class LinearModule(Module):
    def __str__(self):
        return "LinearModule"


    def dx(self,diff,xs,*args,**kwargs):
        var=kwargs['variables'][0][0]
        x=xs[var]

        r=kwargs['matrix'].get('diagonal',None)
        A=kwargs['matrix'].get('interactions',None)
        #print r.axes, k.axes, x.axes
        if not r is None:
            diff[var] +=  r*x
        if not A is None:
            diff[var] += A.dot(x,axes=kwargs['variables'])

class LogisticModule(Module):
    def __str__(self):
        return "LogisticModule"


    def dx(self,diff,xs,*args,**kwargs):
        var=kwargs['variables'][0][0]
        x=xs[var]

        r=kwargs['matrix']['growth']
        k=kwargs['matrix']['capacity']

        #print r.axes, k.axes, x.axes
        diff[var] +=  r*x*(1-x/k)

class CommunityModule(Module):
    '''Only between-species interactions'''

    def __str__(self):
        return "CommunityModule"


    def dx(self,diff,xs,*args,**kwargs):
        A=kwargs['matrix']['interactions']
        variables=kwargs['variables']
        if len(variables)==1:
            #Interactions within variable
            var=kwargs['variables'][0][0]
            x=xs[var]
            diff[var]+=x *  A.dot(x,axes=kwargs['variables'])
        else:
            #Interactions between variables
            var1,var2=[z[0] for z in kwargs['variables']]
            diff[var1]+=xs[var1] *  A.dot(xs[var2],axes=kwargs['variables'][-1:] )



class LotkaVolterraModule(Module):
    def __str__(self):
        return "LotkaVolterraModule"

    def dx(self,diff,xs,*args,**kwargs):
        A=kwargs['matrix']['interactions']
        r=kwargs['matrix']['growth']
        k=kwargs['matrix']['capacity']
        var=kwargs['variables'][0][0]
        x=xs[var]
        #print A.matrix, r.matrix, k.matrix
        diff[var]+=(r/k )* x * (k - x +  A.dot(x,axes=kwargs['variables']) )
        #assert not np.isnan(diff[var].matrix).any()

        #print 'LV', np.mean(np.abs(((r/k)* x * (k - x +
             #A.dot(x,axes=kwargs['variables'])) ).matrix ))

    def E(self,xs,**kwargs):
        A=kwargs['matrix']['interactions']
        k=kwargs['matrix']['capacity']
        r=kwargs['matrix']['growth']
        var=kwargs['variables'][0][0]
        x=xs[var]
        return np.sum( x * ( r - .5 * A.dot(x,axes=kwargs['variables'] )))


class SimpleLotkaVolterraModule(Module):
    def __str__(self):
        return "SimpleLotkaVolterraModule"

    def dx(self,diff,xs,*args,**kwargs):
        A=kwargs['matrix']['interactions']
        r=kwargs['matrix']['growth']
        try:
            d=kwargs['matrix']['diagonal']
        except:
            d=0
        var=kwargs['variables'][0][0]
        #print A
        #print r
        #print d
        x=xs[var]
        diff[var]+= x * (r - d* x +  A.dot(x,axes=kwargs['variables']) )
        #if (x.matrix>10**3).any():
            #xwhere=np.where(x.matrix>10**3)[0][0]
            #print x.matrix[xwhere], r.matrix[xwhere],d.matrix[xwhere],A.matrix[xwhere],A.matrix[:,xwhere]
            #print diff[var][xwhere]
            #raise
        if 'influx' in kwargs['matrix']:
            # print 'yay', kwargs['matrix']['influx'], x,diff[var]
            diff[var]+=kwargs['matrix']['influx']


    def absdx(self,diff,xs,*args,**kwargs):
        A=kwargs['matrix']['interactions']
        r=kwargs['matrix']['growth']
        try:
            d=kwargs['matrix']['diagonal'].abs()
        except:
            d=0
        var=kwargs['variables'][0][0]
        x=xs[var]
        diff[var]+= x * (r.abs() + d* x +  A.abs().dot(x,axes=kwargs['variables']) )
        if 'influx' in kwargs['matrix']:
            diff[var]+=np.abs(kwargs['matrix']['influx'])

class SimpleFunctionalResponseModule(Module):
    def __str__(self):
        return "SimpleFunctionalResponseModule"

    def dx(self,diff,xs,*args,**kwargs):
        A=kwargs['matrix']['interactions']
        r=kwargs['matrix']['growth']
        d=kwargs['matrix']['diagonal']
        var=kwargs['variables'][0][0]
        x=xs[var]
        AA= A.dot(x,axes=kwargs['variables'])
        factor=AA/kwargs['matrix'].get('threshold',1)/np.mean(A.matrix)
        AA/=1+factor.abs()
        #print AA, kwargs['matrix'].get('threshold',1)*np.mean(A.matrix)
        #print r-d*x+AA
        #print np.mean(A.matrix),np.sum(A.matrix==0)*1./A.matrix.size
        diff[var]+= x * (r - d* x + AA )
        #if (x.matrix>10**3).any():
            #xwhere=np.where(x.matrix>10**3)[0][0]
            #print x.matrix[xwhere], r.matrix[xwhere],d.matrix[xwhere],A.matrix[xwhere],A.matrix[:,xwhere]
            #print diff[var][xwhere]
            #raise


class SimpleFunctionalTrophicModule(Module):
    def __str__(self):
        return "SimpleFunctionalTrophicModule"

    def dx(self,diff,xs,*args,**kwargs):
        A=kwargs['matrix']['fluxes']
        eps=kwargs['matrix']['efficiency']
        r=kwargs['matrix']['growth']
        d=kwargs['matrix']['diagonal']

        var=kwargs['variables'][0][0]
        try:
            q=kwargs['matrix']['exponent']
        except:
            q=1
        x=xs[var]

        try:
            C=kwargs['matrix']['nontrophic'].dot(x,axes=kwargs['variables'])
        except:
            C=0
        # print C
        xq=x**q
        factor=1+ (A/kwargs['matrix'].get('threshold',1)).dot(xq,axes=kwargs['variables']).abs() #/max(10**-10,np.abs(np.mean(A.matrix)))

        Af=A/factor
        AG= (eps*Af).dot(xq,axes=kwargs['variables'])
        AL=(Af).T.dot(x,axes=kwargs['variables'])

        diff[var]+= x * (r + C - d* x + AG ) -xq* AL

        if 'influx' in kwargs['matrix'] and kwargs['matrix']['influx'].matrix.any():
            diff[var]+=kwargs['matrix']['influx']





class DemographicModule(Module):
    def __str__(self):
        return "DemographicModule"

    def dx(self,diff,xs,*args,**kwargs):
        D,axes=kwargs['noise']['demography']
        var=kwargs['variables'][0][0]
        x=xs[var]
        diff[var]+=  np.sqrt(x) * D.get(kwargs['t'])



class MacArthurModule(Module):
    def __str__(self):
        return "LotkaVolterraModule"

    def dx(self,diff,xs,*args,**kwargs):
        C=kwargs['matrix']['consumption']
        var=kwargs['variables'][0][0]
        res=kwargs['variables'][1][0]
        x=xs[var]
        if res in xs:
            r=xs[res]
        else:
            r=kwargs['data'][res]
        ax1,ax2=kwargs['variables'][:1],kwargs['variables'][-1:]
        m=kwargs['matrix'].get('mortality',None)
        if not m is None:
            diff[var]+= x*(C.dot(r-x.dot(C,axes= ax1),axes=ax2) - m)
        else:
            diff[var]+= x*C.dot(r-x.dot(C,axes= ax1),axes=ax2)

class DispersalModule(Module):

    def __str__(self):
        return "LotkaVolterraModule"


    def dx(self,diff,xs,*args,**kwargs):
        A=kwargs['matrix']['dispersal']
        var=kwargs['variables'][0][0]
        x=xs[var]
        res= A.dot(x,axes=kwargs['variables'])
        diff[var]+=res
        #print  'D',np.mean(np.abs(res.matrix))


class NoiseModule(Module):

    def __str__(self):
        return "NoiseModule"


    def dx(self,diff,xs,*args,**kwargs):
        for role,noise in kwargs['noise'].iteritems():
            var=noise.axes[0]
            # print noise.get(kwargs['t']),xs[kwargs['variables'][0][0] ]
            if 'scaling' in kwargs:

                diff[var[0]]+=noise.get(kwargs['t'])*np.clip(xs[kwargs['variables'][0][0]],0,None)**kwargs['scaling']
            else:
                diff[var[0]]+=noise.get(kwargs['t'])
