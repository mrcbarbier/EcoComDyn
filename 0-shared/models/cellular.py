# -*- coding: utf-8 -*-

from imports import *
import itertools
from model import BaseModel,BaseNetwork
from graphs import NeighborGraph

class CellularModel(BaseModel):

    dft_prm={
            'n':{
        'axes':['spec'],
        'type':'variable',
        'shape':(20,),
        },
            'domain':{
        'shape':(100,100),
        'boundary':'periodic',
        'structure':'2d',
        'type':'lattice',
        'save':False,
        },
            'growth':{
        'type':'matrix',
        'variables':'n',
        'mean':1., 'std':0,
        'sign':1,
        'role':'growth',
        },
            'capacity':{
        'type':'matrix',
        'variables':'n',
        'mean':1., 'std':0,
        'sign':1,
        'role':'capacity',
        },
            'community':{
        'type':'matrix',
        'variables':['n','n'],
        'mean':5., 'std':.3,
        'symmetry':1.,
        'structure':'bunin',
        'role':'interactions',
        },

            'interactions':{
        'type':'matrix',
        'variables':['n','n'],
        'role':'interactions',
        'formula':'community'
                }
    }

    generators={
        'matrix':BaseNetwork,
        'lattice':NeighborGraph,
        }



    def evol(self,tmax=500,tsample=.1,dt=.1,keep=['end','all'][1],
            print_msg=0,adaptive_step=1,
            **kwargs):
        prm=self.set_params(**kwargs)
        comm=self.data['interactions'].matrix
        growth=self.data['growth'].matrix
        capa=self.data['capacity'].matrix
        lattice=self.data['domain']
        selfint=growth/capa
        nspec=prm['n']['shape'][0]


        #initial vector
        x0=np.random.randint(-1,prm['n']['shape'][0], prm['domain']['shape'] )
        frac=kwargs.get('filling_fraction',None)
        if frac!=None:
            #Initial filling_fraction
            species=set(range(nspec) )
            x0[np.unravel_index(np.random.randint(0,x0.size, int((1-frac)*x0.size) ),x0.shape) ]=-1
            #Check if species are missing
            setx0=set(list(x0[x0>=0]) )
            while setx0<species:
                for missing in species-setx0:
                    x0[np.unravel_index(np.random.randint(0,x0.size,1),x0.shape)]=missing
                setx0=set(list(x0[x0>=0]))
        x=x0.copy()

        allcoords= np.array(list(itertools.product(*[range(i) for i in x0.shape] )),
            dtype='int').reshape(x.shape+(2,) )

        t=0
        basic_capacity=x0.size#/nspec
        nnei=lattice.props['k']
        while t<tmax:
            print t
            #if adaptive_step and t>0 and dt==0.1:
                #dt=0.1/max( np.max(comm*abund),np.max(selfint*abund),np.max(growth*abund))*np.max(abund)


            if   (keep=='all' and t%tsample<dt ) or (keep=='endpoint' and t+dt >=tmax):
                measure=True
                abund=np.zeros(nspec,dtype='float')
                contact=np.zeros(prm['n']['shape']+prm['n']['shape'],dtype='float')
            else:
                measure=False

            filled=allcoords[x>=0]
            np.random.shuffle(filled)
            for coords in filled :
                coords=tuple(coords)
                speci=x[coords]
                assert speci>=0
                nei=lattice.neighbors[coords]
                neisp=x[zip(*nei)]
                free=[nei[n] for n in np.where(neisp<0)[0]]
                neisp=neisp[neisp>=0]
                if measure:
                    #Measure abundances and pair contacts
                    abund[speci]+=1
                    contact[speci,list(neisp)]+=1

                ints=np.array([comm[speci,specj] if speci!=specj else selfint[speci]
                    for specj in neisp])

                #Birth
                if free:
                    birth=growth[speci]-np.sum(ints[ints<0])/nnei
                    if np.random.random() <birth*dt:
                        pos=free[np.random.randint(0,len(free))]
                        x[pos]=speci
                        #print np.array(coords)-pos

                #Death
                death=np.sum(ints[ints>0])/nnei
                if np.random.random()<death*dt:
                    x[coords]=-1

            if measure:
                self.append_results(n=abund,
                    contact=contact,x=x,index=t)
            t+=dt

        return True


BaseModel.klasses['CellularModel']=CellularModel
