# -*- coding: utf-8 -*-
from imports import *
from utility import *
import numpy as np


def lifted_matrix(A,x=None):
    '''Lifted operator from interaction matrix (and abundances)'''
    #D=np.diag(x) #diagonal abundances
    I=np.eye(*A.shape)
    #print A.shape,I.shape
    Jhat=np.kron(A,I)+np.kron(I,A)
    return Jhat

def variance_interactions(Jhat,lift=False):
    '''Interaction matrix for variances (given a press noise).'''
    if lift:
        Jhat=lifted_matrix(Jhat)
    S=int(np.round(np.sqrt(Jhat.shape[0])))
    return np.linalg.inv(Jhat)[::S,::S]

def activity(x,A,r):
    '''Given interactions A and growth r, check nodes and edges for activity
    (i.e. dynamical role)'''

    #Node states
    inter=np.dot(np.abs(A),x)
    diag=np.diag(A)*x
    node_state=np.ones(x.shape) #By default, other dominated (1)
    node_state[np.abs(r)>np.abs(inter)+np.abs(diag)]=-1 #resource dominated (-1)
    node_state[np.abs(diag)>np.abs(inter/2)]=0 #self dominated (0)
    node_state[x<10.**-15]=-2 #dead (-2)

    #Edge activity
    edg=((A.T*x).T)*x
    comp=(np.ones(A.shape)*diag*x).T
    active=zip(*np.where(np.abs(edg)>np.abs(comp)*1.01 ))
    return node_state,active


def interactions_to_fluxes(A):
    """assuming unit populations, derive biomass fluxes from interactions"""

    epsilon=np.array(A.props['efficiency'])
    F=  A.matrix - (epsilon.T * A.matrix.T + A.matrix )/(1+ epsilon.T * epsilon)
    return F


def naive_trophic_score(A,r,z=.01):
    M=np.eye(A.shape[0])-z*A
    try:
        result= np.log(np.dot(np.linalg.inv(M),r ))/np.log(z)
    except Exception as e:
        print e
        result=np.zeros(r.shape)
    return result

def trophic_score(A,pops=None,influx=None,EPS=10**-15,remove_diag=1,ntries=10):
    #mat=A.matrix
    if pops is None:
        pops=np.ones(A.shape[0])

    fluxes=A.copy()
    fluxes[fluxes<0]=0
    pops=pops.reshape(pops.shape+(1,))
    #print fluxes.shape, pops.shape
    transfer=fluxes*pops*pops.T
    if remove_diag:
        np.fill_diagonal(transfer,0)
    #print transfer, fluxes

    pops = pops.squeeze()
    if influx is None:
        #infer basal from having no prey
        basal=( np.max(A-np.diag(np.diag(A)),axis=1) <=0 ).astype('float')
        #print basal
        #print np.dot(A,basal)
        #print np.dot(A, np.dot(A,basal))
        #print np.sum(A[0])
        #raise
    else:
        # pops=pops.squeeze()
        influx=influx.squeeze()
        influx[influx<0]=0
        basal=influx
        if not basal.any():
            return np.zeros(pops.shape)

    norm=np.sum(transfer,axis=1)+basal
    #print norm,np.sum(transfer,axis=1)
    norm[norm==0]=1
    transfer=transfer/norm.reshape((transfer.shape[0],1))
    transfer[transfer>1]=1-EPS
    dangers=np.where(np.sum(fluxes,axis=1)+np.sum(fluxes,axis=0) ==0 )[0]
    transfer[dangers]=0


    I=np.eye(fluxes.shape[0])
    bad= np.where(np.max(np.abs(transfer-I),axis=1)<EPS)[0]  #np.max(transfer,axis=1),np.sum(A!=0,axis=1)#, fluxes[0]#,transfer[0]

    tries =0
    t=[10**15]
    while np.max(np.abs(t))> 10*transfer.shape[0]:
        if tries>ntries:
            t=naive_trophic_score(A, influx, z=.01)
            t[np.isinf(t)]=np.max(t[np.isinf(t)<1 ])
            t[np.isnan(t)]=np.max(t[np.isnan(t)<1 ])
            return t

        if tries > 0:
            print ' !Trophic score: removing species that cause a divergence'
            dangers=np.unique(np.concatenate(np.where(np.linalg.matrix_power(transfer, 100*transfer.shape[0] ) > 10 ** -3)))
            transfer[np.ix_(dangers,dangers )] = 0

        if tries > 5:
            print ' !Trophic_score error? Probably subgraph disconnected from any source'
            # print np.sum(transfer,axis=1),t
            t[t>A.shape[0]]=np.max(t[t<A.shape[0]])
            # t=naive_trophic_score(A, influx, z=.01)
            # t[np.isinf(t)]=np.max(t[np.isinf(t)<1 ])
            # code_debugger()
            return t
            #ERROR!
            print A, pops, influx
            print t
            from datatools import hist,scatter,plt
            hist(t.ravel(),hold=1)
            plt.figure()
            scatter(t.ravel(),basal.ravel())
        try:
            t = np.dot(np.linalg.inv(I - transfer), np.ones(pops.shape))  # .A1
        except Exception as e:
            print ' !Trophic score ERROR:', e
            if tries > 1:
                code_debugger()
                return np.zeros(pops.shape)
            t = [10 ** 15]
        tries+=1
    # code_debugger()
    return t


if __name__=='__main__':

    def test_trophic_score():
        import networkx as nx
        #from utility_graphs import Graph
        g=Graph()
        g.add_edge(0,2)
        g.add_edge(1,5)
        #g.add_edge(5,3)
        g.add_edge(0,1)
        g.add_edge(2,1)
        g.add_edge(1,3)
        #g.add_edge(2,4)
        g.add_edge(3,4)
        g.add_edge(4,5)

        mat=nx.to_numpy_matrix(g)
        #mat+=np.random.random(mat.shape)/10000.
        #mat[4]=0
        vec=np.array([0,0,0,0,.5,1] )
        from graphs import CommNetwork
        A=CommNetwork(mat*.1,symmetry=-0.1)
        res= trophic_score(np.array(A.matrix),capacity=vec)
        print 'NEW',res
        #print 'OLD',old_trophic_score(mat,vec)
        g.plot(hold=1,labels=['{:.3}'.format(r) for r in res])
        g.plot()
    test_trophic_score()