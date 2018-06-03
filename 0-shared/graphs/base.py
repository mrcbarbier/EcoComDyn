# -*- coding: utf-8 -*-
#import numpy as np
from utility import *
from sampler import BaseSampler


class BaseNetwork(object):
    '''Base class for generic networks'''

    _name='network'

    def __init__(self,matrix=None,**kwargs):
        if hasattr(matrix,'shape') and not matrix.shape:
            matrix=matrix.reshape((1,))
        self.matrix=matrix
        self.props=kwargs
        self.name=kwargs.get('name',self._name)

    def __getitem__(self,*args,**kwargs):
        return self.matrix.__getitem__(*args,**kwargs)

    def __repr__(self):
        return "BaseNetwork(matrix={},**{})".format(self.matrix,self.props)
    def __str__(self):
        return 'Network({})'.format(self.matrix)

    def __add__(self,*args,**kwargs):
        return self.matrix.__add__(*args,**kwargs)
    def __mul__(self,*args,**kwargs):
        return self.matrix.__mul__(*args,**kwargs)

    def __setitem__(self,*args,**kwargs):
        return self.matrix.__setitem__(*args,**kwargs)

    def save(self,path,fname=None):
        if fname is None:
            fname=self.name
        import pandas as pd
        idx,mat=mat_to_cols(self.matrix,drop_zeros=0)
        if len(mat)>1000000:
            t=0
            with  open(path+fname,'w') as f:
                f.write('\tweight\n')
                for i,j in zip(idx,mat):
                    f.write('{}\t{}\n'.format(i,j) )
                    t+=1
                    if t%100000==0:
                        print t
            return
        df=pd.DataFrame([{'weight':j} for j in mat ],index=idx  )
        df.to_csv(path+fname,sep='\t')

    def load(self,path,fname=None):
        if fname is None:
            fname=self.name
        import pandas as pd,scipy.sparse as spar
        df=pd.read_csv(path+fname,sep='\t',header=0,index_col=0)
        self.matrix= dict_to_mat({eval(i):j for i,j in df['weight'].iteritems()})
        return self


    def copy(self):
        return self.__class__(self.matrix.copy(),**deepcopy(self.props))


    @classmethod
    def load_matrix(klass,path):
        return np.loadtxt(path)

    @classmethod
    def make_edges(klass,**attrs):

        shape=tuple(int(x) for x in to_iterable(attrs.get('shape',100),tuple) )
        distribution=attrs.get('distribution','random')
        if distribution!='uniform':
            connectivity=attrs.get('connectivity',None)
            if connectivity is None:
                connectivity=float(attrs.get('mean',shape[0]))/shape[0]

        directed=attrs.get('directed',shape!=shape[::-1] )
        #IF CLASSIC NETWORK ENSEMBLES
        structure=attrs.get('structure','')
        if structure in ('barabasi','smallworld'):
            import networkx as nx, networkx.generators.random_graphs as nxgen
            n=int(shape[0])
            m=int(np.round(connectivity * n))
            if structure == 'barabasi':
                g=nxgen.barabasi_albert_graph(n, m/2,)
            elif structure == 'smallworld':
                p=attrs.get('rewiring',.1)
                g=nxgen.watts_strogatz_graph(n, m, p)

            mat= np.array(nx.to_numpy_matrix(g))
            nedges=np.sum(mat!=0)
            if not directed:
                rdm=(np.random.random(mat.shape)< nedges*1./mat.size/2.)
                rdm+=rdm.T
                replace=(np.random.random(mat.shape)> (1+attrs.get('order',1))/2. )
                replace+=replace.T
            else:
                rdm=(np.random.random(mat.shape)< nedges*1./mat.size)
                replace=(np.random.random(mat.shape)> attrs.get('order',1) )

            mat[replace]=rdm[replace]
            #return mat


        #ELSE HANDCRAFTED ENSEMBLES


        elif  'multipartite' in structure or 'modular' in structure:
            mat=np.ones(shape)
            fractions=attrs.get('fractions',(0.5,0.5))
            if 'modular' in structure:
                order=attrs.get('modularity',1)
            else:
                order=attrs.get('partition',1)
                if 'modularity' in attrs:
                    print 'DANGER: may have confused modulariy and partition! (they are opposite)'
            fractions=np.array(fractions)/np.sum(fractions)
            part=np.random.choice(range(len(fractions)),p=fractions,size=mat.shape[0])
            if 'modular' in structure:
                avoid=(np.add.outer(part,-part)!=0 )
            else:
                avoid=(np.add.outer(part,-part)==0 )
            remove=np.logical_and(avoid,np.random.random(avoid.shape)<1-np.sqrt(1-order) )
            #1-sqrt(1-order) ensures that I remove the correct fraction when considering the OR below
            remove=np.logical_or(remove,remove.T)
            mat[remove]=0

            #print order, np.sum(mat[0]>0),fractions, float(np.sum(samepart))/samepart.size

        elif 'nested' in structure:
            mat=np.ones(shape)
            fij=mat
            mat=np.zeros(fij.shape)
            order=attrs.get('nestedness',1)
            if order==0:
                mat=fij
            else:
                mat=fij.copy()
                N=mat.shape[1]
                idxs=[ i
                 for i in zip(*[z.ravel() for z in tuple(np.indices(fij.shape ))])
                 if i[0]>=N-i[1]
                ]

                #idxs=[i if np.random.random()<(1.+order)/2. else (i[1],i[0])+i[2:]
                    #for i in idxs]
                idxs=[i for i in idxs if np.random.random()<(1.+order)/2. ]
                idxs=zip(*idxs)

                mat[idxs]=0#fij[idxs]
                #print order
            #np.fill_diagonal(mat,0)
            onlyone=(fij*fij.T==0)
            mat[onlyone]=fij[onlyone]

        elif 'triangular' in structure:
            mat=np.ones(shape)
            fij=mat
            mat=np.zeros(fij.shape)
            order=attrs.get('triangular',1)
            if order==0:
                mat=fij
            elif order<1:
                mat=fij.copy()
                idxs=[ i
                 for i in zip(*[z.ravel() for z in tuple(np.indices(fij.shape ))])
                 if i[0]>i[1]
                ]

                #idxs=[i if np.random.random()<(1.+order)/2. else (i[1],i[0])+i[2:]
                    #for i in idxs]
                idxs=[i for i in idxs if np.random.random()<(1.+order)/2. ]
                idxs=zip(*idxs)

                mat[idxs]=0#fij[idxs]
                #print order

            else:
                mat=np.triu(fij)
                np.fill_diagonal(mat,0)
            onlyone=(fij*fij.T==0)
            mat[onlyone]=fij[onlyone]

        else:

            mat=np.zeros(shape)
            rest=np.product(shape[1:])
            np.fill_diagonal(mat,1)
            left=np.zeros(shape[0])
            for i in range(shape[0]):
                #number of edges for row i
                if connectivity==1:
                    nb=shape[0]
                elif distribution =='uniform':
                    nb=np.random.randint(attrs.get('min',1),min(rest,attrs.get('max',rest)) )
                else:
                    mean=connectivity*shape[0]
                    if distribution =='normal':
                        std=attrs.get('std',attrs.get('stdrel',0)*mean )
                        nb=max(1,np.random.normal(mean,std,[1] ).astype('int'))
                    elif distribution =='random':
                        #nb=max(1,np.random.exponential(mean,[1]).astype('int') )
                        nb=mean
                left[i]=nb
            tmp=[]
            for i in range(shape[0]):
                def rndrge(val):
                    rge=range(val)
                    np.random.shuffle(rge)
                    return rge

                nb=left[i]
                if nb<=0:
                    continue
                if directed:
                    tup=np.where(mat[i]==0)
                else:
                    tup=np.where( np.logical_and(mat[i]==0,left>0) )
                    #print len(tup[0])
                coords=rndrge(len(tup[0]))[:int(round(nb)) ]
                tup=tuple([t[c] for c in coords] for t in tup )
                mat[ (i*np.ones(len(tup),dtype='int'),)+tup  ]=1
                if not directed:
                    tmp.append(len(tup[0])*2 )
                    mat[tup+(i*np.ones(len(tup),dtype='int'),) ]=1
                    for t in set(tup[0]):
                        #print t
                        left[t]-=1
                left[i]=0
        np.fill_diagonal(mat,0)


        if 'assortativity' in attrs:
            from datatools import scatter
            for i in range(6):
                degree=np.sum(mat!=0,axis=1).astype('float')
                deg=np.sum(mat!=0,axis=1).astype('float')
                degree/=np.max(degree)
                #print degree
                fij=mat
                mat=np.zeros(fij.shape)
                order=attrs.get('assortativity',1)
                links=zip(*np.where(np.triu(fij)!=0 ) )
                links=sorted(links, key=lambda e:degree[e[0]])
                targets=sorted(links,key=lambda e: order*degree[e[1]] + (1-np.abs(order))*np.random.random()  )
                links=zip(*links)
                targets=zip(*targets)
                #print links[0]
                #scatter(degree[list(links[0])], degree[list(targets[1])])
                #print  links[:1]+targets[1:]
                #mat[ (links[1],targets[0]) ]=1
                mat[ (links[0],targets[1]) ]=1
                mat=mat+mat.T
            newdegree=np.sum(mat!=0,axis=1).astype('float')
            #newdegree/=np.max(newdegree)
            #scatter(deg,newdegree)
            nlinks=np.where(np.triu(mat)!=0 )
            #scatter(newdegree[list(nlinks[1])], newdegree[list(nlinks[0])])
            #print np.corrcoef(newdegree[list(nlinks[1])], newdegree[list(nlinks[0])])[0,1]
            #print degree


        if 'clustering' in attrs:
            #print np.sum(mat)
            order=attrs.get('clustering',1)
            if order!=0:
                fij=np.triu(mat)#np.zeros(fij.shape)
                #where=np.where(fij!=0 )
                #links=zip(*where )

                tri=np.dot(fij,fij)
                con3=zip(*np.where( np.logical_and(tri!=0,fij!=0) ))
                con2=zip(*np.where( np.logical_and(tri!=0,fij==0) ))
                con1=zip(*np.where( np.logical_and(tri==0,fij!=0) ))
                nb=int(np.round(np.abs(order)*min(len(con2),len(con1))))
                for n in range(nb):
                    if order>0:
                        tgt=con2.pop(np.random.randint(0,len(con2)) )
                        #con3.append(tgt)
                        fij[tgt]=1
                    else:
                        src=con3.pop(np.random.randint(0,len(con3)) )
                        #con2.append(src)
                        fij[src]=0

                mat=fij+fij.T
            #mat=fij
            #np.random.shuffle(links)

        #print "DIRECTED",directed
        #if directed :
            #return mat
        #mat= np.clip(mat+mat.T,0,1)
        #if not directed:
            #assert (mat == mat.T).all()

        # code_debugger()
        if np.mean(mat!=0)> connectivity:
            mat[np.random.random(mat.shape)> connectivity/np.mean(mat!=0) ]=0
        return mat

    @classmethod
    def symmetrize(klass,mat,attrs):
        if np.std(mat)>0:
            g=(mat-np.mean(mat))/np.std(mat)
        else:
            g=mat
        if 'symmetry' in attrs:
            gamma=attrs['symmetry']
            if gamma is None:
                gamma=0
            assert -1 <= gamma <= 1
            if gamma==0:
                asymm=0
            else:
                asymm=(1- np.sqrt(1-gamma**2) )/gamma
            g= (g+asymm * g.T) / np.sqrt(1+ asymm**2 )
            mat=g*np.std(mat)+np.mean(mat)
        return mat

    @classmethod
    def randomize(klass,base,**attrs):
        #Randomly redistribute links
        directed=attrs.get('directed',False)
        mat=base.copy()
        rdm=np.random.random(mat.shape)
        if not directed:
            rdm=np.triu(rdm)

        rewire=np.where(np.logical_and(rdm> attrs['ordered'],np.eye(mat.shape[0])!=1) )
        target=np.array(zip(*rewire))
        np.random.shuffle(target)
        target=tuple(target.T)
        if directed:
            mat[rewire]=mat[target]
        else:
            def cat_transpose(lst):
                return tuple(np.concatenate(x) for x in zip(*[lst,lst[::-1]]) )
            mat[cat_transpose(rewire) ]=mat[cat_transpose(target)]
        return mat



    @classmethod
    def make(klass,**attrs):
        attrs.setdefault('connectivity',attrs.get('connectance',1) )
        """Generate a matrix"""
        def get_attr(attr,dft):
            if not attr in attrs:
                return dft
            a=attrs[attr]
            if isinstance(a,basestring):
                a=get_attr(a,dft)
            return a

        OBSOLETES=['edgestruct','edgeorder','edgestdrel','edgemean']
        for i in attrs:
            if i in OBSOLETES:
                print 'WARNING: Graph.make -- obsolete attr:',i

        if 'matrix' in attrs:
            return klass( matrix=attrs.pop('matrix'),**attrs)

        shape=to_iterable(attrs.get('shape',100),tuple)

        sampler=attrs.get('sampler',BaseSampler(shape=shape)).sample


        role=attrs.get('role',None)
        mat=None
        mode=attrs.get('mode',None)
        if mode=='network' and not 'from' in attrs:
            edges=klass.make_edges(**attrs)
                 #**{i:attrs.pop(i) for i in tuple(attrs) if 'edge' in i} )
            connectivity=attrs.pop('connectivity',np.sum(edges>0)*1./edges.size)
            mat=1
            structure=''
        else:
            edges=1
            connectivity=1
            structure=attrs.get('structure','random')

        if 'from_data'  in attrs:
            path=Path(attrs['from_data'])
            mat=klass().load(path,fname='').matrix

        elif 'from' in attrs:
            base=np.array(np.loadtxt(attrs['from']))
            if base.shape!=attrs['shape']:
                if np.sum(base.shape)==attrs['shape'][0]:
                    #Multipartite network
                    import scipy.sparse as sparse
                    G=nx.algorithms.bipartite.from_biadjacency_matrix(sparse.lil_matrix(base))
                    base=np.array(nx.to_numpy_matrix(G))
                else:
                    raise Exception('Incorrect shape for loaded network')
            if not 'randomize' in structure:
                if structure:
                    print 'Graph.make: Ignoring structure due to from:',attrs['from']
            mat=base

        elif 'edgefrom' in attrs:
            frm=attrs['edgefrom']
            if frm in attrs['required']:
                edges=attrs['required'][frm].matrix
            else:
                edges= klass.load_matrix(frm)

        if mode == 'network' and role == 'fluxes':
            print "base.py: HACK WHERE FLUXES ARE EDGES * SAMPLER"
            mat *= sampler(**attrs)

        # ================================ GENERIC STRUCTURE ===================================

        if structure == 'random':
            mat=sampler(**attrs)
        if 'randomize' in structure:
            mode=attrs.get('rdmode','basic')
            directed=attrs.get('directed',False)
            base=mat.copy()
            if mode=='basic':
                mat=klass.randomize(base,**attrs)
                print 'COMPARISON AT END OF RANDOMIZATION',np.sum(mat!=base)

            elif mode in ('cavity','cavsimple'):
                def objective(mat):
                    if mode =='cavsimple':
                        return np.var(mat)
                    return np.var(np.mean(mat,axis=1) )#,np.mean(np.var(mat,axis=1)) ))
                tgtchanges=(1-attrs['ordered'])*np.sum(base!=0)#/2
                changes=0
                trials=30
                ref=objective(base)
                mat=base.copy()
                S=mat.shape[0]
                beta=1./ref**2*S**2
                while changes <tgtchanges*.95 and trials>0:
                    src=[]
                    tgt=[]
                    delta=objective(mat)-ref
                    links=zip(*np.where(np.triu(mat)!=0 ) )
                    np.random.shuffle(links)
                    links=links[:S/10 ]
                    cur = links.pop(0)
                    accepted=None
                    while cur:
                        s,t=cur
                        partner=range(mat.shape[1])
                        partner.remove(t)
                        np.random.shuffle(partner)
                        while partner and not accepted:
                            t2=partner.pop(0)
                            if accepted is None:#first link change
                                break
                            wei1,wei2=mat[s,t],mat[s,t2]
                            if mode=='cavsimple':
                                delta_var_mean=(wei2**2-wei1**2) - (wei2-wei1)**2
                            if mode=='cavity':
                                delta_wei=(wei1-wei2)*1./S
                                mean_row1=np.mean(mat[t])
                                mean_row2=np.mean(mat[t2])
                                delta_var_mean = 2* delta_wei* (delta_wei + mean_row2-mean_row1)
                                tmp=mat.copy()
                                tmp[([s,s,t,t2],[t,t2,s,s])]=mat[([s,s,t2,t],[t2,t,s,s])]
                                #print delta_var_mean, objective(mat)-objective(tmp)
                                delta_var_mean=objective(tmp)-objective(mat)

                            accepted=np.random.random()< np.exp(- beta* delta_var_mean*delta )
                            #print delta_var_mean*delta <0,accepted,- beta* delta_var_mean*delta

                        mat[([s,s,t,t2],[t,t2,s,s])]=mat[([s,s,t2,t],[t2,t,s,s])]
                        accepted=False
                        #print delta<objective(mat)-ref
                        changes+=1.
                        if changes>tgtchanges or not links:
                            break
                        cur = links.pop(0)
                    tgt,src=[tuple(zip(z)) for z in (tgt,src)]
                    mat[tgt]=mat[src]
                    trials-=1
                print 'COMPARISON AT END OF RANDOMIZATION',np.sum(mat!=base),objective(mat),ref
            else:
                raise Exception('Not useful yet')

                prevmat=base.copy()
                target=objective(prevmat)
                current =0
                beta=10.
                orderleft=1
                trialsleft=1000
                nonzero=2.*np.sum(base!=0)
                while orderleft>0 and trialsleft>0:
                    mat=klass.randomize(prevmat,directed=directed,ordered=.8)
                    value=np.sum(np.abs(objective(mat)/target-1))
                    rdm=(np.random.random(mat.shape)<np.exp(-beta*(value-current ) ))
                    prevmat[rdm]=mat[rdm]
                    current=np.sum(np.abs(objective(prevmat)/target)-1)
                    orderleft=1-attrs.get('ordered',1)-np.sum(mat!=base)/nonzero
                    trialsleft-=1
                mat=prevmat
                print 'COMPARISON AT END OF RANDOMIZATION',np.sum(mat!=base),objective(mat),target,current

            if 'maxrow' in attrs:
                #Useful for mutualism
                mat = mat/np.max(np.sum(mat,axis=1))*attrs['maxrow']



        # ================================ SCALINGS ===================================

        if structure == 'bunin':
            S=shape[0]
            effS=float(shape[-1])*connectivity ##Important: if not square matrix, scale with number of sources, not targets
            effS=max(1,effS)
            mu=attrs.get('mean',1)
            sigma=attrs.get('std',np.sqrt(attrs.get('var',1)))
            gamma=np.clip(attrs.get('symmetry',1),-1,1)

            if 'corr_KA' in attrs and np.abs(attrs['corr_KA'])>10**-8:
                #correlation between carrying capacity and A
                K=attrs['required'].values()[0].matrix.copy().reshape((S,1) )
                K/= np.mean(K)

                sigma_K=np.std(K)
                if min(sigma,sigma_K)>10**-8:
                    c=-attrs['corr_KA']/sigma_K**2
                    gamma*= (1+ effS * c**2 * sigma_K**2/sigma**2)
                else:
                    c=0
                sigma=np.sqrt(sigma**2- effS*c**2 *sigma_K**2)
                if np.isnan(sigma):
                    sigma=10**-8
                    c=0

                if not -1 <= gamma <=1:
                    print 'DANGER',gamma, effS, sigma,c, sigma_K,attrs.get('symmetry',1)
                gamma=np.clip(gamma,-1,1)
            else:
                sigma_K=0
                c=0
                K=0

            assert -1 <= gamma <= 1
            if gamma==0:
                asymm=0
            else:
                asymm=(1- np.sqrt(1-gamma**2) )/gamma


            aij=sampler(mean=0,std=1./np.sqrt(effS))

            aij= (aij+asymm * aij.T) / np.sqrt(1+ asymm**2 )

            mat=-mu/effS - sigma*aij

            if not c is 0:
                mat+=c*(K-1 )

            np.fill_diagonal(mat,0)

            if 0:#'corr_KA' in attrs and 0:
                #Test that corr_KA is generated correctly
                from datatools import nonan
                #print np.mean(edges),np.mean(aij),mu,sigma,effS
                mm=-mat.copy()*edges
                np.fill_diagonal(mm,np.nan)
                if not edges is 1:
                    mm[edges==0]=np.nan
                m=nonan(mm)
                if len(m)==0:
                    m=np.zeros(shape)
                mu=attrs.get('mean',1)
                sigma=attrs.get('std',np.sqrt(attrs.get('var',1)))
                test=gamma/(1+effS * c**2 * sigma_K**2/sigma**2)

                #print np.sum(edges)
                print 'mu', np.mean(m)*effS,mu,np.sum(-mat*edges)/S
                print 'sig',np.var(m)*effS,sigma**2,np.sum((mat+mu/effS )**2*edges)/S
                print 'gam',np.corrcoef(mm[~np.isnan(mm) ],mm.T[~np.isnan(mm.T) ] )[0,1],test,# (np.mean(nonan(mm*mm.T))-np.mean(m)**2)/(sigma**2/effS),test,np.sum(( ( mat+mu/effS)*(mat.T+mu/effS)*edges))/np.sum((mat+mu/effS )**2*edges)
                print 'cor',np.mean(nonan(mm*K))-np.mean(m)*np.mean(K), attrs.get('corr_KA',0)#,c
                print 'connectivity',np.mean(np.sum(np.abs(mat*edges)>
                    np.max(np.abs(mat),axis=1)*0.05,axis=1) )/mat.shape[0],connectivity



        # ================================ RESOURCE COMPETITION ===================================
        if structure=='resource':
            if role == 'consumption':
                if 'niche' in attrs['required']:

                    mat0=np.abs(sampler(**attrs ))
                    r=attrs['required']['resource']
                    R=r.shape[0]
                    pos=np.linspace(0,1,R)

                    niches=attrs['required']['niche'].matrix
                    wid=attrs['required']['niche_width'].matrix
                    var=(wid/2)**2
                    dist=np.add.outer(niches,-pos)
                    var=var.reshape(var.shape+tuple(1 for i in dist.shape[1:]) )
                    mat=np.clip(np.exp(- dist**2/(2*var) )/np.sqrt(2*np.pi * var),10**-10,100)
                    mat *= mat0
                    mean=np.mean(mat)
                    std=np.std(mat)
                    if 'randomize' in attrs and attrs['randomize']:
                        #print mat[0]
                        for i in range(mat.shape[0]):
                            tmp=mat[i]
                            np.random.shuffle(tmp)
                            mat[i]=tmp
                        #print mat[1]
                else:
                    mat=np.abs(sampler(**attrs ))
                    mean=attrs['mean']
                    std=attrs['std']
                norm=np.sqrt((mean**2+std**2)*connectivity*mat.shape[1])
                mat/=norm
            elif role=='diagonal':
                cons=attrs['required']['consumption']
                mat=np.diag(np.dot(cons.matrix,cons.matrix.T))+0.1#np.diag(cons.dot(cons,axes=[('r','com')] ).matrix)
                #mat=mat/np.mean(mat)

            elif role=='growth':
                cons=attrs['required']['consumption']
                r=attrs['required']['resource']
                m=attrs['required']['mortality']
                R=r.shape[0]
                r=r/(1-np.mean(m.matrix))/(R* np.mean(cons.matrix))
                mat=np.dot(cons.matrix,r.matrix)
                mat-= np.mean(r.matrix)*np.dot(cons.matrix,np.ones(r.shape) )*m.matrix
                #mat=mat -m.matrix# (np.sum(r.matrix)*np.mean(cons.matrix) -1)

            elif role=='interactions':
                cons=attrs['required']['consumption']
                mat=np.dot(cons.matrix,cons.matrix.T)#cons.dot(cons,axes=[('r','com')] ).matrix
                mean=np.mean(np.diag(mat))
                mat=-mat
                #mat/=mean
                np.fill_diagonal(mat,0)


        if structure=='constraint':
            shape=attrs.pop('shape')
            constraint=attrs.pop('constraint')
            mat =  sampler(shape=shape[0],**attrs )
            hori= sampler(shape=shape[0],min=0,max=attrs['horizon'],distribution='uniform')
            bad=[True]
            while len(bad):
                dico={}
                dico.update(attrs)
                dico.update({'x':hori,'y':mat })
                bad = np.where( eval(constraint,{},dico)==False )[0]
                if len(bad):
                    hori[bad]=sampler(shape=len(bad),min=0,max=attrs['horizon'],distribution='uniform')

            # code_debugger()
            mat=np.array(zip(mat,hori))
        # ================================ NICHE BASED ===================================
        attrs.setdefault('required',[])
        niche=None
        if 'niche'  in attrs['required']:
            niche=attrs['required']['niche'].matrix
            if len(niche.shape)==2:
                niche,other_niche=niche[:,0],niche[:,1]
            else:
                other_niche=None
            if len(niche.shape)>2:
                print 'WARNING: Niche matrix has dim>2, not ready for this!'

        if structure=='ranked':
            mode=attrs.get('mode','linear')
            niche=(niche-np.min(niche))/(np.max(niche)-np.min(niche))
            aij=sampler(mean=0,std=1 )
            aij=klass.symmetrize(aij,attrs)
            if role =='interactions':
                dist=np.add.outer(niche,-niche)
                mu=attrs.get('mean',1)
                c_mu=attrs.get('meanslope',0)
                sigma=attrs.get('std',.1)
                c_sigma=attrs.get('stdslope',0)

                if mode =='linear':
                    mat=mu+c_mu*dist + np.sqrt(np.clip(sigma**2 + c_sigma*dist,0,None))*aij
                if mode =='exponential':
                    c_mu=np.exp(c_mu)
                    c_sigma=np.exp(c_sigma)
                    mat=mu*c_mu**dist + np.sqrt(sigma**2 * c_sigma**dist)*aij
            elif role=='growth':
                if 'min' in attrs:
                    K=attrs.get('min',1)
                    nich=niche
                else:
                    K=attrs.get('mean',1)
                    nich=2*niche-1
                c_K=attrs.get('meanslope',0)
                zeta=attrs.get('std',1)
                c_zeta=attrs.get('stdslope',0)
                if mode =='linear':
                    mat=K+c_K*nich + np.sqrt(np.clip(zeta**2+c_zeta*nich,0,None) )*aij
                if mode =='exponential':
                    c_K=np.exp(c_K)
                    c_zeta=np.exp(c_zeta)
                    mat=K*c_K**nich  + np.sqrt(zeta**2 * c_zeta**nich)*aij

        if structure=='nichecomp':
            cons=attrs['required']['consumption'].matrix
            if role =='interactions':
                nvar=attrs['required']['niche_width'].matrix**2
                dist=np.add.outer(niche,-niche)
                sigs=np.add.outer(nvar,nvar)
                totcons=np.multiply.outer(cons,cons)
                mat=-totcons*np.exp(- dist**2/(2*sigs) )/np.sqrt(2*np.pi * sigs)
                np.fill_diagonal(mat,0)

            elif role=='growth':
                m=attrs['required']['mortality'].matrix
                mat=cons*(1-m)
            elif role=='diagonal':
                nvar=attrs['required']['niche_width'].matrix**2
                mat=cons**2/np.sqrt(2*np.pi* nvar)


        # ================================ TROPHIC ===================================
        if 'efficiency' in attrs['required']:
            epsilon=attrs['required']['efficiency'].matrix
            if epsilon.shape!=attrs['shape']:
                epsilon=epsilon.reshape(epsilon.shape+tuple(1
                    for i in range( len(attrs['shape']) - len(epsilon.shape) )  ) )
        else:
            epsilon=attrs.get('efficiency',1)


        if structure =='trophic':
            if role in ['interactions','Amean']:
                fij=attrs['required']['fluxes'].matrix
                if niche is None:
                    pows=np.arange(fij.shape[0])
                else:
                    pows=niche-1
                fij=(fij.T * attrs['mean']*get_attr('rate',1)**pows ).T
                fij[fij.T>fij ]=0
                mat=epsilon*fij-fij.T
            if role =='diagonal':
                mat= np.ones(shape)*get_attr('rate',1)**(niche-1)*get_attr('mean',1)
            if role in ['growth','rmean']:
                mat=sampler(**attrs)
                if  niche is None:
                    pows = np.arange(mat.shape[0])
                    mat[1:]=-get_attr('mortality',0)*get_attr('rate',1) ** pows[1:]
                else:
                    pows=niche[np.where(niche>=np.min(niche)+1)]-1
                    mat[np.where(niche>=np.min(niche)+1)]=-get_attr('mortality',0)*get_attr('rate',1) ** pows
            if role=='fluxes':
                mat = np.zeros(attrs['shape'], dtype='float')
                mat[kth_diag_indices(mat, -1)] = 1

        if structure=='envelope':
            niche=niche.astype('float')
            if role=='interactions':
                fij=attrs['required']['fluxes'].matrix.copy()
                if not other_niche is None:
                    hori=other_niche
                    horidist=np.add.outer(-hori,hori).astype('float')
                    fij[np.abs(horidist)>.5]=0
                dist=np.add.outer(-niche,niche).astype('float')
                rang=attrs['range']
                if hasattr(rang,'keys'):
                    rg=rang
                    rang=rg.get('dist',1)
                    rgtyp=rg.get('type','fixed')
                    rgshape=rg.get('shape','gate')
                    order=rg.get('order',1)
                    sd=rg.get('wid',2*rang )/2
                    exp=rg.get('exp',1)
                else:
                    rgtyp=attrs.get('rangetype','fixed')
                    rgshape=attrs.get('rangeshape','gate')
                    order=attrs.get('rangeorder',1)
                    sd=attrs.get('rangewid',2*rang)/2
                    exp= attrs.get('rangeexp', 1)
                rang=np.array(rang,dtype='float')

                arg = np.argsort(niche)
                if rgtyp=='fraction':
                    oldrang=rang
                    rang= rang*niche.reshape( niche.shape+(1,) )/2.
                    sd=sd*rang/oldrang

                if order<1:
                    moved=(np.random.random(fij.shape[0])>order)
                    if not rang.shape:
                        rang=rang*np.ones(fij.shape[0])
                    rang[moved]= np.random.random(rang[moved].shape)*rang[moved]

                if rgshape=='gate':
                    fij[dist>-rang+sd ]=0
                    fij[dist<-rang-sd ]=0
                elif rgshape=='exp':
                    fij=fij*np.exp(-(-dist/rang))
                    fij[dist>0]=0
                elif rgshape == 'biexp':
                    fij = fij * np.exp(-np.abs(-dist / rang)**exp )
                    # fij[dist > 0] = 0
                elif rgshape=='exp3':
                    fij=fij*np.exp(-(-dist/rang)**3 )
                    fij[dist>0]=0
                elif rgshape=='normal':
                    fij=fij*np.exp(-(-dist-rang)**2/(2*(0.3*sd)**2) )
                mat=epsilon*fij-fij.T
                if attrs.get('normalize',0):
                    mat/=np.mean(np.sum(fij!=0,axis=1) )
                    # print np.mean(mat[mat>0])
                # code_debugger()

        if structure=='niche':
            if role=='capacity':
                aij=np.abs(sampler(**attrs ))
                mat=aij/ np.sum(niche<1) * niche.size
                mat[niche>1]=0
            elif role=='growth':
                mat=np.abs(sampler(**attrs ))
                thresh=1.
                mat=mat*1./ np.sum(niche-np.min(niche)<thresh) * niche.size
                mat[niche-np.min(niche)>=thresh]=-attrs.get('mortality',0.1)
                mat/=np.max(mat)
            elif role=='fluxes' or role=='interactions':
                if role =='fluxes':
                    fij = sampler(**attrs)
                else:
                    fij=attrs['required']['fluxes'].matrix.copy()
                dist=np.add.outer(-niche,niche)
                rgtyp=attrs.get('rangetype','fixed')
                rang=attrs.get('range',(0.5,1.5) )
                if hasattr(rang,'__iter__'):
                    minrang,maxrang=rang
                else:
                    minrang,maxrang=0,rang
                if rgtyp=='fixed':
                    if not maxrang is None:
                        fij[dist<-maxrang]=0
                    if not minrang is None:
                        fij[dist>-minrang]=0
                elif rgtyp=='fraction':
                    if not maxrang is None:
                        fij[(dist.T<-maxrang*niche).T]=0
                    if not minrang is None:
                        fij[(dist.T>-minrang*niche).T]=0
                np.fill_diagonal(fij,0)
                #If bidirectional fluxes, then randomly drop one
                sgn=np.triu(np.random.randint(0,2,fij.shape)*2.-1 )
                sgn=sgn+sgn.T
                fij[ np.logical_and(np.logical_and(fij!=0,fij.T!=0),sgn*dist>0) ]=0
                if role =='fluxes':
                    mat=fij
                else:
                    mat=epsilon*fij-fij.T
                #from datatools import hist,scatter
                #scatter(dist[mat!=0],mat[mat!=0])
                #hist(dist[fij!=0])
                #assert np.mean(mat)<0
            elif role=='efficiency':
                mat=np.abs(sampler(**attrs ))
                dist=np.add.outer(-niche,niche)
                mat*=np.exp(-dist*np.log(attrs.get('factor',1) ) )

                #print attrs.get('factor',1)
                #epsilon=epsilon*np.exp(-1./4.*dist)



        if 'nicheWM' in structure:
            #Williams & Martinez (2000) niche model
            if attrs['role']=='interactions':
                fij=attrs['required']['fluxes'].matrix.copy()
                beta=1/(2*attrs.get('connectance',.2) ) -1
                rge=np.random.beta(1,beta,size=niche.shape)  * niche
                rge/=2 #Half range
                center=np.random.uniform(rge,niche,niche.shape )
                dist=np.add.outer(-niche,niche)
                fij[(dist.T>center+rge- niche).T]=0
                fij[(dist.T<center-rge - niche).T]=0
                np.fill_diagonal(fij,0)
                mat=epsilon*fij-fij.T

        if 0 and attrs['role']=='interactions':
            ## PLOT INFLUXES AS FUNCTION OF NICHES
            nmat=np.array([niche for n in range(niche.shape[0])])
            from datatools import scatter,plot,plt
            scatter(nmat.ravel(),nmat.T.ravel(),c=fij.ravel())


        if 'levels' in structure:
            if attrs['role']=='fluxes':
                mat=sampler(**attrs)
                np.fill_diagonal(mat,0)
                dist=np.add.outer(-niche,niche)
                order=attrs.get('ordered',1)
                canteat =np.logical_or((dist>=-attrs.get('mindist',0)  ), (dist<-1.00001) )
                mat[canteat]*=np.random.choice((0,1),p=(order,1-order),size=np.sum(canteat) )

            if attrs['role']=='interactions':
                fij=attrs['required']['fluxes'].matrix
                #dist=np.add.outer(-niche,niche)
                mat=epsilon*fij-fij.T
            if attrs['role']=='growth':
                mat=np.abs(sampler(**attrs ))
                order=attrs.get('ordered',1)
                nonzero=(niche>0)
                replace=np.random.choice((-attrs['mortality']/np.mean(mat),1),p=(order,1-order),
                    size=np.sum(nonzero) )
                mat[nonzero] *= replace

                #OLD: STRICT REPLACEMENT BY MORTALITY
                #replace=np.random.choice((-attrs['mortality'],np.nan),p=(order,1-order),
                    #size=np.sum(nonzero) )
                #replace[np.isnan(replace)]=mat[nonzero][np.isnan(replace)]
                #mat[nonzero] = replace

        if structure=='princeton':
            T0=273.15
            T=attrs['required']['temperature'].matrix[0]
            niche=attrs['required']['niche'].matrix.copy()
            sdist=10.**niche
            #print sdist

            L=(sdist / .025)** (1./3.3) / 100. # length of individual in m
            spd=(1.*sdist**0.13)/100.*60.*60.*24 # swim speed from Megrey (mday-1)
            ga= 0.60 #fraction of day spent hunting


            from scipy.special import erf

            def const(shape):
                #Normalization constant for thermal envelope
                xs=np.linspace(-10,10,1000)
                return np.max( np.exp(-xs.reshape(xs.shape+(1,))**2) *(1+erf(np.multiply.outer(xs,shape) )),axis=0  )

            def v(V,z,sig=5,shape=-2.7):
                #Effective search rate incl. thermal envelope
                res= np.exp(-(T-z)**2/sig**2) *(1+erf(shape*(T-z)/sig ) )
                return res*V/const(shape)

            def M(s,e=.71,
                    #c=0,
                    c=18.47,
                    E=0.63,k=0.0000862):
                #Basal metabolic rate (METE units)
                return  np.exp(c+e*np.log(s) - E/(k*(T+T0) )   )

            def m(M,s,energy=7000.):
                #Properly scaled metabolic rate
                return M/s * 86400./energy *np.exp(0.03*spd * 100./86400 )

            K= 0.946*5
            r=0.0015*5
            basal=(niche<=np.percentile(niche,10))
            niche[basal]-=1
            #basal=(np.ones(niche.shape)==1)

            if attrs['role']=='influx':
                mat=np.zeros(sdist.shape )
                #mat[basal]=r

            if attrs['role']=='growth':

                mat=-m(M(sdist),sdist)
                #mat[basal]=-r/K
                mat[basal]=r

            if attrs['role']=='threshold':
                LogP = (0.761*np.log10(sdist/1000)) + 10.85 - np.log10(np.exp(0.63/(.0000862*285.65))) # P is kg per ind per yr
                pmax = (10.**LogP) / sdist * 1000 / 365 #g per g per day
                #print attrs['required']['efficiency']
                mat=np.mean(attrs['required']['efficiency'].matrix,axis=1)/(m(
                    M(sdist),sdist)  + pmax)#/100.

            if attrs['role']=='diagonal':
                mat=attrs['required']['search'].matrix#*0
                mat[basal]=r/K

            if attrs['role']=='search':
                mat= ga*np.pi*(L**2)*spd/sdist

            if attrs['role']=='preference':
                dist=np.add.outer(niche,-niche)
                optdist=attrs['mean']
                vardist=attrs['std']**2
                mat= np.exp(-(dist-optdist)**2/2./vardist )

            if attrs['role']=='fluxes':
                Vdist,zdist,sig,shape,phi=[attrs['required'][x].matrix for x in (
                    'search','thopt','thwidth','thskew','preference')]
                vs=v(Vdist,zdist,sig,shape)
                mat=(phi.T*vs).T
                mat[basal]=0

        # ================================ OTHER ECOLOGICAL MODELS ===================================

        if structure =='colonization':
            #Colonization-competition tradeoff
            order=attrs.get('ordered',1)
            if role!='colonization':
                colon=attrs['required']['colonization'].matrix
            if role =='colonization':
                mat=sampler(**attrs)
                ordered=np.sort(mat)
                if order<1:
                    xs=(np.random.random(mat.shape)<order)
                    mat[xs]=ordered[xs]
                else:
                    mat=ordered
            if role =='growth':
                mat=colon- attrs['required']['mortality'].matrix
            if role =='diagonal':
                mat=colon.copy()
            if role =='interactions':
                S=shape[-1]
                effS=S*connectivity
                mu=attrs.get('mean',1)
                sigma=attrs.get('std',np.sqrt(attrs.get('var',1)))

                fij=-(np.add.outer(colon,colon)*mu +sampler(mean=0
                    ,std=sigma/np.sqrt(effS),connectivity=connectivity))
                #fij=(fij+fij.T)/2.
                mat=np.zeros(fij.shape)
                if order<1:
                    idxs=[ i
                     for i in zip(*[z.ravel() for z in tuple(np.indices(fij.shape ))])
                     if i[0]<i[1]
                    ]

                    idxs=[i if np.random.random()<(1.+order)/2. else (i[1],i[0])+i[2:]
                        for i in idxs]
                    idxs=zip(*idxs)

                    mat[idxs]=fij[idxs]

                else:
                    mat=np.triu(fij)
                    np.fill_diagonal(mat,0)

        if structure=='mixture':
            #Special structure for mixture model test
            req=attrs['required']
            niche=req['niche'].matrix
            def req_or_attr(x):
                if attrs.get('attrs_before_req',0):
                    res = attrs.get(x, req.get(x, None))
                else:
                    res= req.get(x,attrs.get(x,None))
                if hasattr(res,'matrix'):
                    res=res.matrix.copy()
                return res
            if role=='growth':
                mat=np.zeros(attrs['shape'])
                rmean=req_or_attr('rmean')
                rstd=req_or_attr('rstd')
                rstdrel=req_or_attr('rstdrel')
                if not rstdrel is None:
                    rstd=rstdrel*rmean
                if not rmean is None:
                    grps= sorted(set(niche))
                    members=[ np.where(niche==i)[0] for i in grps ]
                    for i1 in range(len(grps)):
                        g1= members[i1]
                        locattrs={}
                        locattrs.update(attrs)
                        locattrs['shape']=(len(g1),)+mat.shape[1:]
                        locattrs['mean']=rmean[i1]
                        if not rstd  is None:
                            locattrs['std']=rstd[i1]
                            if 'stdrel' in locattrs:
                                del locattrs['stdrel']
                        mat[g1]=sampler(**locattrs)
                else:
                    mat[niche>0]=0
            if role=='interactions':
                if not req_or_attr('Amean')  is None:
                    mat=np.zeros(attrs['shape'])
                    grps= sorted(set(niche))
                    members=[ np.where(niche==i)[0] for i in grps ]
                    Amean=req_or_attr('Amean')
                    if hasattr(edges,'shape'):
                        nedges=[ [ np.mean([np.sum(edges[m,m2]) for m in m1]) for m2 in members] for m1 in members]
                    else:
                        nedges=[ [ len(m2) for m2 in members] for m1 in members]
                    nedges=np.array(nedges)
                    for n1,n2 in zip(*np.where(nedges==0)):
                        edges[np.random.choice(members[n1]),np.random.choice(members[n2]) ]=1
                        nedges[n1,n2]=1

                    scaling=attrs.get('scaling',None)
                    if scaling == 'tot':
                        Amean/=nedges
                    if not req_or_attr('Astdrel') is None:
                        Astd=Amean*req_or_attr('Astdrel')
                    else:
                        Astd=req_or_attr('Astd')
                    Asym=req_or_attr("Asym")
                    for i1 in range(len(grps)):
                        g1= members[i1]
                        for i2 in range(i1+1)[::-1]:
                            g2=members[i2]
                            handled=0
                            tries =0
                            prev=None
                            while not handled and tries<1:
                                locmat=BaseSampler().sample(mean=0.,std=1.,shape=(len(g1),len(g2)) )
                                if i2==i1:
                                     locmat=klass.symmetrize(locmat,{'symmetry':Asym[i1,i2]} )[:len(g1),:len(g2)]
                                     locmatT=locmat.T
                                     np.fill_diagonal(locmat,0)
                                else:
                                    c = np.sqrt(1 - Asym[i1, i2] ** 2)
                                    locmatT= c*  BaseSampler().sample(mean=0.,std=1.,shape=(len(g2),len(g1)) ) + Asym[i1,i2]*locmat.T

                                mat[np.ix_(g1,g2)]= Amean[i1,i2]+locmat*Astd[i1,i2]
                                if i2!=i1:
                                    mat[np.ix_(g2,g1)]= Amean[i2,i1]+locmatT*Astd[i2,i1]

                            # if 0:
                                #TODO: REMOVE THIS!!!!!!!!!!!!!!!!!!!!!!!!!
                                ########
                                # from datatools import scatter,plot,plt
                                from itertools import product as iprod
                                Aij=mat
                                mwhere=[g1,g2 ]
                                interactions = [[[Aij[i, j] for i, j in iprod(m1, m2) if i!=j] for m2 in mwhere] for m1 in
                                                mwhere]
                                tAmean = np.array([[np.mean(i) if len(i) > 0 else 0 for i in j] for j in interactions])
                                value=np.max( np.abs(tAmean-Amean[ np.ix_([i1,i2],[i1,i2] )] ))
                                if prev is None or value < prev[1]:
                                    prev=mat[np.ix_(g1, g2)],value

                                if value <10**-3:
                                    handled=1
                                tries+=1
                            if not handled:
                                # print '   Mixture group gen: Using best choice',i1,i2,prev[1]
                                mat[np.ix_(g1, g2)]=prev[0]

                            if 0:
                                tAstd = np.array([[np.std(i) if len(i) > 1 else 0 for i in j] for j in interactions])

                                intsym = [[[Aij[i, j] * Aij[j, i] for i, j in iprod(m1, m2)] for m2 in mwhere] for m1 in
                                          mwhere]
                                tAsym = np.array([[np.mean(i) if len(i) > 0 else 0 for i in j] for j in intsym])


                                gamma = (tAsym - tAmean * tAmean.T) / (tAstd * tAstd.T + 10 ** -15)

                                print 'Asym'
                                print Asym[i1,i1],Asym[i1,i2]
                                print Asym[i2,i1],Asym[i2,i2]
                                print gamma,'\n\n'
                                print 'Amean'
                                print Amean[ np.ix_([i1,i2],[i1,i2] )]
                                print tAmean
                                print '\n\n'
                                print 'Astd'
                                print Astd[ np.ix_([i1,i2],[i1,i2] )]
                                print tAstd
                                print '\n\n'

                    if attrs.get('ordered',1) <1:
                        rewire=np.where(np.random.random(mat.shape) > attrs['ordered'] )
                        target=np.array(zip(*rewire))
                        np.random.shuffle(target)
                        mat[rewire]=mat[tuple(target.T)]
                        # xx = np.concatenate(members)
                        # from datatools import plt
                        # plt.imshow(mat[np.ix_(xx, xx)])
                        # plt.show()


                    # print (req_or_attr('Amean').T /get_attr('Amean_rate',1)**np.arange(len(members)) ).T
                    # print (np.array([[np.mean(np.sum((mat*edges)[np.ix_(m1, m2)],axis=1))  for m2 in members] for m1 in members]).T /get_attr('Amean_rate',1)**np.arange(len(members)) ).T

                else:
                    mat=sampler(**attrs)
                    dist=np.add.outer(niche,-niche)
                    order=(np.random.random(mat.shape) < attrs.get('ordered',1) )
                    fakesame=(np.random.random(mat.shape)< 1/(1+np.sum(dist!=0)*1./np.sum(dist==0)) )

                    from numpy import logical_and as AND, logical_or as OR, logical_not as NOT

                    same=OR(AND((dist==0), order), AND(fakesame, NOT(order)))
                    mat[same]=-np.abs(mat[same])*3
                    diff=OR(AND((dist!=0),order), AND(NOT(fakesame),NOT(order) ))
                    mat[diff]=np.abs(mat[diff])

        diag = attrs.get('diagonal', None)
        if not diag is None:
            if diag == 'laplacian':
                np.fill_diagonal(mat,0)
                np.fill_diagonal(mat,-np.sum(mat,axis=0))
            else:
                if hasattr(edges,'shape'):
                    np.fill_diagonal(edges,1)
                if hasattr(diag, 'keys'):
                    mat[np.eye(mat.shape[0]).astype('bool')] = BaseSampler(shape=(mat.shape[0],), **diag).sample()
                else:
                    if hasattr(diag, '__iter__'):
                        mat[np.eye(mat.shape[0]).astype('bool')] = np.diag(diag)
                    else:
                        np.fill_diagonal(mat, diag)

                # np.fill_diagonal(mat,attrs['diagonal'])
            #print mat

        if not mode=='network' and not structure in ('bunin','niche','nicheWM','mixture','ranked'):
            mat=klass.symmetrize(mat,attrs)

        mat=mat*edges
        if 'rowmax' in attrs:
            mat=mat/np.max(np.sum(mat,axis=1) )*attrs['rowmax']
        if 'colmax' in attrs:
            mat=mat/np.max(np.sum(mat,axis=0) )*attrs['colmax']

        return klass( matrix=mat,**attrs)

    def flux(self,vec):
        '''Biomass influx into each species given abundance vector vec'''
        import numpy as np
        #print vec.shape, self.matrix.shape
        mat=np.clip(self.matrix,0,None)
        mat=np.dot(mat,vec)*vec
        mat[mat<0]=0
        return mat


class CommNetwork(object):
    netprops=['kmean'] #Network properties
    edge_data=False #Data stored on edges (beyond weight)?

    def __init__(self,array=None,**kwargs):
        self.props=kwargs
        self.matrix=array
        self.refresh_graph()

    def __getitem__(self,*args,**kwargs):
        return self.matrix.__getitem__(*args,**kwargs)

    def __setitem__(self,*args,**kwargs):
        val= self.matrix.__setitem__(*args,**kwargs)
        self.add_edge(args[0][0],args[0][1],weight=args[1])
        return val

    def add_edge(self,i,j,**kwargs):
        return self.graph.add_edge(i,j,**kwargs)

    def apply(self,vec):
        '''Instantaneous impact on abundance vector vec'''
        import numpy as np
        return np.dot(self.matrix,vec)

    def flux(self,vec):
        '''Biomass influx into each species given abundance vector vec'''
        import numpy as np
        #print vec.shape, self.matrix.shape
        mat=np.clip(self.matrix,0,None)
        mat=np.dot(mat,vec)*vec
        mat[mat<0]=0
        return mat


    def refresh_graph(self):
        import networkx as nx
        self.graph=nx.DiGraph(self.matrix)

    @property
    def shape(self):#,*args,**kwargs):
        return self.matrix.shape#*args,**kwargs)

    def copy(self):
        return self.__class__(self.matrix.copy(),**deepcopy(self.props))

    def label(self):
        mat=self.matrix.copy()
        frzero=np.sum(mat[mat==0])/float(mat.size)
        np.fill_diagonal(mat,0)
        diag=np.diag(self.matrix)
        ratio=np.log10(np.median(np.abs(diag[diag!=0]) )/np.median(np.abs(mat[mat!=0])))
        if ratio <-1:
            st='weak'
        elif ratio <1:
            st='med'
        else:
            st='str'
        if frzero<=.3:
            sp='sparse'
        else:
            sp='dense'

        lbl='_{}int_'.format(st,sp)
        if 'genre' in self.props:
            return self.props['genre']+lbl
        else:
            return lbl


    def save(self,dirpath,idx=None,fname='network'):
        '''If idx is a sequence, save only columns and rows involving
        species in idx'''
        import pandas as pd
        from itertools import product as iprod

        if self.edge_data is False:
            #Matrix
            index,mat=mat_to_cols(self.matrix,drop_zeros=1,idx=idx)
            df=pd.DataFrame([{'weight':j} for j in mat ],index=pd.MultiIndex.from_tuples(index)  )
            df.to_csv(dirpath+fname+'.mat', sep='\t')
            #np.save(dirpath+'network',self.matrix)
        else:
            self.edge_data.to_csv(dirpath+fname+'.mat', sep='\t',index_label=['',''])

        #Properties
        f=open(dirpath+fname+'.dat','w')
        self.props['size']=self.matrix.shape[0]
        if idx is None:
            idx='all'
        self.props['saved']=idx
        for p in sorted(self.props):
            j=self.props[p]
            if isinstance(j,basestring):
                j=repr(j)
            f.write('{}:{}\n'.format(p,j ))

    def from_df(self,df):
        size=self.props.get('size',1)
        matrix=np.zeros((size,size) )
        if not df is None:
            for idx,series in df.iterrows():
                row=dict(series)
                if not 'weight' in row:
                    row['weight']=row.pop(df.columns[0] )
                matrix[idx]=row['weight']
                self.add_edge(idx[0],idx[1],**row)
        self.matrix=matrix

    def load(self,dirpath,fname='network'):
        import pandas as pd

        #Properties
        self.props.clear()
        f=open(dirpath+fname+'.dat','r')
        for l in f:
            i,j=l.strip().split(':')
            self.props[i]=eval(j)

        #Matrix
        if self.edge_data is False:
            fil=pd.read_csv(dirpath+fname+'.mat',sep='\t',index_col=[0,1])
            #self.matrix=np.load(dirpath+'network.npy')
        else:
            fil=pd.read_csv(dirpath+fname+'.mat',sep='\t',
                header=0,index_col=[0,1])
            self.edge_data=fil

        self.from_df(fil)


        return self

    def get_netprop(self,prop,idx=None):
        '''Compute network property (only over indices idx if specified)'''
        import networkx as nx
        #mat=self.matrix.copy()
        #np.fill_diagonal(mat,0)
        #if not idx is None:
            #mat=mat.ix_(idx,idx)
            #graph=nx.DiGraph(mat)
        #else:
            #graph=self.graph
        if idx is None:
            g=self.graph
        else:
            g=self.graph.subgraph(idx)
        if prop == 'kmean':
            return np.mean([g.degree(x) for x in g.nodes()])
        if prop =='betweenness':
            return np.mean(nx.betweenness_centrality(g).values())
        if prop =='connectivity':
            return nx.average_node_connectivity(g)

    def store_netprop(self,prop=None):
        if prop is None:
            return [self.store_netprop(p) for p in self.netprops]
        self.props[prop]=self.get_netprop(prop)


