# -*- coding: utf-8 -*-

#=================== FAST ASSEMBLY =============================


def simple_assembly_iter_cavity(size=25,mu=.0,sigma=1.,sigma_k=.1,gamma=-.5,
        q=None,v=None,h=None,phi=None,do_print =1,normalize=1):
    #ITERATIVELY ASSEMBLE COMMUNITY
    if None in (q,v,h,phi):
        q,v,h,phi=cavity_solve(size,mu,sigma,sigma_k,gamma)
    #A,r,x=gen_rdm(size,mu,sigma,sigma_k,gamma,eq=False)


    aij=(np.random.normal(0,1./np.sqrt(size),(size,size)) )
    aij+=.5*gamma*aij.T
    r=np.abs(np.random.normal(1,sigma_k,size))
    x=np.ones(size)

    u=(1-mu/size)/sigma
    avgN=1./(sigma*h + mu)
    #a=(A+mu/size)/sigma
    a=aij
    lam= (r-1)/(sigma * avgN)

    #print 'while here', avgN,u,np.var(a),np.var(lam), u-gamma*v
    pops=x/np.mean(x)#/avgN
    t=0
    tmax=1000*size
    prev=-1

    while t < tmax:
        i = np.random.randint(size)
        n0=(lam[i] + h + np.dot(a[i],pops) ) /(u-gamma *v)

        pops[i]=max(0,n0)
        if normalize:
            pops/=np.mean(pops)


        assert np.mean(pops)>.9
        t+=1

        if t%(tmax/100)==0:
            if do_print:
                print np.sum(pops==0)/float(size)
                #print np.mean(pops**2),q,np.mean(pops)
            if np.var(pops)-prev <0.0001:
                if do_print:
                    print 'Breaking at', t
                break
            prev= np.var(pops)
    return pops*avgN


def assembly_iter_cavity(size=25,mu=.0,sigma=1.,sigma_k=.1,gamma=-.5,
        q=None,v=None,h=None,phi=None,do_print =1,normalize=1, nstates=12):
    #ITERATIVELY ASSEMBLE COMMUNITY, SPIN GLASS STYLE

    if None in (q,v,h,phi):
        q,v,h,phi=cavity_solve(size,mu,sigma,sigma_k,gamma)


    aij=(np.random.normal(0,1./np.sqrt(size),(size,size)) )
    aij+=.5*gamma*aij.T
    r=np.abs(np.random.normal(1,sigma_k,size))
    x=np.random.random((size,nstates))

    u=(1-mu/size)/sigma
    utilde=u-gamma*v
    avgN=1./(sigma*h + mu)
    #a=(A+mu/size)/sigma
    a=aij
    lam= (r-1)/(sigma * avgN)
    sigma_lambda=sigma_k/sigma/avgN

    #print 'while here', avgN,u,np.var(a),np.var(lam), u-gamma*v
    pops=x/np.mean(x,axis=0)#/avgN
    t=0
    tmax=1000*size
    prev=-1

    beta=10.

    def normal(mn,vr):
        return np.exp(- mn**2/ 2 / vr) / np.sqrt(2*np.pi*vr )

    def calcF(n,ns):
        mn= n *utilde - h
        vr=np.mean(ns**2) + sigma_lambda**2
        return np.log(u * normal(mn,vr) / ((1-phi) * len(ns) ) )

    truepops=pops.copy()
    while t < tmax:
        i = np.random.randint(size)

        DeltaFs=[]
        temp=[]
        for s in range(nstates):
            ns=truepops[:,s]
            n0=(lam[i] + h + np.dot(a[i],ns) ) /(u-gamma *v)
            DeltaF=calcF(n0,ns)-calcF(pops[i,s],ns)
            DeltaFs.append(DeltaF)
            temp.append(n0)
        probas=np.exp(-beta * np.array(DeltaFs) )
        probas/=np.sum(probas)

        for s in range(nstates):
            pops[i,s]=np.random.choice(temp,p=probas)
            truepops[i,s]=max(0,pops[i,s])
        if normalize:
            mnpop=np.mean(truepops,axis=0).reshape((1,nstates) )
            pops/=mnpop
            truepops/=mnpop


        t+=1

        if t%(tmax/100)==0:
            if do_print:
                print np.sum(truepops==0,axis=0)/float(size)
                #print np.mean(pops**2),q,np.mean(pops)
            if np.var(truepops)-prev <0.0001:
                if do_print:
                    print 'Breaking at', t
                break
            prev= np.var(truepops)
    return truepops.ravel()*avgN

