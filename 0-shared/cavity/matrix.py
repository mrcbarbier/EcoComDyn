# -*- coding: utf-8 -*-
#=================== GENERATE MATRIX ===============================

def gen_pool_matrix(size=25,mu=.0,sigma=1.,gamma=-.5):
    """generate random (pool) matrix"""

    if gamma==0:
        asymm=0
    else:
        asymm=(1- np.sqrt(1-gamma**2) )/gamma

    aij=(np.random.normal(0,1./np.sqrt(size),(size,size)) )

    aij= (aij+asymm * aij.T) / np.sqrt(1+ asymm**2 )

    A=-mu/size + sigma*aij
    np.fill_diagonal(A,-1)
    return A



def gen_pool_matrix_trophic(size=25,mu=.0,sigma=1.,gamma=-.5):
    """generate random (pool) trophic matrix (i.e. with conservation) """

    S=size
    C1=2.*mu/S/(1+gamma)
    C2= ( 4.*sigma/S + C1**2*(1+gamma)**2 )/(1+gamma**2)/2

    aij=(np.random.normal(C1,np.sqrt(C2-C1**2),(size,size)) )

    assert not np.isnan(aij).any()


    if 1:
        np.fill_diagonal(aij,0)
        test=np.tril( np.random.randint(0,2,aij.shape) ) #empty one triangle
        aij=np.tril(aij)
        #xs=np.where(test>0)
        #aij[xs]= gamma*aij.T[xs]
        aij+=gamma*aij.T
        print np.mean(aij*aij.T)/( gamma*C2 )
        #print aij[:5,:5], np.mean(aij),np.mean(aij*aij.T)#,  xs
        #aij[aij>np.abs(aij.T)]=gamma*aij.T[aij>np.abs(aij.T)]
        print aij[:5,:5]
        print np.mean(aij)*S,mu,np.var(aij)*S,sigma,np.mean(aij*aij.T)*S,np.mean( ((aij-np.mean(aij)) *(aij.T - np.mean(aij))) ) * size


        #print np.mean(aij*aij.T),np.mean( ((aij-np.mean(aij)) *(aij.T - np.mean(aij))) ),np.var(aij)
        print xs

    A=  sigma*aij/np.std(aij)/np.sqrt(size)
    A+=-mu/(size-1.)-np.mean(A)
    np.fill_diagonal(A,-1)
    from itertools import product
    xs=[(i,j) for i,j in product(range(A.shape[0]),range(A.shape[0]) ) if i!=j]
    xs=[list(z) for z in np.array(xs).T]
    B=A[xs]
    print 'sigma', np.std(B), sigma/np.sqrt(size), 'mu',- np.mean(B),mu/size, 'gamma', np.mean( ((A-np.mean(B)) *(A.T - np.mean(B)))[xs] )/np.var(B),gamma
    return A



def gen_x_cavity(size,q,v,h,phi,mu,sigma,sigma_k,gamma):
    avgN=1./(sigma*h + mu)
    u=(1-mu/size)/sigma
    sigma_lambda=sigma_k/sigma/avgN
    x=[]
    while len(x)<size:
        x0=np.random.normal(0,1)*np.sqrt(q+sigma_lambda)
        x0+=h
        x0/= u-gamma *v
        if x0>0:
            x.append(x0*avgN)
    return np.array(x)

def muAB(q,h,mu,sigma,gamma):
    avgN=1./(sigma*h + mu)
    varN=q* avgN**2

    muB=(1.-mu*avgN)/varN
    muA=1.+ gamma *muB*avgN
    muA/=varN
    return avgN,varN,muA,muB

def gen_cavity(size=25,mu=.0,sigma=1.,sigma_k=.1,gamma=-.5):
    """Directly generate stable matrix using Bunin correlations"""


    q,v,h,phi=cavity_solve(size,mu,sigma,sigma_k,gamma)
    x=gen_x_cavity(size,q,v,h,phi,mu,sigma,sigma_k,gamma)
    #hist(x)
    avgN,varN,muA,muB=muAB(q,h,mu,sigma,gamma)
    print np.mean(x),avgN/phi
    sigma_lambda=sigma_k/sigma/avgN
    r=np.abs(np.random.normal(1,sigma_k,size))

    cov=np.zeros((size**2+size,size**2+size)) #covariance of a_ij, a_kl and r_i

    def ind4(i,j):
        return i//size, i%size, j//size,j%size

    #First, set variances
    it=np.nditer(cov,flags=['multi_index'])
    for v in it:
        i,j=it.multi_index
        if j>i:
            continue
        if i==j:
            if i<size**2:
                #Variance of a_ij
                ii,jj,k,l=ind4(i,j)
                #print ii,jj,k,l
                cov[i,i]=1./size
            else:
                #variance of r_i
                cov[i,i]=sigma_lambda**2


    it=np.nditer(cov,flags=['multi_index'])
    for v in it:
        i,j=it.multi_index
        if i>=j:
            continue
        if j<size**2:
            ii,jj,k,l=ind4(i,j)
            #Covariance of a_ij and a_kl
            varprod=np.sqrt( cov[i,i]*cov[j,j] ) / avgN**2
            #print varprod, avgN, cov[i,i],cov[j,j]
            #Bunin gives correlations, I need covariances, hence multiply by varprod
            if ii==k:
                cov[i,j]=cov[j,i]=-x[jj]*x[l]/size**2/(q+sigma_lambda**2) *varprod
            elif ii==l:
                cov[i,j]=cov[j,i]=-gamma*x[jj]*x[k]/size**2/(q+sigma_lambda**2) *varprod
            elif jj==l:
                cov[i,j]=cov[j,i]=-gamma**2*x[ii]*x[k]/size**2/(q+sigma_lambda**2) *varprod

        elif i<size**2:
            #Covariance of a_ij and r_i
            ii=i/size
            jj=j-size**2
            cov[i,jj]=x[ii]/size/(1+q/sigma_lambda**2)/ avgN
            cov[jj,i]=gamma*x[jj]/size/(1+q/sigma_lambda**2)/ avgN

    xx=np.vstack([x for i in range(x.size)]).T

    mu_eff=(mu-muA*xx*xx.T + muB *(gamma * xx + xx.T ))
    mu_eff=mu_eff.ravel()
    #hist(mu_eff)

    gen=np.random.multivariate_normal(np.concatenate((mu_eff,np.ones(size))),cov)

    aij=gen[:size**2].reshape((size,size))
    aij+=.5*gamma*aij.T

    r=gen[size**2:]


    A=-mu/size + sigma*aij
    #A[np.random.random(A.shape)>connectance ]=0
    #scatter(x,(np.dot(np.linalg.inv(A),r) ))
    x=np.dot(np.linalg.inv(A),r)
    return A,r,x

def test(nsamples=100,size=15,mu=1.,sigma=1.,sigma_k=.2,gamma=-.5,connectance=.5):
    s=0
    measure={
        'meanJv':[],
        'varJv':[],
        'corAJv':[],
        }
    Jvs=[]
    As=[]
    #print np.var(np.random.normal(0,1./np.sqrt(size),(size,size))), 1./size
    #return
    while s  < (nsamples):
        #test.append((np.mean(aij*aij.T),gamma/size))

        #hist(x)
        #print np.min(x)
        #A,r,x=gen_random(size,mu,sigma,sigma_k,gamma)
        A,r,x=gen_cavity(size,mu,sigma,sigma_k,gamma)
        if (x<0).any():
            print 'fail'
            continue
        print 'Sample',s
        Jhat=lifted_matrix(A,0)
        Jv=variance_interactions(Jhat)
        Jvs.append(Jv)
        As.append(A)
        s+=1
    #plot((-1,1),(-1,1),hold=1)
    #scatter(*zip(*test))
    means=[np.mean(Jv) for Jv in Jvs]
    hist(means,hold=1,normed=1,bins=200,log='y')
    plt.figure()
    vrs=[np.var(Jv) for Jv in Jvs]
    print 'Mean var',np.mean(vrs)
    hist(vrs,hold=1,log='xy',normed=1,bins=200)
    plt.figure()
    cors=[ np.mean(A*Jv)-np.mean(A)*np.mean(Jv) for A, Jv in zip(As,Jvs)]
    print 'Mean cor',np.mean(cors)
    hist(cors,normed=1 ,bins=200,log='y')

