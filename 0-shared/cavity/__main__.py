# -*- coding: utf-8 -*-
from __init__ import *




if 1:
    #test()
    sigma=2.8
    gamma=-0.5
    sigma_k=1.2
    mu=.1
    size=200

    if 0:
        xs=np.linspace(-50,50,100)
        dist1=cavity_distri(S=size,sigma=sigma,sigma_k=sigma_k,mu=mu,gamma=gamma)
        plot(xs,[dist1(x) for x in xs],hold=1)

        pops=[p for x in range(10)
            for p in assembly_iter_cavity(size=size,sigma=sigma,sigma_k=sigma_k,mu=mu,gamma=gamma)]
        hist(pops,normed=1,log='y')
    #phi_vs_sigma(sigma=sigma,gamma=gamma,sigma_k=sigma_k,mu=mu,do_scale=1,do_plot=1)
    #phi_vs_sigma(hold=1,sigma=sigma,gamma=gamma,sigma_k=sigma_k,mu=mu)


    #FIRST RESULT FIGURES, (sigma,mu) space
    resolution=128
    size=1000
    article_set=['phi','variability','totN','simpsoninv','prod','chi2']
    for cols in (article_set, ):
        for gamma in (1,0,-1):
            mat=cavity_show_loop(title='',cols=cols,sigma=np.linspace(.2,1.5,resolution),
                mu=np.linspace(-2,4,resolution),sigma_k=0.,gamma=gamma,size=size,NTRIALS=3 )

    plt.show()
    mat=cavity_show_loop(sigma=np.linspace(.2,1.5,resolution),sigma_k=np.linspace(.2,4.5,resolution),mu=0,gamma=0,size=1000 )
    mat=cavity_show_loop(sigma=np.linspace(.2,2.5,resolution),sigma_k=.5,mu=1.5,gamma=np.linspace(-.99,.99,resolution),size=100 )
    mat=cavity_show_loop(sigma=np.linspace(.2,2.5,resolution),sigma_k=np.linspace(.2,1.2,resolution),mu=1.5,gamma=-.99,size=100 )
    mat=cavity_show_loop(sigma=np.linspace(.2,1.5,resolution),sigma_k=np.linspace(.2,1.2,resolution),mu=1.5,gamma=.99,size=100 )
    mat=cavity_show_loop(sigma=np.linspace(.2,2.5,resolution),mu=np.linspace(-1,8,resolution),sigma_k=5,gamma=-.99,size=100 )
    mat=cavity_show_loop(sigma=np.linspace(.2,1.5,resolution),mu=np.linspace(-1,12,resolution),sigma_k=5,gamma=.99,size=100 )
    mat=cavity_show_loop(sigma=np.linspace(.2,1.5,resolution),mu=np.linspace(-1,12,resolution),sigma_k=.1,gamma=.99,size=100 )
    mat=cavity_show_loop(sigma=.5,mu=np.linspace(-1,4,resolution),sigma_k=1,gamma=np.linspace(-.99,.99,resolution),size=100 )

    plt.show()