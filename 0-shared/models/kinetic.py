# -*- coding: utf-8 -*-

from lotvol import *
from graphs import DistriNetwork


class KineticModel(DynamicModel):
    '''Model for Lotka-Volterra-like kinetics with following options:
        - time-dependent particle-wise interaction strength w/ other species

    '''
    Network=DistriNetwork

    def __str__(self):
        return "KineticModel"

    class Particle(object):
        '''A particle is a tuple of values + a creation time'''
        def __init__(self,seed,t,species):
            self.seed=tuple(seed)
            self.t0=t
            self.species=species

        def __getitem__(self,i):
            return self.seed[i]

        def __repr__(self):
            return str(id(self))
                #'Particle(species {}, {}, t0={})'.format(self.species,self.seed,self.t0)


    def get_transitions(self,particles,abundance,t,oldW=None):
        if oldW is None:
            oldW={}
        tgt={}
        W={}
        rate={}
        nei,edges={},{}
        maxrate=0
        A=self.interactions

        for p in particles:
            s=p.species
            if not s in nei:
                edg= [(i,j) for i,j in zip(*A.edges(s)) if abundance[i]]
                nei[s],edges[s]=(tuple(x) for x in zip(*edg))
            tgt[p]=nei[s]

            if p in oldW and len(oldW[p])==len(tgt[p]):
                W[p]=oldW[p]
            else:
                #if p in oldW:
                    #print 'CHANGE',len(oldW[p]),len(tgt[p])
                W[p]=tuple(d.sample(val,t=t-p.t0)#*abundance[n]/omega
                    for n,d,val in zip(nei[s],edges[s],p.seed) )

            #print W[p]
        #print s
        #print p, W[p],tgt[p]
        return tgt, W

    def evol(self,x0,tmax,Mavg=None,tsample=1.0,force_match=False,print_msg=0,
            keep='all',store=False,ftol=0.01,**kwargs):
        '''
        Mavg: typical number of particles (M)
        force_match: only perform reaction on ij edge if suitable partner
        (particle reacting on ji edge)'''

        A=self.interactions
        if Mavg is None:
            #DEFAULT NUMBER OF PARTICLES PER SPECIES
            Mavg=300*x0.shape[-1]
        Particle=self.Particle
        def create(t,spec):
            return Particle(A.generate_seed(spec),t,spec)

        if hasattr(x0[0],'__iter__'):
            #x0 given as dictionary of particles, do not create them
            nspec=len(x0)
            init={s: x0[s] for s in xrange(nspec)}
            omega=Mavg/np.sum(x0)
        else:
            #x0 given as abundances, create particles
            nspec=len(x0)
            sumx0=np.sum(x0)
            init={s:[create(0,s) for n in xrange(int(Mavg*x0[s]/sumx0)) ]
                for s in xrange(nspec) }
            omega=Mavg/sumx0
        #Omega is a "timescale" that is important due to nonlinearity (abundances in transition rates)

        species=init
        def get_x(species):
            return np.array([len(species.get(s,[])) for s in xrange(nspec) ])

        r=self.growth

        t=0.0
        x0=get_x(species)
        x=x0.copy()
        traj=Trajectory(self,x0 )

        if tmax is None:
            tmax =1000#10**6
            converge=1
        else:
            converge=kwargs.get('converge',0)

        lastchange=1
        ftol=max(ftol,1./Mavg )

        #ARE INTERACTIONS TIME DEPENDENT?
        W={}
        if  A.props.get('freq',None) is None and A.props.get('vanish',None) is None:
            time_dependence=False
            orig=[part for s in species for part in species[s]]
            tgt,W=self.get_transitions(orig,x,t)
            #Reuse same W wherever possible
        else:
            time_dependence=True

        created=[]
        while t<tmax and not (t>0 and converge and lastchange < ftol) :


            particles=[part for s in species for part in species[s]]

            if time_dependence:
                #W has to change every timestep
                tgt,W=self.get_transitions(particles,x,t)
            else:
                #Only update for newly created particles
                #tgt,W=self.get_transitions(particles,x,t)
                ptgt,pW=self.get_transitions(created,x,t)
                tgt.update(ptgt)
                W.update(pW)
                created[:]=[]

            dt=tsample

            #Interactions
            if W:
                reactions={} #Particle pools indexed by ij-edge b/w species
                for part in particles:
                    xloc=x[list(tgt[part])]/omega
                    if not xloc.any():
                        continue
                    Wloc=np.abs(W[part])*xloc
                    rate=sum(Wloc)
                    if np.random.random()>rate*dt:
                        continue
                    #print 'Yooo', sum(np.abs(W[part])),rate[part],xloc,np.abs(W[part])*xloc
                    nxt=np.random.choice(tgt[part],p=Wloc/rate )
                    reactions.setdefault((part.species,nxt),[]).append(part)
                for edge in reactions:
                    particles=reactions[edge]
                    if force_match:
                        partners=reactions.get(edge[::-1],[])[:]
                        particles=particles[:len(partners)]
                    spec=edge[0]
                    for part in particles:
                        if W[part][tgt[part].index(edge[1])]<0:
                            #destroy particle
                            species[spec].remove(part)
                            if not species[spec]:
                                del species[spec]
                            if not time_dependence :
                                del W[part]
                                del tgt[part]
                        else:
                            #create particle of same species
                            newp=create(t,spec)
                            species[spec].append(newp)
                            created.append(newp)


            #Growth
            for spec in xrange(nspec):
                species.setdefault(spec,[])
                rate= len(species[spec])* r[spec]# +1
                rate*=dt
                rest=rate%1.0
                rate=np.floor(rate)
                if np.random.random()<rest:
                    rate+=1.
                for n in xrange(int(rate) ):
                    newp=create(t,spec)
                    species[spec].append(newp)
                    created.append(newp)

            t+=dt
            if t %tsample <= dt:
                N=sum( len(x) for x in species.values() )
                if print_msg:
                    print 't=',t, 'N=',N
                if keep=='all':
                    traj.append(get_x(species) ,index=float(t))
                    #dX= traj.matrix[-1]-traj.matrix[-2]
                    #print 'A:',dX/dt
                    #print 'B:',dX/dt/traj.matrix[-2]#/traj.matrix[-2,::-1]*omega

            x,xold=get_x(species),x
            lastchange=np.abs(1.-x.astype('float')[xold>0]/xold[xold>0])
            if len(lastchange):
                lastchange=np.max(lastchange)/tsample


        #END LOOP =======================

        traj.append(get_x(species) ,index=float(t))
        if store:
            self.trajectories.append(traj)
        return traj
