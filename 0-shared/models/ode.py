# -*- coding: utf-8 -*-
from model import *
from modules import *

dftparams={
        'n':{
    'type':'variable',
    'axes': ('com',),
    'death': 10.**-10,
    'shape':(50,),
    'mean':1.,
    'std':1,
    'sign':1,
    },

        'growth':{
    'type':'matrix',
    'variables':['n'],
    'dynamics':'nlog',
    'role':'growth',
    'mean':1., 'std':0,
    'sign':1,
    },
        'capacity':{
    'type':'matrix',
    'variables':[('n','com')],
    'dynamics':'nlog',
    'role':'capacity',
    'mean':1., 'std':0,
    'sign':1,
    },
        'community':{
    'type':'matrix',
    'variables':[('n','com'),('n','com') ],
    'dynamics':'ncom',
    'role':'interactions',
    'mean':.5, 'std':.3,
    'symmetry':1.,
    'structure':'bunin'
    },
        #'dispersal':{
    #'type':'matrix',
    #'variables':[('n','spa'),('n','spa')],
    #'dynamics':'ndisp',
    #'role':'dispersal',
    #'mean':.01, 'std':0.,
    #'diagonal':'laplacian',
    #'sign':1,
    #},

}

dftdynamics={
    'nlv':{'type':'lotvol',
        'variables':[('n','com'),],
        },
    'nlog':{'type':'logistic',
        'variables':[('n','com'),],
        },
    'ncom':{'type':'community',
        'variables':[('n','com'),],
        },
    #'ndisp':{'type':'dispersal',
        #'variables':[('n','spa'),],
        #},

    }

class DynamicModel(BaseModel):
    submodels=('dynamics',) #

    dft_prm=dftparams
    dft_dyn=dftdynamics

    modules={
        'linear':LinearModule,
        'logistic':LogisticModule,
        'community':CommunityModule,
        'lotvol':LotkaVolterraModule,
        'simplelv':SimpleLotkaVolterraModule,
        'simplefr':SimpleFunctionalResponseModule,
        'trophicfr':SimpleFunctionalTrophicModule,
        'macarthur':MacArthurModule,
        'dispersal':DispersalModule,
        'noise':NoiseModule,
        }
    generators={
        'variable':Trajectory,
        'constant':LabeledArray,
        'matrix':BaseNetwork,
        'noise':Noise,
        }

    def __init__(self,
            parameters=None,
            dynamics=None,
            data=None,
            # sampler={},
            **kwargs):

        if parameters is None:
            parameters=deepcopy(dftparams)
        if dynamics is None:
            dynamics=deepcopy(dftdynamics)
        if data is None:
            data={}
        self.data=data
        self.parameters=parameters
        self.set_params(**kwargs) #Update parameters if given in kwargs

        self.dynamics=dynamics

        self.attrs={'parameters':self.parameters,
            'dynamics':self.dynamics}

        if kwargs.get('initialize',True):
            self.initialize()

        if 'results' in kwargs:
            for i,j in kwargs['results'].iteritems():
                if i in self.results and not j is None:
                    self.results[i]=j.copy()

    def set_params(self,**kwargs):
        prm= BaseModel.set_params(self,**kwargs)
        self.update_variables()
        return prm


    def update_variables(self):
        """Derive variables from parameter list."""
        prm=self.parameters
        for v in sorted(prm):
            if not hasattr(prm[v],'keys'):
                raise Exception('Model parameter "{}" badly formatted (probably left a comma at the end).'.format(v) )
        self.variables=OrderedDict([(v,prm[v]) for v in sorted(prm) if prm[v]['type']=='variable' ] )
        if not hasattr(self,'results'):
            self.results={i:None for i in self.variables}
        self.noises=OrderedDict([(v,prm[v]) for v in sorted(prm) if prm[v]['type']=='noise' ] )

    def initialize(self,labels=None):
        """Internal use. Generate model data from parameter."""

        if labels is None:
            prms=self.parameters.items()
        else:
            prms=[(label,self.parameters[label]) for label in labels ]
        done=[]
        tries=0
        while prms:
            label,prm=prms.pop(0)
            tries+=1
            if tries>5000:
                raise Exception("COULD NOT INITIALIZE:{}".format(prms))
            if '_' in label:
                print 'WARNING: Underscore in variable name can cause failure of parameter update.'
            if prm['type']=='variable':
                done.append(label)
                continue
            if 'formula' in prm or 'requires' in prm:
                #One data requires another
                import re
                formula=prm.get('formula','')
                formulacomps=re.findall('([\s\w,-_]*)',formula)
                #print formulacomps
                req=[l for l in self.parameters if (l in formulacomps
                     or l in prm.get('requires',[])) and not l in done]
                if req:
                    #print 'STILL REQ',label,req
                    prms.append( (label,prm) )
                    continue
            done.append(label)
            #print 'DONE',label
            if 'dynamics' in prm:
                if not set( to_iterable(prm['dynamics'])).intersection( set(self.dynamics)):
                    continue
            if label in self.data:
                dat=self.data[label]
            elif 'formula' in prm:
                self.data[label]=dat=self.generate_from_formula(label,prm)
            else:
                self.data[label]=dat=self.generate(label,prm)
            #print label, dat


    @property
    def description(self):
        """Information that allows to reconstruct the model."""
        dic={}
        dic.update(self.attrs)
        dic['klass']=self.__class__.__name__
        return dic

    def get_module(self,dyn):
        """Load correct dynamical module"""
        return self.modules.get(dyn,LogisticModule)(parent=self)

    def get_alive(self,dic,death=-1.):
        al={}
        eff={}
        shape=[]
        for var in dic:
            mat=dic[var].matrix
            alive=tuple(sorted(set(j)) for j in np.where(mat>death))
            shape.append( tuple(len(i) for i in alive) )
            #al.append(np.ix_(*alive))
            al[var]= alive
            eff[var]=LabeledArray(mat[np.ix_(*alive)],axes=self.get_axes(var) )
        return eff,al,shape

    def flatten(self,x,alive=None):
        """Creates a flat array (for purposes of integration) out of the
        alive terms in x"""
        if alive is None:
            return np.concatenate( [x[i].ravel() for i in self.variables ] )
        return np.concatenate( [x[i][np.ix_(*alive[i])] for i in self.variables ] )

    def reshape(self,x,shapes):
        if hasattr(shapes,'keys'):
            shapes=[shapes[v] for v in self.variables]
        sizes=[np.prod(z) for z in shapes]
        cumul=[0]+list(np.cumsum(sizes))
        #print cumul, [cumul[i+1] for i in range(len(shapes))], type(x),type(shapes),shapes
        return [x[cumul[i]:cumul[i+1]].reshape(shapes[i]) for i in range(len(shapes))]


    def update(self,x,xflat,shapes,alive=None):
        #print 'updating',x,xflat
        new=self.reshape(xflat,shapes)
        for i,var in enumerate(self.variables):
            if alive is None:
                x[var].matrix[:]=new[i]
            else:
                #print x[var].shape, alive[var], new[i].shape
                x[var].matrix[np.ix_(*alive[var])]=new[i]

    def get_E(self,x):
        print "TODO: Lyapunov"
        return 0


    def get_dx(self,t,x):
        """External use, for measurements."""
        #xflat=self.flatten(x)
        diff ={var:LabeledArray(np.zeros( self.variables[var]['shape'] ),
                axes= self.get_axes(var)) for var in self.variables}
        for dyn,data in self.make_dynamics().iteritems():
            data['module'].dx(diff,x,**data)
        return diff# self.flatten(diff)

    def make_dynamics(self):
        """Make dicts of dynamical modules and their required properties."""
        dynamics={}
        for dyn, data in self.dynamics.iteritems():
            if data.get('disabled',0):
                continue
            dynamics[dyn]={}
            dynamics[dyn].update(data)
            module= self.get_module(data['type'])
            dynamics[dyn]['module']=module
            dynamics[dyn]['matrix']={}
            dynamics[dyn]['noise']={}
            dynamics[dyn]['roles']={}

        #TODO: be clever and make everything into labeled objects from their creation!!
        for label,prm in self.parameters.iteritems():
            if not 'dynamics' in prm or not set( to_iterable(prm['dynamics'])
                        ).intersection( set(self.dynamics) ):
                continue
            dyn=prm['dynamics']
            role=prm['role']
            data=dynamics[dyn]
            data['roles'][prm['role']]=label
            prm=self.parameters[label]
            axes=[]
            for v in prm['variables']:
                if isinstance(v,basestring):
                    axes+=self.get_axes(v)
                else:
                    axes.append(v)
            typ=prm['type']
            if typ=='matrix':
                matrix=LabeledArray(self.data[label].matrix,axes=axes )
                data['matrix'][role]=matrix
            elif typ=='noise':
                data['noise'][role]=self.data[label]
        return dynamics

    def set_init(self,**kwargs):
        for i,j in kwargs.iteritems():
            self.results[i]=Trajectory(init=j )

    def evol(self,tmax=3.,tsample=None,keep=['end','all'][1],
            print_msg=0,
            remove_dead=0,
            ftol=0.00001,
            dftol=None,
            divergence=10**8,
            **kwargs):

        death=kwargs.get('death',np.min([prm.get('death',1) for prm in self.variables.values() ] ))
        if dftol is None:
            dftol=death*.9

        t=0.
        if tmax is None:
            tmax =5000#10**6
            converge=1
        else:
            converge=kwargs.get('converge',0)
        if self.noises:
            converge=0
        if tsample is None:
            tsample=tmax/20


        self.set_params(**kwargs)
        import scipy.integrate as scint

        #Externally controlled parameters ("scenarios")
        scenarios={}
        for i,j in kwargs.get('scenarios',{}).iteritems():
            if i in self.parameters:
                scenarios[i]=j.copy()

        #Dealing with repeated runs of the model (randomize or keep same conditions)

        x=OrderedDict()
        reseed=kwargs.get('reseed', True )
        extend=kwargs.get('extend',False)
        if extend and (reseed is True):
            reseed=None
        if not reseed:
            reseed=[]
        elif reseed is True:
            reseed=list(self.variables)+list(self.noises)
        for var in self.variables:
            if var in reseed or self.results[var] is None:
                #Normal situation: reseed random factors in the model at every run
                attrs=self.variables[var]
                self.results[var]=self.generate(var,attrs)
            elif extend:
                #Start from end of previous run
                self.results[var]=Trajectory(init=self.results[var][-1] )
            else:
                #Restart from same conditions as previous run
                self.results[var]=Trajectory(init=self.results[var][0] )
            x[var]=LabeledArray(self.results[var][-1].copy(),axes=self.get_axes(var))

        for noise in self.noises:
            if noise in reseed or self.data[noise].matrix is None:
                self.data[noise].generate(0,2*tmax,tsample )

        xeff,alive,effshapes=self.get_alive(x,death=-1)#death)

        #Prepare dynamics: create effective matrices used here
        #TODO:(reduced as dead species are removed)
        diff ={var:LabeledArray(np.zeros( self.variables[var]['shape'] ),
                axes= self.get_axes(var)) for var in self.variables}


        #To store absdx in order to check how close to equilibrium a species is (judged as |dx|/absdx)
        absdiff ={var:LabeledArray(np.zeros( self.variables[var]['shape'] ),
                axes= self.get_axes(var)) for var in self.variables}


        dynamics=self.make_dynamics()

        deltax,xdeltax=10**6,0
        xflat=self.flatten(xeff)
        error=None
        old=np.ones(xflat.shape)
        if tsample=='adaptive':
            deltat=tmax/100.
        else:
            deltat=tsample

        lastfail=0
        eqs=[]
        while t<=tmax-deltat and not (t>0 and converge and np.abs(deltax) < ftol) :
            cumul=[0]+list(np.cumsum([np.prod(z) for z in effshapes])) #cumulative shapes for get_dx
            def get_dx_ode(t,xflat):
                for i,var in enumerate(self.variables):
                    xeff[var].matrix[:]=xflat[cumul[i]:cumul[i+1]]

                for var in diff:
                    diff[var][:]=0
                #Dynamical modules
                for dyn,data in dynamics.iteritems():
                    if not data['roles']:
                        continue
                    data['module'].dx(diff,xeff,t=t,data=self.data, **data)
                flat= self.flatten(diff)
                if dftol>0:
                    flat[xflat<dftol]=np.maximum(0,flat[xflat<dftol])
                # print xflat,flat
                return flat


            if self.noises:
                integ='lvode'
            else:
                integ='dop853'
            integrator=scint.ode(get_dx_ode).set_integrator(integ,nsteps=500000)
            integrator.set_f_params()

            success=False
            tries=0
            while not success:
                if tries>kwargs.get("MAX_TRIES",500):
                    raise
                if tries>0 and deltat>10**-15:
                    lastfail=t
                    if deltat/2 >= kwargs.get('mintsample',0):
                        deltat/=2.
                        if print_msg:
                            print 'Failure, trying tsample', deltat
                    else:
                        error = 'ERROR: tsample < mintsample required for integrator to converge'
                        break
                integrator.set_initial_value(xflat, t)  # .set_f_params(reff)
                result=integrator.integrate(integrator.t+min(deltat,tmax-integrator.t) )
                success=integrator.successful()
                if success and not tries and t>deltat*100 and t-lastfail>deltat*20:
                    deltat*=2
                tries+=1



            if np.isnan(result).any():
                error = 'ERROR: NaN'
                break
            elif (result>divergence).any():
                error= 'ERROR: Divergence'
                break
            elif death>=0 and (result<0).any():
                error='WARNING: Negative (|x|={})'.format(np.max(np.abs(result[result<0] )))
                result[result<0]=death*.9
                #break
            result[result<death]=max(0,min(death*.9,dftol))

            if remove_dead:
                result[result<death]=0
                self.update(x,result,effshapes,alive=alive)
                xeff,alive,effshapes=self.get_alive(x,death=death)
            else:
                self.update(x,result,effshapes,alive=alive)
            xeff,alive,effshapes=self.get_alive(x,death=min(0,death))
            xflat=self.flatten(xeff)
            t=min(tmax,t+deltat)
            #External dynamics
            for i,j in scenarios.iteritems():
                self.parameters[i]=j.get(t,interpolate=False)

            deltax=[z if (not np.isnan(z) and result[i]>death)
                else 0 for i,z in enumerate((result/old -1)) ]




            argdeltax=np.argmax(np.abs(deltax))
            xdeltax=result[argdeltax]
            deltax=deltax[argdeltax]

            if print_msg:
                print t,deltax,xdeltax,argdeltax#,result[argdx],old[argdx]#, x['n'][np.where(xflat>death)[0]]
                #print result
                #print old
                #if remove_dead:
                    #print 'DIED',xeff.shape
            if keep=='all' and t<tmax:
                #raise
                self.append_results(index=t,**x)
            old=result.copy()
            if converge=='force':
                # SEE HOW CLOSE TO EQUILIBRIUM
                def get_absdx_ode(t, xflat):
                    for i, var in enumerate(self.variables):
                        xeff[var].matrix[:] = xflat[cumul[i]:cumul[i + 1]]

                    for var in absdiff:
                        absdiff[var][:] = 0
                    # Dynamical modules
                    for dyn, data in dynamics.iteritems():
                        if not data['roles']:
                            continue
                        if hasattr(data['module'], 'absdx'):
                            data['module'].absdx(absdiff, xeff, t=t, data=self.data, **data)
                    flat = self.flatten(absdiff)
                    return flat
                absdx=get_absdx_ode(t,xflat.copy())
                dx=get_dx_ode(t,xflat.copy())
                movt=np.abs(dx[xflat>death])/absdx[xflat>death]
                # print 'Movt mean',np.mean(movt<1./xflat.shape[0])

                if np.mean(movt<1./xflat.shape[0] )>0.1:
                    if print_msg:
                         print 'Checking if close to equilibrium'
                    movt = np.abs(dx) / absdx
                    if 1:
                        import scipy.optimize as sopt

                        x2=xflat.copy()
                        decaying = np.logical_and(movt > 1. / xflat.shape[0], dx < 0)
                        if sum(decaying)!=0:
                            def get_dx(xx):
                                x=np.zeros(x2.shape)
                                x[decaying==False] = 10. ** np.clip(xx, -20, 15)
                                return get_dx_ode(t,x)[decaying==False]/x[decaying==False]

                            r = sopt.root(get_dx, np.log10(x2[decaying==False]))
                            eqx = np.zeros(x2.shape)+death/2
                            eqx[decaying == False] = 10 ** np.clip(r.x, -20, 15)
                            if r.success and not np.max(get_dx_ode(t,eqx)[decaying])>death:
                                print 'FOUND EQ', np.max(np.abs(get_dx_ode(t,eqx)/eqx))
                                eqs.append(eqx)
                                # if np.sum(eqx>10**-10)>1:
                                #     code_debugger()
                                self.update(x, np.clip(eqx, death, None), effshapes, alive=alive)
                                self.append_results(index=tmax, **x)
                                t=tmax
                                break

        #END LOOP=================================

        if not error is None:
            print "ERR", error
            if 'ERROR' in error:
                return False
        self.append_results(index=t,**x)


        if converge and not eqs:
            if t>=tmax:
                print 'FAILED TO CONVERGE, DX/X:{} (X: {})'.format(deltax,xdeltax)
        # if converge =='force':
        #     for idx,eqx in enumerate(eqs):
        #         self.update(x, np.clip(eqx, death, None), effshapes, alive=alive)
        #         self.append_results(index=t+idx+1,**x)

        return True

BaseModel.klasses['DynamicModel']=DynamicModel
