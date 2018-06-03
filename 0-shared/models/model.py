# -*- coding: utf-8 -*-
from imports import *


class BaseModel(object):
    PARALLEL_LOCK=False
    klasses={}
    dft_prm={}
    submodels=() #Override in children classes if different submodels need different pieces of data

    generators={
        'variable':Trajectory,
        'constant':LabeledArray,
        'noise':Noise,
        }

    def __init__(self,results=None,**kwargs):
        self.parameters=deepcopy(kwargs.pop('parameters',self.dft_prm))
        self.set_params(**kwargs)
        self.data={}
        if results is None:
            results={}
        self.results=deepcopy(results)
        if kwargs.get('initialize',True):
            self.initialize()

    def initialize(self,labels=None):
        """Internal use. Generate model data from parameter."""

        if labels is None:
            prms=self.parameters.items()[:]
        else:
            prms=[(label,self.parameters[labels]) for label in labels ]
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
                req=[l for l in self.parameters if (l in formulacomps
                     or l in prm.get('requires',[])) and not l in done]
                if req:
                    prms.append( (label,prm) )
                    continue
            done.append(label)

            if True in [mod  in prm and not prm[mod] in getattr(self,mod)
                for mod in self.submodels]:
                    continue

            if label in self.data:
                dat=self.data[label]
            elif 'formula' in prm:
                self.data[label]=dat=self.generate_from_formula(label,prm)
            else:
                self.data[label]=dat=self.generate(label,prm)

    def generate(self,label,kwargs):
        """Generate data (matrix,noise,initial conditions) from parameters"""
        flat=self.export_params()
        prm={}
        prm.update(flat)
        prm.update(kwargs)
        typ=prm['type']
        if not 'shape' in prm and 'variables' in prm:
            prm['shape']=self.get_shape(prm)
        if 'requires' in prm:
            prm['required']={self.parameters[d]['role']:self.data[d] for d in prm['requires']}

        if typ=='variable':
            return Trajectory(init=BaseSampler().sample(**prm))

        if typ=='constant':
            if not 'shape' in prm:
                #Take shape from variables
                shape=[]
                for var,ax in prm['axes']:
                    vpm=self.parameters[var]
                    shape.append(vpm['shape'][vpm['axes'].index(ax)]  )
                prm['shape']=shape
            prm['shape']=tuple(self.get_param(s) if isinstance(s,basestring) else s for s in prm['shape'])
            if 'axes' in prm:
                axes=[(label,ax) if not hasattr(ax,'__iter__') else ax for ax in  prm['axes']  ]
            else:
                axes=[None for s in prm['shape']]
            return LabeledArray(BaseSampler().sample(**prm),axes )


        if typ=='noise':
            shape=self.get_shape(prm)
            prm.setdefault('shape',shape)
            axes=[]
            for v in prm['variables']:
                if isinstance(v,basestring):
                    axes+=self.get_axes(v)
                else:
                    axes.append(v)
            prm['axes']=axes
            prm['shape']=shape
            return LabeledNoise(**prm)

        if 'shape' in prm:
            prm['shape']=tuple(self.get_param(s) if isinstance(s,basestring) else s for s in prm['shape'])
        return self.generators[typ].make(model=self,**prm)

    def generate_from_formula(self,label,prm):
        """Generate data (matrix,noise,initial condition) from formula"""
        typ=prm['type']
        #print [(label,self.parameters[label].keys()) for label in self.data]
        data={lab:self.get_labeled_data(lab) for lab in self.data}
        if typ=='matrix':
            shape= self.get_shape(prm)
            res= self.generators[typ](matrix=eval(prm['formula'],{},data).matrix )
            assert res.matrix.shape==shape
            return res
        raise Exception('generate_from_formula not implemented yet for type {}'.format(typ) )


    def renew(self):
        self.results={}
        self.initialize()

    def append_results(self,index=None,**kwargs):
        for i,j in kwargs.iteritems():
            if j is None:
                continue
            res=j
            if not hasattr(res,'__iter__') and not hasattr(res,'matrix'):
                res=np.array([res])
            if not i in self.results:
                self.results[i]=Trajectory(init=res,t=index)
            else:
                self.results[i].append(res,index=index)


    @property
    def description(self):
        '''Information that allows to reconstruct the model.'''
        dic={}
        dic['parameters']=self.parameters
        dic['klass']=self.__class__.__name__
        return dic

    @classmethod
    def load(klass,dirpath,initialize=True,**kwargs):
        '''Open model+results stored in a specific directory.
        Supports batch loading
        If initialize==False, do not launch the modules' initialize method
        (can speed up loading significantly). '''
        import os
        dirpath=Path(dirpath)
        try:
            assert os.path.isdir(dirpath)
        except:
            raise Exception('Directory {} not found'.format(dirpath))


        attrs=loads(open(dirpath+'model.dat','r'))
        use_class=attrs.pop('klass',None)
        if not use_class is None:
            #Change class!
            klass= BaseModel.klasses[use_class]
        elif klass == BaseModel:
            klass=BaseModel.klasses.values()[0]
            print "WARNING: Model not specified, using {}".format(klass.__name__)
        model= klass(initialize=False,**attrs)

        #Try to load stored data (matrix, noise)
        for label,prm in model.parameters.iteritems():

            genname=prm.get('type',label)
            if not genname in model.generators:
                continue
            gen=model.generators[genname](name=label)
            if not hasattr(gen,'load') :
                continue
            try:
                network=gen.load(dirpath,label)
                model.data[label]=network
            except Exception as e:
                if not (label in model.results or prm['type']=='variable') and prm.get('save',1):
                    print 'WARNING: Model.load',e


        models=[]
        if 'results.compressed.csv' in os.listdir(dirpath):
            #Uncompress to temporary files
            txt=''.join([l for l in
                open(dirpath+'results.compressed.csv','r').readlines() if l.strip()]).split('$END')
            from cStringIO import StringIO
            files=[StringIO(f) for f in txt if f.strip()]
            if 0:
                txt=''
                i=0
                for l in open(dirpath+'results.compressed.csv','r'):
                    if not '$END' in l:
                        txt+=l
                    else:
                        f=open(dirpath+'results.{}.csv~'.format(i),'w' )
                        f.write(txt)
                        f.close()
                        txt=''
                        i+=1
        else:
            files =[dirpath+ fname for fname in os.listdir(dirpath)
                if 'results' in fname and '.csv' in fname and not '.compressed.' in fname]

        for fil in files:
            index,matrix=load_table(fil)
            locmodel=model.copy(initialize=initialize)
            for name,mat in matrix.iteritems():
                locmodel.results[name]=Trajectory(matrix=mat,index=index)
            models.append(locmodel)
            if initialize == True:
                locmodel.initialize()
            if 0:
                if isinstance(fil,basestring) and fil[-1]=='~':
                    #Remove temporary file
                    os.remove(dirpath+fname)
        if len(models)==1 and not kwargs.get('force_list',False):
            return models[0]
        return models

    def compress_results(self,dirpath,remove=True):
        """If many results"""
        import os
        f= open(dirpath+'results.compressed.csv','w')
        for fname in os.listdir(dirpath):
            if 'results' in fname and '.csv' in fname and not '.compressed.' in fname:
                for l in open(dirpath+fname,'r'):
                    f.write(l)
                f.write('\n$END\n')
                if remove:
                    os.remove(dirpath+fname)
        f.close()


    def save(self,path,overwrite=True,save_matrix='all',suffix=None):
        '''Save model+results at path.'''
        import os
        path=Path(path)
        path.norm()
        if not self.PARALLEL_LOCK:
            #If parallel lock, do not re-save model parameters and module contents
            path.mkdir()

            for label,data in self.data.iteritems():
                if self.parameters[label].get('save',True) == False:
                    continue
                data.save(path,label)

            def filtered(dic):
                if not hasattr(dic,'keys'):
                    return dic
                return {i:filtered(j) for i, j in dic.items() if i!='matrix'}

            dumps(open(path+'model.dat','w'),filtered(self.description) )

        if suffix is None:
            fname='results.csv'
        elif suffix is True:
            i=0
            fname='results.{}.csv'.format(i)
            while fname in os.listdir(path):
                i+=1
                fname='results.{}.csv'.format(i)
        else:
            fname='results.{}.csv'.format(suffix)
        save_table(self.results,path=path+fname )
        return

    def copy(self,initialize=True):
        copy=self.__class__(initialize=initialize,**deepcopy(self.description))
        for name, res  in self.results.iteritems():
            if res is None:
                continue
            copy.results[name]=res.copy()
        for name, data  in self.data.iteritems():
            if data is None:
                continue
            copy.data[name]=data.copy()
        return copy


    def export_params(self,parameters=None):
        """Flatten parameter structure to the format 'label_name'. """
        if parameters is None:
            parameters=self.parameters
        final={}
        #for label,prm in self.parameters.iteritems():
            #for i,j in prm.iteritems():
                #params['{}_{}'.format(label,i)]=j
        def recurs(path,dico):
            if hasattr(dico,'keys'):
                for i,j in dico.iteritems():
                    recurs('{}_{}'.format(path,i),j )
            else:
                final[path]=dico
        for label in parameters:
            recurs(label,self.parameters[label])

        return final

    def set_params(self,**kwargs):
        """If some parameter 'name' in the template 'label' is changed by
        external arguments to the model of the form 'label_name', update it."""
        for i,j in kwargs.iteritems():
            if i in self.dft_prm:
                self.parameters.setdefault(i,{})
                for k in j:
                    self.parameters[i][k]=j[k]
                continue
            try:
                part=i.split('_')
                dico=self.parameters
                p=part.pop(0)
                while part:
                    if p in dico:
                        dico=dico[p]
                        p=part.pop(0)
                    else:
                        p+='_'+part.pop(0)
                if p in dico or kwargs.get('force',1) and p!=i:
                    if not p in dico:
                        print 'set_params forcing:',i,p,j
                    dico[p]=j
            except:
                pass
        return self.parameters

    def get_param(self,param):
        '''Get parameter value from flattened version'''
        part=param.split('_')
        dico=self.parameters
        p=part.pop(0)
        while part:
            if p in dico:
                dico=dico[p]
                p=part.pop(0)
            else:
                p+='_'+part.pop(0)
        return dico.get(p,None)

    def get_param_label(self,**kwargs):
        '''Query: fetch parameter label from attributes.'''
        for i,j in self.parameters.iteritems():
            handled=True
            for k,l in kwargs.iteritems():
                if kwargs[k]!= self.get_param('{}_'.format(i)+k):
                    handled=False
                    break
            if handled:
                return i
        return None

    def find_data(self,**attrs):
        candidates=sorted(self.data)
        for i, j in attrs.iteritems():
            candidates=filter(lambda e:self.parameters[e].get(i,None)==j,candidates )
        cd= [self.data[c] for c in candidates ]
        if len(cd)==1:
            return cd[0]
        return cd

    def get_labeled_data(self,label=None):
        if label is None:
            return {lab:self.get_labeled_data(lab) for lab in self.data}
        var=[]
        prm=self.parameters[label]
        if not 'variables' in prm and not 'axes' in prm:
            #No labels possible
            return self.data[label]
        if 'axes' in prm:
            axes=[]
            for ax in prm['axes']:
                if len(ax)==1:
                    axes.append((label,ax) )
                else:
                    axes.append(ax)
            return LabeledArray(self.data[label].matrix,axes=tuple(axes))
        for v in prm['variables']:
            if isinstance(v,tuple):
                var.append(v)
            else:
                var+=[(v,ax) for ax in self.parameters[v]['axes'] ]
        return LabeledArray(self.data[label].matrix,axes=tuple(var))

    def get_labeled_result(self,label=None,idx=None):
        if label is None:
            return {lab:self.get_labeled_result(lab,idx) for lab in self.results}
        if idx is None:
            return LabeledArray(self.results[label].matrix,axes=('t',)+
                tuple( (label,ax) for ax in self.parameters[label]['axes']))
        return LabeledArray(self.results[label].matrix[idx],axes=
            tuple( (label,ax) for ax in self.parameters[label]['axes']))




    def get_shape(self,prm):
        if 'shape' in prm:
            return prm['shape']
        shape=()
        for v in prm['variables']:
            if isinstance(v,basestring):
                shape+=self.parameters[v]['shape']
            else:
                shape+=self.parameters[v[0]]['shape'][ self.get_axes(v[0]).index(v) ],
        return shape


    def get_axes(self,v):
        return [(v,ax) if not v in ax else ax for ax in self.parameters[v]['axes']]
