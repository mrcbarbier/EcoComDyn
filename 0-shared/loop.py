# -*- coding: utf-8 -*-

from models import *
from imports import *

def mp_run(array):
    #Multiprocessing function
    path,Model,prm,results,tmax,tsample,keep,replica=array
    model=Model(results=results,**prm)
    model.PARALLEL_LOCK=True
    model.evol(tmax=tmax,tsample=tsample,keep=keep)
    model.save(path,overwrite=0,suffix=replica )
    model.module['prey'].tree=None

class SkipLoop(Exception):
    """Exception thrown out to force to skip a loop key."""
    def __init__(self,*args,**kwargs):
        Exception.__init__(self,*args,**kwargs)

class Looper(object):

    '''
    Options:
        clean_directory
        reseed
        replicas

    Override evol() for more complex scenarios (perturbations, coalescence...)'''

    MAX_TRIALS=10
    Model=DynamicModel

    def __init__(self,model=Model,*args,**kwargs):
        self.Model=model
        self.attrs={}
        self.attrs.update(kwargs)

    def loop_core(self,tsample=50.,tmax=None, path='.',**kwargs):
        import pandas as pd
        for i,j in self.attrs.iteritems():
            kwargs.setdefault(i,j)


        path=Path(path)
        path.norm()
        path.mkdir(rmdir=kwargs.get('clean_directory',False))
        table=[]

        dic={
            'tsample':tsample,
            'tmax':tmax,
            'path':path,
            'sys':kwargs.get('sys',None),
            }


        if 'model' in kwargs:
            model=kwargs['model']
            model.renew()
            model.set_params(**kwargs)
        else:
            model=self.Model(**kwargs)
        if not kwargs.get('reseed',1):
            kwargs['model']=model #transmit the model further down the loop

        dic.update(model.export_params())

        if 'replicas' in kwargs:
            '''Multiprocessing test'''
            from multiprocessing import Pool,freeze_support
            freeze_support()
            systems=range(kwargs['replicas'])

            array=[(path.copy(),self.Model,model.export_params(),model.results,tmax,tsample,kwargs.get('keep','endpoint'),sys)
                     for sys in systems  ]

            pool=Pool()

            pool.map(mp_run,array )
            pool.close()
            pool.join()
            pool.terminate()
            model.compress_results(path)
        else:
            kwargs.setdefault('converge',1)
            kwargs.setdefault('keep','endpoint')
            model,success=self.evol(model=model,tmax=tmax,tsample=tsample,path=path,**kwargs)
            if not success and not 'model' in kwargs:
                kwargs['loop_trials']=kwargs.get('loop_trials',0)+1

                if kwargs['loop_trials']<self.MAX_TRIALS:
                    print '### model diverged, trying again! ###'
                    return self.loop_core(tsample=tsample,tmax=tmax, path=path,**kwargs)
                else:
                    print '### model diverged {} times, setting results to 0! ###'.format(kwargs['loop_trials'])
                    for var in model.results:
                        model.results[var].matrix[:]=0
            if 'loop_trials' in kwargs:
                del kwargs['loop_trials']
            model.save(path,overwrite=1 )

        table.append(dic)
        pd.DataFrame(table).to_csv(path+'files.csv')
        return table

    def evol(self,model=None,tmax=None,tsample=None,path=None,**kwargs):
        success=model.evol(tmax=tmax,tsample=tsample,**kwargs)
        return model,success


    def loop(self,axes=None,path='.',*args,**kwargs):
        '''External loop function'''
        respath=Path(path)

        if axes is None or len(axes)==0:
            return []
        axis=axes[0]
        axes=axes[1:]


        key=axis[0]
        if 'check_key' in kwargs and not kwargs['check_key']==False:
            if kwargs['check_key'] in (True,1):
                kwargs['check_key'] = self.Model(**kwargs).export_params()
            to_check=to_iterable(key)
            for c in to_check:
                if not c in kwargs['check_key']:
                    raise Exception('Key "{}" not found among possible parameters.'.format(c))

        table =[]
        labels=[]
        refkwargs=kwargs
        for x in axis[1]:
            kwargs={}
            kwargs.update(refkwargs) #Prevents changes to kwargs from carrying over

            if hasattr(key,'__iter__'):
                #Case where multiple factors must change jointly (CANNOT BE SYS OR INIT)
                for i,j in zip(key,x):
                    kwargs[i]=j
                folder=Path('{}-{}_etc'.format(key[0],x[0]) )
                base=str(folder)
                inti=2
                while str(folder) in labels:
                    folder=Path(base+'{}'.format(inti))
                    inti+=1
                labels.append(str(folder))

            else:
                kwargs[key]=x
                folder=Path('{}-{}'.format(key,x) )

            if axes or not key in ('init','replica','sys'):
                #Do not print the innermost level of recursion if it is either
                #initial conditions, or replica (for stochastic systems)
                if kwargs.get('print_msg',1):
                    print folder

            try:
                self.key_react(key,kwargs)
            except SkipLoop as e:
                print e
                continue


            if axes:
                table+=self.loop(*args,axes=axes,path=respath+folder, **kwargs)
            else:
                table+=self.loop_core(*args,path=respath+folder, **kwargs)
        return table


    def key_react(self,key,kwargs):
        pass


def loop(axes=None,path='.',*args,**kwargs):
    if os.path.isdir(path):
        print 'WARNING: Path {} already exists.'.format(path)
    return Looper().loop(axes,path,*args,**kwargs)
