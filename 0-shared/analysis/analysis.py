# -*- coding: utf-8 -*-

from measures import *
from plotting import *
from movies import *
from models import BaseModel
import os, gc

def rebuild_filelist(path,verbose=True):
    import pandas as pd


    if verbose:
        print 'Rebuilding files.csv for {}'.format(path)
    final=None
    idx=0
    for root,dirs,files in os.walk(str(path)):
        if str(root)!=str(path):
            if 'files.csv' in files:
                pass
            elif 'model.dat' in files and 'results.csv' in files:
                print 'Missing files.csv, regenerating'
                dic={'path':root }
                #dic.update( eval(open('model.dat','r' ) ) )
                df=pd.DataFrame([dic])
                df.to_csv(Path(root)+'files.csv')
            else:
                continue
            if verbose:
                print '{} Found files.csv in {}'.format(idx, root)
                idx+=1
            if final is None:
                final=pd.read_csv(Path(root)+'files.csv',index_col=0)
            else:
                final=final.append(pd.read_csv(Path(root)+'files.csv',index_col=0),ignore_index=1)
            final.set_value(final.index[-1],'path',str(Path(root).osnorm()))

    if not final is None:
        final.to_csv(path+'files.csv')
    else:
        print 'Error: No file found! {}'.format(path)



def extract_data(path,**kwargs):
    '''Go through arborescence inside path, opening and yielding all the
    encountered model outputs.'''

    path=Path(path)
    import pandas as pd

    PROPERTIES=['path','label','sys','init','tmax','tsample']+kwargs.get('keep_cols',[])

    if not 'files.csv' in os.listdir(path):
        rebuild_filelist(path,verbose=True)
    filelist=read_csv(path+'files.csv')
    measures=[]
    for idx,f in filelist.iterrows():
        if kwargs.get('print_msg',True):
            print idx, f['path']
        measure={}
        measure.update({z:f.get(z,None) for z in set(PROPERTIES) } )
        if measure['sys']==None and '/sys' in measure['path']:
            try:
                measure['sys']=int(measure['path'].split('/sys-')[1].split('/')[0].split('_')[0])
            except Exception as e:
                print 'extract_data:',e
                measure['sys']=None

        #print [v for v in measure if measure[v] is None]


        try:
            models=kwargs.get('Model',BaseModel).load(f['path'],initialize=False)
        except Exception as e:
            print 'extract_data:',e
            continue
        models = to_iterable(models)
        for im, model in enumerate(models):
            measure.update({z:f.get(z,model.get_param(z) ) for z in model.export_params() } )

            if 'read_filter' in kwargs and not kwargs['read_filter'] is None:
                to_read=True
                filt=kwargs['read_filter']
                if not hasattr(filt[0],'__iter__'):
                    filt=[filt]
                for key,operator,val in filt:
                    if isinstance(val,basestring):
                        val='"{}"'.format(val)
                    if not eval('measure.get("{}",None) {} {}'.format(key,operator,val)):
                        if kwargs.get('print_msg',True):
                            print 'Skipping because wrong: {} {} {}'.format(key,operator,val)
                        to_read=False
                        break
                if not to_read:
                    continue
            if len(models)>1:
                measure['replica']=im
            yield idx,f,model,measure

def do_split_arrays(mlist,split_arrays,ref=None):
    newlist=[]
    for nxt in mlist:
        for s in nxt:
            if isinstance(nxt[s],basestring):
                if 'array' in nxt[s] or '(' in nxt[s]:
                    print '   WARNING!!! Incorrectly retupleified', s
        if isinstance(split_arrays,dict):
            candidates=[]
            if 'contains' in split_arrays:
                candidates=[n for n in nxt if split_arrays['contains'] in n ]
            if 'add' in split_arrays:
                candidates+=split_arrays['add']
            if 'remove' in split_arrays:
                candidates=[c for c in candidates if not c in split_arrays['remove'] ]
            if '!contains' in split_arrays:
                candidates = [c for c in candidates if not True in [x in c for x in to_iterable(split_arrays['!contains']) ] ]
            locsplit=candidates
        else:
            locsplit=list(split_arrays)

        locsplit=[s for s in locsplit if not isinstance(nxt[s],basestring) and hasattr(nxt[s],'__iter__')  ]
        if not locsplit:
            continue
        if ref is None:
            maxlen=max(len(nxt[s]) for s in locsplit )
        else:
            maxlen=len(nxt[ref])
        if not newlist:
            bad=[s for s in locsplit if len(nxt[s])< maxlen]
            if bad:
                print 'do_split_arrays: too short:' + ','.join(bad)
        locsplit=[s for s in locsplit if len(nxt[s])==maxlen]
        try:
            vals= zip(*[nxt[s] for s in locsplit   ] )
        except Exception as e:
            for i,s in enumerate(locsplit):
                print '->',i,s,nxt[s]
            raise e
        base=[n for n in nxt if not n in locsplit]
        for vs in vals:
            dic={i:nxt[i] for i in base}
            for i,j in zip(locsplit,vs):
                dic[i]=j
            newlist.append(dic)
    return newlist

def make_measures(path,use_measures=(),maxsize=10*10**9,filespath=None,
                  # aggregate_over=None,
                  fname='measures.csv',**kwargs):
    '''
    Kwargs:
        use_measure: specify other measurement functions than the usual.'''
    if filespath==None:
        filespath=path
    if not '.' in fname:
        fname=fname+'.csv'
    import pandas as pd
    measures=[]
    grouped={}
    use_prev=kwargs.get('use_prev',0)

    prevpaths={}
    killprev=False
    prevmeas=None
    if use_prev:
        if fname in os.listdir(path):
            prevmeas=pd.read_csv(path+fname,index_col=0)
            if os.path.getsize(path+fname) > 5.*10**8:
                print 'Previous {} is large (>500 Mb), will not be kept open at all times.'.format(fname)
                killprev=True
            prevpaths.update( {k:fname for k in prevmeas['path'].values })
        for i in os.listdir(path):
            if  fname+'.' in i:
                print i
                try:
                    red=pd.read_csv(path + i, index_col=0)
                except Exception as e:
                    if os.path.getsize(path+i)< 8:
                        print 'make_measures: empty file'
                    else:
                        print 'make_measures:',e
                for k in red['path'].values:
                    print '     ',k, k in prevpaths
                    if k in prevpaths:
                        other=prevpaths[k]
                        print 'Duplicate of',other
                        if os.path.getmtime( path+other) > os.path.getmtime( path+ i):
                            os.remove(path+i)
                            print "Removed",i
                            continue
                        else:
                            os.remove(path+other )
                            print "Removed", other
                    prevpaths[k]=i

    for content in extract_data(filespath,**kwargs):
        idx,f,model,props=content

        if  use_prev:# and not aggregate_over:
            if f['path'] in prevpaths:
                handled = True
                if isinstance(use_prev, tuple):
                    key, operator, val = use_prev
                    if isinstance(val, basestring):
                        val = '"{}"'.format(val)
                    if not eval('props.get("{}",None) {} {}'.format(key, operator, val)):
                        handled = False
                if handled:
                    print '    Already in {}'.format(prevpaths[f['path']])
                    if prevpaths[f['path']] == fname :
                        # Export this part of prevmeas so it can be merged later
                        if isinstance(prevmeas,basestring):
                            prevmeas = pd.read_csv(path +prevmeas, index_col=0)

                        locidx = idx
                        while fname+'.{}'.format(locidx) in os.listdir(path):
                            locidx+=1
                        prevmeas[prevmeas['path'] == f['path']].to_csv(path + fname + '.{}'.format(locidx), index=False)
                        if killprev:
                            del prevmeas
                            gc.collect()
                            prevmeas=fname
                    continue
            # if not prevmeas is None and not isinstance(prevmeas,basestring):
            #     import psutil
            #     print psutil.virtual_memory()
            #     del prevmeas
            #     prevmeas= fname #CLOSE prevmeas while doing something else
            #     print "NOW",psutil.virtual_memory()
        measure={}
        measure.update(props)
        if use_measures:
            #Externally provided measure functions
            txts=[]
            for meas in use_measures:
                if meas is None:
                    continue
                if isinstance(meas,basestring):
                    txts.append(meas)
                else:
                    if txts:
                        measure_gen(model,measure,typ=txts,**kwargs)
                        txts=[]
                    meas(model,measure,**kwargs)
            if txts:
                measure_gen(model,measure,typ=txts,**kwargs)

        else:
            measure_usual(model,measure)


        mlist=[measure]
        measure.update(props)

        split_arrays=kwargs.get('split_arrays',None)
        if split_arrays:
            #Measure contains arrays that should be made into multiple lines
            #of measures file
            mlist=do_split_arrays(mlist,split_arrays)
            #print mlist
        # else:
        #     def totuple(a):
        #         if not hasattr(a,'__iter__'):
        #             return a
        #         try:
        #             return tuple(totuple(i) for i in a)
        #         except TypeError:
        #             return a
        #     for i,j in measure.iteritems():
        #         if isinstance(j, np.ndarray):
        #             measure[i]='array({})'.format(totuple(j))
        #         else:
        #             measure[i]=totuple(j)

        # if not aggregate_over:

        mlist=[tupleify(x) for x in mlist]

        locidx = idx
        while fname + '.{}'.format(locidx) in os.listdir(path):
            locidx += 1
        pd.DataFrame(mlist).to_csv(path+fname+'.{}'.format(locidx),index=False )

        del mlist
        del measure
        del model
        gc.collect()

        # else:
        #     dico={}
        #     dico.update(props)
        #     if hasattr(aggregate_over,'__iter__'):
        #         identifier=sorted( dico.pop(agg) for agg in aggregate_over )
        #     else:
        #         identifier=dico.pop(aggregate_over)
        #     key=str(tuple(sorted(dico.items())))
        #     grouped.setdefault(key,[identifier,dico,[],[]])
        #     for m in mlist:
        #         grouped[key][2].append(model)
        #         grouped[key][3].append(m)


    used_tmps=[]
    totsize=0
    # if aggregate_over:
    #     for group,key in enumerate(sorted(grouped)):
    #         identifier,dico,models,meas=grouped[key]
    #         for dico in meas:
    #             dico['group']=group
    #             lst.append(pd.DataFrame(dico) , ignore_index=1)
    # else:
    for i in os.listdir(path):
        if fname+'.' in i:
            used_tmps.append(i)
            totsize+=os.path.getsize(path+i)


    if totsize < maxsize:
        lst=[pd.read_csv(path+i) for i in used_tmps if len([l for l in open(path+i,'r')])>1 ]
        if lst:
            measures=pd.concat(lst, ignore_index=1)
        else:
            print 'ERROR IN MAKE_MEASURES: EMPTY LIST'
            code_debugger()

        if kwargs.get('update_prev',False):
            #Update previous dataframe, keeping info that was not measured this time
            #prev.update(measures)
            try:
                prev = open_measures(path,fname)
                prev[[k for k in measures]]=measures
            except Exception as e:
                print e
                prev=measures
            measures=prev
        # code_debugger()
        measures.to_csv(path+fname)
        for i in  used_tmps:
            os.remove(path+i)
    else:
        print( 'make_measures: Size exceeded maxsize ({} > {}), not joining tmp files!'.format(totsize,maxsize) )
    return measures
