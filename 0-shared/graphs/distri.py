# -*- coding: utf-8 -*-
from base import *

### REWRITE THIS WITH NEW DISTRIBUTIONS

class DistriNetwork(CommNetwork):
    '''Extension of CommNetwork to case where edges are probability distributions'''
    edge_data=True

    def __init__(self,edge_data=None,**kwargs):
        import pandas as pd, networkx as nx
        self.props=kwargs
        self.graph=nx.DiGraph()
        if isinstance(edge_data,np.ndarray):
            self.from_matrix(edge_data)
        else:
            self.edge_data=pd.DataFrame(edge_data)
            self.from_df(self.edge_data)
        self.correlations=kwargs.get('correlations',None)

    def __getitem__(self,*args,**kwargs):
        return self.matrix.__getitem__(*args,**kwargs)

    def __setitem__(self,*args,**kwargs):
        d=Distribution(weight=args[1],**self.props)
        #kwargs['weight']=args[1]
        kwargs['distri']=d
        kwargs['weight']=args[1]
        self.add_edge(args[0][0],args[0][1],**kwargs)
        return self.matrix.__setitem__(*args)

    def set(self,i,j,**kwargs):
        props={}
        props.update(self.props)
        props.update(kwargs)
        d=Distribution(**props)
        kwargs['distri']=d
        self.matrix[i,j]=kwargs['weight']
        return self.add_edge(args[0][0],args[0][1],**kwargs)

    def add_edge(self,i,j,**kwargs):
        import pandas as pd
        edgdat={}
        edgdat.update(kwargs)
        edgdat.update(edgdat.pop('distri').attr )
        if (i,j) in self.edge_data.index:
            self.edge_data.loc[i,j]=pd.Series(edgdat)
        else:
            df=pd.DataFrame([edgdat],index=[(i,j)] )
            df.index = pd.MultiIndex.from_tuples(df.index)
            self.edge_data=self.edge_data.append(df )
        return self.graph.add_edge(i,j,**kwargs)


    def from_matrix(self,matrix):
        import pandas as pd
        df=pd.DataFrame.from_records(matrix).stack()
        df.name='weight'
        self.edge_data=pd.DataFrame(None)
        self.from_df(pd.DataFrame(df))

    def from_df(self,df):
        size=self.props.get('size',1)
        matrix=np.zeros((size,size) )
        if not df is None:
            for idx,series in df.iterrows():
                row={i:j for i,j in self.props.iteritems() if i in Distribution.attributes}
                row.update(series)
                if not 'weight' in row:
                    row['weight']=row.pop(df.columns[0] )
                matrix[idx]=row['weight']
                row['distri']=Distribution(**row)
                self.add_edge(idx[0],idx[1],**row)
        self.matrix=matrix


    def edges(self,s):
        lst=[]
        edg=self.graph.edge
        nei=sorted(edg[s])
        for i in nei:
            lst.append( edg[s][i]['distri'])
        return nei,lst

    def copy(self):
        return self.__class__(self.edge_data.copy(),**deepcopy(self.props))

    def generate_seed(self,species):
        '''For particle generation'''
        tup=[]
        edges=self.graph.edges(species)
        if self.correlations is None:
            for e in edges:
                tup.append(np.random.random())
        else:
            import scipy.stats
            if hasattr(self.correlations,'__iter__'):
                print 'Correlation graph not implemented yet'
                cov=0
            else:
                nedg=len(edges)
                cov=np.ones((nedg,nedg))*self.correlations
                np.fill_diagonal(cov,1)
            rdm=np.random.multivariate_normal(np.zeros(nedg),cov)
            tup=[scipy.stats.norm.cdf(x,0,1) for x in rdm ]

        return tuple(tup)
