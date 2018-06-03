# -*- coding: utf-8 -*-

from plotting import *
from utility import *

def run_movie(path,destination='./Movie',death=None,dpi=150,tsample=1,imgformat='jpg',maxpic=800,noticks=False ):
    '''Creates a movie at destination. Images are stored at destination/img/'''
    destination = Path(destination)
    imgfolder=destination+Path('img')
    imgfolder.mkdir()
    from subprocess import call

    import matplotlib.patches as mpatches
    from matplotlib.collections import PatchCollection
    from models import BaseModel

    ## Remove old figs
    try:
        call("rm -f {}.{}".format(imgfolder+'*',imgformat),shell=True )
    except:
        print("Your OS is not supported by my awesome movie plotter. Old figures not removed")

    ## Load in data
    m=BaseModel.load(path,initialize=False)
    var=m.variables.keys()[0]
    traj=m.results[var]

    prm=m.export_params()
    colors=['b','g','r','y','k',]

    if 'death' in m.parameters[var]:
        death=m.parameters[var]['death']


    if death is None:
        death=0

    #Size conversion (data to figure size)
    def convert_dist(ax,pos):
        if not hasattr(pos,'__iter__'):
            pos=(0,pos)
        tarea=( ax.transData.transform((area[2],area[3] ) )-ax.transData.transform((area[0],area[1]) ) )
        return (ax.transData.transform(pos ) - ax.transData.transform( (0,0) ) )/tarea

    graph=Graph(nx.from_numpy_matrix( m.data['community'].matrix) )

    ## plot
    alive=graph.nodes()[:]
    for t in range(0,min(maxpic*tsample,traj.shape[0])):
        if t%tsample!=0:
            continue

        plt.clf()
        fig = plt.figure( edgecolor=[.4,.4,.4]);
        ax = fig.add_subplot(121)


        dead=np.where( traj[t]<death)[0]
        for d in dead:
            if d in alive:
                graph.remove_node(d)
                alive.remove(d)

        graph.plot(hold=1,newfig=0,node_size=traj[t][alive] ,edge_width=.1)
        ax = fig.add_subplot(122)
        plt.xlim(xmax=traj.index[-1] )
        plot(np.array([traj[t2].ravel() for t2 in range(t)]),hold=1,xs=traj.index[:t] )


        if noticks:
            plt.xticks([])
            plt.yticks([])

        # save
        ti = str(t/tsample+100000)
        plt.savefig(imgfolder+'Fig_{}.{}'.format(ti[1:],imgformat),dpi=dpi,bbox_inches='tight')
        print(t)

        plt.close('all')

    try:
        call("ffmpeg -y -r 18  -i {}%05d.{} {}".format(imgfolder+'Fig_',
            imgformat,destination+'movie.mp4'),shell=True)
    except:
        pass