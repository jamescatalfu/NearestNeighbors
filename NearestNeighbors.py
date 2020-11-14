import numpy as np
import matplotlib.pyplot as plt
import random
import timeit
from scipy.spatial import Voronoi, voronoi_plot_2d, KDTree, Delaunay, delaunay_plot_2d
# from sklearn.neighbors import NearestNeighbors



### Random points
def getRandomPoints():
    radius = 20
    dataRangeX = (0, 250)
    dataRangeY = (0, 250)
    ###numOfPoints hits limit at 117; 118+ does not process for some reason...
    numOfPoints = 10

    deltas = set()
    for x in range(-radius, radius + 1):
        for y in range(-radius, radius + 1):
            if x*x + y*y <= radius*radius:
                deltas.add((x,y))

    randPoints = []
    excluded = set()
    i = 0
    while i < numOfPoints:
        x = random.randrange(*dataRangeX)
        y = random.randrange(*dataRangeY)
        if (x,y) in excluded:
            continue
        randPoints.append((x,y))
        i += 1
        excluded.update((x+dx, y+dy) for (dx,dy) in deltas)
    
    return randPoints
### end RandomPoints

# Gets Nearest Neighbors based on circular radius/range from point
def fixedDistanceNN(npPoints):
    
    vor = Voronoi(npPoints)

    fig = voronoi_plot_2d(vor, show_vertices=False)

    ax = fig.add_subplot(1,1,1)
    
    for i in range(0, len(npPoints)):
        ax.annotate(text = i, xy = (npPoints[i][0], npPoints[i][1]))
    
    circ = plt.Circle(npPoints[0], radius=100, fill=False, edgecolor='green')
    ax.add_patch(circ)

    print('Fixed Length Radius: ')
    print(KDTree(npPoints).query_ball_point(npPoints[0], 100))

    plt.savefig('graph_FixedDistance.png')
    
    return


def delaunayGraph(npPoints):
    delaunay = Delaunay(npPoints, qhull_options='Qc')

    fig = delaunay_plot_2d(delaunay)

    ax = fig.add_subplot(1,1,1)
    
    for i in range(0, len(npPoints)):
        ax.annotate(text = i, xy = (npPoints[i][0], npPoints[i][1]))
    
    plt.savefig('graph_Delaunay.png')
    
    indptr = delaunay.vertex_neighbor_vertices[0]
    indices = delaunay.vertex_neighbor_vertices[1]
    # ([indptr], [indices])
    # for index k, neighbors are indices[indptr[k]:indptr[k+1]]
    delaunayDict = dict()
    for k in range(0, len(npPoints)):
        startIndex = indptr[k]
        endIndex = indptr[k + 1]
        # print(k, ": ")
        for i in range(startIndex, endIndex):
            # print(indices[i])
            if k in delaunayDict:
                delaunayDict[k].append(indices[i])
            else:
                delaunayDict[k] = [ indices[i] ]

    return delaunayDict

# Gets neighbors connecting to voronoi cell
def voronoiTessellation(npPoints):
    #print(KDTree(npPoints).query_pairs(r=150))
    delaunayDict = delaunayGraph(npPoints)
    
    #print(delaunayDict)
    
    vor = Voronoi(npPoints, qhull_options='Qc')

    fig = voronoi_plot_2d(vor, show_vertices=False)

    ax = fig.add_subplot(1,1,1)

    
    for i in range(0, len(npPoints)):
        ax.annotate(text = i, xy = (npPoints[i][0], npPoints[i][1]))

    # For each vertex
    for entry in delaunayDict:
        # print(entry, delaunayDict[entry])
        polygon = [ vor.vertices[i] for neighbor in delaunayDict[entry] ]
        plt.fill(*zip(*polygon))

    print('Voronoi Tesselation Neighbors: ' + delaunayDict[0])
    for indx in delaunayDict[0]:
        #print(vor.point_region[indx])
        region = vor.regions[vor.point_region[indx]]
        if not -1 in region:
            polygon = [vor.vertices[i] for i in region]
            plt.fill(*zip(*polygon))
    
    plt.savefig('graph_VoronoiTessellation.png',)
    
    return


def performance():
    
    mySetup = '''
        import numpy as np
        import matplotlib.pyplot as plt
        import random
        from scipy.spatial import Voronoi, voronoi_plot_2d, KDTree, Delaunay, delaunay_plot_2d


    ### Random points
    radius = 20
    dataRangeX = (0, 250)
    dataRangeY = (0, 250)
    ###numOfPoints hits limit at 117; 118+ does not process for some reason...
    numOfPoints = 117

    deltas = set()
    for x in range(-radius, radius + 1):
        for y in range(-radius, radius + 1):
            if x*x + y*y <= radius*radius:
                deltas.add((x,y))

    randPoints = []
    excluded = set()
    i = 0
    while i < numOfPoints:
        x = random.randrange(*dataRangeX)
        y = random.randrange(*dataRangeY)
        if (x,y) in excluded:
            continue
        randPoints.append((x,y))
        i += 1
        excluded.update((x+dx, y+dy) for (dx,dy) in deltas)
    ### end RandomPoints
    '''
    fixedDistanceRadius = '''
        vor = Voronoi(npPoints)

        fig = voronoi_plot_2d(vor, show_vertices=False)

    ax = fig.add_subplot(1,1,1)
    
    for i in range(0, len(npPoints)):
        ax.annotate(text = i, xy = (npPoints[i][0], npPoints[i][1]))
    
    circ = plt.Circle(npPoints[0], radius=100, fill=False, edgecolor='green')
    ax.add_patch(circ)

    KDTree(npPoints).query_ball_point(npPoints[0], 100)
    '''

    #timeit statement
    print('Fixed Distance Performance Analysis: ')
    print (timeit.timeit(setup = mySetup, stmt = fixedDistanceRadius, number = 10000))
    return




# Call functions
npPoints = np.array(getRandomPoints())
fixedDistanceNN(npPoints)
performance()
#voronoiTessellation(npPoints)
#delaunayGraph(npPoints)