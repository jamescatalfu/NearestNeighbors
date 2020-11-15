import numpy as np
import matplotlib.pyplot as plt
import random
import time
from scipy.spatial import Voronoi, voronoi_plot_2d, KDTree, Delaunay, delaunay_plot_2d
# from sklearn.neighbors import NearestNeighbors

##########################################################################
##### Global variables
# Index to test getting neighbors on
indexToGetNeighbors = 0

# Radius to get neighbors (fixed Distance)
radiusParam = 50

# Number of points to get for testing
numOfPointsToGen = 50

# Random points ranges
xDataRangeLower, xDataRangeUpper = 0, 200
yDataRangeLower, yDataRangeUpper = 0, 200

# Range and domain of graph to output
xLowerLimit = -5
xHigherLimit = 250
yLowerLimit = -5
yHigherLimit = 250

##### Global variables
##########################################################################

##########################################################################
### Random points for testing the functions
def getRandomPoints():
    
    coords = [( round(random.random() * xDataRangeUpper), round(random.random() * yDataRangeUpper) )
              for _ in range(numOfPointsToGen)]
    randPoints = coords
    
    # Points further out to fix infinite Voronoi edges (coloring)
    randPoints = np.append(randPoints, [[999,999], [-999,999], [999,-999], [-999,-999]], axis = 0)
    return randPoints
### end RandomPoints
##########################################################################

##########################################################################
# Fixed Distance: Gets Nearest Neighbors based on circular radius/range from point
def fixedDistanceNN(npPoints):
    
    vor = Voronoi(npPoints)
    _ = voronoi_plot_2d(vor, show_vertices=False)
    
    # Annotate the points
    for i in range(0, len(npPoints)):
        plt.annotate(text = i, xy = (npPoints[i][0], npPoints[i][1]))
    
    # Plot the radius used to get the neighbors
    circle1 = plt.Circle(npPoints[0], radiusParam, ec='g', fc='w')
    plt.gcf().gca().add_artist(circle1)

    # Gets neighbors at indicated index for specified range
    # --> Used KDTree
    print('Fixed Length Radius Neighbors: ', 
          KDTree(npPoints).query_ball_point(npPoints[indexToGetNeighbors], radiusParam))

    # Change and print the plot to a square equal axis
    plt.axis('equal')
    plt.xlim([xLowerLimit, xHigherLimit]), plt.ylim([yLowerLimit, yHigherLimit])
    plt.savefig('graph_FixedDistance.png')
    
    return
### end fixedDistance
##########################################################################

##########################################################################
# Helper function: Delaunay graph, for Voronoi cell neighbor coloring
def delaunayGraph(npPoints):
    
    # Use delaunay triangulation graph to get neighbors
    delaunay = Delaunay(npPoints[:-4], qhull_options='Qc')
    _ = delaunay_plot_2d(delaunay)
  
    # Annotate the points  
    for i in range(0, len(npPoints[:-4])):
        plt.annotate(text = i, xy = (npPoints[:-4][i][0], npPoints[:-4][i][1]))

    # Use delaunay triangulation to get the neighbors
    indptr, indices = delaunay.vertex_neighbor_vertices
    # ([indptr], [indices])
    # for index k, neighbors are indices[indptr[k]:indptr[k+1]]
    delaunayDict = dict()
    for k in range(0, len(npPoints[:-4])):
        startIndex = indptr[k]
        endIndex = indptr[k + 1]
        for i in range(startIndex, endIndex):
            if k in delaunayDict:
                delaunayDict[k].append(indices[i])
            else:
                delaunayDict[k] = [ indices[i] ]
                
    # Change and print the plot to a square equal axis
    plt.axis('equal')
    plt.xlim([xLowerLimit, xHigherLimit]), plt.ylim([yLowerLimit, yHigherLimit])
    plt.savefig('graph_Delaunay.png')

    return delaunayDict

# Voronoi Tessellation: Gets neighbors connecting to voronoi cell
def voronoiTessellation(npPoints):
    
    # Use delaunay triangulation graph to get neighbors
    delaunayDict = delaunayGraph(npPoints)
    vor = Voronoi(npPoints, qhull_options='Qc')
    _ = voronoi_plot_2d(vor, show_vertices=False)
    
    # Annotate the points
    for i in range(0, len(npPoints)):
        plt.annotate(text = i, xy = (npPoints[i][0], npPoints[i][1]))
    
    # Color neighbors to 'indexToGetNeighbors'
    for indx in delaunayDict[indexToGetNeighbors]:
        region = vor.regions[ vor.point_region[indx] ]
        if not -1 in region:
            polygon = [vor.vertices[i] for i in region]
            plt.fill(*zip(*polygon))
    
    # Change and print the plot to a square equal axis
    plt.axis('equal')
    plt.xlim([xLowerLimit, xHigherLimit]), plt.ylim([yLowerLimit, yHigherLimit])
    plt.savefig('graph_VoronoiTessellation.png',)
    
    print('Voronoi Tesselation Neighbors: ', delaunayDict[indexToGetNeighbors])
    
    return
### end voronoi tessellation and helper function
##########################################################################

##########################################################################
# Performance: Measures timing performance between fixedDistance and Voronoi functions
def performance():
    npPoints = np.array(getRandomPoints())
    
    #Fixed Distance
    startFD = time.time()
    fixedDistanceNN(npPoints)
    endFD = time.time()
    fixedDistanceTime = (endFD - startFD)
    
    #Voronoi Tessellation
    startVT = time.time()
    voronoiTessellation(npPoints)
    endVT = time.time()
    voronoiTessellationTime = (endVT - startVT)
    
    print("Time to execute FixedDistance : ", fixedDistanceTime, " seconds")
    print("Time to execute VoronoiTess:    ", voronoiTessellationTime, " seconds")
    
### end performance
##########################################################################

##########################################################################
# Main: Call functions
performance()