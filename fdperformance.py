import timeit

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
vor = Voronoi(randPoints)

fig = voronoi_plot_2d(vor, show_vertices=False)

ax = fig.add_subplot(1,1,1)

for i in range(0, len(randPoints)):
    ax.annotate(text = i, xy = (randPoints[i][0], randPoints[i][1]))

circ = plt.Circle(randPoints[0], radius=100, fill=False, edgecolor='green')
ax.add_patch(circ)

KDTree(randPoints).query_ball_point(randPoints[0], 100)
'''

#timeit statement
print('Fixed Distance Performance Analysis: ')
print (timeit.timeit(setup = mySetup, stmt = fixedDistanceRadius, number = 40))