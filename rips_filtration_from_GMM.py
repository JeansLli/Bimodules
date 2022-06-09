#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from sklearn.metrics.pairwise import euclidean_distances
from sklearn.neighbors import KernelDensity
from sklearn.datasets import make_blobs
import matplotlib.pyplot as plt
import time

from function_rips import function_rips,function_rips_sparse
from generating_pointclouds import draw_points_from_sphere
from generating_pointclouds import draw_points_from_annulus

##########################################################
# Draw samples from a sphere, and compute the scc2020 file
##########################################################

from sys import argv

argc=len(argv)

if(argc>=2):
    no_points=int(argv[1])
else:
    no_points=1000

if(argc>=3):
    no_simplices=int(argv[2])
else:
    no_simplices=200000

if(argc>=4):
    file_name=argv[3]
else:
    file_name = 'function_rips_GMM.txt' 

print('number of points:' + str(no_points))
print('number of simplices:' + str(no_simplices))
print('output filename:' + file_name)


start = time.time()

pointcloud, y_true = make_blobs(n_samples=no_points, centers=[[0,0],[1,0],[0,1],[1,1]],
                                cluster_std=0.20, n_features=2, random_state=4)
plt.scatter(pointcloud[:, 0], pointcloud[:, 1])
plt.show()

end = time.time()
runtime = end - start
print('Sampling: ' + str(runtime))

band = 0.2

kde = KernelDensity(kernel='linear', bandwidth=band).fit(pointcloud)
function_vals = np.exp(kde.score_samples(pointcloud))

end = time.time()
runtime = end - start
print('Generating densities: ' + str(runtime))

from scipy.spatial import cKDTree
kd_tree = cKDTree(pointcloud)

end = time.time()
runtime = end - start
print("KD constructir" + str(runtime))

num_pairs=0
r=0.001

while 2*num_pairs < no_points*(no_points-1) and num_pairs < no_simplices:
    pairs = list(kd_tree.query_pairs(r=r))
    num_pairs = len(pairs)
    print("r, num_pairs",r,num_pairs)
    print("count neighbors", kd_tree.count_neighbors(kd_tree,r=r))
    r = 2 * r



end = time.time()
runtime = end - start
print('KD tree: ' + str(runtime))

from scipy.spatial import distance
from math import sqrt

def dist(p,q):
    return sqrt((p[0]-q[0])**2 + (p[1]-q[1])**2)
#def dist(p,q):
#    return sqrt((p[0]-q[0])**2 + (p[1]-q[1])**2+(p[2]-q[2])**2)

distances=np.empty((len(pairs),3))
print('Memory: ' + str(runtime))
for i in range(len(pairs)):
    distances[i,0]=pairs[i][0]
    distances[i,1]=pairs[i][1]
    p=pointcloud[pairs[i][0]]
    q=pointcloud[pairs[i][1]]
    distances[i,2]=dist(p,q)
    
print(distances.shape)

end = time.time()
runtime = end - start
print('Distance entries: ' + str(runtime))

if True:


    function_rips_sparse(no_points, distances, function_vals=function_vals, d=2,
                         file_name=file_name, 
                  give_radius=False, r=0.1, 
                  give_num_simplices=True,
                  num_simplices=no_simplices)
else:
    function_rips(dm, function_vals=function_vals, d=2, 
                  file_name=file_name, 
                  give_radius=False, r=0.1, 
                  give_num_simplices=True,
                  num_simplices=no_simplices)


end = time.time()
runtime = end - start
print('Function rips: ' + str(runtime))

