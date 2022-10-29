#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 15 13:14:45 2021

@author: arolle
"""

import numpy as np
from scipy.special import comb
from itertools import combinations

def function_rips_sparse(num_points, distances, function_vals, d, 
                  file_name, 
                  give_radius=False, r=0.1, 
                  give_num_simplices=False,
                  num_simplices=1000):

    ########################################################################
    # check that exactly one of give_radius and give_num_simplices are true.
    ########################################################################
    
    if give_radius == True and give_num_simplices == True:
        print('exactly one of give_radius and give_num_simplices must be true!')
        return None
    
    if give_radius == False and give_num_simplices == False:
        print('exactly one of give_radius and give_num_simplices must be true!')
        return None

    print(distances)

    # sort the array of distances from smallest distance to largest
    sort_distances = distances[distances[:, 2].argsort()]

    num_distances=len(sort_distances)
    print("Sorted",num_distances)
    

    neighbors=[set({}) for i in range(num_points)]

    

    # (simplex, grade)
    zero_simplices = [(tuple([i]), (0, -function_vals[i])) for i in range(num_points)]
    one_simplices = []  
    two_simplices = []
    three_simplices = []
    
    
    if give_num_simplices == True: # Yes, it's our case
        
        # since all vertices are born at radius zero,
        # we implicitly add all vertices immediately.
        current_num_simplices = num_points
        
        for i in range(num_distances):
            
            x = int(sort_distances[i, 0])
            y = int(sort_distances[i, 1])
            
            if x != y and current_num_simplices < num_simplices:
                
                # compute the grade of the edge
                gr = grade(simplex=(x,y), 
                           function_vals=function_vals, 
                           r=sort_distances[i,2])
                
                # add the edge
                edge = [x,y]
                edge.sort()
                one_simplices.append((tuple(edge), gr))
                current_num_simplices += 1

                common_neighbors = neighbors[x].intersection(neighbors[y])
                neighbors[x].add(y)
                neighbors[y].add(x)
                #print("Common neighbors:",common_neighbors)
                
                # add two-simplices with this edge
                if d > 1:
                    
                    for p in common_neighbors:
                            
                        # compute the grade of the two-simplex
                        gr = grade(simplex=(p,x,y), 
                                   function_vals=function_vals, 
                                   r=sort_distances[i,2])
                        
                        triangle = [p,x,y]
                        triangle.sort()
                        two_simplices.append((tuple(triangle), gr))
                        current_num_simplices += 1
                    
                        # add three-simplices that have
                        # this triangle on the boundary
                        if d > 2:
                            common_neighbors_with_p=common_neighbors.intersection(neighbors[p])
                            
                            # Find points q such that the edges
                            # (q, p), (q, x) and (q, y) were already added.
                            # If q < p, then we already added this 
                            # three-simplex
                            for q in common_neighbors_with_p:
                                if p>q:
                                    continue # avoids double-creation
                                
                                # compute the grade of the edge
                                gr = grade(simplex=(q,p,x,y), 
                                           function_vals=function_vals, 
                                           r=sort_distances[i,2])
                                
                                tetra = [q,p,x,y]
                                tetra.sort()
                                three_simplices.append((tuple(tetra), gr))
                                current_num_simplices += 1
                                        
            if current_num_simplices >= num_simplices:
                break
            
    if give_radius == True:
        
        current_num_simplices = num_points
        
        for i in range(num_distances):
            
            x = int(sort_distances[i, 0])
            y = int(sort_distances[i, 1])
            
            if x != y and sort_distances[i, 2] <= r:
                
                # compute the grade of the edge
                gr = grade(simplex=(x,y), 
                           function_vals=function_vals, 
                           r=sort_distances[i,2])
                
                # add the edge
                edge = [x,y]
                edge.sort()
                one_simplices.append((tuple(edge), gr))
                current_num_simplices += 1

                common_neighbors = neighbors[x].intersection(neighbors[y])
                neighbors[x].add(y)
                neighbors[y].add(x)
                #print("Common neighbors:",common_neighbors)
                
                # add two-simplices with this edge
                if d > 1:
                    
                    for p in common_neighbors:
                            
                        # compute the grade of the two-simplex
                        gr = grade(simplex=(p,x,y), 
                                   function_vals=function_vals, 
                                   r=sort_distances[i,2])
                        
                        triangle = [p,x,y]
                        triangle.sort()
                        two_simplices.append((tuple(triangle), gr))
                        current_num_simplices += 1
                    
                        # add three-simplices that have
                        # this triangle on the boundary
                        if d > 2:
                            common_neighbors_with_p=common_neighbors.intersection(neighbors[p])
                            
                            # Find points q such that the edges
                            # (q, p), (q, x) and (q, y) were already added.
                            # If q < p, then we already added this 
                            # three-simplex
                            for q in common_neighbors_with_p:
                                if p>q:
                                    continue # avoids double-creation
                                
                                # compute the grade of the edge
                                gr = grade(simplex=(q,p,x,y), 
                                           function_vals=function_vals, 
                                           r=sort_distances[i,2])
                                
                                tetra = [q,p,x,y]
                                tetra.sort()
                                three_simplices.append((tuple(tetra), gr))
                                current_num_simplices += 1
                                        
            if sort_distances[i, 2] > r:
                break        
            
            
    ###############################
    # print scc2020 file ##########
    ###############################

    print("3-simplices: ",100.0*len(three_simplices)/current_num_simplices)
    print("2-simplices: ",100.0*len(two_simplices)/current_num_simplices)
    print("1-simplices: ",100.0*len(one_simplices)/current_num_simplices)
    print("0-simplices: ",100.0*len(zero_simplices)/current_num_simplices)

    print("Printing")
    
    # sort simplices
    # by colex order on the grades
    zero_simplices.sort(key=lambda x : (x[1][1], x[1][0]))
    one_simplices.sort(key=lambda x : (x[1][1], x[1][0]))
    two_simplices.sort(key=lambda x : (x[1][1], x[1][0]))
    three_simplices.sort(key=lambda x : (x[1][1], x[1][0]))
    
    # record the indices of the simplices using a dictionary
    # such that given a simplex (as a tuple)
    # the dictionary returns the index of the simplex
    # in the ordered list of simplices
    zero_simplices_indices = {}
    for counter, item in enumerate(zero_simplices):
        zero_simplices_indices[item[0]] = counter
    
    one_simplices_indices = {}
    for counter, item in enumerate(one_simplices):
        one_simplices_indices[item[0]] = counter
        
    two_simplices_indices = {}
    for counter, item in enumerate(two_simplices):
        two_simplices_indices[item[0]] = counter
        
    # start printing file
    f = open(file_name, 'w')
    print('scc2020', file=f)
    
    # number of persistence parameters
    print('2', file=f)
    
    # sizes of generating sets
    if d == 3:
        m1 = len(three_simplices)
        m2 = len(two_simplices)
        m3 = len(one_simplices)
        m4 = num_points
        line = str(m1) + ' ' + str(m2) + ' ' + str(m3) + ' ' + str(m4) + ' 0'
        print(line, file=f)
        
    if d == 2:
        m1 = len(two_simplices)
        m2 = len(one_simplices)
        m3 = num_points
        line = str(m1) + ' ' + str(m2) + ' ' + str(m3) + ' 0'
        print(line, file=f)

    
        
    # second parameter is contravariant
    # comment out because neg value is currently taken
    #line = '--reverse 2'
    #print(line, file=f)
    
    # boundary matrices
    if d == 3:
        for i in range(3):
            
            if i == 0:
                
                current_simplices = three_simplices
                boundary_indices = two_simplices_indices
                
            if i == 1:
                
                current_simplices = two_simplices
                boundary_indices = one_simplices_indices
                
            if i == 2:
                
                current_simplices = one_simplices
                boundary_indices = zero_simplices_indices     
                
            for item in current_simplices:

                
                simplex = item[0]
                gr = item[1]
                
                line = str(gr[0]) + ' ' + str(gr[1]) + ' ;'
                
                col = []
            
                for j in range(len(simplex)):
                    face = [simplex[k] for k in range(len(simplex)) if k != j]
                    row_index = boundary_indices[tuple(face)]
                    col.append(row_index)
                
                for row_index in col:
                    line += ' '
                    line += str(row_index)
                    
                print(line, file=f)

                
    if d == 2:
        for i in range(2):
                
            if i == 0:
                
                current_simplices = two_simplices
                boundary_indices = one_simplices_indices
                
            if i == 1:
                
                current_simplices = one_simplices
                boundary_indices = zero_simplices_indices     
                
            for item in current_simplices:
                
                simplex = item[0]
                gr = item[1]
                
                line = str(gr[0]) + ' ' + str(gr[1]) + ' ;'
                
                col = []
            
                for j in range(len(simplex)):
                    face = [simplex[k] for k in range(len(simplex)) if k != j]
                    row_index = boundary_indices[tuple(face)]
                    col.append(row_index)
                
                for row_index in col:
                    line += ' '
                    line += str(row_index)
                    
                print(line, file=f)
                
    # print grades of vertices
    for item in zero_simplices:
        
        gr = item[1]
        line = str(gr[0]) + ' ' + str(gr[1])
        print(line, file=f)
        
    f.close()



# Compute an scc2020 file representing the chain complex
# F^3 -> F^2 -> F^1 -> F^0
# where F^i is the free bipersistence module 
# generated by the i-simplices of function-Rips.
#
# The base field is Z/2.
#
# Input:
#
# dm : a distance matrix.
#
# function_vals : a list or array of function values for each point.
#
# d : maximum dimension of simplices to consider.
#     Can be 2, or 3.
#     E.g., if d=2, we only compute F^2 -> F^1 -> F^0.
#
# file_name : name of the output file.
#
# give_radius : if True, we compute function-Rips up to a threshold radius.
#
# r : the threshold radius.
#
# give_num_simplices : if True, we compute function-Rips up to a threshold
#                      number of simplices.
#
# num_simplices : the threshold number of simplices.
#
# Note: exactly one of give_radius and give_num_simplices must be true.

def function_rips(dm, function_vals, d, 
                  file_name, 
                  give_radius=False, r=0.1, 
                  give_num_simplices=False,
                  num_simplices=1000):
    
    
    
    ########################################
    # get sorted array of pairwise distances
    ########################################
    
    num_points = dm.shape[0]
    
    # number of pairwise distances is num_points choose 2, plus num_points
    num_distances = comb(num_points, 2, exact=True) + num_points
    
    # create array of pairwise distances
    distances = np.zeros(shape=(num_distances, 3), dtype=np.float64)
    
    index = 0
    for p in range(num_points):
        for q in range(num_points):
            if p <= q:
                
                distances[index, :] = [p, q, dm[p,q]]
                index += 1
                
    
    return function_rips_sparse(num_points, distances, function_vals, d, 
                  file_name, 
                  give_radius, r, 
                  give_num_simplices,
                  num_simplices)
    

def grade(simplex, function_vals, r):

    k = function_vals[simplex[0]]
    for vertex in simplex:
        k=min(k,function_vals[vertex])
    
    return (r, -k)
