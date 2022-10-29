#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 16 17:36:08 2021

@author: arolle
"""

import numpy as np


# Draw points uniformly from the unit sphere in R^3,
# sampling less from a disc around the north pole,
# plus uniform outliers in the cube [-2, 2]^3.
#
# Input:
#
# num_points : number of points to draw in total.
#
# sample_weight : portion of points to draw from the sphere.
#                 We should have 0 <= sample_weight <= 1.
#                 Outliers in the cube [-2, 2]^3 will be drawn
#                 with probability 1 - sample_weight.
#
# r : the radius around the north pole for the disc.
#     We use Euclidean distance, not the geodesic distance.
#
# north_pole_weight : relative weight of the disc.
#                     We should have 0 <= north_pole_weight <= 1.
#                     E.g., if north_pole_weight = 0.5, then
#                     a point in the disc around the north pole is half
#                     as likely to appear as a point outside the disc.

def draw_points_from_sphere(num_points, sample_weight, r, north_pole_weight):
    
    N = [0, 0, 1]
    
    ######################################
    # choose number of samples from sphere
    ######################################
    
    num_samples = 0
    
    # initialize random generator
    rng = np.random.default_rng()
    
    # draw num_points values uniformly from [0, 1)
    values = rng.random(size=num_points, dtype=np.float64)
    
    for value in values:
        if value < sample_weight:
            num_samples += 1
            
    num_outliers = num_points - num_samples
    
    ##########################
    # draw samples from sphere
    ##########################
            
    sample = np.zeros(shape=(num_samples, 3), dtype=np.float64)
    
    current_num_samples = 0
    
    while current_num_samples < num_samples:
        
        # draw a vector from the cube [-1, 1]^3
        rng = np.random.default_rng()
        u = rng.random(size=3, dtype=np.float64)
        v = 2 * u - [1,1,1]
        magnitude = np.linalg.norm(v)
        
        if magnitude < 1 and magnitude != 0:
            
            p = v / magnitude
            
            if np.linalg.norm(v - N) >= r:
                
                sample[current_num_samples, :] = p
                current_num_samples += 1
                
            else:
                
                rng = np.random.default_rng()
                coin = rng.random(size=1, dtype=np.float64)
                
                if coin < north_pole_weight:
                    
                    sample[current_num_samples, :] = p
                    current_num_samples += 1
                    
    ###############
    # draw outliers
    ###############
                    
    outliers = np.zeros(shape=(num_outliers, 3), dtype=np.float64)
    
    for i in range(num_outliers):
        
        rng = np.random.default_rng()
        v = rng.random(size=3, dtype=np.float64)
        outliers[i, :] = 4 * v - [2,2,2]
        
    pointcloud = np.concatenate([sample, outliers], axis=0)
    
    return pointcloud


# Draw points from gaussians whose means live in an annulus in R^2,
# plus uniform outliers.
#
# Input:
#
# num_points : number of points to draw in total.
#
# sample_weight : portion of points to draw from the gaussians.
#                 We should have 0 <= sample_weight <= 1.
#                 Outliers in a square containing the annulus will be drawn
#                 with probability 1 - sample_weight.
#
# r : inner radius of the annulus.
#
# R : outer radius of the annulus (we should have r < R).
#
# num_gaussians : number of gaussians.
#
# var : the covariance matrix for each gaussian will be var times
#       the 2 x 2 identity matrix.

def draw_points_from_annulus(num_points, sample_weight, r, R, 
                             num_gaussians, var):
    
    ############################
    # draw means for gaussians
    ############################
    
    means = []
    
    current_num_means = 0
    
    while current_num_means < num_gaussians:
        
        # draw a vector from the square [-R, R]^2
        rng = np.random.default_rng()
        u = rng.random(size=2, dtype=np.float64)
        v = (2 * R) * u - [R, R]
        magnitude = np.linalg.norm(v)
        
        if magnitude > r and magnitude < R:
            means.append(v)
            current_num_means += 1
    
    #######################################################    
    # choose how many points to draw from each distribution
    #######################################################
    
    # k = number of distributions
    k = num_gaussians + 1
    
    gaussian_weight = sample_weight / num_gaussians
    
    # partial_sums: the j^th entry is the sum weight_0 + ... + weight_{j}    
    partial_sums = np.zeros(shape=k, dtype=np.float64)
    
    for j in range(k):
        
        if j == 0:
            partial_sums[j] = gaussian_weight
            
        if j > 0 and j < k-1:
            partial_sums[j] = partial_sums[j-1] + gaussian_weight
            
        if j == k-1:
            partial_sums[j] = 1.0
    
    # samples_per_distribution:
    # the j^th entry is the number of samples we draw from the j^th distribution 
    samples_per_distribution = np.zeros(shape=k, dtype=np.int64)
    
    # initialize random generator
    rng = np.random.default_rng()
    
    # draw num_samples values uniformly from [0, 1)
    values = rng.random(size=num_points, dtype=np.float64)
    
    for value in values:
        
        for j in range(k + 1):
            
            if value < partial_sums[j]:
                
                samples_per_distribution[j] += 1
                
                break
            
    #############################
    # draw samples from gaussians
    #############################
    
    cov_matrix = var * np.identity(2)
    
    samples = {}
    
    for j in range(k-1):
        
        # initialize random generator
        rng = np.random.default_rng()
        
        samples[j] = rng.multivariate_normal(mean=means[j],
                                             cov=cov_matrix,
                                             size=samples_per_distribution[j])
        
        
    ################################################
    # draw noise samples from a uniform distribution
    ################################################    
        
    rng = np.random.default_rng()
    
    outliers = rng.uniform(low=-R,
                           high=R,
                           size=(samples_per_distribution[k-1],2))
        
    samples[k-1] = outliers
    
    ############################################
    # concatenate samples from all distributions
    ############################################
    
    pointcloud = np.concatenate([samples[j] for j in range(k)], axis=0)
    
    return pointcloud