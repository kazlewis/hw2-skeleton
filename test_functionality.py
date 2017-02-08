#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 31 12:39:23 2017

@author: kazlewis
"""

# Ensure all dependencies are imported from __init__.py (modifed)
from hw2skeleton import *

# Read in active sites
active_sites = read_active_sites("data")

# Test partitioning clustering
p_clusters = cluster_by_partitioning(active_sites)
p_sum_dist_Q = sum_dist_Q(p_clusters)
print("the clusters resulting from K-means partition clustering are ")
print(p_clusters)
print("After the final run, Q was " + str(p_sum_dist_Q))

# Test hierarchical clustering
h_clusters = cluster_hierarchically(active_sites)
h_sum_dist_Q = sum_dist_Q(h_clusters)
print("the clusters resulting from nearest-centroid hierarchical clustering are ")
print(h_clusters)
print("After the final run, Q was " + str(h_sum_dist_Q))

# Test optimal number of clusters for K-means. This will break unless cluster_by_partitioning is adjusted in cluster.py
# to accept an initial K value as input.
Q = math.inf
Q_ratio_prev = 0
K = 2
for i in range(len(active_sites)):
    p_clusters = cluster_by_partitioning(active_sites, K)
    Q_prev = Q
    Q = sum_dist_Q(p_clusters)
    Q_ratio = Q/Q_prev
    if Q_ratio < Q_ratio_prev or len(p_clusters) == 1:  
        break
    else:
        Q_ratio_prev = Q_ratio
    K += 1

print(str(K) + " clusters were found to be optimal with k-means.")

