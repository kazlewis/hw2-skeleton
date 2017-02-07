#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 31 12:39:23 2017

@author: kazlewis
"""

from hw2skeleton import *
# Ensure all dependencies are imported

#test_active_sites = read_active_sites("testdata")
#site1 = test_active_sites[0]
#site2 = test_active_sites[1]
#site3 = test_active_sites[2]

active_sites = read_active_sites("testdata")


#print ('num_residues in site1 is ' + str(len(site1.residues)))

p_clusters = cluster_by_partitioning(active_sites)

print(p_clusters)



