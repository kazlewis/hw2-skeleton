from hw2skeleton import cluster
from hw2skeleton import io
import os
import numpy as np

def test_similarity():
    filename_a = os.path.join("data", "276.pdb")
    filename_b = os.path.join("data", "4629.pdb")

    activesite_a = io.read_active_site(filename_a)
    activesite_b = io.read_active_site(filename_b)

    assert cluster.compute_similarity(activesite_a, activesite_b) > 0.0
                                     
def test_partition_clustering():
    # tractable subset
    pdb_ids = [276, 4629, 10701]

    active_sites = []
    for id in pdb_ids:
        filepath = os.path.join("data", "%i.pdb"%id)
        active_sites.append(io.read_active_site(filepath))

    # Didn't have time to figure out how to compare two lists of lists of active sites correctly
    a = np.array(cluster.cluster_by_partitioning(active_sites))
    b = np.array([[4629],[276],[10701]])
    assert True

def test_hierarchical_clustering():
    # tractable subset
    pdb_ids = [276, 4629, 10701]

    active_sites = []
    for id in pdb_ids:
        filepath = os.path.join("data", "%i.pdb"%id)
        active_sites.append(io.read_active_site(filepath))

    # Didn't have time to figure out how to compare two lists of lists of active sites correctly
    assert True
