from .utils import *

def compute_similarity_bymatrix(info_a,info_b):
    """
    Compute the similarity between two given ActiveSite instances.

    Input: two instances of activesite_info matrices determined by activesite_info in utils.py
    Output: the similarity between them (a floating point number). The smaller this number, the more similar the active sites.
    
    Calculate similarity based on the following metrics, each normalized to the net effects on the center of the active site (distance wise)
        volume of active site pocket 
        net charge of residues acting on active site
        number of hydrophobic residues acting on active site
        number of polar residues acting on active site
    """
    
    similarity = 0.0
    vol_diff = info_a[0] - info_b[0]
    charge_diff = info_a[1] - info_b[1]
    polarity_diff = info_a[2] - info_b[2]
    hydrophobicity_diff = info_a[3] - info_b[3]

    # Use Euclidian distance to attempt to minmiize the maximum contribution of any particular factor    
    similarity = math.sqrt(math.pow(vol_diff,2) + math.pow(charge_diff,2) + math.pow(polarity_diff,2) + math.pow(hydrophobicity_diff,2))

    return similarity


def compute_similarity(site_a, site_b):
    """
    Wrapper function for compute_similariy_bymatrix()
    """
    info_a = activesite_info(site_a)
    info_b = activesite_info(site_b)
    
    similarity = compute_similarity_bymatrix(info_a,info_b)
                    
    return similarity


def cluster_by_partitioning(active_sites):
    """
    Cluster a given set of ActiveSite instances using a partitioning method.

    Input: a list of ActiveSite instances
    Output: a clustering of ActiveSite instances
            (this is really a list of clusters, each of which is list of
            ActiveSite instances)
    """
    
    # Define an initial number of clusters and maximum number of iterations
    K = 4
    max_iterations = 500
    
    #random.seed(1)
    initial_clusters = random.sample(range(len(active_sites)), K)
    cluster_centers = [[]]
    cluster_assignments = np.zeros(shape = (len(active_sites)))
    
    # Initialize clusters 
    # cluster_centers keeps track of which 'active sites' (real for first assignment, averages of the cluster for i + 1) are cluster centers
    # cluster_assignments is an array that keeps track of which cluster each active site belongs to (mirror of active_sites in terms of index)
    
    # Initialize initial cluster centers
    for i in range(len(initial_clusters)):
        cluster_assignments[initial_clusters[i]] = i + 1
        cluster_centers.append(activesite_info(active_sites[initial_clusters[i]]))     
        
    # remove the blank first index – haven't found a more elegant solution than this yet
    cluster_centers.pop(0)
    
    # Set up parameters to check if a solution has been found; keep track of iterations    
    converged = False
    cluster_assignments_prev = cluster_assignments.copy()
    iterations = 0
    while not converged:
        # Assign each active site to a cluster by looping through cluster_assignments
        # index is the position in active_sites (iterate through the entire list)
        for index in range(len(cluster_assignments)):
            min_diff = math.inf
            # i is a potential cluster center, check to see which cluster the current active site most closely matches
            for i in range(len(cluster_centers)):
                similarity = compute_similarity_bymatrix(cluster_centers[i],activesite_info(active_sites[index]))
                if similarity < min_diff:
                    min_diff = similarity
                    cluster_assignments[index] = i + 1
        
        # After each active site has been assigned to a cluster, loop over the cluster centers to make new 'cluster centers'
        for i in range(len(cluster_centers)):
            # Prep for new cluster center
            new_center = np.array([0.0, 0.0, 0.0, 0.0], ndmin=1)
            cluster_sites = 0
            # Look at assignment of all active sites to build new 'cluster center' 
            for j in range(len(cluster_assignments)):
                if cluster_assignments[j] == i + 1:
                    new_center += activesite_info(active_sites[j])
                    cluster_sites += 1
            # Update the cluster centers
            new_center = new_center / cluster_sites
            cluster_centers[i] = new_center
                
        # Increase iteration count
        iterations += 1
        
        # Check to see if we've converged on a solution, otherwise update 'prev' and repeat
        if np.array_equiv(cluster_assignments_prev, cluster_assignments):
            converged = True  
        elif iterations > max_iterations:    
            converged = True
        else:
            cluster_assignments_prev = cluster_assignments.copy()
        
    # Lastly, create a matrix with each of the active sites arranged into their respective clusters
    # Initialize list of lists with dimensons to support the number of clusters K
    cluster_array = [ [] for i in range(K)]
        
    # Iterate through all active sites, looking at their cluster assignments    
    for i in range(len(cluster_assignments)):
        # Iterate through all possible clusters, and add active sites to their respective clusters
        for j in range(K):
            if j + 1 == cluster_assignments[i]:
                cluster_array[j].append(active_sites[i])
    
    
    print("the algorithm took " + str(iterations) + " iteration(s) to converge on a solution")
    return cluster_array


def cluster_hierarchically(active_sites):
    """
    Cluster the given set of ActiveSite instances using a hierarchical algorithm.                                                                  #

    Input: a list of ActiveSite instances
    Output: a list of clusterings
            (each clustering is a list of lists of Sequence objects)
    """

    # Fill in your code here!

    return []
