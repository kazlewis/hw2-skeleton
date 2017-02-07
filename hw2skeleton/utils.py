# Some utility classes to represent a PDB structure
# Also included is a python script to calculate the area of a 3d polygon
# Script taken from http://stackoverflow.com/questions/12642256/python-find-area-of-polygon-from-xyz-coordinates
# Lastly included is a function to generate a matrix of information about a given active site

import numpy as np
import math
import random

class Atom:
    """
    A simple class for an amino acid residue
    """

    def __init__(self, type):
        self.type = type
        self.coords = (0.0, 0.0, 0.0)

    # Overload the __repr__ operator to make printing simpler.
    def __repr__(self):
        return self.type

class Residue:
    """
    A simple class for an amino acid residue
    """

    def __init__(self, type, number):
        self.type = type
        self.number = number
        self.atoms = []

    # Overload the __repr__ operator to make printing simpler.
    def __repr__(self):
        return "{0} {1}".format(self.type, self.number)

class ActiveSite:
    """
    A simple class for an active site
    """

    def __init__(self, name):
        self.name = name
        self.residues = []

    # Overload the __repr__ operator to make printing simpler.
    def __repr__(self):
        return self.name

# Calculating the area of a 3d polygon
# unit normal vector of plane defined by points a, b, and c
def unit_normal(a, b, c):
    x = np.linalg.det([[1,a[1],a[2]],
         [1,b[1],b[2]],
         [1,c[1],c[2]]])
    y = np.linalg.det([[a[0],1,a[2]],
         [b[0],1,b[2]],
         [c[0],1,c[2]]])
    z = np.linalg.det([[a[0],a[1],1],
         [b[0],b[1],1],
         [c[0],c[1],1]])
    magnitude = (x**2 + y**2 + z**2)**.5
    return (x/magnitude, y/magnitude, z/magnitude)

#area of polygon poly
def poly_area(poly):
    if len(poly) < 3: # not a plane - no area
        return 0
    total = [0, 0, 0]
    N = len(poly)
    for i in range(N):
        vi1 = poly[i]
        vi2 = poly[(i+1) % N]
        prod = np.cross(vi1, vi2)
        total[0] += prod[0]
        total[1] += prod[1]
        total[2] += prod[2]
    result = np.dot(total, unit_normal(poly[0], poly[1], poly[2]))
    return abs(result/2)

#def activesite_vol(site):
#    if type(site) != ActiveSite:
#        raise IOError("Not a valid active site.")
#    residue_vol_array = np.array([0.0, 0.0, 0.0], ndmin=2)
#    for index, all_residues in enumerate(site.residues):
#        #print(index)
#        #print(all_residues.atoms)
#        total_residue_coord = np.array([0.0, 0.0, 0.0], ndmin =2)
#        for all_atoms in all_residues.atoms:
#            total_residue_coord += np.array(all_atoms.coords)
#        avg_residue_coord = total_residue_coord/(len(all_residues.atoms))
#        residue_vol_array = np.concatenate((residue_vol_array,avg_residue_coord), axis=0)
#    residue_vol_array = np.delete(residue_vol_array, 0, axis=0)
#    return poly_area(residue_vol_array)

def activesite_info(site):
    if type(site) != ActiveSite:
        raise IOError("Not a valid active site.")
        
    # Create a matrix that stores the following info
    #   approximate volume of the active site    
    #   weighted effects of charged residues on active site centroid    
    #   weighted effects of polar residues on active site centroid
    #   weighted effects of hydrophobic residues on active site centroid
    
    info_matrix = np.array([0.0, 0, 0, 0])
    
    # Get a matrix containing the (average) coordinate of each of the residues within the average site
    residue_coord_array = np.array([0.0, 0.0, 0.0], ndmin=2)
    for index, residues in enumerate(site.residues):
        total_residue_coord = np.array([0.0, 0.0, 0.0], ndmin =2)
        for all_atoms in residues.atoms:
            total_residue_coord += np.array(all_atoms.coords)
        avg_residue_coord = total_residue_coord/(len(residues.atoms))
        residue_coord_array = np.concatenate((residue_coord_array,avg_residue_coord), axis=0)
    residue_coord_array = np.delete(residue_coord_array, 0, axis=0)
    # Convert matrix to volume using script pulled from the internet
    info_matrix[0] = poly_area(residue_coord_array)

    # Determine the centroid (average) of the active site
    total_centroid_coord = np.array([0.0, 0.0, 0.0], ndmin =1)
    for index, residues in enumerate(residue_coord_array):
        total_centroid_coord += residues
    activesite_centroid = total_centroid_coord/(len(residue_coord_array))
    
    # Iterate through the residues and assign charge, polar, and hydrophobic contributions to the centroid of the active site
    activesite_charge = 0.0
    activesite_polarity = 0.0
    activesite_hydrophobicity = 0.0
    for index, residue in enumerate(site.residues):
        x_dist = residue_coord_array[index][0] - activesite_centroid[0]
        y_dist = residue_coord_array[index][1] - activesite_centroid[1]
        z_dist = residue_coord_array[index][2] - activesite_centroid[2]
        distance_to_centroid = math.sqrt(x_dist ** 2 + y_dist ** 2 + z_dist ** 2)
        if residue.type in ("HIS","ARG","LYS"):
            activesite_charge += 1 * distance_to_centroid
        elif residue.type in ("ASP","GLU"):
            activesite_charge -= 1 * distance_to_centroid
        if residue.type in ("HIS","ARG","LYS","SER","THR","ASN","GLN","ASP","GLU"):
            activesite_polarity += 1 * distance_to_centroid
        if residue.type in ("ALA","VAL","ILE","LEU","MET","PHE","TYR","TRP"):
            activesite_hydrophobicity += 1 + distance_to_centroid    
    info_matrix[1] = activesite_charge
    info_matrix[2] = activesite_polarity
    info_matrix[3] = activesite_hydrophobicity
               
    # Normalize values in info_matrix
    m = 0.25 * (info_matrix[0] + info_matrix[1] + info_matrix[2] + info_matrix[3])
    S = 0.25 * (abs(info_matrix[0] - m) + abs(info_matrix[1] - m) + abs(info_matrix[2] - m) + abs(info_matrix[3] - m))
    info_matrix[0] = (info_matrix[0] - m)/S
    info_matrix[1] = (info_matrix[1] - m)/S
    info_matrix[2] = (info_matrix[2] - m)/S
    info_matrix[3] = (info_matrix[3] - m)/S            
               
    return info_matrix




    
#def num_residues(site):
#    if type(site) != ActiveSite:
#        raise IOError("Not a valid active site")
#    return len(site.residues)
#
#def site_charge(site):
#    charge = 0
#    if type(site) != ActiveSite:
#        raise IOError("Not a valid active site")
#    for residue in site.residues:
#        if residue.type in ("HIS","ARG","LYS"):
#            charge += 1
#        elif residue.type in ("ASP","GLU"):
#            charge -= 1
#    return charge
#
#def site_polarity(site):
#    polarity = 0
#    if type(site) != ActiveSite:
#        raise IOError("Not a valid active site")
#    for residue in site.residues:
#        if residue.type in ("HIS","ARG","LYS","SER","THR","ASN","GLN","ASP","GLU"):
#            polarity += 1
#    return polarity
#
#def site_hydrophobicity(site):
#    hydrophobicity = 0
#    if type(site) != ActiveSite:
#        raise IOError("Not a valid active site")
#    for residue in site.residues:
#        if residue.type in ("ALA","VAL","ILE","LEU","MET","PHE","TYR","TRP"):
#            hydrophobicity += 1
#    return hydrophobicity
#
#
#
#
#
#






