import time
t1 = time.time()
import argparse, os, sys
import pandas as pd
import numpy as np
import random
from astropy import units as u
from astropy.coordinates import SkyCoord
import math
from sympy import symbols, Eq, solve
from haversine import haversine, Unit

output_dir = "/home/simon/condor/output/"

parser = argparse.ArgumentParser()
parser.add_argument('-i', type=int, required=True, 
    help="job number")

args = parser.parse_args()
i = args.i

# The original catalog containing 1177 blazars:
df = pd.read_csv("/home/simon/condor/blazar_clean.csv")

# First, a set of angles and radii for the shifting of the coordinates is created:
c = np.random.uniform(0,100, 10000)
ra1 = np.random.choice([-1,1], 10000)*np.random.uniform(0,np.sqrt(c))
dec1 = np.random.choice([-1,1], 10000)*np.sqrt(-ra1**2+c)
ra2 = np.random.choice([-1,1], 10000)*np.random.uniform(0,np.sqrt(c))
dec2 = np.random.choice([-1,1], 10000)*np.sqrt(-ra2**2+c)

matches_list = [] # This list is the output of this .py file and will contain lists consisting of 21 entries (for each combination of L and r), each list represents a new generated catalog. After generating one catalog the matching algorithm will do the matching for all 21 settings, save the result in matches_list and move to generating a new catalog.
    
def generate_catalog():
    # Now to the actual scrambling procedure:
    catalog = ([ra+np.random.choice([np.random.choice(ra1),np.random.choice(ra2)]) for ra in df["RA"]], [dec+np.random.choice([np.random.choice(dec1),np.random.choice(dec2)]) for dec in df["DEC"]]) # Here the original catalog is scrambled, but the scramling process shifts some blazars outside the desired boundaries.
    # Shift all blazars outside the outer boundaries back inside:
    for j in range(1177):
        while np.logical_or(catalog[1][j] > -5, np.logical_or(catalog[1][j] < -85, np.logical_or(catalog[0][j] > 360, catalog[0][j] < 0))):
            catalog[0][j] = df["RA"][j] + np.random.choice([np.random.choice(ra1),np.random.choice(ra2)])
            catalog[1][j] = df["DEC"][j] + np.random.choice([np.random.choice(dec1),np.random.choice(dec2)])
    # Shift all blazars near the galactic plane into the boundaries:
    catalog_gal = SkyCoord(ra=catalog[0]*u.degree, dec=catalog[1]*u.degree, frame="icrs").galactic
    for j in range(1177):
        while np.logical_or(0 < catalog_gal.b.value[j] < 10, -10 < catalog_gal.b.value[j] < 0):
            if np.logical_and(df["DEC"][j] + np.random.choice([np.random.choice(dec1),np.random.choice(dec2)]) < -5, df["DEC"][j] + np.random.choice([np.random.choice(dec1),np.random.choice(dec2)]) > -85):
                catalog[0][j] = df["RA"][j] + np.random.choice([np.random.choice(ra1),np.random.choice(ra2)])
                catalog[1][j] = df["DEC"][j] + np.random.choice([np.random.choice(dec1),np.random.choice(dec2)])
                catalog_gal.b.value[j] = SkyCoord(ra=catalog[0][j]*u.degree, dec=catalog[1][j]*u.degree, frame="icrs").galactic.b.value
    return catalog


# The matching algorithm as a defined function: 
def match_catalog(catalog, L, r):
        if L == 3.5:
            df_centers = pd.read_csv("/home/simon/condor/neutrino_hotspots_35.csv").to_numpy()
        if L == 4.0:
            df_centers = pd.read_csv("/home/simon/condor/neutrino_hotspots_40.csv").to_numpy()
        if L == 4.5:
            df_centers = pd.read_csv("/home/simon/condor/neutrino_hotspots_45.csv").to_numpy()
        match_count = []
        for n in range(df_centers.shape[0]):
            for j in range(len(catalog[0])):
                if haversine([df_centers[n,2], df_centers[n,1]], [catalog[1][j], catalog[0][j]], unit=Unit.DEGREES, normalize=True) < r:
                    match_count.append(n)
        return len(list(set(match_count)))
    
    
# For each setting of L and r, the number of found matches are saved in a list. This list is then appended to the output matches_list.  
matches = []
for L in [3.5,4,4.5]:
    for r in [0.4,0.45,0.5,0.55,0.6,0.65,0.7]:
        matches.append(match_catalog(generate_catalog(), L, r))
matches_list.append(matches)
    
outfname = "{0:s}/output{1:n}.txt".format(output_dir, i)

t2 = time.time()

f = open(outfname, "w")
f.write(str(matches_list)+"\n")
f.write(str("Runtime: "+t2-t1))
f.close()
