import argparse, os, sys, math
import pandas as pd
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from haversine import haversine, Unit

parser = argparse.ArgumentParser()
parser.add_argument('-i', type=int, required=True, help="job number")
i = parser.parse_args().i

# The original catalog containing 1177 blazars:
df = pd.read_csv("/home/simon/condor/blazar_clean.csv")

# First, a set of angles and radii for the shifting of the coordinates is created:
num_of_stifting_positions = 1000000
shifting_distance = np.sqrt(np.random.uniform(0, 100, num_of_stifting_positions))
angle = np.random.uniform(0, 360, num_of_stifting_positions)
ra_shift = np.sin(angle * (np.pi / 180.0)) * shifting_distance
dec_shift = np.cos(angle * (np.pi / 180.0)) * shifting_distance

# Now to the actual scrambling procedure:
time_gen = []
def generate_catalog():
    catalog = ([ra+np.random.choice(ra_shift) for ra in df["RA"]], [dec+np.random.choice(dec_shift) for dec in df["DEC"]]) # The original catalog is scrambled, but the scramling process shifts some blazars outside the desired boundaries.
    # Shift all blazars outside the outer boundaries back inside:
    for j in range(1177):
        while np.logical_or(catalog[1][j] > -5, np.logical_or(catalog[1][j] < -85, np.logical_or(catalog[0][j] > 360, catalog[0][j] < 0))):
            catalog[0][j] = df["RA"][j] + np.random.choice(ra_shift)
            catalog[1][j] = df["DEC"][j] + np.random.choice(dec_shift)
    # Shift all blazars near the galactic plane into the boundaries:
    catalog_gal = SkyCoord(ra=catalog[0]*u.degree, dec=catalog[1]*u.degree, frame="icrs").galactic
    for j in range(1177):
        while np.logical_or(0 < catalog_gal.b.value[j] < 10, -10 < catalog_gal.b.value[j] < 0):
            dec_gal = np.random.choice(dec_shift)
            # One has to be careful not to shift blazars back outside the outer boundaries, therefore check:
            if np.logical_and(df["DEC"][j] + dec_gal < -5, df["DEC"][j] + dec_gal > -85):
                ra_gal = np.random.choice(ra_shift)
                catalog[0][j] = df["RA"][j] + ra_gal
                catalog[1][j] = df["DEC"][j] + dec_gal
                catalog_gal.b.value[j] = SkyCoord(ra=catalog[0][j]*u.degree, dec=catalog[1][j]*u.degree, frame="icrs").galactic.b.value
    return catalog


# The matching algorithm as a function. Matches one catalog for all possible combinations of L and r:
time_match = []
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
    return len(set(match_count))
    

# For each setting of L and r, the number of found matches are saved in a list.
result = ""
num_of_catalogs = 10
for num in range(num_of_catalogs):
    matches = []
    for L in [3.5,4,4.5]:
        for r in [0.4,0.45,0.5,0.55,0.6,0.65,0.7]:
            matches.append(match_catalog(generate_catalog(), L, r))
    result += str(matches).translate({ord('['): None,ord(']'): None})+"\n"
    
output_dir = "/home/simon/condor/output/"
outfname = "{0:s}/output{1:n}.txt".format(output_dir, i)
f = open(outfname, "w")
f.write("L1r1,L1r2,L1r3,L1r4,L1r5,L1r6,L1r7,L2r1,L2r2,L2r3,L2r4,L2r5,L2r6,L2r7,L3r1,L3r2,L3r3,L3r4,L3r5,L3r6,L3r7"+"\n")
f.write(str(result))
f.close()