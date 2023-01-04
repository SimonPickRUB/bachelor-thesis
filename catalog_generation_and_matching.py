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

# First, a list of angles and radii for the shifting of the coordinates is created:
x=symbols('x') 
y=symbols('y') 
solutions = []
for k in range(600): # The -RA and +RA direction have to be treated separately, thus the total number of shifting variables is double the argument of the range function.
    # Shifting in -RA direction
    equation_1 = Eq((x**2 + y**2), random.uniform(0,100))
    equation_2 = Eq((y/x), math.tan(math.radians(random.uniform(0,360))))
    solutions.append(solve((equation_1, equation_2), (x, y))[0])
    # Shifting in +RA direction
    equation_1 = Eq((x**2 + y**2), random.uniform(0,100))
    equation_2 = Eq((y/x), math.tan(math.radians(random.uniform(0,360))))
    solutions.append(solve((equation_1, equation_2), (x, y))[1])

# Now to the actual scrambling procedure:
matches_list = [] # This list is the output of this .py file and will contain lists consisting of 21 entries (for each combination of L and r), each list represents a new generated catalog. After generating one catalog the matching algorithm will do the matching for all 21 settings, save the result in matches_list and move to generating a new catalog.
catalog_amount = 1000 # Put the number of catalogs to be generated here
catalogs={}
for l in range(catalog_amount):
    random_polar=[] # The next 3 lines load a random choice of shifting values for the generator
    for j in range(1177): 
        random_polar.append(random.choice(solutions))
    catalogs[l] = ([ra+random.choice(random_polar)[0] for ra in df["RA"]], [dec+random.choice(random_polar)[1] for dec in df["DEC"]]) # Here the original catalog is scrambled, but the scramling process shifts some blazars outside the desired boundaries.
    # Shift all blazars outside the outer boundaries back inside:
    for j in range(1177):
        while np.logical_or(catalogs[l][1][j] > -5, np.logical_or(catalogs[l][1][j] < -85, np.logical_or(catalogs[l][0][j] > 360, catalogs[l][0][j] < 0))):
            catalogs[l][0][j] = df["RA"][j] + random.choice(random_polar)[0]
            catalogs[l][1][j] = df["DEC"][j] + random.choice(random_polar)[1]
    # Shift all blazars near the galactic plane into the boundaries:
    catalogs_gal = SkyCoord(ra=catalogs[l][0]*u.degree, dec=catalogs[l][1]*u.degree, frame="icrs").galactic
    for j in range(1177):
        while np.logical_or(0 < catalogs_gal.b.value[j] < 10, -10 < catalogs_gal.b.value[j] < 0):
            n = random.choice(random_polar)
            if np.logical_and(df["DEC"][j] + n[1] < -5, df["DEC"][j] + n[1] > -85):
                catalogs[l][0][j] = df["RA"][j] + n[0]
                catalogs[l][1][j] = df["DEC"][j] + n[1]
                catalogs_gal.b.value[j] = SkyCoord(ra=catalogs[l][0][j]*u.degree, dec=catalogs[l][1][j]*u.degree, frame="icrs").galactic.b.value
                
    # The matching algorithm as a defined function: 
    def match_catalog(catalogs, L, r):
        if L == 3.5:
            df_centers = pd.read_csv("/home/simon/condor/neutrino_hotspots_35.csv").to_numpy()
        if L == 4.0:
            df_centers = pd.read_csv("/home/simon/condor/neutrino_hotspots_40.csv").to_numpy()
        if L == 4.5:
            df_centers = pd.read_csv("/home/simon/condor/neutrino_hotspots_45.csv").to_numpy()
        match_count = []
        for n in range(df_centers.shape[0]):
            for j in range(len(catalogs[0])):
                if haversine([df_centers[n,2], df_centers[n,1]], [catalogs[1][j], catalogs[0][j]], unit=Unit.DEGREES, normalize=True) < r:
                    match_count.append(n)
        return len(list(set(match_count)))
    
    
    # For each setting of L and r, the number of found matches are saved in a list. This list is then appended to the output matches_list.
    globals()['matches_%d' % l] = []
    for L in [3.5,4,4.5]:
        for r in [0.4,0.45,0.5,0.55,0.6,0.65,0.7]:
            globals()['matches_%d' % l].append(match_catalog(catalogs[l], L, r))
    matches_list.append(globals()['matches_%d' % l])
    
outfname = "{0:s}/output{1:n}.txt".format(output_dir, i)

f = open(outfname, "w")
f.write(str(matches_list))
f.close()
