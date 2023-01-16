#!/usr/bin/env python
# coding: utf-8

# In[1]:


import matplotlib.pyplot as plt

import numpy as np

import healpy as hp
from healpy.newvisufunc import projview

import math

from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits

import pandas as pd


# ## Visualisation of HEALPix map with healpy

# In[ ]:


Lmap = hp.read_map("/mnt/c/Users/joerg/Desktop/log10p_maps.fits")
projview(Lmap, coord=["E"], projection_type="mollweide", title="p-value map (equatorial)", graticule=True, cmap="gray", unit="-log10(p-value)")


# ## Converting pixels to eq. coordinates and weighing with p-values

# In[ ]:


ra, dec = hp.pix2ang(nside=512, ipix=[i for i in range(3145728)], lonlat=True)
plt.hist2d(ra, dec, bins=[2048,2048], weights=Lmap)
plt.title("p-value histogram (all data)")
plt.xlabel(r"right ascension $\alpha$ / °")
plt.ylabel(r"declination $\delta$ / °")
plt.yticks([-90, -60, -30, 0, 30, 60, 90])
plt.xticks([0, 60, 120, 180, 240, 300, 359.93], ["0", "60", "120", "180", "240", "300", "360"])
plt.show()


# ## Having a look at the data with astropy

# In[ ]:


fit = fits.open("/mnt/c/Users/joerg/Desktop/log10p_maps.fits")
data = fit[1].data.field(0)
print(f"The data of the fits file is stored in {data.shape[0]} rows and {data.shape[1]} columns representing pixels of the HEALPix pixelization scheme.")
print("The value of each individual entry is -log(p-value) of the corresponding pixel on the map.")
print(data)


# #  Cleaning data with pandas DataFrames

# In[ ]:


galactic_coords = SkyCoord(ra=ra*u.degree, dec=dec*u.degree, frame="icrs").galactic
coord_list = list(zip(ra, dec, Lmap, galactic_coords.l.degree, galactic_coords.b.degree))
df = pd.DataFrame(coord_list, columns=["RA","DEC","-log[p-value]","l","b"])
df1 = df[(df["DEC"]<-5) & (df["DEC"]>-85)]
df2 = df1[df1["b"]>10]
df3 = df1[df1["b"]<-10]
df4 = pd.concat([df2, df3], ignore_index=True)
df4


# In[ ]:


plt.hist2d(df4["RA"].values, df4["DEC"].values, bins=[1024,1024], weights=df4["-log[p-value]"])
plt.title("p-value histogram (equatorial)")
plt.xlabel(r"right ascension $\alpha$ / °")
plt.ylabel(r"declination $\delta$ / °")
plt.yticks([-85, -80, -70, -60, -50, -40, -30, -20, -10, -5])
plt.xticks([0, 60, 120, 180, 240, 300, 359.93], ["0", "60", "120", "180", "240", "300", "360"])
plt.show()


# In[ ]:



