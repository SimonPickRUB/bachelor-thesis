#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import pandas as pd
from astropy.io import fits
from astropy import units as u
from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt

file = fits.open("/mnt/c/Users/joerg/Desktop/3561_objects.fits")
data = file[1].data


# ## Blazar Catalog Data Histogram (LowRes)

# In[2]:


#catalog data
plt.hist2d(data.field(2), data.field(3), bins =[500, 500],cmap = "gray")
plt.title("2D Histogram of Blazar Catalog Data (LowRes)")
plt.xlabel(r"right ascension $\alpha$ / °")
plt.ylabel(r"declination $\delta$ / °")

#galactic plane and galactic center
l_interval = np.linspace(0, 360, 360) - np.linspace(243, 244, 360)
galactic_plane_coords = SkyCoord(l=l_interval*u.degree, b=0*u.degree, frame='galactic')
galactic_center_coords = SkyCoord(l=0*u.degree, b=0*u.degree, frame="galactic")
gp = plt.plot(galactic_plane_coords.icrs.ra, galactic_plane_coords.icrs.dec, linewidth=0.4, label="galactic plane")
gc = plt.scatter(galactic_center_coords.icrs.ra, galactic_center_coords.icrs.dec, s=70, marker="*", label="galactic center")
plt.legend(loc=1, fontsize=7, fancybox=False)

plt.savefig("/mnt/c/Users/joerg/Desktop/2D_Histogram_Blazar_LowRes.png", dpi=500)


# ## Blazar Catalog Data Histogram (HighRes)

# In[3]:


#catalog data
plt.hist2d(data.field(2), data.field(3), bins =[2800, 2800], cmap = "gray")
plt.title("2D Histogram of Blazar Catalog Data (HighRes)")
plt.xlabel(r"right ascension $\alpha$ / °")
plt.ylabel(r"declination $\delta$ / °")

#galactic plane and galactic center
l_interval = np.linspace(0, 360, 360) - np.linspace(243, 244, 360)
galactic_plane_coords = SkyCoord(l=l_interval*u.degree, b=0*u.degree, frame='galactic')
galactic_center_coords = SkyCoord(l=0*u.degree, b=0*u.degree, frame="galactic")
gp = plt.plot(galactic_plane_coords.icrs.ra, galactic_plane_coords.icrs.dec, linewidth=0.2, label="galactic plane")
gc = plt.scatter(galactic_center_coords.icrs.ra, galactic_center_coords.icrs.dec, s=30, marker="*", label="galactic center")
plt.legend(loc=1, fontsize=7, fancybox=False)

plt.savefig("/mnt/c/Users/joerg/Desktop/2D_Histogram_Blazar_HighRes.png", dpi=1300)


# # Now: Data Cleaning!
# ## Remove data in equatorial coordinate system first

# In[4]:


print(len(data["DEC"]))
clean_data0 = data[data["DEC"]>-85]
print(len(clean_data0))
clean_data1 = clean_data0[clean_data0["DEC"]<-5]
print(len(clean_data1))
clean_data2 = clean_data1[clean_data1["OBJECT_TYPE"] != "BL Lac Candidate"]
print(len(clean_data2))


# ## Remove data in galactic coordinate system
# Coordinate system transformation from equatorial to galactic with astropy

# In[5]:


# Coordinate transformation
coords = SkyCoord(ra=clean_data2.field(2), dec=clean_data2.field(3), unit="deg").galactic
ra = []
dec = []
for i_ra in coords:
    ra.append((i_ra.l.degree))
for i_dec in coords:
    dec.append(i_dec.b.degree)
coords = np.moveaxis(np.asarray([ra,dec]),-1,0)

# Choose values for the interval limits for the area around the galactic plane
galactic_plane_min = -10
galactic_plane_max = 10
coords_clean = np.copy(coords)[::-1]
bad_b = []
bad_l=[]
for i,coord in enumerate(coords):
    l=coord[0]
    b=coord[1]
    if galactic_plane_min <b<galactic_plane_max:
        coords_clean = np.delete(coords_clean,coords.shape[0]-i-1,axis=0)
        bad_b.append(b)
        bad_l.append(l)
print(f"Number of deleted entries: {coords.shape[0]-coords_clean.shape[0]}")


# ## Blazar Catalog Data Histogram (Cleaned Data)

# In[6]:


coords_clean_galactic = SkyCoord(l=coords_clean[:,0]*u.degree, b=coords_clean[:,1]*u.degree, frame="galactic").icrs
plt.hist2d(coords_clean_galactic.ra.degree, coords_clean_galactic.dec.degree, bins =[1500, 1500], cmap = "gray")
plt.title("2D Histogram of Blazar Catalog Data (Cleaned Data)")
plt.xlabel(r"right ascension $\alpha$ / °")
plt.ylabel(r"declination $\delta$ / °")

l_interval = np.linspace(0, 360, 360) - np.linspace(240, 250, 360)
galactic_plane_coords = SkyCoord(l=l_interval*u.degree, b=0*u.degree, frame='galactic')
gp = plt.plot(galactic_plane_coords.icrs.ra, galactic_plane_coords.icrs.dec, linewidth=0.1, label="galactic plane")

#galactic_plane_upper = SkyCoord(l=l_interval*u.degree, b=-10*u.degree, frame='galactic')
#galactic_plane_lower = SkyCoord(l=l_interval*u.degree, b=10*u.degree, frame='galactic')
#gp1 = plt.plot(galactic_plane_upper.icrs.ra, galactic_plane_upper.icrs.dec, linewidth=0.1, label="galactic plane1")
#gp2 = plt.plot(galactic_plane_lower.icrs.ra, galactic_plane_lower.icrs.dec, linewidth=0.1, label="galactic plane2")

plt.savefig("/mnt/c/Users/joerg/Desktop/2D_Histogram_Blazar_Test.png", dpi=1300)
print(f"The number of files in the final blazar data set is {len(coords_clean)}.")


# # Redshift 
# The previous steps have only cleaned the coordinates of the data set but not the whole data set with all the other information. This will be done next:

# In[7]:


# This cell is working but needs to be fixed to solve following issue:
# clean_data2 gets overwritten and range is only up to 1177, should be 1214
# To fix this issue rewrite this cell using pandas and use the np.copy() function!

bad_icrs = SkyCoord(l=bad_l*u.degree, b=bad_b*u.degree, frame="galactic").icrs
for i in range(1177):
    for j in range(37):
        if round(bad_icrs.dec.degree[j],5) == round(clean_data2["DEC"][i],5):
            clean_data2 = clean_data2[clean_data2["DEC"]!=round(bad_icrs.dec.degree[j],6)]
print(clean_data2.shape)

# Switch to pandas for simpler data management
from astropy.table import Table
hdu = fits.TableHDU(clean_data2)
hdu.writeto("/mnt/c/Users/joerg/Desktop/blazar.fits", overwrite=True)
datatable = Table.read("/mnt/c/Users/joerg/Desktop/blazar.fits", format="fits")
df = datatable.to_pandas()
df


# ## Working with cleaned data: Redshift distribution

# In[8]:


plt.hist(clean_data2.field(4), bins=100, rwidth=0.9)
plt.title("Blazars with no redshift information are shown as z=0")
plt.suptitle("Redshift distribution")
plt.xlabel("Redshift z")
plt.ylabel("Number of blazars per bin")
plt.xticks([0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5])
plt.show()


# In[ ]:




