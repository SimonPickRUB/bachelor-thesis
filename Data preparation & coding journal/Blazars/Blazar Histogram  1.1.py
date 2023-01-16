#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import pandas as pd
from astropy.io import fits
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table
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


# # DATA CLEANING:
# Take-home message: Always use pandas if possible

# In[4]:


# converting fits data to pandas dataframe
tablehdu = fits.TableHDU(data)
tablehdu.writeto("/mnt/c/Users/joerg/Desktop/blazar1.fits", overwrite=True)
datatable = Table.read("/mnt/c/Users/joerg/Desktop/blazar1.fits", format="fits")
df = datatable.to_pandas()

# add galactic coordinates to dataframe
galactic_coords = SkyCoord(ra=data["RA"]*u.degree, dec=data["DEC"]*u.degree, frame="icrs").galactic
dataframe_list = list(zip(galactic_coords.l.degree, galactic_coords.b.degree))
df1 = pd.DataFrame(dataframe_list, columns=["l","b"])
df_complete_dataset = pd.concat([df, df1], axis=1)
df_complete_dataset


# In[5]:


df3 = df_complete_dataset[(df_complete_dataset["DEC"]<-5) & (df_complete_dataset["DEC"]>-85)]
df4 = df3[df3["b"]>10]
df5 = df3[df3["b"]<-10]
df6 = pd.concat([df4, df5], ignore_index=True)
df7 = df6.copy()
for i in range(1207):
    if str(df6["OBJECT_TYPE"].values[i]) == "b'BL Lac Candidate        '":
        df7 = df7.drop(i)
new_index = list(range(1177))
df_clean = df7.reindex(new_index)
df_clean


# ## Blazar Catalog Data Plot (Cleaned Data)

# In[15]:


plt.scatter(df_clean["RA"].values, df_clean["DEC"].values, s=10)
plt.title("2D Histogram of Blazar Catalog Data (Cleaned Data)")
plt.xlabel(r"right ascension $\alpha$ / °")
plt.ylabel(r"declination $\delta$ / °")

l_interval = np.linspace(0, 360, 360) - np.linspace(141, 332, 360)
galactic_plane_coords = SkyCoord(l=l_interval*u.degree, b=0*u.degree, frame='galactic')
gp = plt.plot(galactic_plane_coords.icrs.ra, galactic_plane_coords.icrs.dec, linewidth=0.1, label="galactic plane")

#galactic_plane_upper = SkyCoord(l=l_interval*u.degree, b=-10*u.degree, frame='galactic')
#galactic_plane_lower = SkyCoord(l=l_interval*u.degree, b=10*u.degree, frame='galactic')
#gp1 = plt.plot(galactic_plane_upper.icrs.ra, galactic_plane_upper.icrs.dec, linewidth=0.1, label="galactic plane1")
#gp2 = plt.plot(galactic_plane_lower.icrs.ra, galactic_plane_lower.icrs.dec, linewidth=0.1, label="galactic plane2")


# ## Working with cleaned data: Redshift distribution

# In[8]:


plt.hist(df_clean["REDSHIFT"], bins=100, rwidth=0.9)
plt.title("Blazars with no redshift information are shown as z=0")
plt.suptitle("Redshift distribution")
plt.xlabel("Redshift z")
plt.ylabel("Number of blazars per bin")
plt.xticks([0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5])
plt.show()


# In[ ]:




