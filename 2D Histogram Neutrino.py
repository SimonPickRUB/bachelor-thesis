#!/usr/bin/env python
# coding: utf-8

# In[1]:


import matplotlib.pyplot as plt
import numpy as np
import healpy as hp
from healpy.newvisufunc import projview


# In[2]:


map = hp.read_map("/mnt/c/Users/joerg/Desktop/log10p_maps.fits")
projview(map, coord=["G", "E"], projection_type="mollweide", title="p-value map", graticule=True, cmap="gray", unit="-log10(p-value)")


# In[3]:


from astropy.io import fits
fit = fits.open("/mnt/c/Users/joerg/Desktop/log10p_maps.fits")
data = fit[1].data.field(0)
print(data.shape)
print(data)


# In[4]:


ra, dec = hp.pix2ang(nside=512, ipix=[i for i in range(3145728)], lonlat=True)
print(ra)
print(dec)


# In[7]:


dec_clean = []
for i in range(data.shape[0]*data.shape[1]):
    if -5>dec[i]>-85:
        dec_clean.append(dec[i])
print(np.array(dec_clean))


# In[ ]:




