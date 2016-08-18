import netCDF4
from netCDF4 import Dataset, MFDataset, num2date, date2num, date2index
import numpy as np
import numpy.ma as ma
from datetime import datetime
import pandas as pd
from pandas import Series,DataFrame, Panel
from datetime import date,datetime,timedelta
from mpl_toolkits.basemap import Basemap,cm
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import ocw.data_source.local as local
import ocw.utils as utils
from glob import glob

#
target_lat = 15.5


#file_object = Dataset('merra2_TQ_1980_Central_America.nc')
file_object = Dataset('merra2_TQH_1980_Central_America.nc')
lat = file_object.variables['lat'][:]
lon = file_object.variables['lon'][:]
lev = file_object.variables['lev'][:]
y_index = np.where(lat == target_lat)[0]

#nc_files = glob('merra2_TQ_*.nc')
nc_files = glob('merra2_TQH_*.nc')
nc_files.sort()

nyear = 26
nz=25
nx=63
MSE = ma.zeros([nyear*365, nz,  nx])
times = []
for ifile,file in enumerate(nc_files):
    f = Dataset(file)
    #MSE[365*ifile:365+365*ifile,:] = ma.squeeze(f.variables['T'][:,:,y_index,:]*1004.+f.variables['QV'][:,:,y_index,:]*2.5e6)/1.e+3 
    MSE[365*ifile:365+365*ifile,:] = ma.squeeze(f.variables['T'][:,:,y_index,:]*1004.+f.variables['H'][:,:,y_index,:]*9.8+f.variables['QV'][:,:,y_index,:]*2.5e6)/1.e+3 
    times.extend(num2date(f.variables['time'][:],units='days since 1980-01-01 00:00:00'))

times = np.array(times)


# onset or demise dates
dates = [datetime(1980, 5, 18),
         datetime(1981, 4, 14),
         datetime(1982, 5, 17),
         datetime(1983, 5, 25),
         datetime(1984, 5, 8 ),
         datetime(1985, 6, 9 ),
         datetime(1986, 5, 18),
         datetime(1987, 5, 29),
         datetime(1988, 6, 6 ),
         datetime(1989, 6, 13),
         datetime(1990, 5, 31),
         datetime(1991, 6, 6 ),
         datetime(1992, 5, 15),
         datetime(1993, 5, 9 ),
         datetime(1994, 8, 1 ),
         datetime(1995, 6, 10),
         datetime(1996, 6, 3 ),
         datetime(1997, 5, 31),
         datetime(1998, 5, 24),
         datetime(1999, 6, 4 ),
         datetime(2000, 5, 13),
         datetime(2001, 5, 23),
         datetime(2002, 5, 14),
         datetime(2003, 5, 25),
         datetime(2004, 4, 27),
         datetime(2005, 5, 27)]

# find time indicies for onset or demise         
t_index = [np.where(times == date)[0] for date in dates]
nlag = 23
time_lag = -12 + np.arange(nlag)  # 12 days before ~ 10 days after

# deseasonalize
daily_clim = ma.mean(MSE.reshape([26,365,nz,nx]), axis=0) 
values2 = MSE.reshape([26,365,nz,nx])
for iday in np.arange(365):
    values2[:,iday,:] = values2[:,iday,:] - daily_clim[iday,:]
values2 = values2.reshape([26*365,nz,nx])

clim_var = ma.zeros([nlag, nz, nx])

for ilag, lag in enumerate(time_lag):
    clim_var[ilag,:] = ma.mean(values2[t_index+lag,:], axis=0)

fig = plt.figure()
clevs = np.arange(21)*0.2-2
ax=fig.add_subplot(1,1,1)
plot=ax.contourf(lon,lev,clim_var[12,:], levels=clevs, extend='both')
ax.invert_yaxis()
#ax.set_yscale('log')
ax.set_xlabel('Longitude (degrees)')
ax.set_ylabel('Pressure (mb)')
ax.set_title('MSE anomaly (including gz term) on the onset dates [kJ/kg]')
cbar_ax = fig.add_axes([0.92, 0.15, 0.01, 0.7])
fig.colorbar(plot, cax=cbar_ax)

plt.show()

#fig.savefig('MSE_anomalies_plus_gzterm_onset_dates_latitude_%d' %target_lat, dpi=600,bbox_inches='tight')
