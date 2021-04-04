import numpy.ma as ma
from netCDF4 import Dataset
from numpy import dtype
import numpy as np
from ncdf_read import select_region_TLL_files
from ncdf_read import select_region_TLLL_files
from stat_ncl import clmmonTLL,calcmonanomTLL
import datetime as dt
import matplotlib.pyplot as plt
from textwrap import wrap
plt.rcParams['font.family'] = 'YuGothic'
from mpl_toolkits.axes_grid1 import make_axes_locatable
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import cmocean
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# Python script for displaying composited SST for specific events (In this case, Sep-Nov of pIOD years)
# as a map
# Figure 1

dir_in="../data/"

lon1=30.0
lon2=120.0
lat1=-30.0
lat2=30.0
lev1=0.0
lev2=0.0

start_year=1959
end_year=2017
sec_to_month=60.0*60.0*24.0*30.0

# Targetedv events
target_years=[1961, 1963, 1967, 1972, 1977, 1994, 1997, 2006]
ntarget_years=[1960, 1964, 1974, 1975, 1996, 1998, 2005, 2010]

figname="map_sst.png"

dt1=dt.datetime(start_year,1,1,0,0,0) # 1959-01-01 00:00
dt2=dt.datetime(end_year,12,31,0,0,0) # 2015-01-01 00:00

#read file
fnames_temp=[]

for iyear in range(start_year,end_year+1):
    fnames_temp.append(dir_in+"roms_sst_"+str(iyear)+"_monthly.nc")

varname_temp="temp"
lon_temp,lat_temp,lev_temp,time_temp,temp=select_region_TLLL_files(fnames_temp,varname_temp,dt1,dt2,lev1,lev2,lat1,lat2,lon1,lon2,return_dims=True)
time_temp=np.asarray(time_temp) # Make time list to arra
temp[temp > 9000]=None
temp=np.nanmean(temp,axis=1) # Average along depth-direction
temp_clm=clmmonTLL(temp) # Get climatology
temp_anm=calcmonanomTLL(temp,temp_clm) # Get interannual anomalies
nx=lon_temp.shape[0]
ny=lat_temp.shape[0]
nwgt=72
wgt=np.ones(nwgt)/nwgt
wgt=np.ones(72)
wgt=wgt/np.sum(wgt)
temp_smooth=np.apply_along_axis(lambda m:np.convolve(m,wgt,"same"),0,temp_anm)
temp_anm=temp_anm-temp_smooth

# Slice targeted year
nevent=len(target_years)
temp_comp=np.copy(temp_anm[0:nevent,:,:])
temp_compn=np.copy(temp_anm[0:nevent,:,:])

for i in range(0,nevent):
    # Get time
    target_dt1=dt.datetime(target_years[i],9,1,0,0,0)   # sep 1st of targeted year
    target_dt2=dt.datetime(target_years[i],11,30,0,0,0) # nov 30th of targeted year
    itime=np.where((time_temp >= target_dt1) & (time_temp <= target_dt2))[0]
    temp_comp[i,:,:]=np.nanmean(temp_anm[itime,:,:],axis=0)
for i in range(0,nevent):
    # Get time
    ntarget_dt1=dt.datetime(ntarget_years[i],9,1,0,0,0)   # sep 1st of targeted year
    ntarget_dt2=dt.datetime(ntarget_years[i],11,30,0,0,0) # nov 30th of targeted year
    nitime=np.where((time_temp >= ntarget_dt1) & (time_temp <= ntarget_dt2))[0]
    temp_compn[i,:,:]=np.nanmean(temp_anm[nitime,:,:],axis=0)

# Composited value
temp_view=np.mean(temp_comp[:,:,:,],axis=0)
temp_view.shape
temp_view1=np.mean(temp_comp[:,:,:,],axis=0)
temp_viewn=np.mean(temp_compn[:,:,:,],axis=0)
temp_viewn.shape
temp_viewn1=np.mean(temp_compn[:,:,:,],axis=0)
nx=temp_view.shape[0]
ny=temp_view.shape[1]
pn_temp_view=temp_view+temp_viewn

#Significance test
t95=2.364#95%
temp_t95=abs(temp_view[:,:])*(nevent)**0.50/(np.std(temp_comp[:,:,:],axis=0)+1.0e-12)-t95
ntemp_t95=abs(temp_viewn[:,:])*(nevent)**0.50/(np.std(temp_compn[:,:,:],axis=0)+1.0e-12)-t95
t_diff_95=abs(temp_view+temp_viewn)*(nevent*nevent*(nevent*2-2)/nevent*2)**0.50/(nevent*np.std(temp_comp,axis=0)**2+nevent*np.std(temp_compn,axis=0)**2)**0.50-t95

temp_view1[temp_t95 < 0] = None
temp_viewn1[ntemp_t95 < 0] = None
pn_temp_view[t_diff_95 < 0] = None

#ERSST
from matplotlib import rcParams
rcParams['contour.negative_linestyle']='dashdot'
from mpl_toolkits.axes_grid1 import make_axes_locatable
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import cmocean
import cartopy.crs as ccrs
import cartopy.feature as cfeature

#read data from netcdf, convert invalid values as None and retrieve array size
ncin = Dataset('../data/ersst_all.nc', 'r')
sst = ncin.variables['SST'][1152:1968,:,:]
lons = ncin.variables['LONN89_91'][:]
lats = ncin.variables['LAT'][:]
sst[sst > 100] = None
sizearray = np.shape(sst)
nx=181
ny=89

#calculating climatology
sst_year_mon = sst.reshape(68, 12, ny, nx)
sst_mon_clim = ma.empty((12, ny, nx))
for k in range(12):
    sst_mon_clim[k,:,:] = sst_year_mon[:,k,:,:].mean(axis=0)

#calculating anomaly
nt=816  -72
sst_year_mon = sst.reshape(68, 12, ny, nx)
sst_anom = ma.empty((68, 12, ny, nx))
for k in range(12):
    sst_anom[:,k,:,:] = sst_year_mon[:,k,:,:]-sst_mon_clim[k,:,:]
sst_anom2 = sst_anom.reshape(816, ny, nx)

sst_anom_smooth = ma.empty((nt,ny,nx))
sst_anom_rm = ma.empty((nt,ny,nx))
sst_anom_smooth[0,:,:] = sst_anom2[0:73,:,:].mean(axis=0)
for i in range(1,nt):
    sst_anom_smooth[i,:,:] = sst_anom_smooth[i-1,:,:] - sst_anom2[i-1,:,:]/73 \
    + sst_anom2[i+72,:,:]/73
for i in range(0,nt):
    sst_anom_rm[i,:,:] = (sst_anom2[i+35,:,:]+sst_anom2[i+36,:,:]+sst_anom2[i+37,:,:])/3.0 \
    - sst_anom_smooth[i,:,:]

sst_anom_rm1 = sst_anom_rm.reshape(62,12,ny,nx)

 #Composite analysis
nevent=8
t95=2.364 
pIOD=(1961, 1963, 1967, 1972, 1977, 1994, 1997, 2006)#posi
nIOD=(1958, 1960, 1964, 1974, 1975, 1996, 1998, 2010)#nega

sep = []
oct = []
nov = []
nsep = []
noct = []
nnov = []

for k in range (nevent):
    for i in range (nt):
        it = (pIOD[k]-1953)*12
        if (i) == (it+8):
            sep.append(sst_anom_rm[i,:,:])
        if (i) == (it+9):
            oct.append(sst_anom_rm[i,:,:])
        if (i) == (it+10):
            nov.append(sst_anom_rm[i,:,:])            

for k in range (nevent):
    for i in range (nt):
        it = (nIOD[k]-1953)*12
        if (i) == (it+8):
            nsep.append(sst_anom_rm[i,:,:])
        if (i) == (it+9):
            noct.append(sst_anom_rm[i,:,:])
        if (i) == (it+10):
            nnov.append(sst_anom_rm[i,:,:])

sep=np.array(sep)
oct=np.array(oct)
nov=np.array(nov)
nsep=np.array(nsep)
noct=np.array(noct)
nnov=np.array(nnov)

son=(sep+oct+nov)/3
nson=(nsep+noct+nnov)/3

son_comp = np.nanmean(son, axis=0)
son_comp1 = np.nanmean(son, axis=0)
son_t = abs(son_comp)*(nevent)**0.50/(np.std(son, axis=0)+0.000001)-t95
nson_comp = np.nanmean(nson, axis=0)
nson_comp1 = np.nanmean(nson, axis=0)
nson_t = abs(nson_comp)*(nevent)**0.50/(np.std(nson, axis=0)+0.000001)-t95
pnson_t=abs(son_comp[:,:]+nson_comp[:,:])*(nevent*nevent*(nevent*2-2)/nevent*2)**0.50/(nevent*np.std(nson[:,:,:],axis=0)**2+nevent*np.std(son[:,:,:],axis=0)**2)**0.50-t95

pnson_comp=son_comp+nson_comp
pnson_comp1=son_comp+nson_comp

son_comp1[son_t < 0.0] = None
nson_comp1[nson_t < 0.0] = None
pnson_comp1[pnson_t < 0.0] = None

#generate values for colorbar 
clevs1 = np.linspace(-2.0,2.0,21)

#plot figures
f1 = plt.figure(figsize=(40,20))
plt.subplots_adjust(wspace=0.25, hspace=0.05)
plt.subplots_adjust(left=0.06, right=0.95, bottom=0.1, top=0.95)
interval = np.arange(-3.0,3.0,0.2)

ax = f1.add_subplot(231, projection=ccrs.PlateCarree(central_longitude=70))
plt.rcParams["font.size"] = 30
ax.set_xticks([30, 60, 90, 110], crs=ccrs.PlateCarree())
ax.set_yticks([-30, -15, 0, 15, 30], crs=ccrs.PlateCarree())
lon_formatter = LongitudeFormatter(number_format='.0f',dateline_direction_label=True)
lat_formatter = LatitudeFormatter(number_format='.0f',)
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)
ax.set_extent([30, 110, -30, 30], ccrs.PlateCarree())
ax.outline_patch.set_linewidth(1)
ax.tick_params(axis='both', labelsize=25)
fill = plt.contourf(lons, lats, son_comp1, clevs1, extend='both', transform=ccrs.PlateCarree(), cmap=cmocean.cm.balance, zorder=0)
plt.contour(lons, lats, son_comp, interval, colors=('k'), extend='both', transform=ccrs.PlateCarree(), zorder=1)
ax.add_feature(cfeature.COASTLINE, linewidth=2.0, facecolor='w', edgecolor='k', zorder=2)
plt.title('(a) pIOD (ERSST)', loc='center')
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.1, axes_class=plt.Axes)
cb = plt.colorbar(fill, orientation='vertical', cax=cax)
cb.ax.set_yticklabels(cb.ax.get_yticklabels())
cb.set_ticks([-2.0,-1.5,-1.0,-0.5,0,0.5,1.0,1.5,2.0])
cb.set_label('℃')

ax = f1.add_subplot(232, projection=ccrs.PlateCarree(central_longitude=70))
plt.rcParams["font.size"] = 30
ax.set_xticks([30, 60, 90, 110], crs=ccrs.PlateCarree())
ax.set_yticks([-30, -15, 0, 15, 30], crs=ccrs.PlateCarree())
lon_formatter = LongitudeFormatter(number_format='.0f',dateline_direction_label=True)
lat_formatter = LatitudeFormatter(number_format='.0f',)
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)
ax.set_extent([30, 110, -30, 30], ccrs.PlateCarree())
ax.outline_patch.set_linewidth(1)
ax.tick_params(axis='both', labelsize=25)
fill = plt.contourf(lons, lats, nson_comp1, clevs1, extend='both', transform=ccrs.PlateCarree(), cmap=cmocean.cm.balance, zorder=0)
plt.contour(lons, lats, nson_comp, interval, colors=('k'), extend='both', transform=ccrs.PlateCarree(), zorder=1)
ax.add_feature(cfeature.COASTLINE, linewidth=2.0, facecolor='w', edgecolor='k', zorder=2)
plt.title('(b) nIOD (ERSST)', loc='center')
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.1, axes_class=plt.Axes)
cb = plt.colorbar(fill, orientation='vertical', cax=cax)
cb.ax.set_yticklabels(cb.ax.get_yticklabels())
cb.set_ticks([-2.0,-1.5,-1.0,-0.5,0,0.5,1.0,1.5,2.0])
cb.set_label('℃')

ax = f1.add_subplot(233, projection=ccrs.PlateCarree(central_longitude=70))
plt.rcParams["font.size"] = 30
ax.set_xticks([30, 60, 90, 110], crs=ccrs.PlateCarree())
ax.set_yticks([-30, -15, 0, 15, 30], crs=ccrs.PlateCarree())
lon_formatter = LongitudeFormatter(number_format='.0f',dateline_direction_label=True)
lat_formatter = LatitudeFormatter(number_format='.0f',)
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)
ax.set_extent([30, 110, -30, 30], ccrs.PlateCarree())
ax.outline_patch.set_linewidth(1)
ax.tick_params(axis='both', labelsize=25)
fill = plt.contourf(lons, lats, pnson_comp1, clevs1, extend='both', transform=ccrs.PlateCarree(), cmap=cmocean.cm.balance, zorder=0)
plt.contour(lons, lats, pnson_comp, interval, colors=('k'), extend='both', transform=ccrs.PlateCarree(), zorder=1)
ax.add_feature(cfeature.COASTLINE, linewidth=2.0, facecolor='w', edgecolor='k', zorder=2)
plt.title('(c) pIOD+nIOD (ERSST)', loc='center')
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.1, axes_class=plt.Axes)
cb = plt.colorbar(fill, orientation='vertical', cax=cax)
cb.ax.set_yticklabels(cb.ax.get_yticklabels())
cb.set_ticks([-2.0,-1.5,-1.0,-0.5,0,0.5,1.0,1.5,2.0])
cb.set_label('℃')

ax = f1.add_subplot(234, title="(d) pIOD (ROMS)",projection=ccrs.PlateCarree(central_longitude=70))
ax.set_xticks([30, 60, 90, 110], crs=ccrs.PlateCarree())
ax.set_yticks([-30, -15, 0, 15, 30], crs=ccrs.PlateCarree())
lon_formatter = LongitudeFormatter(number_format='.0f',dateline_direction_label=True)
lat_formatter = LatitudeFormatter(number_format='.0f',)
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)
ax.set_extent([30, 110, -30, 30], ccrs.PlateCarree())
ax.outline_patch.set_linewidth(1)
ax.tick_params(axis='both', labelsize=25)
fill=plt.contourf(lon_2d,lat_2d,temp_view1[:,:],levels=clevs1,extend="both",transform=ccrs.PlateCarree(),cmap=cmocean.cm.balance,zorder=0)
plt.contour(lon_2d,lat_2d,temp_view[:,:],levels=interval,colors="black",transform=ccrs.PlateCarree(),zorder=1)
ax.add_feature(cfeature.COASTLINE, linewidth=2.0, facecolor='w', edgecolor='k', zorder=2)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.1, axes_class=plt.Axes)
cb = plt.colorbar(fill, orientation='vertical', cax=cax)
cb.ax.set_yticklabels(cb.ax.get_yticklabels())
cb.set_label('℃')

ax = f1.add_subplot(235, title="(e) nIOD (ROMS)",projection=ccrs.PlateCarree(central_longitude=70))
ax.set_xticks([30, 60, 90, 110], crs=ccrs.PlateCarree())
ax.set_yticks([-30, -15, 0, 15, 30], crs=ccrs.PlateCarree())
lon_formatter = LongitudeFormatter(number_format='.0f',dateline_direction_label=True)
lat_formatter = LatitudeFormatter(number_format='.0f',)
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)
ax.set_extent([30, 110, -30, 30], ccrs.PlateCarree())
ax.outline_patch.set_linewidth(1)
ax.tick_params(axis='both', labelsize=25)
fill=plt.contourf(lon_2d,lat_2d,ntemp_view1[:,:],levels=clevs1,extend="both",transform=ccrs.PlateCarree(),cmap=cmocean.cm.balance,zorder=0)
plt.contour(lon_2d,lat_2d,ntemp_view[:,:],levels=interval,colors="black",transform=ccrs.PlateCarree(),zorder=1)
ax.add_feature(cfeature.COASTLINE, linewidth=2.0, facecolor='w', edgecolor='k', zorder=2)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.1, axes_class=plt.Axes)
cb = plt.colorbar(fill, orientation='vertical', cax=cax)
cb.ax.set_yticklabels(cb.ax.get_yticklabels())
cb.set_ticks([-2.0,-1.5,-1.0,-0.5,0,0.5,1.0,1.5,2.0])
cb.set_label('℃')

ax = f1.add_subplot(236, title="(f) pIOD+nIOD (ROMS)",projection=ccrs.PlateCarree(central_longitude=70))
ax.set_xticks([30, 60, 90, 110], crs=ccrs.PlateCarree())
ax.set_yticks([-30, -15, 0, 15, 30], crs=ccrs.PlateCarree())
lon_formatter = LongitudeFormatter(number_format='.0f',dateline_direction_label=True)
lat_formatter = LatitudeFormatter(number_format='.0f',)
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)
ax.set_extent([30, 110, -30, 30], ccrs.PlateCarree())
ax.outline_patch.set_linewidth(1)
ax.tick_params(axis='both', labelsize=25)
fill=plt.contourf(lon_2d,lat_2d,pntemp_view1[:,:],levels=clevs1,extend="both",transform=ccrs.PlateCarree(),cmap=cmocean.cm.balance,zorder=0)
plt.contour(lon_2d,lat_2d,pntemp_view[:,:],levels=interval,colors="black",transform=ccrs.PlateCarree(),zorder=1)
ax.add_feature(cfeature.COASTLINE, linewidth=2.0, facecolor='w', edgecolor='k', zorder=2)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.1, axes_class=plt.Axes)
cb = plt.colorbar(fill, orientation='vertical', cax=cax)
cb.ax.set_yticklabels(cb.ax.get_yticklabels())
cb.set_ticks([-2.0,-1.5,-1.0,-0.5,0,0.5,1.0,1.5,2.0])
cb.set_label('℃')

plt.savefig(figname, bbox_inches='tight')
print(figname)

'''
from numpy import dtype
nc = Dataset('/Volumes/Promise-Pegasus_2/nakazato/HEATBUD/comp_map_sst2.nc', 'w')
#nc.createDimension('time1', None)
nc.createDimension('lat', nx)
nc.createDimension('lon', ny)
#nc.createDimension('time', nt)
lat = nc.createVariable('lat', 'f', ('lat'))
lat.long_name = 'latitude'
lat.units = 'degrees north'
lat.axis = 'Y'
lat[:] = lat_temp

lon = nc.createVariable('lon', 'f', ('lon'))
lon.long_name = 'longitude'
lon.units = 'degrees east'
lon.axis = 'X'
lon[:] = lon_temp

view = nc.createVariable('view', 'f', ('lat', 'lon'))
view.long_name = 'temp'
view.units = 'K'
#hgt.missing_value = 1.e+20
view[:,:] = temp_view[:,:]

viewn = nc.createVariable('viewn', 'f', ('lat', 'lon'))
viewn.long_name = 'ntemp'
viewn.units = 'K'
#hgt.missing_value = 1.e+20
viewn[:,:] = temp_viewn[:,:]

t_view = nc.createVariable('t_test', 'f', ('lat', 'lon'))
t_view.long_name = 't_test'
t_view.units = ''
#hgt.missing_value = 1.e+20
t_view[:,:] = temp_t95[:,:]

t_nview = nc.createVariable('nt_test', 'f', ('lat', 'lon'))
t_nview.long_name = 'nt_test'
t_nview.units = ''
#hgt.missing_value = 1.e+20
t_nview[:,:] = ntemp_t95[:,:]

t_pnview = nc.createVariable('pnt_test', 'f', ('lat', 'lon'))
t_pnview.long_name = 'pnt_test'
t_pnview.units = ''
#hgt.missing_value = 1.e+20
t_pnview[:,:] = t_diff_95[:,:]

nc.close()
'''


