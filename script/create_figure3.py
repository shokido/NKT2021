from netCDF4 import Dataset
import numpy as np
from ncdf_read import select_region_TLL_files
from stat_ncl import clmmonTLL,calcmonanomTLL
import datetime as dt
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import cmocean
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# Python script for displaying composited 
# ML_heatbud term
# as a map
# Figure 3

dir_in="/Volumes/Promise-Pegasus_2/kido/ROMS/Outputs/IO_033/ZLEV/HEATBUD/"
lon1=90.0
lon2=110.0
lat1=-10.0
lat2=0.0
lev1=0.0
lev2=0.0

sign="posi"

start_year=1959
end_year=2017

figname="map_heatbud.png"

t95=2.364 #dof=nyear_iod - 1
# Targetedv events
#posi
target_years=[1961, 1963, 1967, 1972, 1977, 1994, 1997, 2006]
#nega
ntarget_years=[1960, 1964, 1974, 1975, 1980, 1996, 1998, 2010]
    
dt1=dt.datetime(start_year,1,1,0,0,0) # 1959-01-01 00:00
dt2=dt.datetime(end_year,12,31,0,0,0) # 2015-01-01 00:00

fnames_temp_rate=[]

for iyear in range(start_year,end_year+1):
    fnames_temp_rate.append(dir_in+"tio_hb_"+str(iyear)+"_monthly.nc")

fnames_temp_xadv=fnames_temp_rate
fnames_temp_yadv=fnames_temp_rate;
fnames_temp_vadv=fnames_temp_rate
fnames_temp_hdiff=fnames_temp_rate;
fnames_temp_vdiff=fnames_temp_rate;
fnames_temp_entr=fnames_temp_rate
fnames_temp_nhf=fnames_temp_rate
fnames_temp_ssrc=fnames_temp_rate

varname_temp_rate="temp_rate"
varname_temp_xadv="temp_xadv"; varname_temp_yadv="temp_yadv"
varname_temp_vadv="temp_vadv"; varname_temp_hdiff="temp_hdiff"
varname_temp_vdiff="temp_vdiff"; varname_temp_entr="temp_entr"
varname_temp_nhf="temp_nhf"; varname_temp_ssrc="temp_ssrc"

lon_temp_rate,lat_temp_rate,time_temp_rate,temp_rate=select_region_TLL_files(fnames_temp_rate,varname_temp_rate,dt1,dt2,lat1,lat2,lon1,lon2,return_dims=True)
temp_rate_clm=clmmonTLL(temp_rate) # Get climatology                            
temp_rate_anm=calcmonanomTLL(temp_rate,temp_rate_clm) # Get interannual anomalies                                                                              
nwgt=3
wgt=np.ones(nwgt)/nwgt
wgt=np.ones(3)
wgt=wgt/np.sum(wgt)
lnwgt=72
lwgt=np.ones(lnwgt)/lnwgt
lwgt=np.ones(72)
lwgt=lwgt/np.sum(lwgt)
lview=np.apply_along_axis(lambda m:np.convolve(m,lwgt,"same"),0,temp_rate_anm)
temp_rate_anm=np.apply_along_axis(lambda m:np.convolve(m,wgt,"same"),0,temp_rate_anm)-lview
nx=lon_temp_rate.shape[0]
ny=lat_temp_rate.shape[0]
nt=temp_rate_anm.shape[0]

lon_temp_xadv,lat_temp_xadv,time_temp_xadv,temp_xadv=select_region_TLL_files(fnames_temp_xadv,varname_temp_xadv,dt1,dt2,lat1,lat2,lon1,lon2,return_dims=True)
temp_xadv_clm=clmmonTLL(temp_xadv) # Get climatology                            
temp_xadv_anm=calcmonanomTLL(temp_xadv,temp_xadv_clm) # Get interannual anomalies
lview=np.apply_along_axis(lambda m:np.convolve(m,lwgt,"same"),0,temp_xadv_anm)
temp_xadv_anm=np.apply_along_axis(lambda m:np.convolve(m,wgt,"same"),0,temp_xadv_anm)-lview

lon_temp_yadv,lat_temp_yadv,time_temp_yadv,temp_yadv=select_region_TLL_files(fnames_temp_yadv,varname_temp_yadv,dt1,dt2,lat1,lat2,lon1,lon2,return_dims=True)
temp_yadv_clm=clmmonTLL(temp_yadv) # Get climatology                            
temp_yadv_anm=calcmonanomTLL(temp_yadv,temp_yadv_clm) # Get interannual anomalies                                                                              
lview=np.apply_along_axis(lambda m:np.convolve(m,lwgt,"same"),0,temp_yadv_anm)
temp_yadv_anm=np.apply_along_axis(lambda m:np.convolve(m,wgt,"same"),0,temp_yadv_anm)-lview

lon_temp_vadv,lat_temp_vadv,time_temp_vadv,temp_vadv=select_region_TLL_files(fnames_temp_vadv,varname_temp_vadv,dt1,dt2,lat1,lat2,lon1,lon2,return_dims=True)
temp_vadv_clm=clmmonTLL(temp_vadv) # Get climatology                            
temp_vadv_anm=calcmonanomTLL(temp_vadv,temp_vadv_clm) # Get interannual anomalies         
lview=np.apply_along_axis(lambda m:np.convolve(m,lwgt,"same"),0,temp_vadv_anm)
temp_vadv_anm=np.apply_along_axis(lambda m:np.convolve(m,wgt,"same"),0,temp_vadv_anm)-lview

lon_temp_hdiff,lat_temp_hdiff,time_temp_hdiff,temp_hdiff=select_region_TLL_files(fnames_temp_hdiff,varname_temp_hdiff,dt1,dt2,lat1,lat2,lon1,lon2,return_dims=True)
temp_hdiff_clm=clmmonTLL(temp_hdiff) # Get climatology                          
temp_hdiff_anm=calcmonanomTLL(temp_hdiff,temp_hdiff_clm) # Get interannual anomalies    
lview=np.apply_along_axis(lambda m:np.convolve(m,lwgt,"same"),0,temp_hdiff_anm)
temp_hdiff_anm=np.apply_along_axis(lambda m:np.convolve(m,wgt,"same"),0,temp_hdiff_anm)-lview

lon_temp_vdiff,lat_temp_vdiff,time_temp_vdiff,temp_vdiff=select_region_TLL_files(fnames_temp_vdiff,varname_temp_vdiff,dt1,dt2,lat1,lat2,lon1,lon2,return_dims=True)
temp_vdiff_clm=clmmonTLL(temp_vdiff) # Get climatology                          
temp_vdiff_anm=calcmonanomTLL(temp_vdiff,temp_vdiff_clm) # Get interannual anomalies         
lview=np.apply_along_axis(lambda m:np.convolve(m,lwgt,"same"),0,temp_vdiff_anm)
temp_vdiff_anm=np.apply_along_axis(lambda m:np.convolve(m,wgt,"same"),0,temp_vdiff_anm)-lview                                                                  

lon_temp_entr,lat_temp_entr,time_temp_entr,temp_entr=select_region_TLL_files(fnames_temp_entr,varname_temp_entr,dt1,dt2,lat1,lat2,lon1,lon2,return_dims=True)
temp_entr_clm=clmmonTLL(temp_entr) # Get climatology                            
temp_entr_anm=calcmonanomTLL(temp_entr,temp_entr_clm) # Get interannual anomalies                                                                              
lview=np.apply_along_axis(lambda m:np.convolve(m,lwgt,"same"),0,temp_entr_anm)
temp_entr_anm=np.apply_along_axis(lambda m:np.convolve(m,wgt,"same"),0,temp_entr_anm)-lview

lon_temp_nhf,lat_temp_nhf,time_temp_nhf,temp_nhf=select_region_TLL_files(fnames_temp_nhf,varname_temp_nhf,dt1,dt2,lat1,lat2,lon1,lon2,return_dims=True)
temp_nhf_clm=clmmonTLL(temp_nhf) # Get climatology                              
temp_nhf_anm=calcmonanomTLL(temp_nhf,temp_nhf_clm) # Get interannual anomalies
lview=np.apply_along_axis(lambda m:np.convolve(m,lwgt,"same"),0,temp_nhf_anm)
temp_nhf_anm=np.apply_along_axis(lambda m:np.convolve(m,wgt,"same"),0,temp_nhf_anm)-lview

lon_temp_ssrc,lat_temp_ssrc,time_temp_ssrc,temp_ssrc=select_region_TLL_files(fnames_temp_ssrc,varname_temp_ssrc,dt1,dt2,lat1,lat2,lon1,lon2,return_dims=True)
temp_ssrc_clm=clmmonTLL(temp_ssrc) # Get climatology                                  
temp_ssrc_anm=calcmonanomTLL(temp_ssrc,temp_ssrc_clm) # Get interannual anomalies    
lview=np.apply_along_axis(lambda m:np.convolve(m,lwgt,"same"),0,temp_ssrc_anm)
temp_ssrc_anm=np.apply_along_axis(lambda m:np.convolve(m,wgt,"same"),0,temp_ssrc_anm)-lview

temp_nhf_anm=temp_nhf_anm+temp_ssrc_anm
temp_anm=np.concatenate((temp_rate_anm.reshape(-1,nt,ny,nx),temp_xadv_anm.reshape(-1,nt,ny,nx),temp_yadv_anm.reshape(-1,nt,ny,nx),temp_vadv_anm.reshape(-1,nt,ny,nx),temp_hdiff_anm.reshape(-1,nt,ny,nx),temp_vdiff_anm.reshape(-1,nt,ny,nx),temp_entr_anm.reshape(-1,nt,ny,nx),temp_nhf_anm.reshape(-1,nt,ny,nx)),axis=0)#(8para,8year,18month)
npara=temp_anm.shape[0]
# Slice targeted year
nevent=len(target_years)

# Get targeted time                                                                   
time_temp_rate=np.asarray(time_temp_rate)
time_temp_xadv=np.asarray(time_temp_xadv)
time_temp_yadv=np.asarray(time_temp_yadv)
time_temp_vadv=np.asarray(time_temp_vadv)
time_temp_hdiff=np.asarray(time_temp_hdiff)
time_temp_vdiff=np.asarray(time_temp_vdiff)
time_temp_entr=np.asarray(time_temp_entr)
time_temp_nhf=np.asarray(time_temp_nhf)
time_temp_ssrc=np.asarray(time_temp_ssrc)

temp_comp=np.copy(temp_anm[:,0:nevent,:,:])
ntemp_comp=np.copy(temp_anm[:,0:nevent,:,:])
for i in range(0,nevent):
    # Get time
    target_dt1=dt.datetime(target_years[i],4,1,0,0,0)   # Jun 1st of targeted year
    target_dt2=dt.datetime(target_years[i],9,30,0,0,0) # Augast 30th of targeted year
    itime=np.where((time_temp_rate >= target_dt1) & (time_temp_rate <= target_dt2))[0]
    temp_comp[:,i,:,:]=np.nanmean(temp_anm[:,itime,:,:],axis=1)

for i in range(0,nevent):
    # Get time
    target_dt1=dt.datetime(ntarget_years[i],4,1,0,0,0)   # Jun 1st of targeted year
    target_dt2=dt.datetime(ntarget_years[i],9,30,0,0,0) # Augast 30th of targeted year
    itime=np.where((time_temp_rate >= target_dt1) & (time_temp_rate <= target_dt2))[0]
    ntemp_comp[:,i,:,:]=np.nanmean(temp_anm[:,itime,:,:],axis=1)

# Composited value
view=np.mean(temp_comp[:,:,:,:,],axis=1)
test=np.mean(temp_comp[:,:,:,:,],axis=1)
nview=np.mean(ntemp_comp[:,:,:,:,],axis=1)
ntest=np.mean(ntemp_comp[:,:,:,:,],axis=1)

#statistical significant
temp_t = abs(view)*(nevent)**0.50/(np.std(temp_comp, axis=1)+1.0e-12)-t95
with np.errstate(invalid='ignore'):
    test[temp_t <  0.0] = None
ntemp_t = abs(nview)*(nevent)**0.50/(np.std(ntemp_comp, axis=1)+1.0e-12)-t95
with np.errstate(invalid='ignore'):
    ntest[ntemp_t <  0.0] = None

pnview=view+nview
pntest=view+nview

for i in range(npara):
    t_diff_95=abs(view[i,:,:]+nview[i,:,:])*(nevent*nevent*(nevent*2-2)/nevent*2)**0.50/(nevent*np.std(ntemp_comp[i,:,:,:],axis=0)**2+nevent*np.std(temp_comp[i,:,:,:],axis=0)**2)**0.50-t95
    pntest[i,:,:][t_diff_95 <  0.0] = None

#plot figures
clevs1 = np.linspace(-0.8,0.8,21)
interval1 = np.arange(-2.0,2.0,0.2)
f1=plt.figure(figsize=(16,8))
f1.suptitle("pIOD", x=0.2, y=0.98, size=30, weight=3)
lon_2d,lat_2d=np.meshgrid(lat_temp_rate,lon_temp_rate)
plt.rcParams["font.size"] = 14

plt.subplots_adjust(wspace=0.4, hspace=0.5)

ax = f1.add_subplot(3,3,1, title="(a) MLT tendency",projection=ccrs.PlateCarree(central_longitude=70))
ax.set_xticks([90,100,110], crs=ccrs.PlateCarree())
ax.set_yticks([-10,-5,0], crs=ccrs.PlateCarree())
lon_formatter = LongitudeFormatter(number_format='.0f',dateline_direction_label=True)
lat_formatter = LatitudeFormatter(number_format='.0f',)
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)
ax.set_extent([90, 110, -10, 0], ccrs.PlateCarree())
ax.outline_patch.set_linewidth(1)
ax.tick_params(axis='both', labelsize=14)
fill=plt.contourf(lon_2d,lat_2d,test[0,:,:],levels=clevs1,extend="both",transform=ccrs.PlateCarree(),cmap=cmocean.cm.balance,zorder=0)
plt.contour(lon_2d,lat_2d,view[0,:,:],linewidth=0.05,levels=interval1,colors="black",transform=ccrs.PlateCarree(),zorder=1)
ax.add_feature(cfeature.COASTLINE, linewidth=2.0, facecolor='oldlace', edgecolor='k', zorder=2)

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.1, axes_class=plt.Axes)
cb = plt.colorbar(fill, orientation='vertical', cax=cax)
cb.ax.set_yticklabels(cb.ax.get_yticklabels())
cb.set_ticks([-0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8])

ax = f1.add_subplot(3,3,2, title="(b) Zonal advection",projection=ccrs.PlateCarree(central_longitude=70))
ax.set_xticks([90,100,110], crs=ccrs.PlateCarree())
ax.set_yticks([-10,-5,0], crs=ccrs.PlateCarree())
lon_formatter = LongitudeFormatter(number_format='.0f',dateline_direction_label=True)
lat_formatter = LatitudeFormatter(number_format='.0f',)
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)
ax.set_extent([90, 110, -10, 0], ccrs.PlateCarree())
ax.outline_patch.set_linewidth(1)
ax.tick_params(axis='both', labelsize=14)
fill=plt.contourf(lon_2d,lat_2d,test[1,:,:],levels=clevs1,extend="both",transform=ccrs.PlateCarree(),cmap=cmocean.cm.balance,zorder=0)
plt.contour(lon_2d,lat_2d,view[1,:,:],linewidth=0.05,levels=interval1,colors="black",transform=ccrs.PlateCarree(),zorder=1)
ax.add_feature(cfeature.COASTLINE, linewidth=2.0, facecolor='oldlace', edgecolor='k', zorder=2)

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.1, axes_class=plt.Axes)
cb = plt.colorbar(fill, orientation='vertical', cax=cax)
cb.ax.set_yticklabels(cb.ax.get_yticklabels())
cb.set_ticks([-0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8])

ax = f1.add_subplot(3,3,3, title="(c) Meridional advection",projection=ccrs.PlateCarree(central_longitude=70))
ax.set_xticks([90,100,110], crs=ccrs.PlateCarree())
ax.set_yticks([-10,-5,0], crs=ccrs.PlateCarree())
lon_formatter = LongitudeFormatter(number_format='.0f',dateline_direction_label=True)
lat_formatter = LatitudeFormatter(number_format='.0f',)
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)
ax.set_extent([90, 110, -10, 0], ccrs.PlateCarree())
ax.outline_patch.set_linewidth(1)
ax.tick_params(axis='both', labelsize=14)
fill=plt.contourf(lon_2d,lat_2d,test[2,:,:],levels=clevs1,extend="both",transform=ccrs.PlateCarree(),cmap=cmocean.cm.balance,zorder=0)
plt.contour(lon_2d,lat_2d,view[2,:,:],linewidth=0.05,levels=interval1,colors="black",transform=ccrs.PlateCarree(),zorder=1)
ax.add_feature(cfeature.COASTLINE, linewidth=2.0, facecolor='oldlace', edgecolor='k', zorder=2)

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.1, axes_class=plt.Axes)
cb = plt.colorbar(fill, orientation='vertical', cax=cax)
cb.ax.set_yticklabels(cb.ax.get_yticklabels())
cb.set_ticks([-0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8])
cb.set_label('$\mathrm{℃/month}$')

ax = f1.add_subplot(3,3,4, title="(d) Vertical advection",projection=ccrs.PlateCarree(central_longitude=70))
ax.set_xticks([90,100,110], crs=ccrs.PlateCarree())
ax.set_yticks([-10,-5,0], crs=ccrs.PlateCarree())
lon_formatter = LongitudeFormatter(number_format='.0f',dateline_direction_label=True)
lat_formatter = LatitudeFormatter(number_format='.0f',)
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)
ax.set_extent([90, 110, -10, 0], ccrs.PlateCarree())
ax.outline_patch.set_linewidth(1)
ax.tick_params(axis='both', labelsize=14)
fill=plt.contourf(lon_2d,lat_2d,test[3,:,:],levels=clevs1,extend="both",transform=ccrs.PlateCarree(),cmap=cmocean.cm.balance,zorder=0)
plt.contour(lon_2d,lat_2d,view[3,:,:],linewidth=0.05,levels=interval1,colors="black",transform=ccrs.PlateCarree(),zorder=1)
ax.add_feature(cfeature.COASTLINE, linewidth=2.0, facecolor='oldlace', edgecolor='k', zorder=2)

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.1, axes_class=plt.Axes)
cb = plt.colorbar(fill, orientation='vertical', cax=cax)
cb.ax.set_yticklabels(cb.ax.get_yticklabels())
cb.set_ticks([-0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8])

ax = f1.add_subplot(3,3,5, title="(e) Horizontal diffusion",projection=ccrs.PlateCarree(central_longitude=70))
ax.set_xticks([90,100,110], crs=ccrs.PlateCarree())
ax.set_yticks([-10,-5,0], crs=ccrs.PlateCarree())
lon_formatter = LongitudeFormatter(number_format='.0f',dateline_direction_label=True)
lat_formatter = LatitudeFormatter(number_format='.0f',)
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)
ax.set_extent([90, 110, -10, 0], ccrs.PlateCarree())
ax.outline_patch.set_linewidth(1)
ax.tick_params(axis='both', labelsize=14)
fill=plt.contourf(lon_2d,lat_2d,test[4,:,:],levels=clevs1,extend="both",transform=ccrs.PlateCarree(),cmap=cmocean.cm.balance,zorder=0)
plt.contour(lon_2d,lat_2d,view[4,:,:],linewidth=0.05,levels=interval1,colors="black",transform=ccrs.PlateCarree(),zorder=1)
ax.add_feature(cfeature.COASTLINE, linewidth=2.0, facecolor='oldlace', edgecolor='k', zorder=2)

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.1, axes_class=plt.Axes)
cb = plt.colorbar(fill, orientation='vertical', cax=cax)
cb.ax.set_yticklabels(cb.ax.get_yticklabels())
cb.set_ticks([-0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8])

ax = f1.add_subplot(3,3,6, title="(f) Vertical diffusion",projection=ccrs.PlateCarree(central_longitude=70))
ax.set_xticks([90,100,110], crs=ccrs.PlateCarree())
ax.set_yticks([-10,-5,0], crs=ccrs.PlateCarree())
lon_formatter = LongitudeFormatter(number_format='.0f',dateline_direction_label=True)
lat_formatter = LatitudeFormatter(number_format='.0f',)
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)
ax.set_extent([90, 110, -10, 0], ccrs.PlateCarree())
ax.outline_patch.set_linewidth(1)
ax.tick_params(axis='both', labelsize=14)
fill=plt.contourf(lon_2d,lat_2d,test[5,:,:],levels=clevs1,extend="both",transform=ccrs.PlateCarree(),cmap=cmocean.cm.balance,zorder=0)
plt.contour(lon_2d,lat_2d,view[5,:,:],linewidth=0.05,levels=interval1,colors="black",transform=ccrs.PlateCarree(),zorder=1)
ax.add_feature(cfeature.COASTLINE, linewidth=2.0, facecolor='oldlace', edgecolor='k', zorder=2)

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.1, axes_class=plt.Axes)
cb = plt.colorbar(fill, orientation='vertical', cax=cax)
cb.ax.set_yticklabels(cb.ax.get_yticklabels())
cb.set_ticks([-0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8])
cb.set_label('$\mathrm{℃/month}$')

ax = f1.add_subplot(3,3,7, title="(g) Entrainment",projection=ccrs.PlateCarree(central_longitude=70))
ax.set_xticks([90,100,110], crs=ccrs.PlateCarree())
ax.set_yticks([-10,-5,0], crs=ccrs.PlateCarree())
lon_formatter = LongitudeFormatter(number_format='.0f',dateline_direction_label=True)
lat_formatter = LatitudeFormatter(number_format='.0f',)
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)
ax.set_extent([90, 110, -10, 0], ccrs.PlateCarree())
ax.outline_patch.set_linewidth(1)
ax.tick_params(axis='both', labelsize=14)
fill=plt.contourf(lon_2d,lat_2d,test[6,:,:],levels=clevs1,extend="both",transform=ccrs.PlateCarree(),cmap=cmocean.cm.balance,zorder=0)
plt.contour(lon_2d,lat_2d,view[6,:,:],linewidth=0.05,levels=interval1,colors="black",transform=ccrs.PlateCarree(),zorder=1)
ax.add_feature(cfeature.COASTLINE, linewidth=2.0, facecolor='oldlace', edgecolor='k', zorder=2)

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.1, axes_class=plt.Axes)
cb = plt.colorbar(fill, orientation='vertical', cax=cax)
cb.ax.set_yticklabels(cb.ax.get_yticklabels())
cb.set_ticks([-0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8])

ax = f1.add_subplot(3,3,8, title="(h) Surface heat flux",projection=ccrs.PlateCarree(central_longitude=70))
ax.set_xticks([90,100,110], crs=ccrs.PlateCarree())
ax.set_yticks([-10,-5,0], crs=ccrs.PlateCarree())
lon_formatter = LongitudeFormatter(number_format='.0f',dateline_direction_label=True)
lat_formatter = LatitudeFormatter(number_format='.0f',)
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)
ax.set_extent([90, 110, -10, 0], ccrs.PlateCarree())
ax.outline_patch.set_linewidth(1)
ax.tick_params(axis='both', labelsize=14)
fill=plt.contourf(lon_2d,lat_2d,test[7,:,:],levels=clevs1,extend="both",transform=ccrs.PlateCarree(),cmap=cmocean.cm.balance,zorder=0)
plt.contour(lon_2d,lat_2d,view[7,:,:],linewidth=0.05,levels=interval1,colors="black",transform=ccrs.PlateCarree(),zorder=1)
ax.add_feature(cfeature.COASTLINE, linewidth=2.0, facecolor='oldlace', edgecolor='k', zorder=2)

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.1, axes_class=plt.Axes)
cb = plt.colorbar(fill, orientation='vertical', cax=cax)
cb.ax.set_yticklabels(cb.ax.get_yticklabels())
cb.set_ticks([-0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8])
cb.set_label('$\mathrm{℃/month}$')

f2=plt.figure(figsize=(16,8))
f2.suptitle("nIOD", x=0.2, y=0.98, size=30, weight=3)
lon_2d,lat_2d=np.meshgrid(lat_temp_rate,lon_temp_rate)
plt.rcParams["font.size"] = 14

plt.subplots_adjust(wspace=0.4, hspace=0.5)

ax = f2.add_subplot(3,3,1, title="(i) MLT tendency",projection=ccrs.PlateCarree(central_longitude=70))
ax.set_xticks([90,100,110], crs=ccrs.PlateCarree())
ax.set_yticks([-10,-5,0], crs=ccrs.PlateCarree())
lon_formatter = LongitudeFormatter(number_format='.0f',dateline_direction_label=True)
lat_formatter = LatitudeFormatter(number_format='.0f',)
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)
ax.set_extent([90, 110, -10, 0], ccrs.PlateCarree())
ax.outline_patch.set_linewidth(1)
ax.tick_params(axis='both', labelsize=14)
fill=plt.contourf(lon_2d,lat_2d,ntest[0,:,:],levels=clevs1,extend="both",transform=ccrs.PlateCarree(),cmap=cmocean.cm.balance,zorder=0)
plt.contour(lon_2d,lat_2d,nview[0,:,:],linewidth=0.05,levels=interval1,colors="black",transform=ccrs.PlateCarree(),zorder=1)
ax.add_feature(cfeature.COASTLINE, linewidth=2.0, facecolor='oldlace', edgecolor='k', zorder=2)

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.1, axes_class=plt.Axes)
cb = plt.colorbar(fill, orientation='vertical', cax=cax)
cb.ax.set_yticklabels(cb.ax.get_yticklabels())
cb.set_ticks([-0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8])

ax = f2.add_subplot(3,3,2, title="(j) Zonal advection",projection=ccrs.PlateCarree(central_longitude=70))
ax.set_xticks([90,100,110], crs=ccrs.PlateCarree())
ax.set_yticks([-10,-5,0], crs=ccrs.PlateCarree())
lon_formatter = LongitudeFormatter(number_format='.0f',dateline_direction_label=True)
lat_formatter = LatitudeFormatter(number_format='.0f',)
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)
ax.set_extent([90, 110, -10, 0], ccrs.PlateCarree())
ax.outline_patch.set_linewidth(1)
ax.tick_params(axis='both', labelsize=14)
fill=plt.contourf(lon_2d,lat_2d,ntest[1,:,:],levels=clevs1,extend="both",transform=ccrs.PlateCarree(),cmap=cmocean.cm.balance,zorder=0)
plt.contour(lon_2d,lat_2d,nview[1,:,:],linewidth=0.05,levels=interval1,colors="black",transform=ccrs.PlateCarree(),zorder=1)
ax.add_feature(cfeature.COASTLINE, linewidth=2.0, facecolor='oldlace', edgecolor='k', zorder=2)

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.1, axes_class=plt.Axes)
cb = plt.colorbar(fill, orientation='vertical', cax=cax)
cb.ax.set_yticklabels(cb.ax.get_yticklabels())
cb.set_ticks([-0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8])

ax = f2.add_subplot(3,3,3, title="(k) Meridional advection",projection=ccrs.PlateCarree(central_longitude=70))
ax.set_xticks([90,100,110], crs=ccrs.PlateCarree())
ax.set_yticks([-10,-5,0], crs=ccrs.PlateCarree())
lon_formatter = LongitudeFormatter(number_format='.0f',dateline_direction_label=True)
lat_formatter = LatitudeFormatter(number_format='.0f',)
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)
ax.set_extent([90, 110, -10, 0], ccrs.PlateCarree())
ax.outline_patch.set_linewidth(1)
ax.tick_params(axis='both', labelsize=14)
fill=plt.contourf(lon_2d,lat_2d,ntest[2,:,:],levels=clevs1,extend="both",transform=ccrs.PlateCarree(),cmap=cmocean.cm.balance,zorder=0)
plt.contour(lon_2d,lat_2d,nview[2,:,:],linewidth=0.05,levels=interval1,colors="black",transform=ccrs.PlateCarree(),zorder=1)
ax.add_feature(cfeature.COASTLINE, linewidth=2.0, facecolor='oldlace', edgecolor='k', zorder=2)

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.1, axes_class=plt.Axes)
cb = plt.colorbar(fill, orientation='vertical', cax=cax)
cb.ax.set_yticklabels(cb.ax.get_yticklabels())
cb.set_ticks([-0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8])
cb.set_label('$\mathrm{℃/month}$')

ax = f2.add_subplot(3,3,4, title="(l) Vertical advection",projection=ccrs.PlateCarree(central_longitude=70))
ax.set_xticks([90,100,110], crs=ccrs.PlateCarree())
ax.set_yticks([-10,-5,0], crs=ccrs.PlateCarree())
lon_formatter = LongitudeFormatter(number_format='.0f',dateline_direction_label=True)
lat_formatter = LatitudeFormatter(number_format='.0f',)
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)
ax.set_extent([90, 110, -10, 0], ccrs.PlateCarree())
ax.outline_patch.set_linewidth(1)
ax.tick_params(axis='both', labelsize=14)
fill=plt.contourf(lon_2d,lat_2d,ntest[3,:,:],levels=clevs1,extend="both",transform=ccrs.PlateCarree(),cmap=cmocean.cm.balance,zorder=0)
plt.contour(lon_2d,lat_2d,nview[3,:,:],linewidth=0.05,levels=interval1,colors="black",transform=ccrs.PlateCarree(),zorder=1)
ax.add_feature(cfeature.COASTLINE, linewidth=2.0, facecolor='oldlace', edgecolor='k', zorder=2)

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.1, axes_class=plt.Axes)
cb = plt.colorbar(fill, orientation='vertical', cax=cax)
cb.ax.set_yticklabels(cb.ax.get_yticklabels())
cb.set_ticks([-0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8])

ax = f2.add_subplot(3,3,5, title="(m) Horizontal diffusion",projection=ccrs.PlateCarree(central_longitude=70))
ax.set_xticks([90,100,110], crs=ccrs.PlateCarree())
ax.set_yticks([-10,-5,0], crs=ccrs.PlateCarree())
lon_formatter = LongitudeFormatter(number_format='.0f',dateline_direction_label=True)
lat_formatter = LatitudeFormatter(number_format='.0f',)
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)
ax.set_extent([90, 110, -10, 0], ccrs.PlateCarree())
ax.outline_patch.set_linewidth(1)
ax.tick_params(axis='both', labelsize=14)
fill=plt.contourf(lon_2d,lat_2d,ntest[4,:,:],levels=clevs1,extend="both",transform=ccrs.PlateCarree(),cmap=cmocean.cm.balance,zorder=0)
plt.contour(lon_2d,lat_2d,nview[4,:,:],linewidth=0.05,levels=interval1,colors="black",transform=ccrs.PlateCarree(),zorder=1)
ax.add_feature(cfeature.COASTLINE, linewidth=2.0, facecolor='oldlace', edgecolor='k', zorder=2)

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.1, axes_class=plt.Axes)
cb = plt.colorbar(fill, orientation='vertical', cax=cax)
cb.ax.set_yticklabels(cb.ax.get_yticklabels())
cb.set_ticks([-0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8])

ax = f2.add_subplot(3,3,6, title="(n) Vertical diffusion",projection=ccrs.PlateCarree(central_longitude=70))
ax.set_xticks([90,100,110], crs=ccrs.PlateCarree())
ax.set_yticks([-10,-5,0], crs=ccrs.PlateCarree())
lon_formatter = LongitudeFormatter(number_format='.0f',dateline_direction_label=True)
lat_formatter = LatitudeFormatter(number_format='.0f',)
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)
ax.set_extent([90, 110, -10, 0], ccrs.PlateCarree())
ax.outline_patch.set_linewidth(1)
ax.tick_params(axis='both', labelsize=14)
fill=plt.contourf(lon_2d,lat_2d,ntest[5,:,:],levels=clevs1,extend="both",transform=ccrs.PlateCarree(),cmap=cmocean.cm.balance,zorder=0)
plt.contour(lon_2d,lat_2d,nview[5,:,:],linewidth=0.05,levels=interval1,colors="black",transform=ccrs.PlateCarree(),zorder=1)
ax.add_feature(cfeature.COASTLINE, linewidth=2.0, facecolor='oldlace', edgecolor='k', zorder=2)

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.1, axes_class=plt.Axes)
cb = plt.colorbar(fill, orientation='vertical', cax=cax)
cb.ax.set_yticklabels(cb.ax.get_yticklabels())
cb.set_ticks([-0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8])
cb.set_label('$\mathrm{℃/month}$')

ax = f2.add_subplot(3,3,7, title="(o) Entrainment",projection=ccrs.PlateCarree(central_longitude=70))
ax.set_xticks([90,100,110], crs=ccrs.PlateCarree())
ax.set_yticks([-10,-5,0], crs=ccrs.PlateCarree())
lon_formatter = LongitudeFormatter(number_format='.0f',dateline_direction_label=True)
lat_formatter = LatitudeFormatter(number_format='.0f',)
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)
ax.set_extent([90, 110, -10, 0], ccrs.PlateCarree())
ax.outline_patch.set_linewidth(1)
ax.tick_params(axis='both', labelsize=14)
fill=plt.contourf(lon_2d,lat_2d,ntest[6,:,:],levels=clevs1,extend="both",transform=ccrs.PlateCarree(),cmap=cmocean.cm.balance,zorder=0)
plt.contour(lon_2d,lat_2d,nview[6,:,:],linewidth=0.05,levels=interval1,colors="black",transform=ccrs.PlateCarree(),zorder=1)
ax.add_feature(cfeature.COASTLINE, linewidth=2.0, facecolor='oldlace', edgecolor='k', zorder=2)

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.1, axes_class=plt.Axes)
cb = plt.colorbar(fill, orientation='vertical', cax=cax)
cb.ax.set_yticklabels(cb.ax.get_yticklabels())
cb.set_ticks([-0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8])

ax = f2.add_subplot(3,3,8, title="(p) Surface heat flux",projection=ccrs.PlateCarree(central_longitude=70))
ax.set_xticks([90,100,110], crs=ccrs.PlateCarree())
ax.set_yticks([-10,-5,0], crs=ccrs.PlateCarree())
lon_formatter = LongitudeFormatter(number_format='.0f',dateline_direction_label=True)
lat_formatter = LatitudeFormatter(number_format='.0f',)
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)
ax.set_extent([90, 110, -10, 0], ccrs.PlateCarree())
ax.outline_patch.set_linewidth(1)
ax.tick_params(axis='both', labelsize=14)
fill=plt.contourf(lon_2d,lat_2d,ntest[7,:,:],levels=clevs1,extend="both",transform=ccrs.PlateCarree(),cmap=cmocean.cm.balance,zorder=0)
plt.contour(lon_2d,lat_2d,nview[7,:,:],linewidth=0.05,levels=interval1,colors="black",transform=ccrs.PlateCarree(),zorder=1)
ax.add_feature(cfeature.COASTLINE, linewidth=2.0, facecolor='oldlace', edgecolor='k', zorder=2)

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.1, axes_class=plt.Axes)
cb = plt.colorbar(fill, orientation='vertical', cax=cax)
cb.ax.set_yticklabels(cb.ax.get_yticklabels())
cb.set_ticks([-0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8])
cb.set_label('$\mathrm{℃/month}$')

f3=plt.figure(figsize=(16,8))
f3.suptitle("pIOD + nIOD", x=0.2, y=0.98, size=30, weight=3)
lon_2d,lat_2d=np.meshgrid(lat_temp_rate,lon_temp_rate)
plt.rcParams["font.size"] = 14

plt.subplots_adjust(wspace=0.4, hspace=0.5)

ax = f3.add_subplot(3,3,1, title="(q) MLT tendency",projection=ccrs.PlateCarree(central_longitude=70))
ax.set_xticks([90,100,110], crs=ccrs.PlateCarree())
ax.set_yticks([-10,-5,0], crs=ccrs.PlateCarree())
lon_formatter = LongitudeFormatter(number_format='.0f',dateline_direction_label=True)
lat_formatter = LatitudeFormatter(number_format='.0f',)
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)
ax.set_extent([90, 110, -10, 0], ccrs.PlateCarree())
ax.outline_patch.set_linewidth(1)
ax.tick_params(axis='both', labelsize=14)
fill=plt.contourf(lon_2d,lat_2d,pntest[0,:,:],levels=clevs1,extend="both",transform=ccrs.PlateCarree(),cmap=cmocean.cm.balance,zorder=0)
plt.contour(lon_2d,lat_2d,pnview[0,:,:],linewidth=0.05,levels=interval1,colors="black",transform=ccrs.PlateCarree(),zorder=1)
ax.add_feature(cfeature.COASTLINE, linewidth=2.0, facecolor='oldlace', edgecolor='k', zorder=2)

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.1, axes_class=plt.Axes)
cb = plt.colorbar(fill, orientation='vertical', cax=cax)
cb.ax.set_yticklabels(cb.ax.get_yticklabels())
cb.set_ticks([-0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8])

ax = f3.add_subplot(3,3,2, title="(r) Zonal advection",projection=ccrs.PlateCarree(central_longitude=70))
ax.set_xticks([90,100,110], crs=ccrs.PlateCarree())
ax.set_yticks([-10,-5,0], crs=ccrs.PlateCarree())
lon_formatter = LongitudeFormatter(number_format='.0f',dateline_direction_label=True)
lat_formatter = LatitudeFormatter(number_format='.0f',)
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)
ax.set_extent([90, 110, -10, 0], ccrs.PlateCarree())
ax.outline_patch.set_linewidth(1)
ax.tick_params(axis='both', labelsize=14)
fill=plt.contourf(lon_2d,lat_2d,pntest[1,:,:],levels=clevs1,extend="both",transform=ccrs.PlateCarree(),cmap=cmocean.cm.balance,zorder=0)
plt.contour(lon_2d,lat_2d,pnview[1,:,:],linewidth=0.05,levels=interval1,colors="black",transform=ccrs.PlateCarree(),zorder=1)
ax.add_feature(cfeature.COASTLINE, linewidth=2.0, facecolor='oldlace', edgecolor='k', zorder=2)

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.1, axes_class=plt.Axes)
cb = plt.colorbar(fill, orientation='vertical', cax=cax)
cb.ax.set_yticklabels(cb.ax.get_yticklabels())
cb.set_ticks([-0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8])

ax = f3.add_subplot(3,3,3, title="(s) Meridional advection",projection=ccrs.PlateCarree(central_longitude=70))
ax.set_xticks([90,100,110], crs=ccrs.PlateCarree())
ax.set_yticks([-10,-5,0], crs=ccrs.PlateCarree())
lon_formatter = LongitudeFormatter(number_format='.0f',dateline_direction_label=True)
lat_formatter = LatitudeFormatter(number_format='.0f',)
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)
ax.set_extent([90, 110, -10, 0], ccrs.PlateCarree())
ax.outline_patch.set_linewidth(1)
ax.tick_params(axis='both', labelsize=14)
fill=plt.contourf(lon_2d,lat_2d,pntest[2,:,:],levels=clevs1,extend="both",transform=ccrs.PlateCarree(),cmap=cmocean.cm.balance,zorder=0)
plt.contour(lon_2d,lat_2d,pnview[2,:,:],linewidth=0.05,levels=interval1,colors="black",transform=ccrs.PlateCarree(),zorder=1)
ax.add_feature(cfeature.COASTLINE, linewidth=2.0, facecolor='oldlace', edgecolor='k', zorder=2)

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.1, axes_class=plt.Axes)
cb = plt.colorbar(fill, orientation='vertical', cax=cax)
cb.ax.set_yticklabels(cb.ax.get_yticklabels())
cb.set_ticks([-0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8])
cb.set_label('$\mathrm{℃/month}$')

ax = f3.add_subplot(3,3,4, title="(t) Vertical advection",projection=ccrs.PlateCarree(central_longitude=70))
ax.set_xticks([90,100,110], crs=ccrs.PlateCarree())
ax.set_yticks([-10,-5,0], crs=ccrs.PlateCarree())
lon_formatter = LongitudeFormatter(number_format='.0f',dateline_direction_label=True)
lat_formatter = LatitudeFormatter(number_format='.0f',)
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)
ax.set_extent([90, 110, -10, 0], ccrs.PlateCarree())
ax.outline_patch.set_linewidth(1)
ax.tick_params(axis='both', labelsize=14)
fill=plt.contourf(lon_2d,lat_2d,pntest[3,:,:],levels=clevs1,extend="both",transform=ccrs.PlateCarree(),cmap=cmocean.cm.balance,zorder=0)
plt.contour(lon_2d,lat_2d,pnview[3,:,:],linewidth=0.05,levels=interval1,colors="black",transform=ccrs.PlateCarree(),zorder=1)
ax.add_feature(cfeature.COASTLINE, linewidth=2.0, facecolor='oldlace', edgecolor='k', zorder=2)

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.1, axes_class=plt.Axes)
cb = plt.colorbar(fill, orientation='vertical', cax=cax)
cb.ax.set_yticklabels(cb.ax.get_yticklabels())
cb.set_ticks([-0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8])

ax = f3.add_subplot(3,3,5, title="(u) Horizontal diffusion",projection=ccrs.PlateCarree(central_longitude=70))
ax.set_xticks([90,100,110], crs=ccrs.PlateCarree())
ax.set_yticks([-10,-5,0], crs=ccrs.PlateCarree())
lon_formatter = LongitudeFormatter(number_format='.0f',dateline_direction_label=True)
lat_formatter = LatitudeFormatter(number_format='.0f',)
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)
ax.set_extent([90, 110, -10, 0], ccrs.PlateCarree())
ax.outline_patch.set_linewidth(1)
ax.tick_params(axis='both', labelsize=14)
fill=plt.contourf(lon_2d,lat_2d,pntest[4,:,:],levels=clevs1,extend="both",transform=ccrs.PlateCarree(),cmap=cmocean.cm.balance,zorder=0)
plt.contour(lon_2d,lat_2d,pnview[4,:,:],linewidth=0.05,levels=interval1,colors="black",transform=ccrs.PlateCarree(),zorder=1)
ax.add_feature(cfeature.COASTLINE, linewidth=2.0, facecolor='oldlace', edgecolor='k', zorder=2)

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.1, axes_class=plt.Axes)
cb = plt.colorbar(fill, orientation='vertical', cax=cax)
cb.ax.set_yticklabels(cb.ax.get_yticklabels())
cb.set_ticks([-0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8])

ax = f3.add_subplot(3,3,6, title="(v) Vertical diffusion",projection=ccrs.PlateCarree(central_longitude=70))
ax.set_xticks([90,100,110], crs=ccrs.PlateCarree())
ax.set_yticks([-10,-5,0], crs=ccrs.PlateCarree())
lon_formatter = LongitudeFormatter(number_format='.0f',dateline_direction_label=True)
lat_formatter = LatitudeFormatter(number_format='.0f',)
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)
ax.set_extent([90, 110, -10, 0], ccrs.PlateCarree())
ax.outline_patch.set_linewidth(1)
ax.tick_params(axis='both', labelsize=14)
fill=plt.contourf(lon_2d,lat_2d,pntest[5,:,:],levels=clevs1,extend="both",transform=ccrs.PlateCarree(),cmap=cmocean.cm.balance,zorder=0)
plt.contour(lon_2d,lat_2d,pnview[5,:,:],linewidth=0.05,levels=interval1,colors="black",transform=ccrs.PlateCarree(),zorder=1)
ax.add_feature(cfeature.COASTLINE, linewidth=2.0, facecolor='oldlace', edgecolor='k', zorder=2)

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.1, axes_class=plt.Axes)
cb = plt.colorbar(fill, orientation='vertical', cax=cax)
cb.ax.set_yticklabels(cb.ax.get_yticklabels())
cb.set_ticks([-0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8])
cb.set_label('$\mathrm{℃/month}$')

ax = f3.add_subplot(3,3,7, title="(w) Entrainment",projection=ccrs.PlateCarree(central_longitude=70))
ax.set_xticks([90,100,110], crs=ccrs.PlateCarree())
ax.set_yticks([-10,-5,0], crs=ccrs.PlateCarree())
lon_formatter = LongitudeFormatter(number_format='.0f',dateline_direction_label=True)
lat_formatter = LatitudeFormatter(number_format='.0f',)
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)
ax.set_extent([90, 110, -10, 0], ccrs.PlateCarree())
ax.outline_patch.set_linewidth(1)
ax.tick_params(axis='both', labelsize=14)
fill=plt.contourf(lon_2d,lat_2d,pntest[6,:,:],levels=clevs1,extend="both",transform=ccrs.PlateCarree(),cmap=cmocean.cm.balance,zorder=0)
plt.contour(lon_2d,lat_2d,pnview[6,:,:],linewidth=0.05,levels=interval1,colors="black",transform=ccrs.PlateCarree(),zorder=1)
ax.add_feature(cfeature.COASTLINE, linewidth=2.0, facecolor='oldlace', edgecolor='k', zorder=2)

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.1, axes_class=plt.Axes)
cb = plt.colorbar(fill, orientation='vertical', cax=cax)
cb.ax.set_yticklabels(cb.ax.get_yticklabels())
cb.set_ticks([-0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8])

ax = f3.add_subplot(3,3,8, title="(x) Surface heat flux",projection=ccrs.PlateCarree(central_longitude=70))
ax.set_xticks([90,100,110], crs=ccrs.PlateCarree())
ax.set_yticks([-10,-5,0], crs=ccrs.PlateCarree())
lon_formatter = LongitudeFormatter(number_format='.0f',dateline_direction_label=True)
lat_formatter = LatitudeFormatter(number_format='.0f',)
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)
ax.set_extent([90, 110, -10, 0], ccrs.PlateCarree())
ax.outline_patch.set_linewidth(1)
ax.tick_params(axis='both', labelsize=14)
fill=plt.contourf(lon_2d,lat_2d,pntest[7,:,:],levels=clevs1,extend="both",transform=ccrs.PlateCarree(),cmap=cmocean.cm.balance,zorder=0)
plt.contour(lon_2d,lat_2d,pnview[7,:,:],linewidth=0.05,levels=interval1,colors="black",transform=ccrs.PlateCarree(),zorder=1)
ax.add_feature(cfeature.COASTLINE, linewidth=2.0, facecolor='oldlace', edgecolor='k', zorder=2)

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.1, axes_class=plt.Axes)
cb = plt.colorbar(fill, orientation='vertical', cax=cax)
cb.ax.set_yticklabels(cb.ax.get_yticklabels())
cb.set_ticks([-0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8])
cb.set_label('$\mathrm{℃/month}$')

plt.savefig(figname)
print(figname)
'''
from numpy import dtype
nc = Dataset('/Volumes/Promise-Pegasus_2/nakazato/HEATBUD/comp_map_all_posi_0123.nc', 'w')
#nc.createDimension('time1', None)
nc.createDimension('lat', ny)
nc.createDimension('lon', nx)
nc.createDimension('para', 8)
nc.createDimension('nevent', 8)

lat = nc.createVariable('lat', 'f', ('lat'))
lat.long_name = 'latitude'
lat.units = 'degrees north'
lat.axis = 'Y'
lat[:] = lat_temp_rate

lon = nc.createVariable('lon', 'f', ('lon'))
lon.long_name = 'longitude'
lon.units = 'degrees east'
lon.axis = 'X'
lon[:] = lon_temp_rate

para = nc.createVariable('para', 'f', ('para'))
para.long_name = 'parameter'
para.units = ''
para.axis = 'p'
para[:] = np.arange(8)

nevent = nc.createVariable('nevent', 'f', ('nevent'))
nevent.long_name = 'IOD year'
nevent.units = 'year'
nevent.axis = 'yr'
nevent[:] = np.arange(8)


comp = nc.createVariable('comp', 'f', ('para', 'nevent', 'lat', 'lon'))
comp.long_name = 'compsite'
comp.units = '℃/month'
comp[:,:,:,:] = temp_comp[:,:,:,:]*60*60*24*30

view = nc.createVariable('view', 'f', ('para',  'lat', 'lon'))
view.long_name = 'compsite'
view.units = '℃/month'
view[:,:,:] = temp_view[:,:,:]*60*60*24*30

ttest = nc.createVariable('ttest', 'f', ('para',  'lat', 'lon'))
ttest.long_name = 'ttest'
ttest.units = '℃/month'
ttest[:,:,:] = temp_t[:,:,:]*60*60*24*30

nc.close()
print("comp_map_all_posi_0123.nc")

temp_sum_view=temp_xadv_view+temp_yadv_view+temp_vadv_view+temp_hdiff_view+temp_vdiff_view+temp_entr_view+temp_nhf_view+temp_ssrc_view


#plot
clevs=np.linspace(-2.0*pow(10,-7),2.0*pow(10,-7),num=21)
clevs1=np.linspace(-2.0*pow(10,-9),2.0*pow(10,-9),num=21)

lon_2d,lat_2d=np.meshgrid(lon_temp_rate,lat_temp_rate)
fig=plt.figure(figsize=(13.0,6.0))
plt.subplots_adjust(wspace=0.1, hspace=0.6)
plt.rcParams["font.size"] = 10

#ax=plt.subplot(431,projection=ccrs.PlateCarree(central_longitude=75))
plt.subplot(4,3,1)
plt.contourf(lon_2d,lat_2d,temp_rate_view,levels=clevs,extend="both",cmap="bwr")#,transform=ccrs.PlateCarree(),cmap="bwr")
#ax.add_feature(cfeature.COASTLINE, linewidth=2.0)
plt.title("SST tendency "+sign)
plt.colorbar()
#plt.gcf().coastlines()

plt.subplot(4,3,2)
plt.contourf(lon_2d,lat_2d,temp_sum_view,levels=clevs,extend="both",cmap="bwr")
plt.title("Sum of RHS "+sign)
plt.colorbar()

plt.subplot(4,3,3)
plt.contourf(lon_2d,lat_2d,temp_xadv_view+temp_yadv_view,levels=clevs,extend="both",cmap="bwr")
plt.title("Zonal + Meridional advection "+sign)
plt.colorbar()

plt.subplot(4,3,4)
plt.contourf(lon_2d,lat_2d,temp_xadv_view,levels=clevs,extend="both",cmap="bwr")
plt.title("Zonal advection "+sign)
plt.colorbar()

plt.subplot(4,3,5)
plt.contourf(lon_2d,lat_2d,temp_yadv_view,levels=clevs,extend="both",cmap="bwr")
plt.title("Meridional advection "+sign)
plt.colorbar()

plt.subplot(4,3,6)
plt.contourf(lon_2d,lat_2d,temp_vadv_view,levels=clevs,extend="both",cmap="bwr")
plt.title("Vertical advection "+sign)
plt.colorbar()

plt.subplot(4,3,7)
plt.contourf(lon_2d,lat_2d,temp_hdiff_view,levels=clevs1,extend="both",cmap="bwr")
plt.title("Horizontal diff. "+sign)
plt.colorbar()

plt.subplot(4,3,8)
plt.contourf(lon_2d,lat_2d,temp_vdiff_view,levels=clevs,extend="both",cmap="bwr")
plt.title("Vertical diff. "+sign)
plt.colorbar()


plt.subplot(4,3,9)
plt.contourf(lon_2d,lat_2d,temp_entr_view,levels=clevs,extend="both",cmap="bwr")
plt.title("Vertical entr. "+sign)
plt.colorbar()

plt.subplot(4,3,10)
plt.contourf(lon_2d,lat_2d,temp_nhf_view,levels=clevs,extend="both",cmap="bwr")
plt.title("NHF "+sign)
plt.colorbar()

plt.subplot(4,3,11)
plt.contourf(lon_2d,lat_2d,temp_ssrc_view,levels=clevs,extend="both",cmap="bwr")
plt.title("Penet "+sign)
plt.colorbar()

plt.savefig(figname)
print(figname)
'''

