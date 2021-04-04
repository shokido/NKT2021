# -*- coding: utf-8 -*-
import numpy.ma as ma
import numpy as np
from ncdf_read import select_region_TLL_files
from stat_ncl import clmmonTLL,calcmonanomTLL
import datetime as dt
import matplotlib
import matplotlib.pyplot as plt
from textwrap import wrap
matplotlib.rcParams['font.family'] = 'YuGothic'

# Python script for displaying area-averaged
# ML-heat budget term
# integration from April to September
# Figure 2

dir_in="/Volumes/Promise-Pegasus_2/kido/ROMS/Outputs/IO_033/ZLEV/HEATBUD/"

#east box
lon1=90.0
lon2=110.0
lat1=-10.0
lat2=0.0

start_year=1959
end_year=2017
sec_to_month=60.0*60.0*24.0*30.0

# Targetedv events                                                            
target_years=[1961, 1963, 1967, 1972, 1977, 1994, 1997, 2006]#posi
ntarget_years=[1960, 1964, 1974, 1975, 1996, 1998, 2005, 2010]#nega

figname="integration_heatbud.png"

dt1=dt.datetime(start_year,1,1,0,0,0) # 1959-01-01 00:00
dt2=dt.datetime(end_year,12,31,0,0,0) # 2017-12-31 00:00

#read file
fnames_temp_rate=[]

for iyear in range(start_year,end_year+1):
    fnames_temp_rate.append(dir_in+"roms_hb_"+str(iyear)+"_monthly.nc")

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
nt=temp_rate.shape[0]
ny=temp_rate.shape[1]
nx=temp_rate.shape[2]
temp_rate_clm=clmmonTLL(temp_rate) # Get climatology
temp_rate_anm=calcmonanomTLL(temp_rate,temp_rate_clm) # Get interannual anomalies
#3month smooth
nwgt=3
wgt=np.ones(nwgt)/nwgt
wgt=np.ones(3)
wgt=wgt/np.sum(wgt)
#72month smooth
lnwgt=72
lwgt=np.ones(lnwgt)/lnwgt
lwgt=np.ones(72)
lwgt=lwgt/np.sum(lwgt)
lview=np.apply_along_axis(lambda m:np.convolve(m,lwgt,"same"),0,temp_rate_anm)#長周期成分(72ヶ月）
temp_rate_anm=np.apply_along_axis(lambda m:np.convolve(m,wgt,"same"),0,temp_rate_anm)-lview#3ヶ月移動平均から長周期成分を除去

lon_temp_xadv,lat_temp_xadv,time_temp_xadv,temp_xadv=select_region_TLL_files(fnames_temp_xadv,varname_temp_xadv,dt1,dt2,lat1,lat2,lon1,lon2,return_dims=True)
temp_xadv_clm=clmmonTLL(temp_xadv) # Get climatology
temp_xadv_anm=calcmonanomTLL(temp_xadv,temp_xadv_clm) # Get interannual anomalies
lview=np.apply_along_axis(lambda m:np.convolve(m,lwgt,"same"),0,temp_xadv_anm)#長周期成分
temp_xadv_anm=np.apply_along_axis(lambda m:np.convolve(m,wgt,"same"),0,temp_xadv_anm)-lview#3ヶ月移動平均から長周期成分を除去

lon_temp_yadv,lat_temp_yadv,time_temp_yadv,temp_yadv=select_region_TLL_files(fnames_temp_yadv,varname_temp_yadv,dt1,dt2,lat1,lat2,lon1,lon2,return_dims=True)
temp_yadv_clm=clmmonTLL(temp_yadv) # Get climatology
temp_yadv_anm=calcmonanomTLL(temp_yadv,temp_yadv_clm) # Get interannual anomalies
lview=np.apply_along_axis(lambda m:np.convolve(m,lwgt,"same"),0,temp_yadv_anm)#長周期成分
temp_yadv_anm=np.apply_along_axis(lambda m:np.convolve(m,wgt,"same"),0,temp_yadv_anm)-lview#3ヶ月移動平均から長周期成分を除去

lon_temp_vadv,lat_temp_vadv,time_temp_vadv,temp_vadv=select_region_TLL_files(fnames_temp_vadv,varname_temp_vadv,dt1,dt2,lat1,lat2,lon1,lon2,return_dims=True)
temp_vadv_clm=clmmonTLL(temp_vadv) # Get climatology
temp_vadv_anm=calcmonanomTLL(temp_vadv,temp_vadv_clm) # Get interannual anomalies
lview=np.apply_along_axis(lambda m:np.convolve(m,lwgt,"same"),0,temp_vadv_anm)#長周期成分
temp_vadv_anm=np.apply_along_axis(lambda m:np.convolve(m,wgt,"same"),0,temp_vadv_anm)-lview#3ヶ月移動平均から長周期成分を除去

lon_temp_hdiff,lat_temp_hdiff,time_temp_hdiff,temp_hdiff=select_region_TLL_files(fnames_temp_hdiff,varname_temp_hdiff,dt1,dt2,lat1,lat2,lon1,lon2,return_dims=True)
temp_hdiff_clm=clmmonTLL(temp_hdiff) # Get climatology
temp_hdiff_anm=calcmonanomTLL(temp_hdiff,temp_hdiff_clm) # Get interannual anomalies
lview=np.apply_along_axis(lambda m:np.convolve(m,lwgt,"same"),0,temp_hdiff_anm)#長周期成分
temp_hdiff_anm=np.apply_along_axis(lambda m:np.convolve(m,wgt,"same"),0,temp_hdiff_anm)-lview#3ヶ月移動平均から長周期成分を除去

lon_temp_vdiff,lat_temp_vdiff,time_temp_vdiff,temp_vdiff=select_region_TLL_files(fnames_temp_vdiff,varname_temp_vdiff,dt1,dt2,lat1,lat2,lon1,lon2,return_dims=True)
temp_vdiff_clm=clmmonTLL(temp_vdiff) # Get climatology
temp_vdiff_anm=calcmonanomTLL(temp_vdiff,temp_vdiff_clm) # Get interannual anomalies
lview=np.apply_along_axis(lambda m:np.convolve(m,lwgt,"same"),0,temp_vdiff_anm)#長周期成分
temp_vdiff_anm=np.apply_along_axis(lambda m:np.convolve(m,wgt,"same"),0,temp_vdiff_anm)-lview#3ヶ月移動平均から長周期成分を除去

lon_temp_entr,lat_temp_entr,time_temp_entr,temp_entr=select_region_TLL_files(fnames_temp_entr,varname_temp_entr,dt1,dt2,lat1,lat2,lon1,lon2,return_dims=True)
temp_entr_clm=clmmonTLL(temp_entr) # Get climatology
temp_entr_anm=calcmonanomTLL(temp_entr,temp_entr_clm) # Get interannual anomalies
lview=np.apply_along_axis(lambda m:np.convolve(m,lwgt,"same"),0,temp_entr_anm)#長周期成分
temp_entr_anm=np.apply_along_axis(lambda m:np.convolve(m,wgt,"same"),0,temp_entr_anm)-lview#3ヶ月移動平均から長周期成分を除去

lon_temp_nhf,lat_temp_nhf,time_temp_nhf,temp_nhf=select_region_TLL_files(fnames_temp_nhf,varname_temp_nhf,dt1,dt2,lat1,lat2,lon1,lon2,return_dims=True)
temp_nhf_clm=clmmonTLL(temp_nhf) # Get climatology
temp_nhf_anm=calcmonanomTLL(temp_nhf,temp_nhf_clm) # Get interannual anomalies
lview=np.apply_along_axis(lambda m:np.convolve(m,lwgt,"same"),0,temp_nhf_anm)#長周期成分
temp_nhf_anm=np.apply_along_axis(lambda m:np.convolve(m,wgt,"same"),0,temp_nhf_anm)-lview#3ヶ月移動平均から長周期成分を除去

lon_temp_ssrc,lat_temp_ssrc,time_temp_ssrc,temp_ssrc=select_region_TLL_files(fnames_temp_ssrc,varname_temp_ssrc,dt1,dt2,lat1,lat2,lon1,lon2,return_dims=True)
temp_ssrc_clm=clmmonTLL(temp_ssrc) # Get climatology
temp_ssrc_anm=calcmonanomTLL(temp_ssrc,temp_ssrc_clm) # Get interannual anomalies
lview=np.apply_along_axis(lambda m:np.convolve(m,lwgt,"same"),0,temp_ssrc_anm)#長周期成分
temp_ssrc_anm=np.apply_along_axis(lambda m:np.convolve(m,wgt,"same"),0,temp_ssrc_anm)-lview#3ヶ月移動平均から長周期成分を除去

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


# Slice targeted year                                                           
nevent=len(target_years)

#posi_comp
temp_rate_comp1=ma.empty((nevent,18,30,60))#(IODyear, develop month, lat, lon)
for i in range(0,nevent):
    itime1=np.where((time_temp_rate >= dt.datetime(target_years[i],4,1,0,0,0)) & (time_temp_rate <= dt.datetime(target_years[i]+1,9,30,0,0,0)))
    itime2=list(itime1)
    itime=np.squeeze(itime2)
    temp_rate_comp1[i,:,:,:]=temp_rate_anm[itime,:,:]
temp_rate_comp=np.mean(np.nanmean(temp_rate_comp1,axis=3),axis=2)#area average
temp_rate_view=np.mean(temp_rate_comp,axis=0)#year average

temp_xadv_comp1=ma.empty((nevent,18,30,60))
for i in range(0,nevent):
    itime1=np.where((time_temp_xadv >= dt.datetime(target_years[i],4,1,0,0,0)) & (time_temp_xadv <= dt.datetime(target_years[i]+1,9,30,0,0,0)))
    itime2=list(itime1)
    itime=np.squeeze(itime2)
    temp_xadv_comp1[i,:,:,:]=temp_xadv_anm[itime,:,:]
temp_xadv_comp=np.mean(np.nanmean(temp_xadv_comp1,axis=3),axis=2)
temp_xadv_view=np.mean(temp_xadv_comp,axis=0)

temp_yadv_comp1=ma.empty((nevent,18,30,60))
for i in range(0,nevent):
    itime1=np.where((time_temp_yadv >= dt.datetime(target_years[i],4,1,0,0,0)) & (time_temp_yadv <= dt.datetime(target_years[i]+1,9,30,0,0,0)))
    itime2=list(itime1)
    itime=np.squeeze(itime2)
    temp_yadv_comp1[i,:,:,:]=temp_yadv_anm[itime,:,:]
temp_yadv_comp=np.mean(np.nanmean(temp_yadv_comp1,axis=3),axis=2)
temp_yadv_view=np.mean(temp_yadv_comp,axis=0)

temp_hadv_view=temp_xadv_view+temp_yadv_view

temp_vadv_comp1=ma.empty((nevent,18,30,60))
for i in range(0,nevent):
    itime1=np.where((time_temp_vadv >= dt.datetime(target_years[i],4,1,0,0,0)) & (time_temp_vadv <= dt.datetime(target_years[i]+1,9,30,0,0,0)))
    itime2=list(itime1)
    itime=np.squeeze(itime2)
    temp_vadv_comp1[i,:,:,:]=temp_vadv_anm[itime,:,:]
temp_vadv_comp=np.mean(np.nanmean(temp_vadv_comp1,axis=3),axis=2)
temp_vadv_view=np.mean(temp_vadv_comp,axis=0)

temp_hdiff_comp1=ma.empty((nevent,18,30,60))
for i in range(0,nevent):
    itime1=np.where((time_temp_hdiff >= dt.datetime(target_years[i],4,1,0,0,0)) & (time_temp_hdiff <= dt.datetime(target_years[i]+1,9,30,0,0,0)))
    itime2=list(itime1)
    itime=np.squeeze(itime2)
    temp_hdiff_comp1[i,:,:,:]=temp_hdiff_anm[itime,:,:]
temp_hdiff_comp=np.mean(np.nanmean(temp_hdiff_comp1,axis=3),axis=2)
temp_hdiff_view=np.mean(temp_hdiff_comp,axis=0)

temp_vdiff_comp1=ma.empty((nevent,18,30,60))
for i in range(0,nevent):
    itime1=np.where((time_temp_vdiff >= dt.datetime(target_years[i],4,1,0,0,0)) & (time_temp_vdiff <= dt.datetime(target_years[i]+1,9,30,0,0,0)))
    itime2=list(itime1)
    itime=np.squeeze(itime2)
    temp_vdiff_comp1[i,:,:,:]=temp_vdiff_anm[itime,:,:]
temp_vdiff_comp=np.mean(np.nanmean(temp_vdiff_comp1,axis=3),axis=2)
temp_vdiff_view=np.mean(temp_vdiff_comp,axis=0)

temp_entr_comp1=ma.empty((nevent,18,30,60))
for i in range(0,nevent):
    itime1=np.where((time_temp_entr >= dt.datetime(target_years[i],4,1,0,0,0)) & (time_temp_entr <= dt.datetime(target_years[i]+1,9,30,0,0,0)))
    itime2=list(itime1)
    itime=np.squeeze(itime2)
    temp_entr_comp1[i,:,:,:]=temp_entr_anm[itime,:,:]
temp_entr_comp=np.mean(np.nanmean(temp_entr_comp1,axis=3),axis=2)
temp_entr_view=np.mean(temp_entr_comp,axis=0)

temp_nhf_comp1=ma.empty((nevent,18,30,60))
for i in range(0,nevent):
    itime1=np.where((time_temp_nhf >= dt.datetime(target_years[i],4,1,0,0,0)) & (time_temp_nhf <= dt.datetime(target_years[i]+1,9,30,0,0,0)))
    itime2=list(itime1)
    itime=np.squeeze(itime2)
    temp_nhf_comp1[i,:,:,:]=temp_nhf_anm[itime,:,:]
temp_nhf_comp=np.mean(np.nanmean(temp_nhf_comp1,axis=3),axis=2)
temp_nhf_view=np.mean(temp_nhf_comp,axis=0)

temp_ssrc_comp1=ma.empty((nevent,18,30,60))
for i in range(0,nevent):
    itime1=np.where((time_temp_ssrc >= dt.datetime(target_years[i],4,1,0,0,0)) & (time_temp_ssrc <= dt.datetime(target_years[i]+1,9,30,0,0,0)))
    itime2=list(itime1)
    itime=np.squeeze(itime2)
    temp_ssrc_comp1[i,:,:,:]=temp_ssrc_anm[itime,:,:]
temp_ssrc_comp=np.mean(np.nanmean(temp_ssrc_comp1,axis=3),axis=2)
temp_ssrc_view=np.mean(temp_ssrc_comp,axis=0)

#nega_comp                                                                                                                       
temp_rate_comp2=ma.empty((nevent,18,30,60))
for i in range(0,nevent):
    nitime1=np.where((time_temp_rate >= dt.datetime(ntarget_years[i],4,1,0,0,0)) & (time_temp_rate <= dt.datetime(ntarget_years[i]+1,9,30,0,0,0)))
    nitime2=list(nitime1)
    nitime=np.squeeze(nitime2)
    temp_rate_comp2[i,:,:,:]=temp_rate_anm[nitime,:,:]
temp_rate_compn=np.mean(np.nanmean(temp_rate_comp2,axis=3),axis=2)
temp_rate_viewn=np.mean(temp_rate_compn,axis=0)

temp_xadv_comp2=ma.empty((nevent,18,30,60))
for i in range(0,nevent):
    nitime1=np.where((time_temp_xadv >= dt.datetime(ntarget_years[i],4,1,0,0,0)) & (time_temp_xadv <= dt.datetime(ntarget_years[i]+1,9,30,0,0,0)))
    nitime2=list(nitime1)
    nitime=np.squeeze(nitime2)
    temp_xadv_comp2[i,:,:,:]=temp_xadv_anm[nitime,:,:]
temp_xadv_compn=np.mean(np.nanmean(temp_xadv_comp2,axis=3),axis=2)
temp_xadv_viewn=np.mean(temp_xadv_compn,axis=0)

temp_yadv_comp2=ma.empty((nevent,18,30,60))
for i in range(0,nevent):
    nitime1=np.where((time_temp_yadv >= dt.datetime(ntarget_years[i],4,1,0,0,0)) & (time_temp_yadv <= dt.datetime(ntarget_years[i]+1,9,30,0,0,0)))
    nitime2=list(nitime1)
    nitime=np.squeeze(nitime2)
    temp_yadv_comp2[i,:,:,:]=temp_yadv_anm[nitime,:,:]
temp_yadv_compn=np.mean(np.nanmean(temp_yadv_comp2,axis=3),axis=2)
temp_yadv_viewn=np.mean(temp_yadv_compn,axis=0)

temp_hadv_viewn=temp_xadv_viewn+temp_yadv_viewn

temp_vadv_comp2=ma.empty((nevent,18,30,60))
for i in range(0,nevent):
    nitime1=np.where((time_temp_vadv >= dt.datetime(ntarget_years[i],4,1,0,0,0)) & (time_temp_vadv <= dt.datetime(ntarget_years[i]+1,9,30,0,0,0)))
    nitime2=list(nitime1)
    nitime=np.squeeze(nitime2)
    temp_vadv_comp2[i,:,:,:]=temp_vadv_anm[nitime,:,:]
temp_vadv_compn=np.mean(np.nanmean(temp_vadv_comp2,axis=3),axis=2)
temp_vadv_viewn=np.mean(temp_vadv_compn,axis=0)

temp_hdiff_comp2=ma.empty((nevent,18,30,60))
for i in range(0,nevent):
    nitime1=np.where((time_temp_hdiff >= dt.datetime(ntarget_years[i],4,1,0,0,0)) & (time_temp_hdiff <= dt.datetime(ntarget_years[i]+1,9,30,0,0,0)))
    nitime2=list(nitime1)
    nitime=np.squeeze(nitime2)
    temp_hdiff_comp2[i,:,:,:]=temp_hdiff_anm[nitime,:,:]
temp_hdiff_compn=np.mean(np.nanmean(temp_hdiff_comp2,axis=3),axis=2)
temp_hdiff_viewn=np.mean(temp_hdiff_compn,axis=0)

temp_vdiff_comp2=ma.empty((nevent,18,30,60))
for i in range(0,nevent):
    nitime1=np.where((time_temp_vdiff >= dt.datetime(ntarget_years[i],4,1,0,0,0)) & (time_temp_vdiff <= dt.datetime(ntarget_years[i]+1,9,30,0,0,0)))
    nitime2=list(nitime1)
    nitime=np.squeeze(nitime2)
    temp_vdiff_comp2[i,:,:,:]=temp_vdiff_anm[nitime,:,:]
temp_vdiff_compn=np.mean(np.nanmean(temp_vdiff_comp2,axis=3),axis=2)
temp_vdiff_viewn=np.mean(temp_vdiff_compn,axis=0)

temp_entr_comp2=ma.empty((nevent,18,30,60))
for i in range(0,nevent):
    nitime1=np.where((time_temp_entr >= dt.datetime(ntarget_years[i],4,1,0,0,0)) & (time_temp_entr <= dt.datetime(ntarget_years[i]+1,9,30,0,0,0)))
    nitime2=list(nitime1)
    nitime=np.squeeze(nitime2)
    temp_entr_comp2[i,:,:,:]=temp_entr_anm[nitime,:,:]
temp_entr_compn=np.mean(np.nanmean(temp_entr_comp2,axis=3),axis=2)
temp_entr_viewn=np.mean(temp_entr_compn,axis=0)

temp_nhf_comp2=ma.empty((nevent,18,30,60))
for i in range(0,nevent):
    nitime1=np.where((time_temp_nhf >= dt.datetime(ntarget_years[i],4,1,0,0,0)) & (time_temp_nhf <= dt.datetime(ntarget_years[i]+1,9,30,0,0,0)))
    nitime2=list(nitime1)
    nitime=np.squeeze(nitime2)
    temp_nhf_comp2[i,:,:,:]=temp_nhf_anm[nitime,:,:]
temp_nhf_compn=np.mean(np.nanmean(temp_nhf_comp2,axis=3),axis=2)
temp_nhf_viewn=np.mean(temp_nhf_compn,axis=0)

temp_ssrc_comp2=ma.empty((nevent,18,30,60))
for i in range(0,nevent):
    nitime1=np.where((time_temp_ssrc >= dt.datetime(ntarget_years[i],4,1,0,0,0)) & (time_temp_ssrc <= dt.datetime(ntarget_years[i]+1,9,30,0,0,0)))
    nitime2=list(nitime1)
    nitime=np.squeeze(nitime2)
    temp_ssrc_comp2[i,:,:,:]=temp_ssrc_anm[nitime,:,:]
temp_ssrc_compn=np.mean(np.nanmean(temp_ssrc_comp2,axis=3),axis=2)
temp_ssrc_viewn=np.mean(temp_ssrc_compn,axis=0)

#parameter concatenation
temp_comp=np.concatenate((temp_rate_comp.reshape(-1,18),temp_xadv_comp.reshape(-1,18),temp_yadv_comp.reshape(-1,18),temp_vadv_comp.reshape(-1,18),temp_hdiff_comp.reshape(-1,18),temp_vdiff_comp.reshape(-1,18),temp_entr_comp.reshape(-1,18),temp_nhf_comp.reshape(-1,18)+temp_ssrc_comp.reshape(-1,18)),axis=0)#(8para,8year,18month)
temp_comp=temp_comp.reshape(8,8,18)


temp_compn=np.concatenate((temp_rate_compn.reshape(-1,18),temp_xadv_compn.reshape(-1,18),temp_yadv_compn.reshape(-1,18),temp_vadv_compn.reshape(-1,18),temp_hdiff_compn.reshape(-1,18),temp_vdiff_compn.reshape(-1,18),temp_entr_compn.reshape(-1,18),temp_nhf_compn.reshape(-1,18)+temp_ssrc_compn.reshape(-1,18)),axis=0)
temp_compn=temp_compn.reshape(8,8,18)

temp_view=np.concatenate((temp_rate_view.reshape(-1,18),temp_xadv_view.reshape(-1,18),temp_yadv_view.reshape(-1,18),temp_vadv_view.reshape(-1,18),temp_hdiff_view.reshape(-1,18),temp_vdiff_view.reshape(-1,18),temp_entr_view.reshape(-1,18),temp_nhf_view.reshape(-1,18)+temp_ssrc_view.reshape(-1,18)),axis=0)#(8para,18month)

temp_viewn=np.concatenate((temp_rate_viewn.reshape(-1,18),temp_xadv_viewn.reshape(-1,18),temp_yadv_viewn.reshape(-1,18),temp_vadv_viewn.reshape(-1,18),temp_hdiff_viewn.reshape(-1,18),temp_vdiff_viewn.reshape(-1,18),temp_entr_viewn.reshape(-1,18),temp_nhf_viewn.reshape(-1,18)+temp_ssrc_viewn.reshape(-1,18)),axis=0)
npara=len(temp_view)


#integration from April to September
growth_comp=ma.empty((npara,8,6))
ngrowth_comp=ma.empty((npara,8,6))
growth_view=ma.empty((npara,6))#use in t_test
ngrowth_view=ma.empty((npara,6))#use in t_test
growth_view1=ma.empty((npara,6))
ngrowth_view1=ma.empty((npara,6))
nday=np.array([30,31,30,31,31,30])
for i in range(npara):
    for k in range(6):
        growth_comp[i,:,k]=np.sum(temp_comp[i,:,0:k+1],axis=1)
        ngrowth_comp[i,:,k]=np.sum(temp_compn[i,:,0:k+1],axis=1)
    
for i in range(npara):
    for k in range(6):
        growth_view[i,k]=np.sum(temp_view[i,0:k+1])
        ngrowth_view[i,k]=np.sum(temp_viewn[i,0:k+1])
    growth_view1[i,:]=growth_view[i,:]*sec_to_month
    ngrowth_view1[i,:]=ngrowth_view[i,:]*sec_to_month

#phase(pIOD, nIOD, p+n) concatenation
growth_view1=np.concatenate((growth_view1[:,5].reshape(-1,8),ngrowth_view1[:,5].reshape(-1,8),(growth_view1[:,5]+ngrowth_view1[:,5]).reshape(-1,8)),axis=0)

#Significance test
t95=2.364#95%
t90=1.894#90%
temp_t95=abs(growth_view)*(nevent)**0.50/(np.std(growth_comp,axis=1)+1.0e-20)-t95
ntemp_t95=abs(ngrowth_view)*(nevent)**0.50/(np.std(ngrowth_comp,axis=1)+1.0e-20)-t95
temp_t90=abs(growth_view)*(nevent)**0.50/(np.std(growth_comp,axis=1)+1.0e-20)-t90
ntemp_t90=abs(ngrowth_view)*(nevent)**0.50/(np.std(ngrowth_comp,axis=1)+1.0e-20)-t90
t_diff_95=abs(growth_view+ngrowth_view)*(nevent*nevent*(nevent*2-2)/nevent*2)**0.50/(nevent*np.std(ngrowth_comp,axis=1)**2+nevent*np.std(growth_comp,axis=1)**2)**0.50-t95
t_diff_90=abs(growth_view+ngrowth_view)*(nevent*nevent*(nevent*2-2)/nevent*2)**0.50/(nevent*np.std(ngrowth_comp,axis=1)**2+nevent*np.std(growth_comp,axis=1)**2)**0.50-t90

#plot
label=np.array(["MLT tendency","zonal advection","meridional advection","vertical advection","horizontal diffusion","vertical diffusion","entrainment","surface heat flux"])
cl=np.array(["black","red","green","blue","purple","orange","grey","brown","violet"])
label2 = ["pIOD","nIOD","pIOD+nIOD"]
x = np.arange(len(label2))
width=0.1
fig = plt.figure(figsize=(18.0,7.0))
plt.rcParams["font.size"] = 30
ax = fig.add_subplot(1,1,1)
for i in range(npara):
    ax.bar(x+width*(i-4), growth_view1[:,i], width, label=label[i], color=cl[i],zorder=1)
k=5
for i in range(8):
    if temp_t95[i,k] > 0.0:
        plt.plot(0+width*(i-4+1/2),growth_view1[0,i]-0.03*growth_view[i,k]/abs(growth_view[i,k]),marker="o",color="black",markerfacecolor="white",markersize=10)
    if (temp_t95[i,k] < 0.0)&(temp_t90[i,k] > 0.0):
        plt.plot(0+width*(i-4+1/2),growth_view1[0,i]-0.03*growth_view[i,k]/abs(growth_view[i,k]),marker="^",color="black",markerfacecolor="white",markersize=10)
    if ntemp_t95[i,k] > 0.0:
        plt.plot(1+width*(i-4+1/2),growth_view1[1,i]-0.03*ngrowth_view[i,k]/abs(ngrowth_view[i,k]),marker="o",color="black",markerfacecolor="white",markersize=10)
    if (ntemp_t95[i,k] < 0.0)&(ntemp_t90[i,k] > 0.0):
        plt.plot(1+width*(i-4+1/2),growth_view1[1,i]-0.03*ngrowth_view[i,k]/abs(ngrowth_view[i,k]),marker="^",color="black",markerfacecolor="white",markersize=10)
    if t_diff_95[i,k] > 0.0:
        plt.plot(2+width*(i-4+1/2),growth_view1[2,i]-0.03*(growth_view[i,k]+ngrowth_view[i,k])/abs(growth_view[i,k]+ngrowth_view[i,k]),marker="o",color="black",markerfacecolor="white",markersize=10)
    if (t_diff_95[i,k] < 0.0)&(t_diff_90[i,k] > 0.0):
        plt.plot(2+width*(i-4+1/2),growth_view1[2,i]-0.03*(growth_view[i,k]+ngrowth_view[i,k])/abs(growth_view[i,k]+ngrowth_view[i,k]),marker="^",color="black",markerfacecolor="white",markersize=10)

ax.set_xticks(x)
ax.set_xticklabels(label2,fontsize=40)
plt.ylabel('℃')
plt.ylim(-1,1)
plt.grid(True)
ax.legend(bbox_to_anchor=(1.0, 1.0), loc='upper left', borderaxespad=0,fontsize=30)

plt.savefig(figname, bbox_inches='tight')
print(figname)
