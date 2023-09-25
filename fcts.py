#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  7 17:05:53 2023

@author: M. El Aabar Ibaoune
"""


import os
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.gridspec import GridSpec
#from sklearn.preprocessing import MinMaxScaler, scale
import numpy as np
import  netCDF4 as nc4
from mpl_toolkits.basemap import Basemap, shiftgrid
import matplotlib.cm as cm
from matplotlib.colors import LinearSegmentedColormap as LSCmap
from datetime import datetime
#import Nio
import pandas as pd
import xarray as xr





def read_xl(path_to_dates_file) : 
	tab = pd.read_excel(path_to_dates_file).to_numpy()
	ymj = tab[:,0:3]
	ymj = (np.asarray(ymj, dtype = 'int')).astype(str)
	y = ymj[:,0]; m = ymj[:,1]; j = ymj[:,2]
	ymj = [datetime(int(a), int(b), int(c)) for a, b, c in zip(y, m, j)] #[a + '-' + b + '-'+ c for a, b, c in zip(y, m, j)]
	WT_idx = tab[:,3]
	return(ymj, WT_idx)

def select_data_among_season(field, dates, season):
	data_of_season = []
	dic_seas = {"DJF":[12,1,2],"MAM":[3,4,5],"JJA":[6,7,8],"SON":[9,10,11]}
	for idate,date in enumerate(dates):
		if int(date.month) in dic_seas[season]:
			#print('> ',season,' eqv:   mnth ',int(date[1]))
			data_of_season.append(np.squeeze(field[idate]))
	return(data_of_season)

def get_var_from_nc(var, filename, read_by_chunks=False):

	ncfile = nc4.Dataset(filename)
	#print(ncfile.variables.keys())
	lons = ncfile['longitude'][:].data # Extraire les longitudes 
	lats = ncfile['latitude'][:].data # Extraire les latitudes 
	temps = ncfile['time'][:].data # Extraire le temps 
	dtime = nc4.num2date(ncfile.variables['time'],ncfile.variables['time'].units)
	#tt = datetime.strptime(str(dtime[0]),'%Y-%m-%d %H:%M:%S')
	#print(tt.year, tt.month, tt.day, tt.hour, tt.second)
	temps_splitted = []
	[temps_splitted.append(datetime.strptime(str(dtime[elem]),'%Y-%m-%d %H:%M:%S')) for elem in range(len(temps))] 
	
	if read_by_chunks:
		# To chunks: proceed 2000 by 2000
		var_all_chunks  = []
		dim_time = len(temps)
		length_check = 4000
		number_chunks  = dim_time//length_check
		#rest = dim_time%length_check
		print('   +> extracting data by chnks for', number_chunks, ' of length ', length_check)
		for k in range(number_chunks):
			bnd1 = k*length_check
			bnd2 = (k+1)*length_check
			print('   +> extracting data between : ', bnd1, ' and ', bnd2) 

	#        name_var = (ncfile["t"][bnd1:bnd2, :, :].data)-273.15
			var_chunk = ncfile[var][bnd1:bnd2, :, :].data/100
			var_all_chunks.extend(var_chunk)
		data_var = var_all_chunks
			
	else:   
		print (' ')
		data_var = ncfile[var][:]
	
	return(lons, lats, temps_splitted, data_var)


def select_data_among_dates(field, dates, temps):
	
	data_of_precipitant_days = []
	date_of_precipitant_days = []
	for date in dates:
		for it in range (len(field)):
			if temps[it].year == date.year : 
				if temps[it].month == date.month :
					if temps[it].day == date.day :
						#print('date : ', date, 'found !')
						data_of_precipitant_days.append(np.squeeze(field[it]))
						date_of_precipitant_days.append(date)
						break
	return(data_of_precipitant_days, date_of_precipitant_days)


def get_coords(file_name):
	print('Getting coords from ', file_name,'...')
	xr_obj = xr.open_dataset(file_name)
	lons = xr_obj.longitude.values
	lats = xr_obj.latitude.values
	time = xr_obj.time.values
	return(lons, lats, time)

def comptime_to_str(date):
	return(str(date.year) +'-'+str(date.month).zfill(2)+'-'+str(date.day).zfill(2) )

def select_dates_of_WT(WTs, dates, WT_nbr):
	print('> Select dates of WT number :',WT_nbr)
	ls_out = []
	for iwt, wt in enumerate(WTs):
		if wt==WT_nbr:
			ls_out.append(dates[iwt])
	return(ls_out)

def get_idx_of_equal_time(dates1, dates2):
	print('> Select communs dates between two lists ...')
	ls = []
	for i in range(len(dates1)):
		if dates1[i] in dates2:
			ls.append(i)
	return(ls)

def get_data_for_idxs(var, file_name, idxes):
	print('Getting ', var, 'from ', file_name,'...')
	ncfile = Nio.open_file(file_name, format="netcdf")
	fd = ncfile.variables[var][-1][1][:]
	for idx in idxes[:-1]:
		#print('idx : ',idx)
		fdi = ncfile.variables[var][idx][0][:][:]
		#print(fdi.shape)
		fd  = np.dstack([fd, fdi])
	ncfile.close()
	return(fd)
	
def get_data_xr(file_name,idxes):
	print('Getting coords from ', file_name,'...')
	xr_obj = xr.open_dataset(file_name)
	fd  = xr_obj.msl.values[0]
	return(fd)

def get_data_nc4(var, filename):

	ncfile = nc4.Dataset(filename)
	#print(ncfile.variables.keys())
	fd = ncfile[var][0] # Extraire les longitudes
	for idx in idxes[:-1]:
		#print('idx : ',idx)
		fdi = ncfile.variables[var][idx][1][:][:]
		#print(fdi.shape)
		fd  = np.dstack([fd, fdi])
	ncfile.close()
	return(fd)
	return( data_var)
	
	

def Basemap_one_fig(lon,lat,data,**kwargs):
	# Create a map instance using basemap classe. 
	# Remarque: Basemap class constructor have many arguments and all are optional,
	# Dafault values are taken if none of the argumnets was not choosen.
	#map = Basemap(projection='ortho', llcrnrx=-3000000, llcrnry=1000000, urcrnrx=3000000, urcrnry=6000000, \
		#lat_0=10, lon_0=0, resolution='c') 
	fig, ax = plt.subplots(dpi=400,figsize=(15, 5))
	#fig, axes = plt.subplots(figsize=(12,12),dpi=300,nrows=nber_of_plos_horz, ncols=nber_of_plos_vrtcl)
	
	mp = Basemap(ax=ax,projection='merc',llcrnrlat=20,urcrnrlat=40,\
	llcrnrlon=-20,urcrnrlon=10,resolution='c')
		
	# Set axis 
	#lon, lat = np.meshgrid(lon, lat)
	# compute native map projection coordinates of lat/lon grid.
	#x, y = mp(lon, lat)

	# Draw coastlines 
	mp.drawcoastlines(linewidth=0.2)
	#mp.drawcountries(linewidth=0.2)

	# draw lat/lon grid lines every 30 degrees.
	#mp.drawmeridians(lat[0:-1:20])#p.arange(-180,180,20))
	#mp.drawparallels(lon[0:-1:20])#np.arange(-90,90,20))
	
	# Data to plot :
	#sc = mp.pcolor(x, y, data, cmap = 'jet', vmin= 270, vmax=300)
	
	colorbar_bounds = (np.min(data),np.max(data))
	if 'colorbar_bounds' in kwargs:
		colorbar_bounds = kwargs['colorbar_bounds']
	
	cmap = 'RdBu_r'
	if 'cmap' in kwargs:
		cmap = kwargs['cmap']
	
	if 'PlotUnderMinWiteColor' in kwargs:
		if kwargs['PlotUnderMinWiteColor'] == True:
			cmap.set_under('w')
	
	marker = 'o'
	if 'marker' in kwargs:
		marker = kwargs['marker']
	size = 30
	if 'size' in kwargs:
		size = kwargs['size']
		
	scatter = False
	if 'scatter' in kwargs:
		size = kwargs['scatter']
			
	
		
	parallels = [ int(x) for x in  lat[0:-1:15]] #np.arange(0.,81,10.)
	# labels = [left,right,top,bottom]
	mp.drawparallels(parallels,labels=[True,False,True,True],fontsize=6)
	meridians = [ int(x) for x in lon[0:-1:40]] #np.arange(10.,351.,20.)
	mp.drawmeridians(meridians,labels=[True,True,False,True],fontsize=6)
		
		
		
	lons, lats = np.meshgrid(lon, lat)
	xx, yy = mp(lons, lats)
	# Draw coastlines 
	mp.drawcoastlines(linewidth=0.1)
	#mp.drawcountries(linewidth=0.25)
	# Data to plot :
	#sc = mp.pcolor(xx, yy, data, cmap = 'jet')
	
		
	#levels=np.arange(-0.1, 0.1, .005) #16; precip(0.03, 0.5, .05) 
	#levels=[-5,-4,-2,-1,-0.5,0.5,1,2,3, 4, 5]
	#levels=np.arange(np.min(data), np.max(data),1)
	#cmap='hsv'#plt.cm.get_cmap('jet', 12) #seismic
	if 'levels' in kwargs.keys():
		levels = kwargs['levels']
	if 'cmap' in kwargs.keys():
		cmap = kwargs['cmap']
		if cmap=='cmapMed':
			cmap = plt.cm.get_cmap('cmapMed')
		if cmap=='cmapMed_centred_on_white':
			cmap = plt.cm.get_cmap('cmapMed_centred_on_white')

	sc = mp.contourf(xx, yy, data)#,levels= levels)# cmap =cmap)#,levels= levels,extend="both")
	#sc.cmap.set_under('k')
	#sc.set_clim(190, 315)
		
	fig.colorbar(sc,orientation='horizontal',fraction=0.06,pad=0.04,ax=ax)
		
	#if 'precip' in kwargs['namefig']: Colorabar nonlinear 
		#fig.colorbar(sc,ticks=[0.5,1,2,3,4,6,8,10,12,15,20], orientation='horizontal',pad=0.04,ax=ax)
		
	sc.levels = np.arange(start=np.min(sc.levels), stop=np.max(sc.levels)-2, step=4)
	#sc = mp.contour(xx, yy, data, levels=sc.levels, linewidths=0.3, colors='k')
		
	#sc.levels = np.arange(start=np.min(sc.levels), stop=np.max(sc.levels)-2, step=2)
	
	ax.clabel(sc,  inline=True, fmt='%.2f', fontsize=4, colors='k')#,levels=sc.levels,)
		
		
	#plt.colorbar(sc,cax=ax)
		
	#sc = mp.contourf(xx, yy, data, cmap = 'jet')#,levels=8)
	
	# Colorbar 
	
		
	#norm = mpl.colors.Normalize(vmin=0,vmax=2)
	#sm = plt.cm.ScalarMappable(cmap="jet")
	#sm.set_array([])
	
	#plt.colorbar(sm)
	#mp.colorbar(cmap=cm.jet)
	
	title = ''
	if 'title' in kwargs.keys():
		title = kwargs['title']
	ax.set_title(title)
	
	namefig = 'Fig.png'
	if 'namefig' in kwargs.keys():
		namefig = kwargs['namefig']
		if namefig[-3:] != 'png':
			namefig=namefig+'.png'
	
	outdir = os.getcwd()
	if 'outdir' in kwargs.keys():
		outdir = kwargs['outdir']
			
	to_save = os.path.join(outdir,namefig)
	print ('  Saving ',to_save)
	plt.savefig(to_save)
	#plt.show()
