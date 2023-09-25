#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  7 17:05:53 2023

@author: M. El Aabar Ibaoune
"""


import os
import sys

import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.gridspec import GridSpec
#from sklearn.preprocessing import MinMaxScaler, scale
import numpy as np

from matplotlib.colors import LinearSegmentedColormap as LSCmap
from mpl_toolkits.basemap import Basemap

from datetime import datetime
import fcts as fc




def read_fld_of_cwt(pltcfgfile):
	# Input files :

	file_xl = pltcfgfile["read_fld_of_cwt"]["path_to_excel_file"]
	out_npz = pltcfgfile["read_fld_of_cwt"]["out_dir_to_npz_files"]
	if not os.path.exists(out_npz):
		os.makedirs(out_npz)
	years = [*range(int(pltcfgfile["read_fld_of_cwt"]["period"][0]),int(pltcfgfile["read_fld_of_cwt"]["period"][1]),1)]#[*range(1959, 2023, 1)]
	# Basic options 
	seasons = pltcfgfile["read_fld_of_cwt"]["seasons"]#'SON' #'DJF', 'MAM', 'JJA', 'SON', "all_months"
	over_seasons = pltcfgfile["read_fld_of_cwt"]["over_seasons"]#True
	# Param options:

	varibales = pltcfgfile["read_fld_of_cwt"]["varibales"]#['MSLP']
	vars_names_in_netcdf = pltcfgfile["read_fld_of_cwt"]["vars_names_in_netcdf"]#['msl']


	times_WT, WTs = fc.read_xl(file_xl)
	WTs = [int(x+0.5) for x in WTs ] # Probelm linked to Pandas (to investigate after)
	# compute the prct
	pure_wts = [1,2,3,4,5,6,7,8,30,40]
	all_wts  = [1,2,3,4,5,6,7,8,30,32,33,34,35,36,37,38,39,40,42,43,44,45,46,47,48,49]
	dic_wt_names = {1:'NE',2:'E',3:'SE',4:'S',5:'SW',6:'W',7:'NW',8:'N',30:'C',32:'CNE',33:'CE',34:'CSE',35:'CS',36:'CSW',37:'CW',38:'CNW',39:'CN',40:'A',42:'ANE',43:'AE',44:'ASE',45:'AS',46:'ASW',47:'AW',48:'ANW',49:'AN'}
	def get_all_pure_wts(ls_pure,wt_glob):
		ls_out_idx = []
		ls_out = []
		for ielm, elem in enumerate(wt_glob):
			if elem in ls_pure:
				ls_out_idx.append(ielm)
				ls_out.append(elem)
		return(ls_out_idx, ls_out)
	idxs_all_pure_wts, all_pure_wts = get_all_pure_wts(pure_wts,WTs)

	len_tot_pure = len(all_pure_wts)
	len_tot_WT = len(WTs)

	ls_yrs = []
	
	for season in seasons:

		for WT_nbr in all_wts : 
			ls_yrs = []
			freqWT=0
			WTname = dic_wt_names[WT_nbr]
			for ivar, var in  enumerate(varibales) : 
				# process data by years 
				for y in years:
					file_nc = 'data/'+var+'_ByYears/'+var+str(y)+'.nc' #'data/Precip_ByYears/TP1950.nc' 
					print(">> Reading : ", var,"from ",file_nc)
					lons, lats, time, data_var = fc.get_var_from_nc(vars_names_in_netcdf[ivar], file_nc)
					data_var = data_var[:,0,:,:]
			
					#print("   + Select data according to WT = ", WT_nbr)
					dates_wtNbr = fc.select_dates_of_WT(WTs, times_WT, WT_nbr)
					#dates_wt = [tm.split('-') for tm in ymd_wt]
					print("   + Select data according to WT = ", WT_nbr)#, 'over date : ', time)
					data_of_WT, dates_wt_same_order_as_data = fc.select_data_among_dates(data_var, dates_wtNbr, time)
					if len(data_of_WT) == 0:
						print('     /!\ No data found for WT_nbr ',WT_nbr)#,' over the year ',y)
						continue
					if season != 'all_months':
						data_season = fc.select_data_among_season(data_of_WT, dates_wt_same_order_as_data, season)
						if len(data_season) == 0:
							print('     /!\ No data found for WT_nbr ',WT_nbr,' over the season ',season)
							continue
							
						fld = np.array(data_season)
						fld_mean = np.mean(fld,axis=0)
						freqWT = freqWT+len(data_season)#if var == 'MSLP':
							#fld_mean = ((fld_mean*0.222441365418949)+98519.2012793173)/100
						outfile = var+'_WTnbr_'+str(WT_nbr)+'_'+season
						ls_yrs.append(fld_mean)
					else:
						fld = np.array(data_of_WT)
						fld_mean = np.mean(fld,axis=0)
						freqWT = freqWT+len(data_of_WT)#if var == 'MSLP':
							#fld_mean = ((fld_mean*0.222441365418949)+98519.2012793173)/100
						
						outfile = var+'_WTnbr_'+str(WT_nbr)#+'_'+WTname
						ls_yrs.append(fld_mean)#
				
				#ls_y = [x for x in ls_yrs if str(x) != 'nan']
				if len(ls_yrs) == 0:
					print(' ---> No data found for WT_nbr ',WT_nbr,' over the whole Period ! ')
					response = input("       Do you want to continue (y/n) ? ")
					flag = True
					while flag == True:
						if response in ['Yes','yes','y','Y']:
							flag = False
							continue
						elif response in ['No','N','no','n']:
							flag = False
							sys.exit()
						else:
							print(response,' not known ! ',' Please answer by y/n ! ') 
				else:
					data = np.mean(np.array(ls_yrs), axis=0)
					out_dir = os.path.join(out_npz, outfile)
					print('>> saving ',out_dir,'  ' )
					np.savez(out_dir, fld=data, lons=lons, lats=lats, WTnbr=WT_nbr, freqWT=freqWT, freqTot=len_tot_WT, WTname=WTname)
									


def plot_fld_of_cwt(pltcfgfile):
	
	maps = pltcfgfile["plot_fld_of_cwt"]["maps"]#  (4,4) # (5,2) | (4,4)
	nmaps=maps[0]*maps[1]

	period = str(pltcfgfile["plot_fld_of_cwt"]["period"][0])+'_'+str(pltcfgfile["plot_fld_of_cwd"]["period"][0])

	#domaine = '30W20E-30N60N'
		
	seasons = pltcfgfile["plot_fld_of_cwt"]["seasons"]#'JJA' #'DJF', 'MAM', 'JJA', 'SON', "all_months"
	type_of_plot = pltcfgfile["plot_fld_of_cwt"]["type_of_plot"]#'hybrid' # pure, hybrid

	plot_anomal = pltcfgfile["plot_fld_of_cwt"]["plot_anomal"]
	
	name_Figs_outdir = pltcfgfile["plot_fld_of_cwt"]["name_Figs_outdir"]+type_of_plot 
	
	nd_nd = str(maps[0])+'_'+ str(maps[1])

	pa_to_hPa = pltcfgfile["plot_fld_of_cwt"]["pa_to_hPa"] #True
	proj_stere = pltcfgfile["plot_fld_of_cwt"]["proj_stere"] #False
	
	bounds = pltcfgfile['plot_fld_of_cwt']['HorizMeans']['cmap_bounds']
	bounds_diff = pltcfgfile['plot_fld_of_cwt']['HorizMeans']['cmap_bounds_diff']
	cmap = pltcfgfile['plot_fld_of_cwt']['HorizMeans']['cmaps']
	cmap_diff = pltcfgfile['plotWhenItRains']['HorizMeans']['cmaps_diff']

	if proj_stere:
		proj="ster"
	else:
		proj="cycl"

	pure_wts = [1,2,3,4,5,6,7,8,30,40]
	pure_wts_names = ['NE','E','SE','S','SW','W','NW','N','C','A']

	hybrid_wts = [32,33,34,35,36,37,38,39,42,43,44,45,46,47,48,49]
	hybrid_wts_names = ['CNE','CE','CSE','CS','CSW','CW','CNW','CN','ANE','AE','ASE','AS','ASW','AW','ANW','AN']

	all_wts  = [1,2,3,4,5,6,7,8,30,32,33,34,35,36,37,38,39,40,42,43,44,45,46,47,48,49]
	all_wt_names = ['CNE','CE','CSE','CS','CSW','CW','CNW','CN','A','ANE','AE','ASE','AS','ASW','AW','ANW','AN']

	dic_wt_names = {1:'NE',2:'E',3:'SE',4:'S',5:'SW',6:'W',7:'NW',8:'N',30:'C',32:'CNE',33:'CE',34:'CSE',35:'CS',36:'CSW',37:'CW',38:'CNW',39:'CN',40:'A',42:'ANE',43:'AE',44:'ASE',45:'AS',46:'ASW',47:'AW',48:'ANW',49:'AN'}

	fig, axs = plt.subplots(maps[0], maps[1],dpi=400,figsize=(25, 15)) # figsize=(15, 15) | figsize=(25, 15)
	fig.subplots_adjust(hspace = .01, wspace=.02) #hspace = .2, wspace=.002| hspace = .01, wspace=.02
	axs = axs.ravel()

	title= 'map n° '
	title= ''

	if type_of_plot == 'hybrid':
		wts = hybrid_wts
		wt_names = hybrid_wts_names
	else:
		wts = pure_wts
		wt_names = pure_wts_names   


	for season in seasons:
		
		file_avg = 'NPZ_avg/MSLP_AVG_'+season+'.npz' 
		input_dir = 'NPZ_allWTs/'+season # NPZ_allWTs
		
		for i in range(nmaps):
			
			if season == "all_months":
				npzfile_name =  os.path.join(input_dir,"MSLP_WTnbr_"+str(wts[i])+".npz")
			else:
				npzfile_name =  os.path.join(input_dir,"MSLP_WTnbr_"+str(wts[i])+'_'+season+".npz")
			
			continu = False
			try:
				npzfile = np.load(npzfile_name)
				if plot_anomal:
					npzfile_avg = np.load(file_avg)
			except OSError:
				print(npzfile_name, 'or',file_avg,'  Not found ! ')
				response = input("       Do you want to continue (y/n) ? ")
				flag = True
				while flag == True:
					if response in ['Yes','yes','y','Y']:
						flag = False
						continu = True
					elif response in ['No','N','no','n']:
						flag = False
						sys.exit()
					else:
						print(response,' not known ! ',' Please answer by y/n ! ') 
			#print(npzfile.files)
			if continu:
				continue
			print('working on '+npzfile_name+'  ....')
			data = npzfile['fld']
			if plot_anomal:
				data_avg = npzfile_avg['fld']
			#print('data : ',data)
			if pa_to_hPa:
				data = data/100 # Covert from Pa to hPa
				if plot_anomal:
					data_avg = data_avg/100
			if plot_anomal:
				data = data-data_avg
			lat = npzfile['lats']
			lon = npzfile['lons']
			#keys = npzfile['pos']
			freq_node = npzfile['freqWT']
			freq_tot = npzfile['freqTot']
			prct_node = round(freq_node/freq_tot *100,2)
			npzfile.close()
			#lon = np.arange(start=-180, stop=180, step=0.25); lat = np.arange(start=-90.25, stop=90, step=0.25)
			
			mx = np.max(data)
			mn = np.min(data)
			mean = np.mean(data)
			
			ax_node = axs[i]
			#ax_node = axs[keys[i][0],keys[i][1]]
			
			#axs[keys[i][0],keys[i][1] ].set_title(title)#+'\n'+ str(i) +\
												  #' (max,min) = ('\
													 # + str(mx)[:5]+','+str(mn)[:5]+')',\
													  #fontsize=5)
			
			#mp = Basemap(ax=axs[keys[i][0],keys[i][1]])
			
			#mp = Basemap(ax=axs[keys[i][0],keys[i][1]], llcrnrlat=5, urcrnrlat=65, llcrnrlon=-45, urcrnrlon=25,resolution='l')
			mp = Basemap(ax=ax_node, llcrnrlat=np.min(lat),urcrnrlat=np.max(lat),llcrnrlon=np.min(lon),urcrnrlon=np.max(lon),resolution='l')
			
			if proj_stere:
				mp = Basemap(ax=ax_node, width=12000000,height=8000000, llcrnrlat=5, urcrnrlat=65, llcrnrlon=-45, urcrnrlon=25, resolution='l',projection='stere',lat_0=5,lon_0=15)
				#mp.drawcoastlines()
				#mp.drawparallels(np.arange(-80.,81.,20.))
				#mp.drawmeridians(np.arange(-180.,181.,40.))
				# draw parallels and meridians.
				# label parallels on right and top
				# meridians on bottom and left
				
			parallels = np.arange(0.,81,20.)
			# labels = [left,right,top,bottom]
			mp.drawparallels(parallels,labels=[False,True,True,False],fontsize=6)
			meridians = np.arange(10.,351.,40.)
			mp.drawmeridians(meridians,labels=[True,False,False,True],fontsize=6)
			
			lon, lat = np.meshgrid(lon, lat)
			xx, yy = mp(lon, lat)
			# Draw coastlines 
			mp.drawcoastlines(linewidth=0.5)
			#mp.drawcountries(linewidth=0.25)
			avail_cmap = [cm  for cm in plt.colormaps() if cm not in dir(plt.cm)]
			#sc = mp.pcolor(xx, yy, data, cmap = 'jet')
			cmapMed_centred_on_white_for_diff = LSCmap.from_list('cmapMed_centred_on_white_for_diff', 
												[(0 , 'blue'),
												(0.2, 'green'),
												(0.3, 'lime'),
												(0.48, 'white'),
												(0.5, 'white'),
												(0.52, 'white'),
												(0.7, 'red'),
												(0.8, 'brown'),
												(1, 'black')])
				  
			if not 'cmapMed_centred_on_white_for_diff' in  avail_cmap : plt.cm.register_cmap(name='cmapMed_centred_on_white_for_diff', cmap=cmapMed_centred_on_white_for_diff)
	
			cmap=plt.cm.get_cmap('jet', 16)#'cmapMed_centred_on_white_for_diff'#plt.cm.get_cmap('jet', 16)
			cmap = cmap
			if plot_anomal:
				cmap = cmap_diff
				bounds = bounds_diff
			levels=np.arange(bounds[0], bounds[1], bounds[2])#np.arange(-12,13,1) # (980,1042,2)
			sc = mp.contour(xx, yy, data, colors='k', linewidths=0.2,linestyles='dashed')
			ax_node.clabel(sc,  inline=True, fmt='%1.0f', fontsize=8, colors='k')
			sc = mp.contour(xx, yy, data, levels=sc.levels, linewidths=0.2, colors='k',linestyles='dashed')
			#plt.colorbar(sc,cax=axs[keys[i][0],keys[i][1]])
			#sc.levels = np.arange(start=np.min(sc.levels), stop=np.max(sc.levels)-2, step=2)
			sc = mp.contourf(xx, yy, data, cmap = cmap,levels=levels)
			#sc = mp.contourf(xx, yy, data, cmap =cmap,levels=15)
		#	sc = mp.contourf(xx, yy, data, cmap =cmap,levels=18)
			
			
			ax_node.set_title(wt_names[i] + '  '+ str(prct_node)+' % ',fontsize=12,fontweight="bold")
			
			#title_ = title + '  max = '+ str(mx)[:5]+'  min = ' +str(mn)[:5]
			
				#map.contour(xx,yy,data)
			#plt.title(title_)
		anomal_title = ''
		anomal_fig_name = ''
		if plot_anomal:
			anomal_title = 'Anomal (fld minus avg)'
			anomal_fig_name = 'ANOMAL'
			name_Figs_outdir = pltcfgfile["plot_fld_of_cwd"]["name_Figs_outdir"]+'_Anomal_'+type_of_plot 
	
				
		fig.suptitle(anomal_title+type_of_plot+' weather types , season : '+ season+', period : '+period , fontsize=18,fontweight="bold")  
		fig.colorbar(sc, ax=axs.ravel().tolist(),pad=0.04)
		outdir = name_Figs_outdir+'/'+season
		if not os.path.exists(outdir):
			os.makedirs(outdir)
		namefig = 'cwt_'+ nd_nd +'_proj_'+proj+'_season_'+season+'_'+anomal_fig_name+'.png'
		outFile = os.path.join(outdir,namefig)
		print('Saving ',outFile, ' ' )
		plt.savefig(outFile, bbox_inches='tight')
		plt.show()
