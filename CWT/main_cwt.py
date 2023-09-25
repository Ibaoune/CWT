#!/usr/bin/env python3
# -*- coding: utf-8 -*-



__doc__ = """
Created on Tue Aug 15 12:57:46 2023
"""
__author__ = "M. El Aabbaribaoune (@_um6p)"

import os
import sys
import json
import core_cwt

def main(pltcfgfile):
    
    # Load the configuration
    with open(pltcfgfile) as file:
        plotcfg = json.load(file)
        
    inputdir = plotcfg['inputdir']
    
    if not os.path.exists(inputdir):
        print(inputdir," does not exists. Exit.")
        sys.exit()
        
    os.chdir(inputdir)
    print("Working inside ",inputdir)
    
            
    if "read_fld_of_cwt" in plotcfg:    
        if plotcfg["read_flds_and_save_to_npz"]:     
            core_cwt.read_fld_of_cwt(plotcfg)
        else:
            print("plotcfg[read_fld_and_save_to_npz] is : ",plotcfg["read_fld_and_save_to_npz"])
    else:
        print('read_fld_of_cwt NOT in plotcfg')
            
    if "plot_fld_of_cwt" in plotcfg:
        
        if plotcfg["plot_flds_of_cwt"]: 
            core_cwt.plot_fld_of_cwt(plotcfg)
        else:
            print("plotcfg[plot_flds_of_cwt] is : ",plotcfg["plot_flds_of_cwt"])
    else:
        print('plot_flds_of_cwt NOT in plotcfg')
        
           
            
if __name__ == "__main__":
    
    pltcfgfile = 'cwt_namelist.json'
    # Read cfg file from user input if given. It overrides default
    #if len(sys.argv) > 1: pltcfgfile = sys.argv[1]
    main(pltcfgfile)  
    
    
