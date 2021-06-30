# -*- coding: utf-8 -*-
"""
Initial Software, Myriam Guillevic and Aurore Guillevic,
Copyright Empa and Inria, 2019 - 2021.

This file is part of ALPINAC.

ALPINAC is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

ALPINAC is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with ALPINAC.  If not, see <https://www.gnu.org/licenses/>.

"""
#!/usr/bin/env python
# Implementation of algorithm from http://stackoverflow.com/a/22640362/6029703
import numpy as np

#import statistics as stat
from scipy.stats import t as Student_t
#import pylab
import math
#import h5py
#https://www.h5py.org/
from scipy.special import gamma  #for Voigt function # wofz, factorial
#from scipy.special import 
import operator #for make_buckets: to sort list
from lmfit import Minimizer
from lmfit import *

from periodic_table import formula_to_mass

import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D #to uncomment

#https://matplotlib.org/
#make graph with two y-axis:
#https://matplotlib.org/examples/api/two_scales.html

from runtype import RunType

#******************
#define colors
#https://xkcd.com/color/rgb/
color1_names = ['xkcd:blue', 'xkcd:green', 'xkcd:orange', 'xkcd:pink', 'xkcd:rose', 
                'xkcd:kelly green', 'xkcd:light orange', 'xkcd:goldenrod', 'xkcd:forest green', 'xkcd:burnt orange', 
                'xkcd:apple green', 'xkcd:reddish purple', 'xkcd:terra cotta', 'xkcd:hot pink', 'xkcd:peach', 
                'xkcd:mauve', 'xkcd:dark red', 'xkcd:taupe', 'xkcd:cerulean', 'xkcd:eggplant',
                'xkcd:blue', 'xkcd:green', 'xkcd:orange', 'xkcd:pink', 'xkcd:rose', 
                'xkcd:kelly green', 'xkcd:light orange', 'xkcd:goldenrod', 'xkcd:forest green', 'xkcd:burnt orange', 
                'xkcd:apple green', 'xkcd:reddish purple', 'xkcd:terra cotta', 'xkcd:hot pink', 'xkcd:peach', 
                'xkcd:mauve', 'xkcd:dark red', 'xkcd:taupe', 'xkcd:cerulean', 'xkcd:eggplant',
                'xkcd:blue', 'xkcd:green', 'xkcd:orange', 'xkcd:pink', 'xkcd:rose', 
                'xkcd:kelly green', 'xkcd:light orange', 'xkcd:goldenrod', 'xkcd:forest green', 'xkcd:burnt orange', 
                'xkcd:apple green', 'xkcd:reddish purple', 'xkcd:terra cotta', 'xkcd:hot pink', 'xkcd:peach', 
                'xkcd:mauve', 'xkcd:dark red', 'xkcd:taupe', 'xkcd:cerulean', 'xkcd:eggplant',
                'xkcd:blue', 'xkcd:green', 'xkcd:orange', 'xkcd:pink', 'xkcd:rose', 
                'xkcd:kelly green', 'xkcd:light orange', 'xkcd:goldenrod', 'xkcd:forest green', 'xkcd:burnt orange', 
                'xkcd:apple green', 'xkcd:reddish purple', 'xkcd:terra cotta', 'xkcd:hot pink', 'xkcd:peach', 
                'xkcd:mauve', 'xkcd:dark red', 'xkcd:taupe', 'xkcd:cerulean', 'xkcd:eggplant']
color2_names = ['xkcd:sky blue', 'xkcd:teal', 'xkcd:peach', 'xkcd:mauve', 'xkcd:dark red', 'xkcd:taupe', 'xkcd:cerulean', 'xkcd:eggplant']



def find_time(rt_suspect, rt_suspect_U, s_per_bin, time_axis):
    rt_time_idx_start=int(round((rt_suspect-rt_suspect_U)/s_per_bin))
    rt_time_idx_stop =int(round((rt_suspect+rt_suspect_U)/s_per_bin))
    rt_time_idx =[idx for idx in range(rt_time_idx_start, rt_time_idx_stop)]
    time_timeseries = [time_axis[idx] for idx in range(rt_time_idx_start,rt_time_idx_stop)]
    return rt_time_idx, time_timeseries

    
def find_closest_centre(x_exact, array_x):
    """
    Array_x has len>1
    Each array_x[i] is a list of x_values belonging to the same centre of x (mass of TOF-index)
    We want to find which array_x[i] has the centre closest to x_exact
    and return the index of this list only
    """
    idx_best =0
    min_diff = abs(np.median([array_x[0][i][1] for i in range(len(array_x[0]))]) - x_exact)
    diff = abs(np.median([array_x[1][i][1] for i in range(len(array_x[1]))]) - x_exact)
    #print('Mass diff 0: ' + str(min_diff) + '; mass diff 1: ' + str(diff))
    
    
    #as buckets are ordered by increasing masses, 
    #median of masses of bucket should also be ordered by increasing masses
    idx_array_x = 1
    #do len test first!
    while idx_array_x < len(array_x) and diff < min_diff:
        idx_best = idx_array_x #+ 1
        min_diff = diff
        idx_array_x += 1
        #print(idx_array_m)
        if idx_array_x < len(array_x): # need to add this otherwise 'list index out of range' for last line
            diff = abs(np.median([array_x[idx_array_x][i][1] for i in range(len(array_x[idx_array_x]))]) - x_exact)
            #print('Mass diff ' +str(idx_array_x) +  ': ' + str(diff))
 
    return idx_best




def peak_detect(x, y, PW, PT, n_per_bin, graph, mode):
    """
    This function detects zero to multiple peaks in x, y values. x can be time or mass.
    A peak is detected based on slope calculation and detection parameters:
    minimum slope value, slope going from positive to negative value.
    The function returns:
    peak_position_final: list of rough estimate(s) of peak centre position(s) 
    peak_value_final: list of rough estimate(s) of intensity at peak centre position(s)
    x: np.array of floats (not integer! even if they are actually indices values)
    y: np.array of floats
    PW: expected peak width, unit same as x
    PT: peak threshold, minimum abs(slope) to consider a peak is present
    n_per_bin: number of point per unit of x. n_per_bin * PW gives number of points making a peak.
    For the RT domain, n_per_bin is s_per_bin (usually 0.166, 6 pts per seconds)
    For the time of flight domain, n_per_bin is one.
    
    return:
    peak_position_final, np.array of floats
    peak_value_final, np.array of floats
    """
    
    #print('Detection thereshold: ' + str(PT))
    #n has to be an odd number (1,3,5...)
    n = int(PW/n_per_bin/2)
    if (n % 2) == 0: #n is an even number
        n+=1
    
    if n < 5: n = 5
    #if n > 17: n = 17
    
    if len(y)-n >4: #20191002

        slope_list = np.zeros((len(y)-n,4))
        
        #20190905 the x are always evenly spaced so
        #den can be calculated just once
        den = np.sum(x[0:n]**2) - (np.sum(x[0:n]))**2/n
        
        for i in range(len(y)-n): #window for slope calculation, i is starting point
            x_sum = np.sum(x[i:i+n]) #last stated index is included
            #print(x[i:i+n])
            #print(x_sum)
            #print(den)
            if den != 0:
                y_sum = np.sum(y[i:i+n])
                num = np.sum(x[i:i+n] * y[i:i+n]) - x_sum*y_sum/n            
                slope = num/den
            else : slope = 0
            slope_list[i:] = [x[i+int(n/2)], y[i+int(n/2)], slope, 0.0]
            #20190906 mygu the above looks correct but I still don't understand why
            #slope_list[i:] = [x[i+int(n/2)+1], y[i+int(n/2)+1], slope, 0.0]
        #print(slope_list)
        
        for idx in range(len(slope_list)-1):      
            #Detection of peaks:
            if slope_list[idx,2] > 0 and slope_list[idx+1,2] <= 0: #and slope_value3 < 0
                #this is a peak
                slope_list[idx, 3] = 1.0
            #detection of valleys:
            elif slope_list[idx,2] < 0 and slope_list[idx+1,2] >= 0: # and slope_value3 > 0
                #this is a valley
                slope_list[idx, 3] = -1.0
        #print(slope_list)
        #print(len(slope_list))
        #201909096 mygu I think it does not work the following way:
        #slope_list[:, 3] = np.where(slope_list[:,2] > 0 and slope_list[idx+1,2] <= 0, 1.0, 0.0)
        #slope_list[:, 3] = np.where(slope_list[idx,2] < 0 and slope_list[idx+1,2] >= 0, -1.0, 0.0)
        
        #find indexes of slope list where a peak has been detected:
        peak_idx = np.argwhere(slope_list[:,3]>0)
        #same for the valleys:
        valley_idx = np.argwhere(slope_list[:,3]<0)
        #print(peak_idx)
        #print(valley_idx)
        
        idx_p = 0
        while len(peak_idx)>0 and idx_p < len(peak_idx):
            idx_peak = int(peak_idx[idx_p])
    
            if idx_p ==0:
                idx_start = 0
                if len(peak_idx)==1:
                    idx_stop  = len(slope_list)-1
                elif len(peak_idx)>1:
                    idx_stop  = int(peak_idx[idx_p+1]+1)
            elif idx_p > 0: #and idx_p == len(peak_idx)-1:
                idx_start = min(int(peak_idx[idx_p-1]+1), len(slope_list)-2)
                if idx_p == len(peak_idx)-1:
                    idx_stop  = len(slope_list)-1
                else:
                    #idx_start = int(peak_idx[idx_p-1]+1)
                    idx_stop  = int(peak_idx[idx_p+1]+1) #the last index is out of the loop (?)
            
            #print(idx_start)
            #print(idx_stop)
            
            slope_max = max(slope_list[idx_start:idx_peak+1, 2])
            slope_min = min(slope_list[idx_peak:idx_stop, 2])
            slope_extrem = max(slope_max, abs(slope_min))
            
            if slope_extrem >= PT:
                #calculate more precise x position
                #assign y max value
                if idx_peak+1 < len(slope_list):
                    x_value_zero_slope = slope_list[idx_peak,0] - slope_list[idx_peak,2] \
                    * (slope_list[idx_peak+1,0] - slope_list[idx_peak,0]) \
                    / (slope_list[idx_peak+1,2] - slope_list[idx_peak,2])
                    slope_list[idx_peak,0] = x_value_zero_slope
                    #look for max value of y around peak position:
                slope_list[idx_peak,1] = max(slope_list[max(0, idx_peak-1):min(idx_peak+1, len(slope_list)-1),1])
                #if the peak is on the last index of the list, do nothing    
                idx_p += 1  
                
            else:
                slope_list[peak_idx[idx_p],3] = 0
                peak_idx = np.argwhere(slope_list[:,3]>0) #le should be -1
                #print(peak_idx)
    
        #check if data starts by tail of a peak:
        if len(valley_idx) > 0 and len(peak_idx) > 0 and valley_idx[0] < peak_idx [0]:
            #data start by tail of a peak peaking before the considered window
            #we arbitrarily assign a peak to the first point in the window
            slope_min = min(slope_list[0:int(valley_idx[0])+1, 2])
            if abs(slope_min) > PT:
                slope_list[0,3] = 1
                
        #take all valid peaks:
        #peak_position_final = slope_list[slope_list[:,3]>0, 0]
        #peak_value_final = slope_list[slope_list[:,3]>0, 1]
        
        peak_detected = slope_list[slope_list[:,3]>0, 0:2] #last index not included it seems
        
        if graph == True:
            
            #PLOT DATA AND SLOPE
            fig, ax1 = plt.subplots(figsize=(8,4))
            fig.suptitle('Peak detection', fontsize=12)
            ax1.plot(x, y, '--bo') #linestyle='--', marker='o', color='b', , markersize = 10
            ax1.tick_params('y', colors='b')
            #ax1.set_xlim(1300, 1800)
            #ax1.set_ylim(0, 3000)

            #plot slope on second y-axis: 
            ax2 = ax1.twinx()
            ax2.plot(slope_list[:,0], slope_list[:,2], 'g-')
            ax2.plot(slope_list[:,0], np.zeros(len(slope_list)), 'k--')
            ax2.set_ylabel('slope', color='g')
            #ax2.set_ylim(-1000, 1000)
            #ax2.set_ylim(-5, 5)
            ax2.tick_params('y', colors='g')
            
            #legend of axis:
            if mode == 'tofidx_peak_detect':
                ax1.set_xlabel('Time of flight index [no unit]')
                ax1.set_ylabel('Intensity [V]', color='b')
            elif mode == 'rt_peak_detect':
                ax1.set_xlabel('Chromatographic retention time [s]')
                ax1.set_ylabel('Intensity [V * ToF index]', color='b')
            elif mode == 'histo_tofidx': 
                ax1.set_xlabel('Time of flight index [no unit]')
                ax1.set_ylabel('Number of peak centres', color='b')
            elif mode == 'histo_rt': 
                ax1.set_xlabel('Chromatographic retention time [s]')
                ax1.set_ylabel('Number of peak centres', color='b')
            else:
                ax1.set_xlabel('x axis')
                ax1.set_ylabel('y axis', color='b')
                
                
            
        
        if len(peak_detected) > 0 and graph == True:  
            ax1.plot(peak_detected[:,0], peak_detected[:,1], 'o', color = 'xkcd:peach', markersize = 12)
        
        if graph == True:
            #need to add same condition as when creating the figure
            #plt.savefig('peak_detect_'+ str(int(round(peak_position_final[0])))+'.png', bbox_inches='tight', transparent = True, dpi = 300)
            
            #plt.savefig('peak_detect_test293.png', bbox_inches='tight', transparent = True, dpi = 300)
            fig.tight_layout()
            plt.show()
    else: peak_detected = []
    
    return peak_detected

def make_mass_cal_dict_of_frag(dict_spikes, spike_duration, dict_k, max_filament_value_pos, cl_run):
    
    """
    

    Parameters
    ----------
    dict_spikes : TYPE dictionnary
        DESCRIPTION. Each key is the centre rt of a spike.
        Each key contains a list of chemial formula to use for that specific spike.

    Returns
    -------
    dict_frag: TYPE dictionnary
    Each key is a fragment formula. Each key contain a list 
    with rt centre where to search for this mass.

    """
    track_masses_mode_continuous = 0
    track_masses_mode_spike = 1
    mode_continuous_status = False
    dict_frag = {}
    dict_frag_mode = {}
    
    #[frag] -> rt_idx_start, rt_idx_stop, rt_idx_centre
    masscal_list_exact = []
    rt_idx_centre_list = []
        
    i_key = 0
    for key_spike in dict_spikes:
        
        if key_spike == '0.0':
            #these are the fragments to be treated as continuous extraction.
            if mode_continuous_status == False:
                mode_continuous_status = True
            
            #
            nb_pt_per_time_avrg = int(round(dict_k['mc_avrg_min']*60.0/cl_run.s_per_bin))
            nb_pt_per_time_spacing = int(round(dict_k['mc_spacing_min']*60.0/cl_run.s_per_bin))
            
            rt_start = 0.0 #seconds, start of mass calibration
            rt_idx_start = int(rt_start/cl_run.s_per_bin)
            
            #position of last index with minimum xx min average
            #and filament on
            idx_last = min(max_filament_value_pos*cl_run.hdf5_dim2, cl_run.len_time_axis - np.remainder(cl_run.len_time_axis, nb_pt_per_time_avrg) - nb_pt_per_time_avrg)
            rt_idx_list = np.linspace(rt_idx_start, idx_last, int(round((idx_last-rt_idx_start)/nb_pt_per_time_spacing)) +1, dtype = int)

            if len(rt_idx_centre_list) == 0:
                #save centre of each time slice:
                rt_idx_centre_list = rt_idx_list + int(round(nb_pt_per_time_avrg/2))
            

            for candidate_frag in dict_spikes[key_spike]:
                if candidate_frag not in dict_frag:
                    dict_frag[candidate_frag] = []
                    dict_frag_mode[candidate_frag] = track_masses_mode_continuous
                for idx_rt in range(len(rt_idx_list)):
                    if idx_rt == len(rt_idx_list)-1:
                        rt_idx_stop = cl_run.len_time_axis
                    else:
                        rt_idx_stop = rt_idx_list[idx_rt] + nb_pt_per_time_avrg
                    dict_frag[candidate_frag].append([rt_idx_list[idx_rt], 
                                                      rt_idx_stop, 
                                                      rt_idx_list[idx_rt] + int(round(nb_pt_per_time_avrg/2))])                 
                    


        else:
            for candidate_frag in dict_spikes[key_spike]:
                if candidate_frag not in dict_frag:
                    dict_frag[candidate_frag] = []
                    dict_frag_mode[candidate_frag] = track_masses_mode_spike
                
                dict_frag[candidate_frag].append([int(round((float(key_spike)-spike_duration/2)/cl_run.s_per_bin)), 
                                                  int(round((float(key_spike)+spike_duration/2)/cl_run.s_per_bin)), 
                                                  int(round((float(key_spike)/cl_run.s_per_bin)))])
                
        i_key += 1
    masscal_list_exact = [formula_to_mass(key) for key in dict_frag]
    #print([key for key in dict_frag])
    """
    dict_centre_list = {}
    for key in dict_frag:
        for time_slice in key:
            #need to adjust value of centre, should be one
            if time_slice[2] not in dict_centre_list:
                dict_centre_list[time_slice[2]] = 1
            else:
                dict_centre_list[time_slice[2]] += 1
    """
    if mode_continuous_status == True:
        #there is at least one continuous mass.
        #Assign a rt centre from it to the spikes
        #Do not do this if there are only spikes
        for key in dict_frag:
            #print(key)
            if dict_frag_mode[key] == track_masses_mode_spike:
                
                for i_slice in range(len(dict_frag[key])):
                    rt_centre = dict_frag[key][i_slice][2]
                    if rt_centre not in rt_idx_centre_list:
                        min_interval = min([abs(rt - rt_centre) for rt in rt_idx_centre_list])
                        idx_min = [i for i in range(len(rt_idx_centre_list)) if abs(rt_idx_centre_list[i] - rt_centre) == min_interval]
                        dict_frag[key][i_slice][2] = rt_idx_centre_list[idx_min][0]
    else:
       
        rt_idx_centre_list = [int(round((float(key_spike)/cl_run.s_per_bin))) for key_spike in dict_spikes]
        
    
 
    return dict_frag, dict_frag_mode, masscal_list_exact, rt_idx_centre_list
    

def calibrate_mass(rt_tofidx_peak, tofidx_peak, tofidx_peak_u, mass_cal_parameters):
    mass_calib_params = {}
    mass_calib_params['mode'] = 2
    
    #find preceeding and following rt from mass_cal_parameters:
    #do linear interpolation respective to rt of mass_cal_parameters
    if rt_tofidx_peak <= min([line[0] for line in mass_cal_parameters]):
        idx_mass_cal = 0
        mass_calib_params['p1'] = mass_cal_parameters[idx_mass_cal][1]
        mass_calib_params['p2'] = mass_cal_parameters[idx_mass_cal][2]
        mass_calib_params['p3'] = mass_cal_parameters[idx_mass_cal][3]
    elif rt_tofidx_peak >= max([line[0] for line in mass_cal_parameters]):
        idx_mass_cal = len(mass_cal_parameters)-1
        mass_calib_params['p1'] = mass_cal_parameters[idx_mass_cal][1]
        mass_calib_params['p2'] = mass_cal_parameters[idx_mass_cal][2]
        mass_calib_params['p3'] = mass_cal_parameters[idx_mass_cal][3]
    else: 
        i = 0
        #line = tofidx_spectrum_list[i]
        while i < len(mass_cal_parameters) and  mass_cal_parameters[i][0] < rt_tofidx_peak:
            idx_mass_cal = i #take the index of the cal time right before the measured rt position
            i+=1
            
        linear_interpol_weight = (rt_tofidx_peak - mass_cal_parameters[idx_mass_cal][0])  / (mass_cal_parameters[idx_mass_cal+1][0] - mass_cal_parameters[idx_mass_cal][0])
        mass_calib_params['p1'] = mass_cal_parameters[idx_mass_cal][1] * (1-linear_interpol_weight) + mass_cal_parameters[idx_mass_cal+1][1] * linear_interpol_weight
        mass_calib_params['p2'] = mass_cal_parameters[idx_mass_cal][2] * (1-linear_interpol_weight) + mass_cal_parameters[idx_mass_cal+1][2] * linear_interpol_weight
        mass_calib_params['p3'] = mass_cal_parameters[idx_mass_cal][3] * (1-linear_interpol_weight) + mass_cal_parameters[idx_mass_cal+1][3] * linear_interpol_weight   
    mass_centre   = f_tofidx_to_mass(tofidx_peak, mass_calib_params)
    mass_centre_u = f_tofidx_to_mass(tofidx_peak + tofidx_peak_u, mass_calib_params) - mass_centre
        
    return mass_centre, mass_centre_u

    
def do_mass_cal(rt_meas_list, tofidx_meas_list, mass_cal_parameters, mass_cal_mode):
    """
    Do mass calibration:
        for each time, calculate corresponding mass calibration parameters at this time.
        Then apply this set of parameters to the given measured time of flight.
    """
    
    mass_meas_list =[None]*len(rt_meas_list)

    rt_cal_list = [line[0] for line in mass_cal_parameters]
    p1_cal_list = [line[1] for line in mass_cal_parameters]
    p2_cal_list = [line[2] for line in mass_cal_parameters]
    p3_cal_list = [line[3] for line in mass_cal_parameters]
    
    #calculate values of parameters for each time in rt_list:
    p1 = np.interp(rt_meas_list, rt_cal_list,  p1_cal_list)
    p2 = np.interp(rt_meas_list, rt_cal_list,  p2_cal_list)
    p3 = np.interp(rt_meas_list, rt_cal_list,  p3_cal_list)
    
    for i in range(len(rt_meas_list)):
        mass_calib_params = {}
        mass_calib_params['mode'] = mass_cal_mode
        
        if mass_cal_mode ==2:
            mass_calib_params['p1'] = p1[i]
            mass_calib_params['p2'] = p2[i]
            mass_calib_params['p3'] = p3[i]
            
            mass_meas_list[i] = f_tofidx_to_mass(tofidx_meas_list[i], mass_calib_params)
        
    return mass_meas_list

def interpol_mass_cal(rt_meas, mass_cal_parameters, mass_cal_mode):
    """
    Do interpolation of mass calibration parameters:
    for each time, calculate corresponding mass calibration parameters at this time.
    
    rt_meas: float (not list of float)
    """

    mass_cal_param_interp = {}
    mass_cal_param_interp['mode'] = mass_cal_mode
    
    if mass_cal_mode == 2:
    
        rt_cal_list = [line[0] for line in mass_cal_parameters]
        p1_cal_list = [line[1] for line in mass_cal_parameters]
        p2_cal_list = [line[2] for line in mass_cal_parameters]
        p3_cal_list = [line[3] for line in mass_cal_parameters]
        
        #calculate values of parameters for each time in rt_list:
        mass_cal_param_interp['p1'] = np.interp(rt_meas, rt_cal_list,  p1_cal_list)
        mass_cal_param_interp['p2'] = np.interp(rt_meas, rt_cal_list,  p2_cal_list)
        mass_cal_param_interp['p3'] = np.interp(rt_meas, rt_cal_list,  p3_cal_list)
    
        
    return mass_cal_param_interp

def do_mass_cal_np(rt_meas_list, tofidx_meas_list, mass_cal_parameters, mass_cal_mode):
    """
    save as above but returns a numpy array
    """
    
    mass_meas_list =np.empty(len(rt_meas_list))

    rt_cal_list = [line[0] for line in mass_cal_parameters]
    p1_cal_list = [line[1] for line in mass_cal_parameters]
    p2_cal_list = [line[2] for line in mass_cal_parameters]
    p3_cal_list = [line[3] for line in mass_cal_parameters]
    
    #calculate values of parameters for each time in rt_list:
    p1 = np.interp(rt_meas_list, rt_cal_list,  p1_cal_list)
    p2 = np.interp(rt_meas_list, rt_cal_list,  p2_cal_list)
    p3 = np.interp(rt_meas_list, rt_cal_list,  p3_cal_list)
    
    for i in range(len(rt_meas_list)):
        mass_calib_params = {}
        mass_calib_params['mode'] = mass_cal_mode
        
        if mass_cal_mode ==2:
            mass_calib_params['p1'] = p1[i]
            mass_calib_params['p2'] = p2[i]
            mass_calib_params['p3'] = p3[i]
            
        mass_meas_list[i] = f_tofidx_to_mass(tofidx_meas_list[i], mass_calib_params)
        
    return mass_meas_list    


def f_tofidx_to_mass(tofidx, mass_calib_params):
    """
    calculate mass [m/z] based on time of flight index
    Input values:
    tofidx: numpy array containing floats, 
    centre of tof peak, result of tof peak integration
    mass_calib_params: dictionnary containing:
        -the mode of the equation, string
        -the parameters of the tofidx_to_mass equation, floats; how many there are depends on the mode
    Return:
    mass [m/z], numpy array containing floats, corresponding to tofidx
    Note 20190724: mode 3, 4, 5 are not (yet?) supported, mode 2 with default values is used instead
    """
    #print(tofidx)
    if mass_calib_params['mode'] == 0:
        mass = ((tofidx - mass_calib_params['p2']) / mass_calib_params['p1'])**2

    elif mass_calib_params['mode'] == 1:
        mass = (mass_calib_params['p1'] / (tofidx - mass_calib_params['p2']))**2

    elif mass_calib_params['mode'] == 2:
        mass = ((tofidx - mass_calib_params['p2']) / mass_calib_params['p1'])**(1/mass_calib_params['p3'])

    elif mass_calib_params['mode'] == 3:
        """
        here tofidx should be a float
        """
        
        if type(tofidx) == 'numpy.float64':
            len_tofidx = 1
        else:
            len_tofidx = len(tofidx)
        #mass = [None]*len_tofidx
        #coeff = [None]*len_tofidx
        #root = [None]*len_tofidx
        #for i in range(len_tofidx):
        coeff = np.array([
                mass_calib_params['p3'],
                float(0.0),
                2*mass_calib_params['p3']*mass_calib_params['p4'],
                mass_calib_params['p1'],
                mass_calib_params['p2'] + mass_calib_params['p3']*mass_calib_params['p4']**2 - tofidx
                ])
        root = np.roots(coeff)
        #print('root: '+ str(root))
        #for mass_root in root[i]:
        for mass_root in root:
            if mass_root.real > 0 and mass_root.imag == 0 :
                mass = mass_root.real**2
                #print('mass: ' + str(mass[i]))
        
    elif mass_calib_params['mode'] == 4:
        """
        here tofidx should be a float
        """
        #for some reason there are two constants, p2 and p5
        #coeffs should be ordered by decreasing n, with y = p_i * x^n
        """
        if type(tofidx) == 'numpy.float64':
            len_tofidx = 1
        else:
            len_tofidx = len(tofidx)
        """
        len_tofidx = 1
        #mass = [None]*len_tofidx
        #coeff = [None]*len_tofidx
        #root = [None]*len_tofidx
        #for i in range(len_tofidx):
        coeff = np.empty((len_tofidx, 5), dtype = float)
        #coeff = np.empty((len_tofidx), dtype = float)
        root = np.empty(len_tofidx, dtype = float)
        mass = np.zeros(len_tofidx, dtype = float)
        for i in range(len_tofidx):
            coeff[i,:] = [
                    mass_calib_params['p3'],
                    float(0.0),
                    mass_calib_params['p4'], 
                    mass_calib_params['p1'], 
                    mass_calib_params['p2'] + mass_calib_params['p5'] - tofidx[i]
                    ]
            #coeff[:, 4] -= tofidx
            #root[i] = np.roots(coeff[i, :])
            root_0 = np.roots(coeff[i, :])
            #print('root: '+ str(root))
            #for i in range(len_tofidx):
            
            for mass_root in root_0:
                if mass_root.real > 0 and mass_root.imag == 0 :
                    mass[i] = mass_root.real**2
                    #print('mass: ' + str(mass[i]))
                            
    #elif mass_calib_params['mode'] == 5:            
    #else: mass = ((tofidx +3729.0) / 2828.0)**(1/0.4999)
    return mass

def f_mass_to_tofidx(mass, mass_calib_params):
    """Return tof index from mass
    mass is a numpy array containing floats or a single float
    mass_calib_params is a dictionnary containing the calibration mode
    and the corresponding calibration parameters
    """
    if mass_calib_params['mode'] == 0:
        tofidx = mass_calib_params['p1'] * np.sqrt(mass) + mass_calib_params['p2']

    elif mass_calib_params['mode'] == 1:
        tofidx = mass_calib_params('p1') * 1/np.sqrt(mass) + mass_calib_params['p2']

    elif mass_calib_params['mode'] == 2:
        tofidx = mass_calib_params['p1'] * mass**(mass_calib_params['p3']) + mass_calib_params['p2']

    elif mass_calib_params['mode'] == 3:
        tofidx = mass_calib_params['p1'] * np.sqrt(mass) + mass_calib_params['p2']
        + mass_calib_params['p3'] * (mass - mass_calib_params['p4'])**2

    elif mass_calib_params['mode'] == 4: #for files from 2018, it seems this is mode 4
        #this mode does not seem optimal: why having two constants p2 and p5?
        tofidx = mass_calib_params['p1'] * np.sqrt(mass) + mass_calib_params['p2']
        + mass_calib_params['p3'] * mass**2 + mass_calib_params['p4'] * mass + mass_calib_params['p5']
    return tofidx

def cost_f_tof_to_mass(mass_cal_params, tofidx, mass):
    """
    x, np.array of float: the centre of peak in tof domain
    y, np.array of float: the theoretical, exact masses with +1 charge
    area, np.array of float: the area of the tof peaks (signal intensity in V * tof_index)
    """ 
    model = f_tofidx_to_mass(tofidx, mass_cal_params)
    #weight = area /np.sum(area)
    cost = (model - mass) #* weight
    return cost

def fit_tofidx_to_mass_mode(tofidx_meas, mass_exact, params_default, mass_diff_ppm_threshold, graph):
    
    #20191030: use mode 0 to calculate default parameters
    #this should give the general square root shape.
    #x: ToF indexes (result of peak integration), np array of floats.
    #y: exact masses, np array of floats.   
    #transform with M = sqrt(m)
    #to then solve as a linear function
    
    
    M = np.sqrt(mass_exact)
    p1_guess, p2_guess, linear_fit, std_dev_linear_fit = fit_linear_np(M, tofidx_meas)
    #print(str(p1_guess) +'; ' + str(p2_guess))
    #p3_guess_m2 = float(0.5)
    
    count_delete = 1
    mass_cal_params = {}
    
    
    tofidx_meas_archive = tofidx_meas.copy()
    mass_exact_archive = mass_exact.copy()

    #no optimisation needed for mode 0
    if params_default['mode'] == 0:  
        mass_cal_params['mode'] = 0
        mass_cal_params['p1'] = p1_guess
        mass_cal_params['p2'] = p2_guess
        mass_cal_params['p3'] = 0 #to make plotting possible without modification


    #optimisation needed for all other modes
    
    elif params_default['mode'] == 2:
        no_params = 3
        guess_params = Parameters()
        guess_params.add('p3', value = float(0.5), min = float(0.498), max = float(0.502))
        guess_params.add('p4', value = 0, vary = False)
        guess_params.add('p5', value = 0, vary = False)
        

    elif params_default['mode'] == 3:
        #use np.root! no, that's to calculate m
        no_params = 4
        guess_params = Parameters()
        guess_params.add('p3', value = float(0.0))
        guess_params.add('p4', value = float(0.0))
        guess_params.add('p5', value = 0, vary = False)
  
    elif params_default['mode'] == 4:
        #use np.root! no, that's to calculate m
        no_params = 4
        guess_params = Parameters()
        guess_params.add('p3', value = float(0.0))
        guess_params.add('p4', value = float(0.0))
        guess_params.add('p5', value = 0, vary = False)

    
    if params_default['mode'] == 2 \
    or params_default['mode'] == 3 \
    or params_default['mode'] == 4:
        guess_params.add('mode', value = params_default['mode'], vary = False)
        guess_params.add('p1', value = p1_guess)
        guess_params.add('p2', value = p2_guess)
        
        while len(tofidx_meas) > no_params and count_delete > 0:
            
            #print(guess_params)

            minner = Minimizer(cost_f_tof_to_mass, guess_params, fcn_args=(tofidx_meas, mass_exact)) #, method='leastsq'
            tof_to_mass_optimised = minner.minimize(method = 'leastsq')
         
            #resulting optimised parameters are callable with:
            #peak_optimised.params
            guess_params['p1'].value = tof_to_mass_optimised.params['p1'].value
            guess_params['p2'].value = tof_to_mass_optimised.params['p2'].value
            guess_params['p3'].value = tof_to_mass_optimised.params['p3'].value
            guess_params['p4'].value = tof_to_mass_optimised.params['p4'].value
            guess_params['p5'].value = tof_to_mass_optimised.params['p5'].value
            
            model = f_tofidx_to_mass(tofidx_meas, guess_params)
            mass_diff_ppm = abs((model-mass_exact)/mass_exact*10**6)
            #print(mass_diff_ppm)
            
            #here, if one data point is really wrong, this will disturb the entire calibration
            #and no points will fit withn 20 ppm.
            #bad data should be eliminated one by one.
            
            max_mass_diff_ppm = max(mass_diff_ppm)
            
            if max_mass_diff_ppm > mass_diff_ppm_threshold:                
            
                #criteria = abs(mass_diff_ppm) <= mass_diff_ppm_threshold
                criteria = mass_diff_ppm < max_mass_diff_ppm
                
                tofidx_meas = tofidx_meas[criteria]
                mass_exact = mass_exact[criteria]
                #data_idx_list = data_idx_list[criteria]
                
            count_delete = len(mass_diff_ppm) - len(tofidx_meas)
            #print('count delete: ' + str(count_delete))
            #print('len data: ' + str(len(tofidx_meas)))

    if len(tofidx_meas) > no_params:
        mass_cal_params['mode'] = tof_to_mass_optimised.params['mode'].value
        mass_cal_params['p1'] = tof_to_mass_optimised.params['p1'].value
        mass_cal_params['p2'] = tof_to_mass_optimised.params['p2'].value
        mass_cal_params['p3'] = tof_to_mass_optimised.params['p3'].value
            
        if params_default['mode'] == 3 or params_default['mode'] == 4:
            mass_cal_params['p4'] = tof_to_mass_optimised.params['p4'].value
            
        if params_default['mode'] == 4:
            mass_cal_params['p5'] = tof_to_mass_optimised.params['p5'].value
    
        if graph == True:
            print('Optimised parameters:')
            tof_to_mass_optimised.params.pretty_print()
            """
            fig, ax1 = plt.subplots(figsize=(8,4))
            fig.suptitle('mass vs tof index fit', fontsize=12)
            ax1.set_xlabel('mass [m/z]')
            ax1.set_ylabel('ToF index', color='b')
            ax1.tick_params('y', colors='b')
            ax1.plot(np.sqrt(exact_mass_rt_list), measured_tofidx_rt_list, 'b.')
            #ax1.plot(exact_mass_rt_list, apex_guess, 'o', color='xkcd:peach', markersize = 13)
            fig.tight_layout()
            plt.show()
            """
    
        #is already extracted:
        #model = f_tofidx_to_mass(tofidx_meas, mass_cal_params)
        #mass_diff_ppm = (model-mass_exact)/mass_exact*10**6
        #re-calculate mass_diff_ppm over the entire data range:
        model = f_tofidx_to_mass(tofidx_meas_archive, guess_params)
        mass_diff_ppm = (model-mass_exact_archive)/mass_exact_archive*10**6
        
    else: 

        mass_cal_params['mode'] = 2
        mass_cal_params['p1'] = p1_guess
        mass_cal_params['p2'] = p2_guess
        mass_cal_params['p3'] = 0.5 #to make plotting possible without modification

        #mass_cal_params['p1'] = []
        #mass_cal_params['p2'] = []
        #mass_cal_params['p3'] = []
        mass_cal_params['p4'] = []
        mass_cal_params['p5'] = []
        mass_diff_ppm = []
        #data_idx_list = []
        model = f_tofidx_to_mass(tofidx_meas_archive, mass_cal_params)
        mass_diff_ppm = (model-mass_exact_archive)/mass_exact_archive*10**6

    


    return mass_cal_params, mass_diff_ppm #, data_idx_list #, tofidx_meas, mass_exact
   
"""
def fit_tofidx_to_mass_mode(tofidx_meas, mass_exact, data_idx_list, area, params_default, mass_diff_ppm_threshold, graph):
    
    #20191030: use mode 0 to calculate default parameters
    #this should give the general square root shape.
    #x: ToF indexes (result of peak integration), np array of floats.
    #y: exact masses, np array of floats.   
    #transform with M = sqrt(m)
    #to then solve as a linear function
    M = np.sqrt(mass_exact)
    p1_guess, p2_guess, linear_fit, std_dev_linear_fit = fit_linear_np(M, tofidx_meas)
    #print(str(p1_guess) +'; ' + str(p2_guess))
    #p3_guess_m2 = float(0.5)
    
    count_delete = 1
    mass_cal_params = {}

    #no optimisation needed for mode 0
    if params_default['mode'] == 0:  
        mass_cal_params['mode'] = 0
        mass_cal_params['p1'] = p1_guess
        mass_cal_params['p2'] = p2_guess
        mass_cal_params['p3'] = 0 #to make plotting possible without modification


    #optimisation needed for all other modes
    
    elif params_default['mode'] == 2:
        no_params = 3
        guess_params = Parameters()
        guess_params.add('p3', value = float(0.5), min = float(0.4995), max = float(0.5005))
        guess_params.add('p4', value = 0, vary = False)
        guess_params.add('p5', value = 0, vary = False)
        

    elif params_default['mode'] == 3:
        #use np.root! no, that's to calculate m
        no_params = 4
        guess_params = Parameters()
        guess_params.add('p3', value = float(0.0))
        guess_params.add('p4', value = float(0.0))
        guess_params.add('p5', value = 0, vary = False)
  
    elif params_default['mode'] == 4:
        #use np.root! no, that's to calculate m
        no_params = 4
        guess_params = Parameters()
        guess_params.add('p3', value = float(0.0))
        guess_params.add('p4', value = float(0.0))
        guess_params.add('p5', value = 0, vary = False)

    
    if params_default['mode'] == 2 \
    or params_default['mode'] == 3 \
    or params_default['mode'] == 4:
        guess_params.add('mode', value = params_default['mode'], vary = False)
        guess_params.add('p1', value = p1_guess)
        guess_params.add('p2', value = p2_guess)
        
        while len(tofidx_meas) > no_params and count_delete > 0:
            
            #print(guess_params)

            minner = Minimizer(cost_f_tof_to_mass, guess_params, fcn_args=(tofidx_meas, mass_exact, area)) #, method='leastsq'
            tof_to_mass_optimised = minner.minimize(method = 'leastsq')
         
            #resulting optimised parameters are callable with:
            #peak_optimised.params
            guess_params['p1'].value = tof_to_mass_optimised.params['p1'].value
            guess_params['p2'].value = tof_to_mass_optimised.params['p2'].value
            guess_params['p3'].value = tof_to_mass_optimised.params['p3'].value
            guess_params['p4'].value = tof_to_mass_optimised.params['p4'].value
            guess_params['p5'].value = tof_to_mass_optimised.params['p5'].value
            
            model = f_tofidx_to_mass(tofidx_meas, guess_params)
            mass_diff_ppm = abs((model-mass_exact)/mass_exact*10**6)
            #print(mass_diff_ppm)
            
            #here, if one data point is really wrong, this will disturb the entire calibration
            #and no points will fit withn 20 ppm.
            #bad data should be eliminated one by one.
            
            max_mass_diff_ppm = max(mass_diff_ppm)
            
            if max_mass_diff_ppm > mass_diff_ppm_threshold:                
            
                #criteria = abs(mass_diff_ppm) <= mass_diff_ppm_threshold
                criteria = mass_diff_ppm < max_mass_diff_ppm
                
                tofidx_meas = tofidx_meas[criteria]
                mass_exact = mass_exact[criteria]
                data_idx_list = data_idx_list[criteria]
                
            count_delete = len(mass_diff_ppm) - len(tofidx_meas)
            #print('count delete: ' + str(count_delete))
            #print('len data: ' + str(len(tofidx_meas)))

    if len(tofidx_meas) > no_params:
        mass_cal_params['mode'] = tof_to_mass_optimised.params['mode'].value
        mass_cal_params['p1'] = tof_to_mass_optimised.params['p1'].value
        mass_cal_params['p2'] = tof_to_mass_optimised.params['p2'].value
        mass_cal_params['p3'] = tof_to_mass_optimised.params['p3'].value
            
        if params_default['mode'] == 3 or params_default['mode'] == 4:
            mass_cal_params['p4'] = tof_to_mass_optimised.params['p4'].value
            
        if params_default['mode'] == 4:
            mass_cal_params['p5'] = tof_to_mass_optimised.params['p5'].value
    
        if graph == True:
            print('Optimised parameters:')
            tof_to_mass_optimised.params.pretty_print()
    
        #is already extracted:
        #model = f_tofidx_to_mass(tofidx_meas, mass_cal_params)
        #mass_diff_ppm = (model-mass_exact)/mass_exact*10**6
        
    else: 
        mass_cal_params['p1'] = []
        mass_cal_params['p2'] = []
        mass_cal_params['p3'] = []
        mass_cal_params['p4'] = []
        mass_cal_params['p5'] = []
        mass_diff_ppm = []
        data_idx_list = []


    return mass_cal_params, mass_diff_ppm, data_idx_list #, tofidx_meas, mass_exact
   
    
def fit_tofidx_to_mass(x, y, area, p1_guess, p2_guess, p3_guess, graph):
    #all the parameters should be given as a dictionnary
    #to merge all fitting functions
    params = Parameters()
    params.add('p1', value = p1_guess)
    params.add('p2', value = p2_guess)
    params.add('p3', value = p3_guess, min = float(0.48), max = float(0.52))
    #print(params['p1'])
    #print(params['p2'])
    #print(params['p3'])
    
    #print(x)
    #print(y)
    #test = cost_f_tof_to_mass(params, x, y)
    #print(test)

    # do fit, here with leastsq model
    #diff_model_data, cost = 
    minner = Minimizer(cost_f_tof_to_mass, params, fcn_args=(x, y, area)) #, method='leastsq'
    tof_to_mass_optimised = minner.minimize(method = 'leastsq')
    
    #print results of optimisation:
    if graph == True:
        print('Optimised parameters:')
        tof_to_mass_optimised.params.pretty_print()
    
    #resulting optimised parameters are callable with:
    #peak_optimised.params
    tof_to_mass_optimised_params=tof_to_mass_optimised.params
    
    mass_cal_params = {}
    mass_cal_params['mode'] = 2
    mass_cal_params['p1'] = tof_to_mass_optimised_params['p1'].value
    mass_cal_params['p2'] = tof_to_mass_optimised_params['p2'].value
    mass_cal_params['p3'] = tof_to_mass_optimised_params['p3'].value
    
    y_model = f_tofidx_to_mass(x, mass_cal_params)
    mass_diff_ppm = (y_model-y)/y*10**6
    return mass_cal_params['p1'], mass_cal_params['p2'], mass_cal_params['p3'], mass_diff_ppm
"""
#**************************************************************************
def fit_linear(x, y):
    n= len(y)
    x_sum = sum(x)
    y_sum = sum(y)
    x_square_sum = sum([x[i]**2 for i in range(n)])
    x_y_sum = sum([x[i] * y[i] for i in range(n)])
    num_linear_fit_slope = n * x_y_sum - x_sum*y_sum
    num_linear_fit_b = x_square_sum * y_sum - x_sum * x_y_sum
    den_linear_fit = n * x_square_sum - x_sum**2

    if den_linear_fit != 0:
        linear_fit_slope = num_linear_fit_slope/den_linear_fit
        linear_fit_b = num_linear_fit_b/den_linear_fit
    else: 
        linear_fit_slope = 0
        linear_fit_b = sum(y)/n #default: straight line
        
    linear_fit = [x[i] * linear_fit_slope + linear_fit_b for i in range(n)]
    std_dev_linear_fit = np.sqrt(sum([(linear_fit[i]- y[i])**2 for i in range(n)]) / float(n-1))
    return linear_fit_slope, linear_fit_b, linear_fit, std_dev_linear_fit


#**************************************************************************
def fit_linear_np(x, y):
    """linear fit with np.array as inputs
    x, y: np arrays of float
    """
    n= len(y)
    x_sum = np.sum(x)
    den = np.sum(x**2) - x_sum**2/n
    if den != 0:
        y_sum = np.sum(y)
        x_y_sum = np.sum(x * y)
        num = x_y_sum - x_sum*y_sum/n      
        slope = num/den
        num_b = (np.sum(x**2)  * y_sum - x_sum * x_y_sum) / n
        b = num_b/den
    else : 
        slope = 0
        b = np.sum(y)/n
     
    linear_fit = x * slope + b 
    std_dev_linear_fit = np.sqrt(np.sum((linear_fit- y)**2) / float(n-1))
    return slope, b, linear_fit, std_dev_linear_fit

#*************************************************************************
    
def cost_f_iso_profile(params, iso_profiles, cl_comp):
    #y is the target, measured mass profile.
    
    #update mass-profile based on mass-offset:
    cl_comp.do_meas_mass_profile(params['mass_offset'])
    
    model = np.zeros(len(cl_comp.meas_mass_profile_x))
    for idx_peak in range(params['nb_peaks'].value):
        model += params['k'+str(idx_peak)] * iso_profiles[idx_peak]

    cost = (model - cl_comp.meas_mass_profile)  #* weight_list
    return cost
    
    
def fit_iso_profile(cl_comp, delta_mass, run_type:RunType, iso_idx_list, iso_idx_list_k_guesses, iso_profiles, iso_fit_mode, graph):

    fit_mode_fixed_x_axis = 0
    fit_mode_variable_x_axis = 1    
    #count_delete_sum = 1
    #while len(iso_idx_list) > 0 and count_delete_sum > 0:
    #graph = True
    params = Parameters()
    params.add('nb_peaks', value =len(iso_idx_list),  vary = False)
    #params.add('iso_idx_list', value =len(iso_idx_list),  vary = False)
    
    if run_type.is_tofwerk_tof() or iso_fit_mode == fit_mode_variable_x_axis:
        params.add('mass_offset', value = delta_mass , min = delta_mass - 2* cl_comp.delta_mass_u, max = delta_mass + 2* cl_comp.delta_mass_u)
    elif run_type.is_NIST() or run_type.is_no_mass_cal() or run_type.is_unit_mass() or iso_fit_mode == fit_mode_fixed_x_axis:
        params.add('mass_offset', value = delta_mass, vary=False)
        
        
    for idx_peak in range(len(iso_idx_list)):
        #parameters per peak:
        params.add('k'+str(idx_peak), value = iso_idx_list_k_guesses[idx_peak], min = cl_comp.k_guesses_min) #, max = cl_comp.k_guesses_max, vary = iso_idx_list_k_guesses[idx_peak] > 0
        #params.add('iso_idx'+str(idx_peak), value = iso_idx_list[idx_peak], vary = False)
    #***end of make parameters***
    
    minner = Minimizer(cost_f_iso_profile, params, fcn_args=(iso_profiles, cl_comp)) #, method='leastsq'
    peak_optimised = minner.minimize(method = 'leastsq')
    
    delta_mass_opt = peak_optimised.params['mass_offset'].value
    #k_opt = np.zeros(len(iso_idx_list))
    k_opt = [None] *len(iso_idx_list)
    profile_opt = np.zeros(len(cl_comp.meas_mass_profile_x))
    for idx_peak in range(len(iso_idx_list)):
        k_opt[idx_peak] = peak_optimised.params['k'+str(idx_peak)].value
        profile_opt+= k_opt[idx_peak] * iso_profiles[idx_peak]
        
    #test: per iso profile, are all peaks below the LOD? if yes, remove index from list. or not?
    #alternative: fix k to zero and vary to false
    
    if graph == True:
        text_size = 12
        fig, ax1 = plt.subplots(figsize=(5,4))
        #fig.suptitle('Meas mass profile vs isotope profiles', fontsize=12)
        ax1.set_xlabel('Mass [m/z]', fontsize=text_size)
        ax1.set_ylabel('Intensity [area]', fontsize=text_size) #, color='b'
        #ax1.tick_params('y', colors='b')
        ax1.set_yscale('log')
        #ax2 = ax1.twinx()
        #ax2.plot(cl_comp.meas_mass_profile_x, (peak_optimised.residual*(-1)), 'c.')
        #ax2.set_ylabel('residuals', color='c')
        #ax2.set_ylim(-max(cl_comp.meas_mass_profile), max(cl_comp.meas_mass_profile))
        #ax2.tick_params('y', colors='c')    


        for i_iso in range(len(iso_idx_list)):
            ax1.plot(cl_comp.meas_mass_profile_x, iso_profiles[i_iso] * k_opt[i_iso], '-', color = str(color1_names[i_iso]))
        #ax1.plot(mu_guess, apex_guess, 'o', color='xkcd:peach', markersize = 12)
        ax1.plot(cl_comp.meas_mass_profile_x, cl_comp.meas_mass_profile, 'k--')
        
            

        #fig.tight_layout()
        plt.show()
        """
        for i_iso in range(len(iso_idx_list)):
            text_size = 12
            fig, ax1 = plt.subplots(figsize=(5,2))
            #fig.suptitle('Meas mass profile vs isotope profiles', fontsize=12)
            ax1.set_xlabel('Mass [m/z]', fontsize=text_size)
            ax1.set_ylabel('Intensity [area]', fontsize=text_size) #, color='b'
            #ax1.tick_params('y', colors='b')
            ax1.set_yscale('log')
            ax1.plot(cl_comp.meas_mass_profile_x, iso_profiles[i_iso] * k_opt[i_iso], '-', color = str(color1_names[i_iso]))

            plt.show()
        """
        
    
    return k_opt, delta_mass_opt
    
#*************************************************************************    

def f_pseudoVoigt(x, A, mu, sigma, alpha):
    """
    https://lmfit.github.io/lmfit-py/builtin_models.html#pseudovoigtmodel
    weighted sum of a Gaussian and Lorentzian distribution function 
    that share values for amplitude (A), center (mu) and full width at half maximum fwhm 
    (and so have constrained values of sigma and height (maximum peak height). 
    A parameter fraction (alpha) controls the relative weight of the Gaussian and Lorentzian components
    x: np array
    A, mu, sigma, alpha : floats (cf. Gaussian function)
    Note: this function does not include the baseline. 
    If any, it has to be added separately.
    """
    sigma_L= sigma * np.sqrt(2 * np.log(2))
    return (1-alpha) * f_Gauss(x, A, mu, sigma) + alpha * A/np.pi * sigma_L / ((x-mu)**2 + sigma_L**2) 

def cost_f_pseudoVoigt(params, x, y):
    """
    Calculate the cost function to minimise
    This cost function is adapted to the number of expected peaks
    It returns the quantity to be minimised using later the minimise function from limfit
    """
    #model=params['baseline']
    model = np.zeros(len(x)) + params['baseline']
    for idx_peak in range(params['nb_peaks'].value):
        model += f_pseudoVoigt(x, 
                               params['A'+str(idx_peak)], 
                               params['mu'+str(idx_peak)], 
                               params['sigma'+str(idx_peak)], 
                               params['alpha'+str(idx_peak)])
    cost = (model - y) #* weight_list
    return cost

def fit_pseudoVoigt(x, y, mu_guess, apex_guess, alpha_guess, n_per_bin, area_threshold, sigma_tofidx_slope, sigma_tofidx_b, graph, mode):
    """This function computes a (or several) pseudo-Voigt fit(s) of the data. 
    Pseudo-Voigt profile is parameterized by mu, sigma, A
    :param x: np array of floating point numbers, x-coordinates
    :param y: np array of floating point numbers, y-coordinates
    :param mu_guess: np array of guess of mu parameter, float
    :param apex_guess: np array of value of y when x=mu, float, max of peak
    :param n_per_bin: number of points x in between two integer values (sampling rate for x)
    :param area_threshold: minimum value of area, float
    
    returns
    mu_opt, optimised centre of peak
    mu_opt_stderr, 
    A_opt, optimised area
    sigma_opt, 
    alpha_opt, 
    baseline_opt, optimised baseline value.
    """

    #sum_of_signal = sum(y)
    #find indexes to order by increasing values:
    #argort is a numpy method that can be applied to a numpy array.
    y_increasing_indexes = y.argsort()
    y_median_test = np.median(y[y_increasing_indexes[0:int(len(y)/3)]])
    #assuming baseline == 0 is actually quite ok, 20200215
    #baseline_guess = float(0.0)
    
    #assuming same baseline for all peaks in a small window:
    #x should cover a mass range large enough around unit mass (or same in tof_idx)
    #so that minimum can be found where no mass is expected
    #i.e. unit mass +/- 0.3 m/z or more   
    baseline_guess = max(0, y_median_test)


    #plot first guess peak position and apex:
    if graph == True:
        fig, ax1 = plt.subplots(figsize=(8,4))
        fig.suptitle('Pseudo-Voigt fit over ToF index domain', fontsize=12)
        ax1.set_xlabel('ToF index [index]')
        ax1.set_ylabel('Intensity [V]', color='b')
        ax1.tick_params('y', colors='b')
        ax1.plot(mu_guess, apex_guess, 'o', color='xkcd:peach', markersize = 12)
        ax1.plot([min(x), max(x)], [y_median_test, y_median_test], 'r--')
        ax1.plot(x, y, '--bo') #'b.'
        

    #initialise while loop:
    count_delete_sum = 1
    start = 0
    sigma_guess = sigma_tofidx_b + sigma_tofidx_slope*mu_guess    
    alpha_guess = [alpha_guess]*len(mu_guess)
    #A_guess = np.divide((apex_guess - baseline_guess), 
    #                    (1.0-alpha_guess[0])/(sigma_guess * np.sqrt(2*np.pi)) + alpha_guess[0] / (math.pi * sigma_guess * np.sqrt(2 * np.log(2))))
    A_guess = np.divide(apex_guess - baseline_guess, (1.0-alpha_guess[0])/(sigma_guess * np.sqrt(2*np.pi)) + alpha_guess[0] / (math.pi * sigma_guess * np.sqrt(2 * np.log(2))))

    if mode == 'MassCal':
        sigma_guess_vary = True
        alpha_guess_vary = True
    else:
        sigma_guess_vary = False #change to True for improved precision (slower)
        alpha_guess_vary = False            

    while len(mu_guess) > 0 and count_delete_sum > 0:
        #***make parameters***
        params = Parameters()
        params.add('baseline', value=max(1, baseline_guess), min = max(0, min(y)-3*min(y)), max = max(3, baseline_guess*3, y_median_test*3)) #
        #params.add('baseline', value=max(1, baseline_guess), min = max(0, baseline_guess/3), max = max(3, baseline_guess*3)) #
        #params.add('baseline', value=baseline_guess, min = max(0.0, min(y)-3.0*min(y)), max = baseline_guess+y_mean/10.0) #
        #baseline cannot be above max I
        params.add('nb_peaks', value =len(mu_guess),  vary = False)

        mu_guess_min = mu_guess-n_per_bin*2.0
        mu_guess_max = mu_guess+n_per_bin*2.0

        A_guess_min = A_guess*0.1
        A_guess_max = A_guess * 10.0

        if mode == 'MassCal':
            sigma_guess_min = sigma_guess-0.4
            sigma_guess_max = sigma_guess*2.0
            
        else:
            sigma_guess_min = sigma_guess-0.2
            sigma_guess_max = sigma_guess+0.2

        alpha_guess_min = [float(0.2)]*len(mu_guess)
        alpha_guess_max = [float(1.0)]*len(mu_guess)

        for idx_peak in range(len(mu_guess)):
            #parameters per peak:
            params.add('mu'+str(idx_peak), value = mu_guess[idx_peak], min = mu_guess_min[idx_peak], max = mu_guess_max[idx_peak])
            params.add('sigma'+str(idx_peak), value = sigma_guess[idx_peak], min = sigma_guess_min[idx_peak], max = min(3.0, sigma_guess_max[idx_peak]), vary = sigma_guess_vary) #, vary = True
            params.add('alpha'+str(idx_peak), value = alpha_guess[idx_peak] , min = alpha_guess_min[idx_peak], max = alpha_guess_max[idx_peak], vary = alpha_guess_vary)
            params.add('A'+str(idx_peak), value = max(float(1.0), A_guess[idx_peak]) , min = max(1, A_guess_min[idx_peak]), max = max(2, A_guess_max[idx_peak]))
        #***end of make parameters***
        
        minner = Minimizer(cost_f_pseudoVoigt, params, fcn_args=(x, y)) #, method='leastsq'
        peak_optimised = minner.minimize(method = 'leastsq')

        #retrieve optimised parameters
        #peak_optimised_params = peak_optimised.params

        """
        if start == 0 and graph == True:
            #print start parameters just once
            start = 1 
            print('****************************')
            print('Guessed parameters: ')
            params.pretty_print()
        """
      
        A_opt = np.zeros(len(mu_guess))
        mu_opt = np.zeros(len(mu_guess))
        sigma_opt = np.zeros(len(mu_guess))
        alpha_opt = np.zeros(len(mu_guess))
        max_opt = np.zeros(len(mu_guess))
        baseline_opt = peak_optimised.params['baseline'].value
        for idx_peak in range(len(mu_guess)):
            A_opt[idx_peak] = peak_optimised.params['A'+str(idx_peak)].value
            mu_opt[idx_peak] = peak_optimised.params['mu'+str(idx_peak)].value
            sigma_opt[idx_peak] = peak_optimised.params['sigma'+str(idx_peak)].value
            alpha_opt[idx_peak] = peak_optimised.params['alpha'+str(idx_peak)].value
            max_opt[idx_peak] = (1.0 - alpha_opt[idx_peak]) * A_opt[idx_peak] / \
            (sigma_opt[idx_peak] * math.sqrt(2.0 * math.pi)) \
            + alpha_opt[idx_peak] * A_opt[idx_peak]/(math.pi * sigma_opt[idx_peak] * np.sqrt(2 * np.log(2))) + baseline_opt

        if mode == 'MassCal':
            area_test = np.where(max_opt >= 1.3* y_median_test)
        else:
            #area_test = np.where(A_opt >= area_threshold) #this is a boolean
            area_test = np.where(max_opt >= 2.6* y_median_test)
        #print(area_test)
        A_guess = A_opt[area_test]
        #print(A_guess)
        count_delete_sum = len(mu_guess) - len(A_guess)
        #print('nb delete: ' + str(count_delete_sum))
        mu_guess = mu_opt[area_test]
        sigma_guess = sigma_opt[area_test]
        alpha_guess = alpha_opt[area_test]
        baseline_guess = baseline_opt

    #***The optimisation is finished***

    #In case peaks have been all eliminated at the previous loop, 
    #all values needs to be emptied:
    if len(mu_guess) == 0:
        mu_opt = []
        #mu_opt_stderr = []
        A_opt= []
        sigma_opt = []
        alpha_opt = []
        #baseline_opt = []

    else:
    
        if graph == True:
            #print('Optimised parameters:')
            #peak_optimised.params.pretty_print()
            ax2 = ax1.twinx()
            ax2.plot(x, peak_optimised.residual*(-1), 'c.')
            ax2.set_ylabel('residuals', color='c')
            ax2.tick_params('y', colors='c')    
            #to plot the fit with good resolution to visualise something:
            x_plot = np.linspace(min(x), max(x), num = int((max(x) - min(x)) /(n_per_bin/20.0)))
            y_opt_sum_plot = [baseline_opt]*len(x_plot) #add baseline only once!
            #y_opt_sum_plot = [0.0]*len(x_plot)
            max_opt = [None]*len(mu_guess)
        #Extract results of fit:    
        #initialisation
        #mu_opt_stderr = [None]*len(mu_guess)
    
        for idx_peak in range(len(mu_guess)):
            #A_opt[idx_peak] = peak_optimised.params['A'+str(idx_peak)].value
            #mu_opt[idx_peak] = peak_optimised.params['mu'+str(idx_peak)].value
            #mu_opt_stderr[idx_peak] = peak_optimised.params['mu'+str(idx_peak)].stderr
            #sigma_opt[idx_peak] = peak_optimised.params['sigma'+str(idx_peak)].value
            #alpha_opt[idx_peak] = peak_optimised.params['alpha'+str(idx_peak)].value
            if graph == True:
                #max_opt[idx_peak] = f_pseudoVoigt(mu_opt[idx_peak], A_opt[idx_peak], mu_opt[idx_peak], sigma_opt[idx_peak], alpha_opt[idx_peak]) + baseline_opt
                #max_opt[idx_peak] = (1.0 - alpha_opt[idx_peak]) * A_opt[idx_peak] / \
                #(sigma_opt[idx_peak] * math.sqrt(2.0 * math.pi)) \
                #+ alpha_opt[idx_peak] * A_opt[idx_peak]/(math.pi * sigma_opt[idx_peak] * np.sqrt(2 * np.log(2))) #+ baseline_opt
                y_opt_plot = f_pseudoVoigt(x_plot, A_opt[idx_peak], mu_opt[idx_peak], sigma_opt[idx_peak], alpha_opt[idx_peak])  
                #ax1.plot(x_plot, y_opt_plot + baseline_opt, 'g--')
                ax1.plot(x_plot, y_opt_plot, 'g--')
                y_opt_sum_plot += y_opt_plot
        
        if graph == True:
            #y_opt_sum_plot = np.sum(, axis = 0)
            ax1.plot(x_plot, y_opt_sum_plot, 'k-')
            ax1.plot(mu_opt, max_opt, '*', color='xkcd:deep red', markersize = 12)
            
            #plt.savefig('Gauss_fit_on_mass_smallA.png', bbox_inches='tight', transparent = True, dpi = 300)
    if graph == True:
        fig.tight_layout()
        plt.show()

    return mu_opt, A_opt, sigma_opt, alpha_opt, baseline_opt #, mu_opt_stderr sum_of_signal #


def f_chromato(x, H, mu, sigma, alpha):
    """Equation from Pap and Papai, Journal of Chromatography, 2001
    Title: Application of a new mathematical function for describing chromatographic peaks.
    x: np.array, retention time (index (int) or seconds (float)?)
    H: float, peak height (value at centre of peak)
    mu: float, position of centre of peak, same unit as x
    sigma: float, standard deviation of the peak.
    alpha: float, parameter for the asymetry of the peak in RT domain, 0<a< 2.
    Note: this function does not include the baseline. 
    The baseline (if any) has to be added separately and just once.
    Returns:
    chromato, np.array, has 1 dimension of lenght of x.
    """
    if alpha > 0.01:
        chromato = np.zeros(len(x))
        v = float(2.0)*alpha*(x-mu) / (sigma*(float(4.0)-alpha**int(2)))
        indexes_valid = np.flatnonzero(v>float(-0.99))
        v_valid = float(2.0)*alpha*(x[indexes_valid]-mu) / (sigma*(float(4.0)-alpha**int(2)))
        chromato[indexes_valid] =  H * np.exp( (float(4.0)/alpha**2 -float(1.0)) * (np.log(float(1.0) + v_valid) - v_valid) )
        
        #For some unknown reason the version hereafter cause sometimes error with log:
        #np.where(condition, value if True, value if False)
        #print(float(1.0) + v)
        #print('f_chromato, v: ' + str(v))
        #chromato = np.where(v>float(-0.99), 
        #                    H * np.exp( (float(4.0)/alpha**2 -float(1.0)) * (np.log(float(1.0) + v) - v) ), 
        #                    float(0.0))
    else: #Gaussian
        A = H * sigma * np.sqrt(2* np.pi)
        chromato = f_Gauss(x, A, mu, sigma)
        
            
    return chromato

def f_chromato_nparray(x, A, mu, sigma, alpha):
    """Equation from Pap and Papai, Journal of Chromatography, 2001
    Title: Application of a new mathematical function for describing chromatographic peaks.
    x: np.array, retention time (index (int) or seconds (float)?)
    A: float, peak height (value at centre of peak)
    mu: float, position of centre of peak, same unit as x
    sigma: float, standard deviation of the peak.
    alpha: float, parameter for the asymetry of the peak in RT domain, 0<a< 2.
    Note: this function does not include the baseline. 
    The baseline (if any) has to be added separately and just once.
    Returns:
    chromato, np.array, has 2 dimensions, the len of mu and then the len of x.
    #20191120 now adapted to np.array as input for peak parameters.
    """
    #v = float(2.0)*alpha*(x-mu) / (sigma*(float(4.0)-alpha**2.0))
    v_i = float(2.0)*alpha / (sigma*(float(4.0)-alpha**2.0))
    chromato = np.zeros(len(mu), len(x))
    for i_p in range(len(mu)):
        #if alpha < threshold, maybe switch to Gauss function.
        v = v_i[i_p]*(x-mu[i_p])
        #np.where(condition, value if True, value if False)
        chromato[i_p] = np.where(v>float(-0.99), A[i_p] * np.exp( (float(4.0)/alpha[i_p]**2 -float(1.0)) * (np.log(float(1.0) + v) - v) ), float(0.0))
    return chromato

def f_chromato_area(A, sigma, alpha):
    """This function calculate the area corresponding to a peak description
    according to Pap and Papai 2001.
    It takes into accunt tailing of the peak (alpha factor), see f_chromato.
    Warning! In practice, calculation of the area fails for alpha < 0.2.
    We use Gaussian equation instead in this case.
    If alpha<0.2, the shape can be well approximated by the Gaussian function anyway.
    The error raised by Python is:
    "RuntimeWarning: invalid value encountered in double_scalars"
    Do not support array!
    A, sigma and alpha must be floats.
    """    

    if alpha > 0.2:
        #use equation with tailing
        w = float(4.0) / alpha**2
        area = gamma(w) * A * alpha * sigma / (2.0 * np.exp((w-1) * (np.log(w-1)-1)))
    else:
        #use Gaussian approximation to avoid computing error
        area = A* (sigma * np.sqrt(2*np.pi)) #height = area/(sigma * np.sqrt(2*np.pi))
    #print(area)
    return area


def f_chromato_area_nparray(A, sigma, alpha):
    """This function calculate the area corresponding to a peak description
    according to Pap and Papai 2001.
    It takes into accunt tailing of the peak (alpha factor), see f_chromato.
    Warning! In practice, calculation of the area fails for alpha < 0.2.
    We use Gaussian equation instead in this case.
    If alpha<0.2, the shape can be well approximated by the Gaussian function anyway.
    The error raised by Python is:
    "RuntimeWarning: invalid value encountered in double_scalars"
    """
    w = float(4.0) / alpha**2
    #area = np.where(alpha > 0.2, gamma(w).dot(A).dot(alpha).dot(sigma) / (2.0 * np.exp((w-1) * (np.log(w-1)-1))), A* (sigma * np.sqrt(2*np.pi)))
    area = np.where(alpha > 0.2, np.divide(np.dot(np.dot(np.dot(gamma(w), A), alpha), sigma), (2.0 * np.exp((w-1) * (np.log(w-1)-1)))), np.dot(A, (sigma * np.sqrt(2*np.pi))))
    
# =============================================================================
#     if alpha > 0.2:
#         #use equation with tailing
#         w = float(4.0) / alpha**2
#         area = gamma(w) * A * alpha * sigma / (2.0 * np.exp((w-1) * (np.log(w-1)-1)))
#     else:
#         #use Gaussian approximation to avoid computing error
#         area = A* (sigma * np.sqrt(2*np.pi)) #height = area/(sigma * np.sqrt(2*np.pi))
#     
# =============================================================================
    #print(area)
    return area

def cost_f_chromato(params, x, y):
    """
    Calculate the cost function to minimise
    This cost function is adapted to the number of expected peaks
    It returns the quantity to be minimised using later the minimise function from limfit
    """
    #nb_chromato_params = int(4) #4 parameters to define the chromato function
    model=np.zeros(len(x)) #params['baseline'] #+str(idx_peak)
    model += params['baseline']
    #this cost function adapts itself to number of found peaks
    for idx_peak in range(params['nb_peaks'].value):
        model += f_chromato(x, params['A'+str(idx_peak)].value, 
                            params['mu'+str(idx_peak)].value, 
                            params['sigma'+str(idx_peak)].value, 
                            params['alpha'+str(idx_peak)].value)
    cost = model - y
    #print(cost)
    return cost

def compute_rt_sigma_guess(sigma_guess_slope, sigma_guess_exp_k, mu_guess):
    """
    Compute the guess value for peak broadness in time domain.
    This is an exponential function of the retention time.
    Parameters have been empirically determined, by fitting results of 
    sigma (from fitted peaks when using a quite large tolerance around sigma),
    vs rt.
    The parameters sigma_guess_slope, sigma_guess_exp_k 
    depend on the chromatographic method used! 
    Using a different temperature ramp will modify them.
    
    mu_guess: float or np.array, rt where sigma_guess is computed
    sigma_guess_slope: float
    sigma_guess_exp_k: float
    """
    rt_sigma_guess = sigma_guess_slope*np.exp(sigma_guess_exp_k*mu_guess)
    
    if mu_guess <= 1000:
        rt_sigma_guess = 2.5
    
    return rt_sigma_guess
    
def fit_chromato(x, y, mu_guess, apex_guess, baseline_guess_0, n_per_bin, conf_int, broadness, sigma_guess_value, rt_window_min, rt_window_max, tof_series, mass_series, ratio_to_linear_fit, graph):
    """This function computes one (or several) chromatographic fit(s) of the data,
    i.e Gaussian with tailing.
    The chromato functin is parameterized by mu, sigma, A, alpha.
    Compared to Gaussian, there is one more parameter: alpha.
    :param x: list of floating point numbers, x-coordinates
    :param y: list of floating point numbers, y-coordinates
    :param mu_guess: list of guess of mu parameter (position of centre on x axis), float
    :param apex_guess: list of value of y when x=mu (apex of peak, on y axis), float
    :param alpha: tailing of the chromatographic peak. 0 > alpha > 2. 
    So far, it seems difficult to constrain value of alpha, for example finding a correlation with total area of peak or elution time.
    :param n_per_bin: number of points x in between two integer values (sampling rate for x)
    (on the rt domain, this is s_per_bin)
    :param area_threshold: minimum value of optimised area, float
    :param ratio_to_linear_fit: minimum fit improvement for Gauss/Chromato compared to linear fit.
    This was chosen compared to simply using a fit residual threshold, because there is 
    a lot of noise close to the LOD.
    
    returns optimised values for:
    mu_opt
    A_opt
    sigma_opt
    alpha_opt
    area_opt, calculated from optimised parameters.
    area_opt_u, std dev of optimised area.
    
    In addition, for each fitted chromatographic peak, calculate:
    - the average centre of tof, using tof_series
    - the average mass centre, using mass_series
    - LOD
    - SN (signal/noise, this is apex_opt/LOD)
    """
    len_x = len(x)
    len_mu_guess = len(mu_guess)
    
    #plot first guess peak position and apex:
    if graph == True:
        fig, ax1 = plt.subplots(figsize=(8,4))
        #fig.suptitle('Chromatographic fit over RT', fontsize=12)
        fig.suptitle('Chromatographic fit over RT at m/z = ' + '{:.3f}'.format(np.mean(mass_series)), fontsize=12)
        ax1.set_xlabel('Retention time [s]')
        ax1.set_ylabel('Intensity [TOF-index * V]', color='b')
        ax1.tick_params('y', colors='b')
        ax1.plot(mu_guess, apex_guess, 'o', color='xkcd:peach', markersize = 13)
        ax1.plot(x, y, 'b.')
        

    #assuming same baseline for all peaks in a small window:
    #x should cover a RT range large enough around the expected peak centre in RT
    #so that minimum can be found where no mass is expected
    #i.e.  guess RT + / - 4* sigma of a peak

    baseline_guess = max(0, baseline_guess_0) #, min(y)*0.8
    #we have now points that corresponds to signal above baseline
    #so the minimum is the closest to the real baseline
    
    #initialise while loop:
    count_delete_sum = 1
    start = 0
   
    #sigma_guess = compute_rt_sigma_guess(sigma_guess_slope, sigma_guess_exp_k, mu_guess)
    sigma_guess = np.zeros(len_mu_guess) + sigma_guess_value
    sigma_guess_min = sigma_guess*0.6
    sigma_guess_max = sigma_guess * 1.6    
    H_guess = apex_guess-baseline_guess
    alpha_guess = [0.1]*len_mu_guess
    
    while len_mu_guess > 0 and count_delete_sum > 0:
        #***make parameters***
        params = Parameters()
        #params.add('baseline', value=max(1, baseline_guess), min = max(0, min(y)-3*min(y)), max = max(3, baseline_guess*3, min(y)+y_mean/10)) #
        params.add('baseline', value=max(float(1.0), baseline_guess), min = max(float(0.0), baseline_guess/3), max = max(float(3.0), baseline_guess*10)) #
        #baseline cannot be above max I
        params.add('nb_peaks', value =len_mu_guess,  vary = False)
                
        mu_guess_min = mu_guess-n_per_bin*3 #np.array ok
        mu_guess_max = mu_guess+n_per_bin*3 #np.array ok
        
        H_guess_min = H_guess/10.0
        H_guess_max = H_guess*10
        
        #would be good to try to define alpha_max as function of area and/or retention time
        alpha_guess_min = [0.001]*len_mu_guess
        alpha_guess_max = [1.6]*len_mu_guess
                
        for idx_peak in range(len_mu_guess):  
            params.add('mu'+str(idx_peak), value = mu_guess[idx_peak], min = mu_guess_min[idx_peak], max = mu_guess_max[idx_peak]) #abcisse cannot be negative
            params.add('sigma'+str(idx_peak), value = sigma_guess[idx_peak], min = max(float(0.5), sigma_guess_min[idx_peak]), max = max(float(0.8), sigma_guess_max[idx_peak]))
            params.add('A'+str(idx_peak), value = max(float(2.0), H_guess[idx_peak]) , min = max(float(1.0), H_guess_min[idx_peak]), max = max(float(3.0), H_guess_max[idx_peak]))
            params.add('alpha'+str(idx_peak), value = alpha_guess[idx_peak], min = alpha_guess_min[idx_peak], max = alpha_guess_max[idx_peak])                  
        
        #fit, here with leastsq model
        minner = Minimizer(cost_f_chromato, params, fcn_args=(x, y)) #, method='leastsq'
        peak_optimised = minner.minimize(method = 'leastsq')
        
        #retrieve optimised parameters: peak_optimised.params    
        
        #Statistical test: is each fitted peak statistically significant?
        baseline_opt = peak_optimised.params['baseline'].value
        t_test_baseline = 0.0 #np.zeros(len(mu_guess))
        H_opt= np.zeros(len_mu_guess)
        mu_opt = np.zeros(len_mu_guess)
        sigma_opt = np.zeros(len_mu_guess)
        alpha_opt = np.zeros(len_mu_guess) #chromato

        apex_opt_u = np.zeros(len_mu_guess)
        y_opt = np.zeros((len_mu_guess,len_x))
        y_opt_sum = np.zeros(len_x)

        t_test_value = np.zeros(len_mu_guess)
        t_test_apex = np.zeros(len_mu_guess)
        
        np_pts = np.zeros(len_mu_guess)
        #np_pts_min = np.zeros(len(mu_guess))
        
        np_pts_left = np.zeros(len_mu_guess)
        np_pts_right = np.zeros(len_mu_guess)
        np_pts_left_min = np.zeros(len_mu_guess)
        np_pts_right_min = np.zeros(len_mu_guess)
        
        #***t-test***
        for idx_peak in range(len_mu_guess):
            H_opt[idx_peak] = peak_optimised.params['A'+str(idx_peak)].value
            mu_opt[idx_peak] = peak_optimised.params['mu'+str(idx_peak)].value
            sigma_opt[idx_peak] = peak_optimised.params['sigma'+str(idx_peak)].value
            alpha_opt[idx_peak] = peak_optimised.params['alpha'+str(idx_peak)].value
            #calculate fit of specific peak:
            y_opt[idx_peak] = f_chromato(x, H_opt[idx_peak], mu_opt[idx_peak], sigma_opt[idx_peak], alpha_opt[idx_peak]) #+ baseline_opt #add baseline only once!

            
        #extract all y_opt before doing the statistical tests!
        #print('len mu_guess: ' + str(len(mu_guess)))
        if len_mu_guess > 0:
            y_opt_sum_no_baseline = np.sum(y_opt, axis = 0) #+ baseline_opt
            y_opt_sum = y_opt_sum_no_baseline + baseline_opt
            
            #baseline variance weighted by inverse of y distance of fit to baseline:
            #make sure weights are non-zeros:
            b_w = 1/ (np.where(y_opt_sum_no_baseline>1, y_opt_sum_no_baseline, 1))**2
            y_opt_baseline_u = np.sqrt(np.sum(b_w*(y- y_opt_sum)**2)/np.sum(b_w)*len(y_opt_sum)/(len(y_opt_sum)-1))
            #print('y_opt_baseline_u: ' + str(y_opt_baseline_u))
            t_test_baseline = Student_t.ppf(conf_int, len(y)-1) * y_opt_baseline_u
            

    
        for idx_peak in range(len_mu_guess):        
            #here add calculation of std between fit and (measured data minus all other fit)
            #we remove from the measured data the modelled signal due to all other peaks:
            extracted_y = y - y_opt_sum + y_opt[idx_peak] # baseline_opt is part of y_opt_sum
            #indexes of x around mu_opt +/- broadness:
            #mu_opt_x_idx = np.flatnonzero(abs(x - mu_opt[idx_peak]) <= (broadness * sigma_opt[idx_peak]))
            #enough points around the apex?
            mu_opt_x_idx = np.flatnonzero(abs(x - mu_opt[idx_peak]) <= broadness * sigma_opt[idx_peak])
            np_pts[idx_peak] = len(mu_opt_x_idx)            
            np_pts_left[idx_peak] = len(np.flatnonzero((mu_opt[idx_peak] - x <= broadness * sigma_opt[idx_peak]) & (mu_opt[idx_peak] - x >= 0)))
            np_pts_right[idx_peak] = len(np.flatnonzero((x - mu_opt[idx_peak] <= broadness * sigma_opt[idx_peak]) & (x - mu_opt[idx_peak] >= 0)))
            
            
            if mu_opt[idx_peak] - rt_window_min > 0: #mu_opt is not outside left of the scan window
                np_pts_left_min[idx_peak]  = int(round(max(1, min(broadness * sigma_opt[idx_peak], mu_opt[idx_peak] - rt_window_min)/n_per_bin * 0.4)))
            else: np_pts_left_min[idx_peak] = 0
            if rt_window_max - mu_opt[idx_peak] > 0: #mu_opt is not outside right of the scan window
                np_pts_right_min[idx_peak] = int(round(max(1, min(broadness * sigma_opt[idx_peak], rt_window_max - mu_opt[idx_peak])/n_per_bin * 0.4)))
            else: np_pts_right_min[idx_peak] = 0

            sum_over_idx_x = sum(y_opt[idx_peak][idx_x] for idx_x in mu_opt_x_idx)
            if sum_over_idx_x > 0 and np_pts_left[idx_peak] >= np_pts_left_min[idx_peak] and np_pts_right[idx_peak] >= np_pts_right_min[idx_peak]:
                #calculate standard deviation of fit to y, weighted by y value:
                #print('test' + str(sum(y_opt[idx_peak][idx_x] for idx_x in mu_opt_x_idx)))
                
                apex_opt_u[idx_peak] = math.sqrt(sum([(y_opt[idx_peak][idx_x] - extracted_y[idx_x])**2 * y_opt[idx_peak][idx_x] for idx_x in mu_opt_x_idx]) \
                          / sum(y_opt[idx_peak][idx_x] for idx_x in mu_opt_x_idx) ) #/ (len(mu_opt_x_idx)-1)
            else: apex_opt_u[idx_peak] = H_opt[idx_peak] #peak of just 1 points, to eliminate
            t_test_value[idx_peak] = Student_t.ppf(conf_int, max(np_pts[idx_peak]-1, 1)) #added min value of 1 to avoid raised error if nan
            t_test_apex[idx_peak] = t_test_value[idx_peak] * apex_opt_u[idx_peak]
            

        if start == 0 and graph == True:
            #print start parameters just once
            start = 1
            """
            print('****************************')
            #print('Guessed parameters: ')
            #params.pretty_print()
            #print('****************************')
            print('Baseline guess: ' + '{:.1f}'.format(baseline_opt) + '; u: ' + '{:.1f}'.format(t_test_baseline))
            print('mu'+ '\t'  + 'pt l'+ '\t'  + 'pt r' + '\t'+ 'min l' + '\t'+ 'min r' + '\t' + 't_H' + '\t' +'H_min' + '\t' + 'H_opt' + '\t' +'std') #+ '\t' +'std f'
            for idx_peak in range(len_mu_guess):
                print('{:.1f}'.format(mu_opt[idx_peak])  +'\t'+ 
                      str(np_pts_left[idx_peak]) +'\t'+
                      str(np_pts_right[idx_peak]) +'\t'+
                      str(np_pts_left_min[idx_peak]) + '\t' +
                      str(np_pts_right_min[idx_peak]) + '\t' +
                      '{:.2f}'.format(t_test_value[idx_peak])+  '\t' +
                      '{:.1f}'.format(max(t_test_apex[idx_peak] + t_test_baseline, t_test_baseline*2)) +  '\t' +
                      '{:.1f}'.format(H_opt[idx_peak]) + '\t' +
                      '{:.1f}'.format(apex_opt_u[idx_peak]) + '\t'
                      )
            #print('y_opt_baseline_u: ' + str(y_opt_baseline_u))
            #print('t_test_baseline: ' + '{:.2f}'.format(t_test_baseline))
            """
            
            #print('****************************')
        #here remove peaks where height - std dev > threshold
        #And update guess values to be used for next loop

        H_test = np.flatnonzero((H_opt - t_test_apex - t_test_baseline >= 0) &
                                   (H_opt - 2*t_test_baseline >= 0)) 
        
        H_guess = H_opt[H_test]
        count_delete_sum = len_mu_guess - len(H_guess)
        
        mu_guess = mu_opt[H_test]
        sigma_guess = sigma_opt[H_test]
        alpha_guess = alpha_opt[H_test]
        baseline_guess = baseline_opt
        
        len_mu_guess = len(mu_guess)
        
        sigma_guess_min = sigma_guess_min[H_test]
        sigma_guess_max = sigma_guess_max[H_test]
        
    #***The optimisation is finished***

    #In case peaks have been all eliminated at the previous loop, 
    #all values needs to be emptied:
    if len_mu_guess == 0:
        mu_opt = []
        mu_opt_stderr = []
        area_opt = []
        apex_opt_u = []
        sigma_opt = []
        alpha_opt = []
        LOD = []
        SN_opt = []
        tof_centre_opt = []
        tof_centre_opt_u = []
        mass_centre_opt = []
        mass_centre_opt_u = []
    
    else:
        mu_opt_stderr = np.zeros(len_mu_guess)
        area_opt = np.zeros(len_mu_guess)
        SN_opt = np.zeros(len_mu_guess)
        tof_centre_opt = np.zeros(len_mu_guess)
        tof_centre_opt_u = np.zeros(len_mu_guess)
        mass_centre_opt = np.zeros(len_mu_guess)
        mass_centre_opt_u = np.zeros(len_mu_guess)
        
        min_pts_for_LOD = np.min(sigma_guess*3/n_per_bin)
        LOD_H = np.zeros(len_mu_guess)
        LOD = np.zeros(len_mu_guess)
            
        

    for idx_peak in range(len_mu_guess):
        mu_opt_stderr[idx_peak] = peak_optimised.params['mu'+str(idx_peak)].stderr
        area_opt[idx_peak] = f_chromato_area(H_opt[idx_peak], sigma_opt[idx_peak], alpha_opt[idx_peak])
        #we need to add a factor 2 here to get the position of the minimum peak
        #with the left tailing touching the std dev of y.
        #because we want 5% (or so) cases of taking noise as a real peak.
        LOD_H[idx_peak] = max(1, 2.0 * Student_t.ppf(conf_int, min_pts_for_LOD-1) * y_opt_baseline_u)
        SN_opt[idx_peak] = H_opt[idx_peak]/LOD_H[idx_peak]
        
        #now calculate LOD as equivalent area with same peak shape as integrated peak
        LOD[idx_peak] = f_chromato_area(LOD_H[idx_peak], sigma_opt[idx_peak], alpha_opt[idx_peak])
        
        #calculate centre of tof and centre of mass for each individual peak:
        #weight values using y_opt calculated for each peak
        y_opt_peak_sum = np.sum(y_opt[idx_peak])
        tof_centre_opt[idx_peak] = np.sum(tof_series * y_opt[idx_peak]) / y_opt_peak_sum
        tof_centre_opt_u[idx_peak] = np.sqrt(np.sum((tof_series - tof_centre_opt[idx_peak])**2 * y_opt[idx_peak]) / (y_opt_peak_sum * (len_x-1)/len_x))
        mass_centre_opt[idx_peak] = np.sum(mass_series * y_opt[idx_peak]) / y_opt_peak_sum
        mass_centre_opt_u[idx_peak] = np.sqrt(np.sum((mass_series - mass_centre_opt[idx_peak])**2 * y_opt[idx_peak]) / (y_opt_peak_sum * (len_x-1)/len_x))

        #calculation from Daymond 2002:
        #this assumes Gaussian shape
        #mu_opt_stderr = var / area * (1 + 2* sqrt(2) * baseline / height)
        #mu_opt_stderr[idx_peak] = sigma_opt[idx_peak]**2/A_opt[idx_peak]

    plot_linear_fit = False
    if len_mu_guess>0:

        if graph == True:
            #prepare to plot general graph with fit:
            """
            print('Optimised parameters:')
            #peak_optimised.params.pretty_print()
            #print('****************************')
            print('mu'+ '\t'  + 'nb pt' + '\t' + 't_a' + '\t' +'a_min' + '\t' +'Area' + '\t' +'std')
            for idx_peak in range(len_mu_guess):
                print('{:.1f}'.format(mu_opt[idx_peak])  +'\t'+ 
                      str(np_pts[idx_peak]) + '\t' +
                      '{:.2f}'.format(t_test_value[idx_peak])+ '\t' +
                      '{:.1f}'.format(max(t_test_apex[idx_peak] + t_test_baseline, t_test_baseline*2)) +  '\t' + 
                      '{:.1f}'.format(area_opt[idx_peak]) + '\t' + 
                      '{:.1f}'.format(apex_opt_u[idx_peak]) + '\t' #+ '{:.3f}'.format(area_opt_u_frac[idx_peak])
                      )
            #print('t_test_baseline: ' + '{:.2f}'.format(t_test_baseline))
            #print(np.mean(peak_optimised.residual))
            #print('****************************')
            """
            
            ax2 = ax1.twinx()
            ax2.plot(x, peak_optimised.residual, 'c.')
            ax2.set_ylabel('residuals', color='c')
            ax2.tick_params('y', colors='c')
    
            #to plot the fit with good resolution to visualise something:
            x_plot = np.linspace(min(x), max(x), num = int((max(x) - min(x)) /(n_per_bin/4.0))) 
            y_opt_sum_plot = np.zeros(len(x_plot))
            for idx_peak in range(len_mu_guess):
                y_opt_plot=f_chromato(x_plot, H_opt[idx_peak], mu_opt[idx_peak], sigma_opt[idx_peak], alpha_opt[idx_peak]) #+ baseline_opt #add baseline only once!
                ax1.plot(x_plot, y_opt_plot, 'g--')
                y_opt_sum_plot += y_opt_plot
            ax1.plot(x_plot, y_opt_sum_plot + baseline_opt, 'k-')
            ax1.plot(mu_opt, H_opt, '*', color='xkcd:deep red', markersize = 11)
            ax1.plot([min(x), max(x)], [baseline_opt, baseline_opt], '--', color = 'xkcd:slate')


        #if len(mu_opt)>0:
        #fit a linear function, to then compare fit with chromato, to test presence of real peak
        #fit with linear function y = a* x + b
        linear_fit_slope, linear_fit_b, linear_fit, std_dev_linear_fit = fit_linear_np(x, y)

        #or better calculate standard deviation and divide by total area
        std_dev_f = np.sqrt(np.sum((y_opt_sum - y)**2) / float(len(y)-1))
        std_dev_f_rel= std_dev_f /np.sum(area_opt)
        #noise_signal = sum([abs((final[i]- y[i]) / (final[i]- baseline_opt))  for i in range(len(y))])/len(y)
        
        if graph == True:
            #print('Optimised areas:')
            #for area in area_opt:
            #    print('{:.2f}'.format(area))

            #print('Std dev of linear fit: ' + '{:.4f}'.format(std_dev_linear_fit))
            print('Std dev of fit: ' + '{:.4f}'.format(std_dev_f) + '; rel. to area: ' + '{:.4f}'.format(std_dev_f_rel))
            print('std dev peak fit / srd dev linear: ' + '{:.4f}'.format(std_dev_f / std_dev_linear_fit ) )

        #test if Gaussian peak really improve the fit compared to linear fit
        #if not, this is not a peak, just noise (typically for N2, O2)
        #in case of real peak std_dev_Gauss/std_dev_linear_fit < 0.5
        
        if std_dev_f/std_dev_linear_fit > ratio_to_linear_fit:
            
            plot_linear_fit = True
            mu_opt = []
            mu_opt_stderr = []
            area_opt = []
            apex_opt_u = []
            sigma_opt = []
            alpha_opt = []
            LOD = []
            SN_opt = []
            tof_centre_opt = []
            tof_centre_opt_u = []
            mass_centre_opt = []
            mass_centre_opt_u = []
            if graph == True:
                ax1.plot(x, linear_fit, 'r--')


            

        

        
    if graph == True:
        #plt.savefig('Gauss_fit_on_mass_smallA.png', bbox_inches='tight', transparent = True, dpi = 300)
        #fig.tight_layout()
        if len_mu_guess > 0 or plot_linear_fit:
            plt.show()
        else:
            plt.close()
    
    return mu_opt, mu_opt_stderr, area_opt, apex_opt_u, sigma_opt, alpha_opt, LOD, SN_opt, tof_centre_opt, tof_centre_opt_u, mass_centre_opt, mass_centre_opt_u

#*************************************************************************
    
def f_Gauss(x, A, mu, sigma):
    """
    Calculate Gaussian function of known parameters A, mu, sigma, for a given x
    x: np array containing floats: abcisses where we want the Gaussian value
    A: parameter of the Gaussian function, float
    mu: parameter of the Gaussian function (centre of peak), float
    sigma: parameter of the Gaussian function, float
    """
    #do not include baseline here otherwise it will be added 
    #each time the gaussian function is summed up!
    return A / (sigma * np.sqrt(2 * np.pi)) * np.exp(-(x-mu)**2/(2*sigma)**2)
        
def cost_f_Gauss(params, x, y):
    """
    Calculate the cost function to minimise
    This cost function is adapted to the number of expected peaks
    It returns the quantity to be minimised using later the minimise function from limfit
    """
    #if model == 'Gaussian':
    nb_Gauss_params = int(3) #3 parameters to define the Gauss function
    model=params['baseline'] #+str(idx_peak)
    #this cost function should adapt itself to number of found peaks
    for idx_peak in range(int((len(params)-1)/nb_Gauss_params)): #need to be integer value; -1 to remove 1 line for baseline
        mu = params['mu'+str(idx_peak)]
        sigma = params['sigma'+str(idx_peak)]
        A = params['A'+str(idx_peak)]
        model += f_Gauss(x, A, mu, sigma)
    
    #using a weight = math.sqrt(data) to give more weight to detection of large peaks:
    #diff_model_data = model - data
    #if there are more points on the side of the peak 
    #as on the peak, multiplying by y may heuu no
    
    #Cost functions that do not work:
    #-------------------------
    #cost = (model - y) * y -> underweight points at baseline
    #baseline is not captured correctly -> area not correct
    #-------------------------
    #definition of an uncertaint value
    #uncertainty at baseline should represent detection noise
    #but not more
    weight_list = []
    U_max = 1.2
    U_min = 1
    for val in y:
        #weight = 1.0
        if val != 0:
            #we want U to increase towards baseline value but not above a threshold U of 1:
            U = min(U_max, abs(params['baseline']/val))
            #however we don't want too small U values when y is large (will disturb baseline estimate)
            if U < U_min: U = U_min
        else: U = U_max
        weight_list.append(1/U)
    #elif model == 'tofidx_to_mass':
        
        
    cost = (model - y) #* weight_list
    return cost #**2 / len(x) 

def fit_Gauss_make_params(x, y, mu_guess, apex_guess, baseline_guess, n_per_bin):
    #make parameters for Gauss fit
    #x_mean = np.sum(x)/len(x)
    
    #y_mean=np.sum(y)/len(y)
    params = Parameters()
    #params.add('baseline', value=max(1, baseline_guess), min = max(0, min(y)-3*min(y)), max = max(3, baseline_guess*3, min(y)+y_mean/10)) #
    params.add('baseline', value=max(1, baseline_guess), min = max(0, baseline_guess/3), max = max(3, baseline_guess*3)) #
    #baseline cannot be above max I
    #print('baseline guess: ' + str(baseline_guess))
    
    mu_guess_min = mu_guess-n_per_bin*3
    mu_guess_max = mu_guess+n_per_bin*3
    
    #sigma_guess = max(0.5, 1.0 + 0.00357488*(mu_guess-1900))
    sigma_guess = 0.0416569*np.exp(0.0017*mu_guess)
    sigma_guess_min = sigma_guess*float(0.7)
    sigma_guess_max = sigma_guess+1.5
    
    A_guess = (apex_guess-baseline_guess) * sigma_guess * 2.0
    A_guess_min = A_guess/2.0
    A_guess_max = A_guess*10
            
    for idx_peak in range(len(mu_guess)):
        #parameters per peak:
        #len(params) is a multiple of 4
        params.add('mu'+str(idx_peak), value = mu_guess[idx_peak], min = mu_guess_min[idx_peak], max = mu_guess_max[idx_peak]) #abcisse cannot be negative
        #params.add('sigma'+str(idx_peak), value = peak_value_final[idx_peak]/2, min = 0.0) #assumed proportionality between max intensity and sigma
        #params.add('sigma'+str(idx_peak), value = 0.65+1/48000*mu_guess[idx_peak], min = 0.65, max = (0.65+1/48000*mu_guess[idx_peak])*2)
        #print(str(mu_guess[idx_peak]) + ', min: ' + str(mu_guess[idx_peak]-n_per_bin*3) + ', max: ' + str(mu_guess[idx_peak]+n_per_bin*3))
        #sigma_guess = 0.6
        #params.add('sigma'+str(idx_peak), value = sigma_guess, min = 0.3, max = 1.0)

        params.add('sigma'+str(idx_peak), value = sigma_guess[idx_peak], min = sigma_guess_min[idx_peak], max = sigma_guess_max[idx_peak])
        
        #A_guess = apex_guess[idx_peak]/(sigma_guess*math.sqrt(2*math.pi))
        #A_guess = (apex_guess[idx_peak]-baseline_guess) * sigma_guess * 2.0
        params.add('A'+str(idx_peak), value = max(2, A_guess[idx_peak]) , min = max(1, A_guess_min[idx_peak]), max = max(3, A_guess_max[idx_peak]))
        #params.add('A'+str(idx_peak), value = max(2, A_guess[idx_peak]) , min = max(1, A_guess/2.0), max = max(3, A_guess*10)) #, max = A_guess*3.0
        #params.add('A'+str(idx_peak), value = apex_guess[idx_peak]/(PW*math.sqrt(2*math.pi)) , min = 0.0)
        #print(str(max(2, A_guess)) + ', min: ' + str(max(1, A_guess/2.0)) + ', max: ' + str(max(3, A_guess*10)))
       
    #print('Guessed parameters: ')
    #params.pretty_print()
    #fit, here with leastsq model
    minner = Minimizer(cost_f_Gauss, params, fcn_args=(x, y)) #, method='leastsq'
    peak_optimised = minner.minimize(method = 'leastsq')
    
    #resulting optimised parameters are callable with:
    #peak_optimised.params
    return params, peak_optimised

def fit_Gauss(x, y, mu_guess, apex_guess, baseline_guess_0, n_per_bin, area_threshold, ratio_to_linear_fit, graph):
    """This function computes a (or several) Gaussian fit(s) of the data. 
    Gaussian is parameterized by mu, sigma, A
    :param x: list of floating point numbers, x-coordinates
    :param y: list of floating point numbers, y-coordinates
    :param mu_guess: list of guess of mu parameter (mean) for Gaussian, float
    :param apex_guess: list of value of y when x=mu, float, max of Gaussian (at peak)
    :param n_per_bin: number of points x in between two integer values (sampling rate for x)
    :param area_threshold: minimum value of area (of Gaussian), float
    :param sigma_threshold: minimum value of sigma, float
    
    returns 
    """

    #assuming same baseline for all peaks in a small window:
    #x should cover a mass range large enough around unit mass (or same in tof_idx)
    #so that minimum can be found where no mass is expected
    #i.e. unit mass +/- 0.3 m/z or more
    #baseline_guess = min(y)
    #y_sort = [y[i] for i in range(len(y))]
    #y_sort.sort()
    #baseline_guess = np.median([y_sort[i] for i in range(int(len(y_sort)/3))]) 
    baseline_guess = max(0, min(y)*0.8, baseline_guess_0)
    #we have now points that corresponds to signal above baseline
    #so the minimum is the closest to the real baseline

    params, peak_optimised = fit_Gauss_make_params(x, y, mu_guess, apex_guess, baseline_guess, n_per_bin)
    #peak_optimised_params = peak_optimised.params

    if graph == True:
        print('Guessed parameters: ')
        params.pretty_print()
        fig, ax1 = plt.subplots(figsize=(8,4))
        fig.suptitle('Gauss fit over RT', fontsize=12)
        ax1.set_xlabel('Retention time [s]')
        ax1.set_ylabel('Intensity [TOF-index * V]', color='b')
        ax1.tick_params('y', colors='b')
        ax1.plot(x, y, 'b.')
        ax1.plot(mu_guess, apex_guess, 'o', color='xkcd:peach')
   
    
    A_opt = np.zeros(len(mu_guess))
    for idx_peak in range(len(mu_guess)):
        A_opt[idx_peak] = peak_optimised.params['A'+str(idx_peak)].value

    count_delete = (np.where(A_opt < area_threshold, 1, 0))
    count_delete_sum = int(np.sum(count_delete))
        
    while len(mu_guess)-count_delete_sum > 0 and count_delete_sum > 0:
        mu_guess = mu_guess[A_opt > area_threshold]
        apex_guess = apex_guess[A_opt > area_threshold]

        params, peak_optimised=fit_Gauss_make_params(x, y, mu_guess, apex_guess, baseline_guess, n_per_bin)
        peak_optimised.params = peak_optimised.params
        A_opt = np.zeros(len(mu_guess))
        for idx_peak in range(len(mu_guess)):
            A_opt[idx_peak] = peak_optimised.params['A'+str(idx_peak)].value
        count_delete = np.where(A_opt < area_threshold, 1, 0)
        count_delete_sum = int(np.sum(count_delete))


    """
    count_delete =0
    for idx_peak in range(len(mu_guess)-1, -1, -1):
        #delete too small peaks (may be just noise)
        #delete summits outside of RT window (not sure, what if we are on the tail of another peak?)
        if peak_optimised.params['A'+str(idx_peak)].value < area_threshold:
        #or peak_optimised.params['mu'+str(idx_peak)].value < min(x) \
        #or peak_optimised.params['mu'+str(idx_peak)].value > max(x):
            count_delete += 1
            mu_guess.remove(mu_guess[idx_peak])
            apex_guess.remove(apex_guess[idx_peak])
    while len(mu_guess)>0 and count_delete > 0:
        params, peak_optimised=fit_Gauss_make_params(x, y, mu_guess, apex_guess, baseline_guess, n_per_bin)
        peak_optimised_params = peak_optimised.params
        count_delete =0
        for idx_peak in range(len(mu_guess)-1, -1, -1):
            if peak_optimised.params['A'+str(idx_peak)].value < area_threshold:
            #or peak_optimised.params['mu'+str(idx_peak)].value < min(x) \
            #or peak_optimised.params['mu'+str(idx_peak)].value > max(x):     
                count_delete += 1
                mu_guess.remove(mu_guess[idx_peak])
                apex_guess.remove(apex_guess[idx_peak])
    """
 
    
    #fit a linear function, to then compare fit with Gauss, to test presence of real peak
    #fit with linear function y = a* x + b
    linear_fit_slope, linear_fit_b, linear_fit, std_dev_linear_fit = fit_linear_np(x, y)

    mu_opt = []
    mu_opt_stderr = []
    A_opt= []
    sigma_opt = []
    baseline_opt = []
    max_opt = []
    
    baseline_opt.append(peak_optimised.params['baseline'].value)
    final = baseline_opt


    if graph == True:
        print('Optimised parameters:')
        peak_optimised.params.pretty_print()
        ax2 = ax1.twinx()
        ax2.plot(x, peak_optimised.residual, 'c.')
        ax2.set_ylabel('residuals', color='c')
        ax2.tick_params('y', colors='c')

        #to plot the fit with good resolution to visualise something:
        x_plot = []
        x_plot.append(min(x))
        x_plot_step = min(x)
        step = n_per_bin/4.0
        while x_plot_step <= max(x):
            x_plot.append(x_plot_step+step)
            x_plot_step += step
        final_plot = baseline_opt

    for idx_peak in range(len(mu_guess)- count_delete_sum):
        A_opt.append(peak_optimised.params['A'+str(idx_peak)].value)
        mu_opt.append(peak_optimised.params['mu'+str(idx_peak)].value)
        mu_opt_stderr.append(peak_optimised.params['mu'+str(idx_peak)].stderr)
        sigma_opt.append(peak_optimised.params['sigma'+str(idx_peak)].value)
        max_opt.append(A_opt[idx_peak]/(sigma_opt[idx_peak]*math.sqrt(2*math.pi))+baseline_opt) #when x = mu
        y_opt = f_Gauss(x, A_opt[idx_peak], mu_opt[idx_peak], sigma_opt[idx_peak]) + baseline_opt #add baseline only once!
        final +=y_opt - baseline_opt
        if graph == True:
            y_opt_plot=f_Gauss(x_plot, A_opt[idx_peak], mu_opt[idx_peak], sigma_opt[idx_peak]) + baseline_opt #add baseline only once!
            ax1.plot(x_plot, y_opt_plot, 'g--')
            final_plot += y_opt_plot - baseline_opt
    
    if len(mu_opt)>0 and graph == True: 
        ax1.plot(x_plot, final_plot, 'k-')
        ax1.plot(mu_opt,max_opt, '*', color='xkcd:deep red')
    
    #if count_delete>0:
    if len(mu_opt)>0:
        #or better calculate standard deviation and divide by total area
        std_dev_Gauss = np.sqrt(np.sum((final- y)**2) / float(len(y)-1))
        std_dev_Gauss_rel= std_dev_Gauss /np.sum(A_opt)
        #noise_signal = sum([abs((final[i]- y[i]) / (final[i]- baseline_opt))  for i in range(len(y))])/len(y)
        
        if graph == True:
            #print('Std dev of linear fit: ' + '{:.4f}'.format(std_dev_linear_fit))
            print('Std dev of fit: ' + '{:.4f}'.format(std_dev_Gauss) + '; rel. to area: ' + '{:.4f}'.format(std_dev_Gauss_rel))
            print('std dev Gauss / srd dev linear: ' + '{:.4f}'.format(std_dev_Gauss / std_dev_linear_fit ) )

        #test if Gaussian peak really improve the fit compared to linear fit
        #if not, this is not a peak, just noise (typically for N2, O2)
        #in case of real peak std_dev_Gauss/std_dev_linear_fit < 0.5
        if std_dev_Gauss/std_dev_linear_fit > ratio_to_linear_fit: 
            mu_opt = []
            mu_opt_stderr = []
            A_opt= []
            sigma_opt = []
            baseline_opt = []
            max_opt = []
            final = []
            if graph == True:
                ax1.plot(x, linear_fit, 'r--')
    if graph == True:
        #plt.savefig('Gauss_fit_on_mass_smallA.png', bbox_inches='tight', transparent = True, dpi = 300)
        fig.tight_layout()
        plt.show()
    
    return mu_opt, mu_opt_stderr, max_opt, A_opt, sigma_opt, baseline_opt, final, peak_optimised.residual

#******************************************************************************

def find_match(array_in, value, idx_value_in_tuple):
    """
    Aurore Guillevic 20190721
    ex: idx_value_in_tuple = 1 to look for mass in mass_cal_list
    Returns an array with all the array_in[i][j] if
    array_in[i][j][idx_in_tuple] == value
    if no match, return None
    """
    ans = [None]*len(array_in)
    for i in range(len(array_in)):
        j=0
        while j<len(array_in[i]) and array_in[i][j][idx_value_in_tuple] != value:
            j+=1
        if j<len(array_in[i]): #we have found a match because we went out of the while loop
            ans[i] = array_in[i][j]
    #Some ans[i] may be None
    return ans

def find_match_without_none(array_in, value, idx_value_in_tuple):
    """
    Aurore Guillevic 20190721
    ex: idx_value_in_tuple = 1 to look for mass in mass_cal_list
    Returns an array with all the array_in[i][j] if
    array_in[i][j][idx_in_tuple] == value
    if no match, return None
    """
    ans = []*len(array_in)
    i1=0
    for i in range(len(array_in)):
        j=0
        while j<len(array_in[i]) and array_in[i][j][idx_value_in_tuple] != value:
            j+=1
        if j<len(array_in[i]): #we have found a match because we went out of the while loop
            ans[i1] = array_in[i][j]
            i1 += 1
    #at this stage we may have empty cells at the end, we want to remove them from ans
    if i1<len(array_in):
        ans = ans[:i1]
    return ans


def make_histogram(results_list, delta_tofidx_discrete, h_index):
    """Make histograms of data.
    """      
        
    # sort results according to the center of mass peak (second item in tuple)
    results_list.sort(key=operator.itemgetter(h_index))        
            
    # histogram of nb of peak per TOF-index (or mass) bin

    #create empty histogram of right size, add +5 for slope detection algo later
    #histo_mass = [0 for i in range( int((results_list[len(results_list)-1][1] - results_list[0][1])/delta_mass_discrete) +20)]
    histo = [0 for i in range( int((results_list[len(results_list)-1][h_index] - results_list[0][h_index])/delta_tofidx_discrete) +20)]
    i = 0
    j = 10 #artificially leave 10 first values as zero to allow slope detection algo later on
    histo_x_a = results_list[0][h_index]
    #mass_b = mass_a + delta_mass_discrete
    histo_x_b = histo_x_a + delta_tofidx_discrete
    while i < len(results_list): # and j < len(histo_mass)
        while i < len(results_list) and results_list[i][h_index] < histo_x_b:
            histo[j] += 1
            i += 1
        j += 1
        histo_x_a = histo_x_b
        #mass_b += delta_mass_discrete
        histo_x_b += delta_tofidx_discrete
    
    #below we have i-10 because above we have j=10; we start filling the histo at j=10
    histo_x_coord = [results_list[0][h_index] + delta_tofidx_discrete*(i-10) for i in range(len(histo))]
    
    return histo_x_coord, histo


def find_bins_np(x, delta_bin, PWh, PTh, h_nb_per_bin, graph, mode):
    
    #*** Make histogram***
    # sort results according to the center of mass peak (second item in tuple)
    #results_list.sort(key=operator.itemgetter(h_index))
    
    #np.sort(a = results_list_np, axis = -1, order = results_list_np[h_index]) #work only with structured numpy array
    
    #results_list_np =  results_list_np[results_list_np[:,1].argsort()] #sort by 2nd column, that's tofidx value
    hist_len = max(2, int((x[-1]+ PWh - x[0] + PWh +1) / delta_bin))
    #print(hist_len)
    
    #create bin edges with extra data before and after, to make sure 
    #slope is correctly calculated at edges
    hist_bin_edges = np.linspace(x[0]- PWh, x[-1]+ PWh, hist_len)
    #print(hist_bin_edges)
    
    #numpy.histogram(a, bins=10, range=None, normed=None, weights=None, density=None)[source]
    np_hist, np_hist_bin_edges = np.histogram(x[:], bins = hist_bin_edges)
    #print(np_hist_bin_edges)
    np_hist_bin_centres = np_hist_bin_edges[0:len(np_hist_bin_edges)-1] + (np_hist_bin_edges[1] - np_hist_bin_edges[0])/2
    # histogram of nb of peak per TOF-index (or mass) bin


    #we now have an histogram of nb of peaks per mass bin
    #we apply again the slope detection algo
    histo_bin_centres = peak_detect(np_hist_bin_centres, np_hist, 
                                 PWh, PTh, h_nb_per_bin, 
                                 graph = graph, mode = mode)
    
    if len(histo_bin_centres) >0:
        bin_centres = histo_bin_centres[:,0]
    else: bin_centres = []

    return bin_centres

def make_bins_np(data_x, data_y, histo_centre_list, std_dev_p_max, nb_bin_min, h_index, h_index_x, s_per_bin, graph, mode, plot_mode):
    """
    Distribute data into the identified groups (histo_centre_list).
    Filter: eliminate data too far away from the centre.
    These data, not assigned to a centre, are then assigned an idx_bin of -1.
    
    Inputs:
    results_list_np: np.array
    histo_centre_list: list of centre of bins of points, np.array
    std_dev_p_max: float, maximum standard deviation of position of points within a bin
    nb_bin_min: float, minimum number of points per bin. If smaller, bin index is set to -1 (unused).
    h_index: integer, dimention of results_list_np to use to make the histogram
    
    Outputs:
    bin_value_list: Python list, has sise of data_x, contain the corresponding bin value.
    no_bins_final: number of different valid bins (not -1).
    bin_list: list of index values for all various bins.
    """
    #idx_bin = 6 #where is stored the index of the bin
    len_x = len(data_x)
    no_bins = len(histo_centre_list) #how many bins do we have
    no_bins_final = no_bins
    bin_list = []
    bin_value_list = np.zeros(len_x, dtype = int)
    #***Prepare graph***
    if graph == True:
        fig, ax1 = plt.subplots(figsize=(10,5))
        ax1.ticklabel_format(axis = 'both', style = 'plain', useMathText=False)
        ax1.grid(True, color='xkcd:grey', linestyle='--', linewidth=0.5)

        if plot_mode == 'rt_domain': #group fragments by co-elution
            fig.suptitle('Fragments grouped by co-elution', fontsize=14)
            plt.xlabel('Retention time [s]', fontsize=14)
            plt.ylabel('Mass of fragment [m/z]', fontsize=14)
            #rt_mean = np.mean(results_list_np[:, 0])
            ax1.plot(data_y, data_x, '.', color = 'xkcd:slate')
        
        elif plot_mode == 'mass_domain': #mass calibration: group by tof-index
            fig.suptitle('TOF-indexes grouped by same fragment', fontsize=12)
            plt.xlabel('Retention time [s]')
            plt.ylabel('Centre of TOF peak [TOF index]')
            ax1.plot(data_x*s_per_bin, data_y, '.', color = 'xkcd:slate')  
    
    #***Assign bin index according to closest centre of bins***
    #20200611
    #for each bin, eliminate already all points that are away of the centre by more than n * std_dev_p_max
    
    if no_bins ==1: #one identified peak only, makes only one bucket
        #bin_value_list[:] = int(0) #assign bin index to bin index column
        bin_value_list[:] = np.where(abs(data_y - histo_centre_list) <= std_dev_p_max * 4.0, int(0), int(-1))
            
    elif no_bins > 1:
        i = 0
        for histo_idx_peak in range(no_bins-1):
            histo_peak_0 = histo_centre_list[histo_idx_peak] # float numbers, this one is the first peak
            histo_peak_1 = histo_centre_list[histo_idx_peak+1] # float numbers, this one is the second peak
            i0 = i
            while i < len_x and abs(data_y[i] - histo_peak_0) < abs(data_y[i] - histo_peak_1):
                i += 1
            #bin_value_list[i0:i] = histo_idx_peak #assign idx_bin
            bin_value_list[i0:i] = np.where(abs(data_y[i0:i] - histo_peak_0) <= std_dev_p_max * 4.0, histo_idx_peak, -1)
        if i < len_x: #results after the last peak centre
            bin_value_list[i:] = histo_idx_peak + 1
            bin_value_list[i:] = np.where(abs(data_y[i:] - histo_peak_1) <= std_dev_p_max * 4.0, histo_idx_peak+1, -1)
    #***********************************************************
    #print(bin_value_list)
    #now we want to eliminate the points that are too far from the assigned centres
    #those points will be assigned an idx_bin of -1
    
        
    if mode == 'no_slope':
        no_peak_min = max(nb_bin_min, 1)
    elif mode == 'with_slope':
        no_peak_min = max(nb_bin_min, 2) #need two points at least for a linear fit

    for i in range(no_bins):

        count_delete = 1
        limit = std_dev_p_max
        bin_indexes = np.flatnonzero(bin_value_list == i)
        #print('len bin indexes: ' + str(len(bin_indexes)))
        indexes = bin_indexes
        y = data_y[bin_indexes]
        no_pts = len(y)
                
        while no_pts >= no_peak_min and (count_delete > 0  or limit > std_dev_p_max):
            no_pts_old = no_pts
            
            if mode == 'no_slope':
                linear_fit = np.median(y)
                linear_fit_b = linear_fit
                linear_fit_slope = 0.0
                
            elif mode == 'with_slope':
                x = data_x[indexes]
                linear_fit_slope, linear_fit_b, linear_fit, std_dev_linear_fit = fit_linear_np(x, y)
            
            max_diff = np.max(y - linear_fit)
            limit = max(std_dev_p_max, max_diff)
            #print('max diff: ' + str(max_diff) + '; limit: ' + str(limit))
            
            #assign idx_bin of -1 if (y - linear_fit) > limit:
            #take all data again, in case throuing out a bad point brings ok data back to the line
            #bin_value_list[indexes] = np.where(abs(y - linear_fit) < limit, i, -1)
            bin_value_list[bin_indexes] = np.where(abs(data_y[bin_indexes] - (data_x[bin_indexes]*linear_fit_slope+linear_fit_b)) < limit, i, -1)
            indexes = np.flatnonzero(bin_value_list == i)
            #update y
            y = data_y[indexes]
            no_pts = len(y)
            count_delete = no_pts_old - no_pts
            
            #print('count delete: ' + str(count_delete))
        if no_pts < no_peak_min:
            bin_value_list[bin_value_list == i] = -1
            no_bins_final -= 1
        else: bin_list.append(i)
        
    
    #extract valid data only:
    #bucket_peak_final = results_list_np[results_list_np[:,idx_bin] != -1, :]
    
    if graph == True:
        #plot each bucket with a different color
        if plot_mode == 'rt_domain':
            color_idx = 0
            for i in range(no_bins):
                indexes = np.flatnonzero(bin_value_list == i)
                ax1.plot(data_y[indexes], data_x[indexes], '.', color = str(color1_names[color_idx]))
                color_idx += 1
                if color_idx > len(color1_names)-1: color_idx = 0

        elif plot_mode == 'mass_domain': #if h_index == 1 and h_index_x == 0: #mass calibration: group by tof-index
            color_idx = 0
            for i in range(no_bins):
                indexes = np.flatnonzero(bin_value_list == i)
                ax1.plot(data_x[indexes]*s_per_bin, data_y[indexes], '.', color = str(color1_names[color_idx]))          
                color_idx += 1
                if color_idx > len(color1_names)-1: color_idx = 0

            #ax1.set_xlim(1600, 1800)
            #plt.savefig('grouped_peaks_vs_RT_'+ '{:.2f}'.format(rt_mean) + '.png', bbox_inches='tight', transparent = True, dpi = 300)
        #fig.tight_layout()
        
        plt.show()
    #******
              
    return bin_value_list, no_bins_final, bin_list #results_list_np #bucket_peak_final, bucket_peak_final_single

   
def make_buckets(results_list, histo_centre_list, std_dev_p_max, nb_bin_min, h_index, h_index_x, s_per_bin, graph, mode):
    """
    Distribute data into the identified groups (histo_centre_list).
    Filter: eliminate data too far away from the centre.
    """
    #bucket_peak_rm = []
    if graph == True:
        #normal plot routine
        fig, ax1 = plt.subplots(figsize=(10,5))
        ax1.ticklabel_format(axis = 'both', style = 'plain', useMathText=False)
        ax1.grid(True, color='xkcd:grey', linestyle='--', linewidth=0.5)

        if h_index == 0 and h_index_x == 1: #group fragments by co-elution
            fig.suptitle('Fragments grouped by co-elution', fontsize=14)
            plt.xlabel('Retention time [s]', fontsize=14)
            plt.ylabel('Mass of fragment [m/z]', fontsize=14)
            rt_mean = (min([line[h_index] for line in results_list]) + max([line[h_index] for line in results_list]))/2

            for line in results_list:
                #mass = ((line[h_index_x] +3729.0) / 2828.0)**(1/0.4999)
                ax1.plot(line[h_index], line[6], '.', color = 'xkcd:slate')
        
        elif h_index == 1 and h_index_x == 0: #mass calibration: group by tof-index
            fig.suptitle('TOF-indexes grouped by same fragment', fontsize=12)
            plt.xlabel('Retention time [RT index or s]')
            plt.ylabel('Centre of TOF peak [TOF index]')
            for line in results_list:
                ax1.plot(line[h_index_x], line[h_index], '.', color = 'xkcd:slate')  


        else: #if h_index == 1 and h_index_x == 0:
            fig.suptitle('Masses grouped by same fragment', fontsize=12)
            plt.xlabel('Retention time [RT index or s]')
            plt.ylabel('Centre of mass peak [m/z]')
            for line in results_list:
                ax1.plot(line[h_index_x], line[h_index], '.', color = 'xkcd:slate')  



    #print('make buckets')
    
    if len(histo_centre_list) ==1: #one identified peak only, makes only one bucket
        bucket_peak = []
        bucket_peak.append(results_list)
    
    elif len(histo_centre_list) > 1:
        bucket_peak = []
        i = 0
        for histo_idx_peak in range(len(histo_centre_list)-1):
            histo_peak_0 = histo_centre_list[histo_idx_peak] # float numbers, this one is the first peak
            histo_peak_1 = histo_centre_list[histo_idx_peak+1] # float numbers, this one is the second peak
            i0 = i
            while i < len(results_list) and abs(results_list[i][h_index] - histo_peak_0) < abs(results_list[i][h_index] - histo_peak_1):
                i += 1
            bucket_peak.append(results_list[i0:i])
        if i < len(results_list): #results after the last peak centre
            bucket_peak.append(results_list[i:])

    #now we want to eliminate the points that are too far from the assigned centres
    #the new list bucket_mass_peak_rm contains the masses that are close enough to the assigned centre
    #bucket_peak_rm = bucket_peak.copy()
    if mode == 'no_slope':
        no_peak_min = max(nb_bin_min, 1)
    elif mode == 'with_slope':
        no_peak_min = max(nb_bin_min, 2) #need two points at least for a linear fit

    for i in range(len(bucket_peak)):
        #print('Bucket number: ' + str(i) + '; len of entire bucket: ' + str(len(bucket_peak[i])))
        
        #if len(bucket_peak[i])<max(nb_bin_min, 2): #do not evaluate data
        #    #do not throw away all, 
        #    #when looking at fragment close to LOD, this may be one or two real fragments.
        #    bucket_peak_rm.append([bucket_peak[i][j] for j in range(len(bucket_peak[i]))])

        #elif len(bucket_peak[i])>=max(nb_bin_min, 2): #need two points to calculate linear fit
        #if len(bucket_peak_rm[i]) >= max(nb_bin_min, 2): 
        count_delete = 1
        limit = std_dev_p_max
        len_bucket_peak_i = len(bucket_peak[i])
        #bucket_peak_rm.append([bucket_peak[i][j] for j in range(len(bucket_peak[i]))])
        
        while len_bucket_peak_i >= no_peak_min and (count_delete > 0  or limit > std_dev_p_max):
            
            y = np.array([float(bucket_peak[i][j][h_index]) for j in range(len(bucket_peak[i]))])
            
            if mode == 'no_slope':
                linear_fit_slope = 0
                linear_fit_b = np.median(y)
                std_dev_linear_fit = np.sqrt(np.sum((y - linear_fit_b)**2) / (len(y)-1))
                linear_fit = [linear_fit_b]*len(y)
                
            elif mode == 'with_slope':
                x = np.array([float(bucket_peak[i][j][h_index_x]) for j in range(len(bucket_peak[i]))])
                linear_fit_slope, linear_fit_b, linear_fit, std_dev_linear_fit = fit_linear_np(x, y)
            
            
            max_diff = max([abs(y[j] - linear_fit[j]) for j in range(len(y))])
            #eliminate mass values that are too far away from avrg_p:
            #for idx in range(len(bucket_mass_peak[i]),-1,-1):
            #    if bmp_ij[1] < (avrg_p - 3 *std_dev_p) or bmp_ij[1] > (avrg_p + 3 *std_dev_p):
            #        bucket_mass_peak[i].remove(bmp_ij)
            #        #delete_count += 1
            len_bucket_old = len(bucket_peak[i])
            #new_limit = max(2.0, std_dev_linear_fit)
            limit = max(std_dev_p_max, max_diff)
            #print('max diff: ' + str(max_diff) + '; limit: ' + str(limit))
            #maybe the best woud be to remove the items one by one
            #bucket_peak_rm[i] = [bucket_peak_rm[i][j] for j in range(len(bucket_peak_rm[i])) if abs(bucket_peak_rm[i][j][1] - median_p) < std_dev_p_tofidx ]
            bucket_peak[i] = [bucket_peak[i][j] for j in range(len(bucket_peak[i])) if abs(bucket_peak[i][j][h_index] - linear_fit[j]) < limit ]
            
            count_delete = len_bucket_old - len(bucket_peak[i])
            len_bucket_peak_i = len(bucket_peak[i])
            #print('count delete: ' + str(count_delete))
            #print(len(bucket_mass_peak_rm[i]))
    #each [i] in bucket_mass_peak contains (rt_idx, centre of mass peak) belonging to same mass over an entire run
    #mass_centre_optimised_list = []
    
    bucket_peak_final = []
    #bucket_peak_final_single = []
    for i in range(len(bucket_peak)):
        if len(bucket_peak[i]) >= nb_bin_min:
            bucket_peak_final.append(bucket_peak[i])
        #else:
        #    bucket_peak_final_single.append(bucket_peak_rm[i])
    
    
    if graph == True:
        color_idx = 0
        #plot each bucket with a different color
        if h_index == 0 and h_index_x == 1:
            
            for idx in range(len(bucket_peak_final)):
                #if len(bucket_peak[idx]) >= nb_bin_min:
                for line in bucket_peak_final[idx]:
                    #mass = ((line[h_index_x] +3729.0) / 2828.0)**(1/0.4999)
                    ax1.plot(line[h_index], line[h_index_x], '.', color = str(color1_names[color_idx])) # 'xkcd:blue'
                color_idx += 1
                if color_idx > len(color1_names):
                    color_idx = 0


        else: #if h_index == 1 and h_index_x == 0:
            for idx in range(len(bucket_peak_final)):
                #if len(bucket_peak[idx]) >= nb_bin_min:
                for line in bucket_peak_final[idx]:
                    ax1.plot(line[h_index_x], line[h_index], '.', color = str(color1_names[color_idx]))
                color_idx += 1
                if color_idx > len(color1_names):
                    color_idx = 0
            
            #ax1.set_xlim(1600, 1800)
            #plt.savefig('grouped_peaks_vs_RT_'+ '{:.2f}'.format(rt_mean) + '.png', bbox_inches='tight', transparent = True, dpi = 300)
        #fig.tight_layout()
        plt.show()
        #******
   
    #at this stage we still may have some bucket_peak_rm with less than 3 peaks
            
    return bucket_peak_final #bucket_peak_final, bucket_peak_final_single

def split_into_windows(rt_idx_start_bin, rt_idx_stop_bin, dict_k_rt_window, dict_k_rt_window_edge, s_per_bin, len_time_axis):
    """
    
    Parameters
    ----------
    rt_idx_start_bin : int
        Index of start time.
    rt_idx_stop_bin : int
        Index of stop time.
    dict_k_rt_window : float
        maximum duration of a window.
    dict_k_rt_window_edge : float
        duration of an 'edge' window to add before and after when spliting in two.
    s_per_bin : float
        corresponding duration of one rt index.
    len_time_axis : float
        max index of the time axis.

    Returns
    -------
    None.

    """
    if rt_idx_stop_bin - rt_idx_start_bin > 1.2 * (dict_k_rt_window + 2 * dict_k_rt_window_edge)/s_per_bin:
          #split windows:
          nb_windows = int(round((rt_idx_stop_bin - rt_idx_start_bin) / (dict_k_rt_window/s_per_bin)))
          #print(nb_windows)
          #exact number of time indexes per window:
          window_idx_duration = int(round(rt_idx_stop_bin - rt_idx_start_bin)/nb_windows)
          #print(window_idx_duration)
          
          edge_idx_duration = math.ceil(dict_k_rt_window_edge/s_per_bin)
          
          #window edges:
          #the first window may have no edge to the left
          start_edge_idx = max(0, rt_idx_start_bin-edge_idx_duration)
          #last window may have no edge to the right
          stop_edge_idx = min(len_time_axis-1, rt_idx_stop_bin + edge_idx_duration)
          
          windows_start_stop = [None]*nb_windows
          #windows_start_stop[0] = [start_edge_idx, start_edge_idx + window_idx_duration + 2*edge_idx_duration]
          for i_w in range(nb_windows):
              windows_start_stop[i_w] = [start_edge_idx  + window_idx_duration* i_w , start_edge_idx + window_idx_duration* (i_w+1) + 2* edge_idx_duration]
              #print(str(windows_start_stop[i_w][0] + edge_idx_duration) + ' ' + str(windows_start_stop[i_w][1]- edge_idx_duration))
          windows_start_stop[-1][1] = stop_edge_idx
          
    else: #make one window
          windows_start_stop = [[rt_idx_start_bin, rt_idx_stop_bin]]
          
    return windows_start_stop

def tof_peak_extraction(data, list_exact_masses, mode):
    """
    mode: mass calibration
    mode: screening
    list_exact_masses: is used in mc mode only; otherwise is calculated from the data
    """
    return None
    
def make_colors():
    #******************
    #define colors
    #https://xkcd.com/color/rgb/
    color1_names = ['xkcd:blue', 'xkcd:green', 'xkcd:orange', 'xkcd:pink', 'xkcd:rose', 
                    'xkcd:grey blue', 'xkcd:light orange', 'xkcd:goldenrod', 'xkcd:forest green', 'xkcd:burnt orange', 
                    'xkcd:apple green', 'xkcd:reddish purple', 'xkcd:dull blue', 'xkcd:hot pink', 'xkcd:peach', 
                    'xkcd:mauve', 'xkcd:dark red', 'xkcd:taupe', 'xkcd:cerulean', 'xkcd:eggplant',
                    'xkcd:blue', 'xkcd:green', 'xkcd:orange', 'xkcd:pink', 'xkcd:rose', 
                    'xkcd:grey blue', 'xkcd:light orange', 'xkcd:goldenrod', 'xkcd:forest green', 'xkcd:burnt orange', 
                    'xkcd:apple green', 'xkcd:reddish purple', 'xkcd:dull blue', 'xkcd:hot pink', 'xkcd:peach', 
                    'xkcd:mauve', 'xkcd:dark red', 'xkcd:taupe', 'xkcd:cerulean', 'xkcd:eggplant',
                    'xkcd:blue', 'xkcd:green', 'xkcd:orange', 'xkcd:pink', 'xkcd:rose', 
                    'xkcd:grey blue', 'xkcd:light orange', 'xkcd:goldenrod', 'xkcd:forest green', 'xkcd:burnt orange', 
                    'xkcd:apple green', 'xkcd:reddish purple', 'xkcd:dull blue', 'xkcd:hot pink', 'xkcd:peach', 
                    'xkcd:mauve', 'xkcd:dark red', 'xkcd:taupe', 'xkcd:cerulean', 'xkcd:eggplant',
                    'xkcd:blue', 'xkcd:green', 'xkcd:orange', 'xkcd:pink', 'xkcd:rose', 
                    'xkcd:grey blue', 'xkcd:light orange', 'xkcd:goldenrod', 'xkcd:forest green', 'xkcd:burnt orange', 
                    'xkcd:apple green', 'xkcd:reddish purple', 'xkcd:dull blue', 'xkcd:hot pink', 'xkcd:peach', 
                    'xkcd:mauve', 'xkcd:dark red', 'xkcd:taupe', 'xkcd:cerulean', 'xkcd:eggplant',
                    'xkcd:blue', 'xkcd:green', 'xkcd:orange', 'xkcd:pink', 'xkcd:rose', 
                    'xkcd:grey blue', 'xkcd:light orange', 'xkcd:goldenrod', 'xkcd:forest green', 'xkcd:burnt orange', 
                    'xkcd:apple green', 'xkcd:reddish purple', 'xkcd:dull blue', 'xkcd:hot pink', 'xkcd:peach', 
                    'xkcd:mauve', 'xkcd:dark red', 'xkcd:taupe', 'xkcd:cerulean', 'xkcd:eggplant',
                    'xkcd:blue', 'xkcd:green', 'xkcd:orange', 'xkcd:pink', 'xkcd:rose', 
                    'xkcd:grey blue', 'xkcd:light orange', 'xkcd:goldenrod', 'xkcd:forest green', 'xkcd:burnt orange', 
                    'xkcd:apple green', 'xkcd:reddish purple', 'xkcd:dull blue', 'xkcd:hot pink', 'xkcd:peach', 
                    'xkcd:mauve', 'xkcd:dark red', 'xkcd:taupe', 'xkcd:cerulean', 'xkcd:eggplant',
                    'xkcd:blue', 'xkcd:green', 'xkcd:orange', 'xkcd:pink', 'xkcd:rose', 
                    'xkcd:grey blue', 'xkcd:light orange', 'xkcd:goldenrod', 'xkcd:forest green', 'xkcd:burnt orange', 
                    'xkcd:apple green', 'xkcd:reddish purple', 'xkcd:dull blue', 'xkcd:hot pink', 'xkcd:peach', 
                    'xkcd:mauve', 'xkcd:dark red', 'xkcd:taupe', 'xkcd:cerulean', 'xkcd:eggplant',
                    'xkcd:blue', 'xkcd:green', 'xkcd:orange', 'xkcd:pink', 'xkcd:rose', 
                    'xkcd:grey blue', 'xkcd:light orange', 'xkcd:goldenrod', 'xkcd:forest green', 'xkcd:burnt orange', 
                    'xkcd:apple green', 'xkcd:reddish purple', 'xkcd:dull blue', 'xkcd:hot pink', 'xkcd:peach', 
                    'xkcd:mauve', 'xkcd:dark red', 'xkcd:taupe', 'xkcd:cerulean', 'xkcd:eggplant',
                    'xkcd:blue', 'xkcd:green', 'xkcd:orange', 'xkcd:pink', 'xkcd:rose', 
                    'xkcd:grey blue', 'xkcd:light orange', 'xkcd:goldenrod', 'xkcd:forest green', 'xkcd:burnt orange', 
                    'xkcd:apple green', 'xkcd:reddish purple', 'xkcd:dull blue', 'xkcd:hot pink', 'xkcd:peach', 
                    'xkcd:mauve', 'xkcd:dark red', 'xkcd:taupe', 'xkcd:cerulean', 'xkcd:eggplant',
                    'xkcd:blue', 'xkcd:green', 'xkcd:orange', 'xkcd:pink', 'xkcd:rose', 
                    'xkcd:grey blue', 'xkcd:light orange', 'xkcd:goldenrod', 'xkcd:forest green', 'xkcd:burnt orange', 
                    'xkcd:apple green', 'xkcd:reddish purple', 'xkcd:dull blue', 'xkcd:hot pink', 'xkcd:peach', 
                    'xkcd:mauve', 'xkcd:dark red', 'xkcd:taupe', 'xkcd:cerulean', 'xkcd:eggplant',
                    'xkcd:blue', 'xkcd:green', 'xkcd:orange', 'xkcd:pink', 'xkcd:rose', 
                    'xkcd:grey blue', 'xkcd:light orange', 'xkcd:goldenrod', 'xkcd:forest green', 'xkcd:burnt orange', 
                    'xkcd:apple green', 'xkcd:reddish purple', 'xkcd:dull blue', 'xkcd:hot pink', 'xkcd:peach', 
                    'xkcd:mauve', 'xkcd:dark red', 'xkcd:taupe', 'xkcd:cerulean', 'xkcd:eggplant',
                    'xkcd:blue', 'xkcd:green', 'xkcd:orange', 'xkcd:pink', 'xkcd:rose', 
                    'xkcd:grey blue', 'xkcd:light orange', 'xkcd:goldenrod', 'xkcd:forest green', 'xkcd:burnt orange', 
                    'xkcd:apple green', 'xkcd:reddish purple', 'xkcd:dull blue', 'xkcd:hot pink', 'xkcd:peach', 
                    'xkcd:mauve', 'xkcd:dark red', 'xkcd:taupe', 'xkcd:cerulean', 'xkcd:eggplant']
    color2_names = ['xkcd:sky blue', 'xkcd:teal', 'xkcd:peach', 'xkcd:mauve', 'xkcd:dark red', 'xkcd:taupe', 'xkcd:cerulean', 'xkcd:eggplant']

    #grey scale
    #color3_names = ['xkcd:grey', 'xkcd:slate', 'xkcd:black', 'xkcd:grey', 'xkcd:slate', 'xkcd:black', 'xkcd:grey', 'xkcd:slate', 'xkcd:black']

    color3_names = ['xkcd:dusty pink', 'xkcd:grey', 'xkcd:wine red', 'xkcd:grey purple', 'xkcd:dark rose', 
                    'xkcd:dark salmon', 'xkcd:slate', 'xkcd:black', 'xkcd:eggplant', 'xkcd:taupe', 'xkcd:mauve']
    
    return color1_names, color2_names, color3_names
