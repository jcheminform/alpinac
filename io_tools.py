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


Created on Tue Aug 18 10:14:09 2020
"""

import pathlib
from pathlib import Path
from jcamp import JCAMP_reader #to read NIST spectra
from periodic_table import formula_to_mass
import periodic_table as chem_data



#format of fragment file to read from:
ic_rt = 0           #rt, seconds
ic_m = 1            #mass, average of peak
ic_m_u = 2          #uncertainty of mass, ppm (measurement noise only) 201001 to add cal u as well
ic_m_cal_u = 3
ic_area = 4         #area of peak
ic_area_u = 5       #uncertainty of area
ic_m_sigma = 6
ic_alpha = 7
ic_comp_bin = 8     #index of group of co-eluting fragments
ic_LOD = 9          #LOD value
ic_SN = 10           #signal to noise ratio
ic_ionisation = 11



# indices of mass and comp_bin in each tuple of output:
idx_rt = 0
idx_m = 1
idx_m_u = 2
idx_m_cal_u = 3
idx_area = 4
idx_area_u = 5
idx_peak_width = 6
idx_peak_alpha = 7
idx_comp_bin = 8
idx_LOD = 9
idx_ionisation = 10
idx_SN = None

#format of formula file to read from:
iform_percent_norm_I = 0
iform_rt = 1
iform_total_I = 2
iform_total_I_percent = 3
iform_mass_meas = 4
iform_mass_exact = 5
iform_mass_u_ppm = 6
iform_mass_diff_ppm = 7
iform_formula = 8
iform_DBE = 9
iform_I = 10
iform_I_f = 11
iform_ionisation = 12


def make_dict_of_paths():

    #dic_labels = {"number":0,
    Esperanza_data_path = 'E:/TOF-Data'
    mygu_data_path = 'D:/TOF-Data'
    
    is_esperanza = False
    if is_esperanza: 
        data_path = Esperanza_data_path
    else:
        data_path = mygu_data_path    
    #define all fixed paths:
    dict_path = {'path_current': Path.cwd(),
                 'path_hdf5': Path(data_path)}


    
    dict_path['path_hdf5_2018'] = dict_path['path_hdf5'] / '2018'
    dict_path['path_hdf5_2019'] = dict_path['path_hdf5'] / '2019'
    dict_path['path_hdf5_2020'] = dict_path['path_hdf5'] / '2020'
    dict_path['path_hdf5_2021'] = dict_path['path_hdf5'] / '2021'
    dict_path['path_CIEI'] = Path('D:/Tofwerk_CI_EI/Medusa-GC-CIEI-TOF')
    
    #dict_path['bafu_campaign'] = dict_path['path_hdf5'] / 'Bafu_campaign'
    
    
    
    dict_path['path_TOF_computer'] = Path('//152.88.82.57/d/DATA')
    
    dict_path['path_masscal_output'] = dict_path['path_current'] / 'data' / 'TOF_mass_calibration' / 'mass_calibration'
    dict_path['path_masscal_output_image'] = dict_path['path_current'] / 'data' / 'TOF_mass_calibration' / 'mass_calibration_images'

    dict_path['path_GCWerks_hdf5'] = dict_path['path_current'] / 'data' / 'TOF_export_GCWerks' / 'TOF_hdf5'
    dict_path['path_GCWerks_output'] = dict_path['path_current'] / 'data' / 'TOF_export_GCWerks' / 'GCWerks_export'
    dict_path['peak_list'] = dict_path['path_current'] / 'data' / 'peak_list'
    
    #dict_path['path_GCWerks_input'] = dict_path['path_current'] / 'data' / 'TOF_export_GCWerks' / 'GCWerks_input'

    #dict_path['path_nontarget_hdf5'] = dict_path['path_current'] / 'data' / 'nontarget_screening' / 'TOF_hdf5'
    dict_path['path_nontarget_hdf5'] = dict_path['path_hdf5']
    dict_path['path_nontarget_output'] = dict_path['path_current'] / 'data' / 'nontarget_screening' / 'nontarget_output'
    dict_path['path_nontarget_periodtab'] = dict_path['path_current'] / 'data' / 'nontarget_screening' / 'periodic_table'
    
    dict_path['NIST_EI'] = dict_path['path_current'] / 'data' / 'nontarget_screening' / 'NIST_EI'
    
    dict_path['fragments'] = dict_path['path_current'] / 'data' / 'nontarget_screening' / 'fragments'
    dict_path['fragments_train'] = dict_path['path_current'] / 'data' / 'nontarget_screening' / 'fragments' / 'ALPINAC_train'
    dict_path['fragments_val'] = dict_path['path_current'] / 'data' / 'nontarget_screening' / 'fragments' / 'ALPINAC_val'
    dict_path['fragments_no_mc'] = dict_path['path_current'] / 'data' / 'nontarget_screening' / 'fragments_no_mc'
    dict_path['formulas'] = dict_path['path_current'] / 'data' / 'nontarget_screening' / 'formulas' 
    dict_path['mass_spectra'] = dict_path['path_current'] / 'data' / 'nontarget_screening' / 'mass_spectra'
    dict_path['timeseries'] = dict_path['path_current'] / 'data' / 'nontarget_screening' / 'timeseries'

    dict_path['GCWerks_peakID_file'] = dict_path['path_current'] / 'data' / 'peak_list' / 'GCWerks_peakID_file'
    dict_path['Tofware_peak_list'] = dict_path['path_current'] / 'data' / 'peak_list' / 'Tofware_peak_list'
    
    return dict_path



  
    
def get_data_from_frag_file(filename: Path) -> dict:
    """Read data (mass, etc) and return a structured dict.
    
    Read data from a file obtained after calibration, adjust mass: add
    the electron mass to compensate the EI (electron ionisation), return
    a dict indexed by bin number.

    INPUT:
    - ``filename`` -- Path object, can be opened with .open('r')
    - ``electron_mass`` -- float, mass of the electron (0.00054858)
    
    OUTPUT: a dict indexed by bin number (integer value of column `ic_comp_bin`)
    Data are in this order:
    RT	mass	mass_u_ppm	mass_cal_u_ppm	area	area_u	peak_width	peak_alpha	compound_bin	LOD
    """
    # open filename with Path syntax
    f = filename.open('r')
    # 2. Read file
    frag_data_batch = {}
    for line in f:
        raw = line.split()
        if raw[0] == 'RT': # header
            continue
        row = [float(raw[i]) for i in range(ic_alpha+1)] + [int(raw[ic_comp_bin])] + [float(raw[ic_LOD])] + [raw[ic_ionisation]]
        row[ic_m_u] = row[ic_m_u]/10**6 * row[ic_m] #convert u mass in ppm in m/z
        row[ic_m_cal_u] = row[ic_m_cal_u]/10**6 * row[ic_m] #convert mass_cal_u from ppm to m/z
        
        row[ic_m] += chem_data.electron_mass #add mass of one electron to measured, ionised masses.
        #print(row)
        comp_bin = row[idx_comp_bin]
        #print(comp_bin)
        if comp_bin in frag_data_batch:
            frag_data_batch[comp_bin].append(tuple(row))
        else:
            frag_data_batch[comp_bin] = [tuple(row)]
    f.close()
    #print(frag_data_batch)
    
    
    return frag_data_batch

def get_data_from_frag_unit_mass_file(filename: Path) -> dict:
    """Read data (mass, etc) and return a structured dict.
    
    Read data from a file obtained after calibration, adjust mass: add
    the electron mass to compensate the EI (electron ionisation), return
    a dict indexed by bin number.

    INPUT:
    - ``filename`` -- Path object, can be opened with .open('r')
    - ``electron_mass`` -- float, mass of the electron (0.00054858)
    
    OUTPUT: a dict indexed by bin number (integer value of column `ic_comp_bin`)
    """
    # open filename with Path syntax
    f = filename.open('r')
    # 2. Read file
    frag_data_batch = {}
    for line in f:
        raw = line.split()
        if raw[0] == 'RT': # header
            continue
        row = [float(raw[i]) for i in range(ic_alpha)] + [int(raw[ic_comp_bin])] + [float(raw[ic_LOD])] + [raw[ic_ionisation]]
        row[ic_m_u] = row[ic_m_u] #u mass should be 0.24999
        row[ic_m] += chem_data.electron_mass #add mass of one electron to measured, ionised masses.
        comp_bin = row[idx_comp_bin]
        if comp_bin in frag_data_batch:
            frag_data_batch[comp_bin].append(tuple(row))
        else:
            frag_data_batch[comp_bin] = [tuple(row)]
    f.close()
    
    return frag_data_batch

class NISTMetadata():
    def __init__(self):
        self.NIST_chem_name = None
        self.NIST_target_formula = None
        self.NIST_mol_mass_exact = None
        self.NIST_CAS = None
    

def get_data_from_NIST_file(filename: Path) -> dict:
    """load unit mass resolution EI mass spectra from NIST database

    Read data from NIST file. This is unit mass resolution. Format the 
    data so that it matches the format of the function 
    `get_data_from_frag_file`.

    INPUT:
    - ``filename`` -- Path object, can be opened with .open('r')
    
    OUTPUT: a dictionary indexed by the bin number (usually 0 for NIST)
    """


    dict_target = JCAMP_reader(filename)

    frag_data_batch = {0:[]}
    frag_metadata = NISTMetadata()
    
    frag_metadata.NIST_chem_name = dict_target['title']
    frag_metadata.NIST_target_formula = ''.join(dict_target['molform'].split())
    frag_metadata.NIST_mol_mass_exact = formula_to_mass(frag_metadata.NIST_target_formula)
    if 'cas registry no' in dict_target:
        frag_metadata.NIST_CAS = dict_target['cas registry no']
    else:
        frag_metadata.NIST_CAS = 'no_CAS_' + dict_target['title']
    #print('CAS number: ' + str(self.NIST_CAS))
    
    #meas_mass = dict_target['x']
    #meas_I = dict_target['y']
    #meas_rt=[1]*len(meas_mass)
    #meas_mass_u=[0.25]*len(meas_mass)
    LOD = min(dict_target['y'])*0.1

    for i in range(len(dict_target['x'])):
        frag_data_batch[0].append(tuple([
                float(0.0), #rt
                float(dict_target['x'][i]),#*float(self.NIST_mol_mass_exact)/float(546), #meas mass
                float(0.24999), #meas mass u [m/z]
                float(0.0), #m_cal_u
                float(dict_target['y'][i]), #area
                float(0.05)*dict_target['y'][i], #area u
                float(0.0), #ic_width
                float(0.0), #ic_alpha
                int(0), #comp_bin
                LOD,
                float(0),
                'EI'
                ]))
    return frag_data_batch, frag_metadata

# =============================================================================
# #format of formula file to read from:
# iform_percent_norm_I = 0
# iform_rt = 1
# iform_total_I = 2
# iform_total_I_percent = 3
# iform_mass_meas = 4
# iform_mass_exact = 5
# iform_mass_u_ppm = 6
# iform_mass_diff_ppm = 7
# iform_formula = 8
# iform_DBE = 9
# iform_I = 10 #intensity computed for this formula only
# iform_I_f = 11
# iform_ionisation = 12
# =============================================================================

def get_data_from_formula_file(filename: Path) -> dict:
    """Read data (mass, etc) and return a structured dict.
    
    Read data from a file obtained after identification.

    INPUT:
    - ``filename`` -- Path object, can be opened with .open('r')

    
    OUTPUT: a dict indexed by bin number (integer value of column `ic_comp_bin`)
    important will be colun 11, fraction that this formula represent for total signal at this mass peak.
    
    """
    # open filename with Path syntax
    f = filename.open('r')
    # 2. Read file
    frag_data_batch = {}
    extract_data = False
    for line in f:
        raw = line.split()
        print(raw[0])
        #while extract_data == False:
        if extract_data == False and raw[0] != '%': # header
            continue
        else:
            extract_data = True
            #print("extract_data = True")
                
        if extract_data == True and raw[0] != '%' :
            #print("extract")
            #print(raw)
            row = [float(raw[i]) for i in range(iform_formula)] + \
            [raw[iform_formula]] + \
            [float(raw[iform_DBE])] + \
            [float(raw[iform_I])] + \
            [float(raw[iform_I_f])] + \
            [raw[iform_ionisation]]
            #print(row)
            comp_bin = 0
            #print(comp_bin)
            if comp_bin in frag_data_batch:
                frag_data_batch[comp_bin].append(tuple(row))
            else:
                frag_data_batch[comp_bin] = [tuple(row)]
    f.close()
    #print(frag_data_batch)
    
    
    return frag_data_batch

def export_mc(file,
              rt_idx_centre_list, 
              params_optimised, 
              sigma_tofidx_slope, 
              sigma_tofidx_b, 
              sigma_m_slope, 
              sigma_m_b, 
              alpha_pseudo_voigt,
              masscal_list_exact_export, 
              mass_cal_u, 
              masscal_list_exact,
              s_per_bin):
    

    with file.open('w') as outfile:
        outfile.write(str(len(params_optimised)) + '\t' + 'Number of points, mass cal. vs RT' +   '\n')
        outfile.write(str(len(mass_cal_u)) + '\t' + 'Number of points, mass cal. u vs mass' + '\n')
        outfile.write(str(params_optimised[0]['mode']) + '\t' + 'Mass cal. mode'  + '\n')
        outfile.write(str(sigma_tofidx_slope) + '\t' + 'Peak shape, sigma_slope vs TOF index' +  '\n')
        outfile.write(str(sigma_tofidx_b)  + '\t' + 'Peak shape, sigma_b vs TOF index' + '\n')
        outfile.write(str(sigma_m_slope) + '\t' + 'Peak shape, sigma_slope vs mass' +  '\n')
        outfile.write(str(sigma_m_b) + '\t' +'Peak shape, sigma_b vs mass'  +  '\n')
        outfile.write(str(alpha_pseudo_voigt) + '\t' +'alpha, for tofidx Pseudo-Voigt fit'  +  '\n')
        outfile.write(str(s_per_bin) + '\t' + 's_per_bin' +  '\n')
        for idx_rt in range(len(params_optimised)):
            if params_optimised[0]['mode'] == 1:
                outfile.write(str(rt_idx_centre_list[idx_rt]) + '\t' + 
                              str(params_optimised[idx_rt]['p1']) + '\t' + 
                              str(params_optimised[idx_rt]['p2']) + '\t' + 
                              str(params_optimised[idx_rt]['p3']) + '\n')
            elif params_optimised[0]['mode'] == 2:
                outfile.write(str(rt_idx_centre_list[idx_rt]) + '\t' + 
                              str(params_optimised[idx_rt]['p1']) + '\t' + 
                              str(params_optimised[idx_rt]['p2']) + '\t' + 
                              str(params_optimised[idx_rt]['p3']) + '\n')
            elif params_optimised[0]['mode'] == 3:
                outfile.write(str(rt_idx_centre_list[idx_rt]) + '\t' + 
                              str(params_optimised[idx_rt]['p1']) + '\t' + 
                              str(params_optimised[idx_rt]['p2']) + '\t' + 
                              str(params_optimised[idx_rt]['p3']) + '\t' + 
                              str(params_optimised[idx_rt]['p4']) + '\n')
            elif params_optimised[0]['mode'] == 4:
                outfile.write(str(rt_idx_centre_list[idx_rt]) + '\t' + 
                              str(params_optimised[idx_rt]['p1']) + '\t' + 
                              str(params_optimised[idx_rt]['p2']) + '\t' +
                              str(params_optimised[idx_rt]['p3']) + '\t' + 
                              str(params_optimised[idx_rt]['p4']) + '\t' + 
                              str(params_optimised[idx_rt]['p5']) + '\n')
        for idx_m in range(len(mass_cal_u)):
            outfile.write('{:.7f}'.format(masscal_list_exact_export[idx_m]) + '\t' + '{:.2f}'.format(mass_cal_u[idx_m]) + '\n')
            
        #20200319 new: write detailed list of masses used for mc per time bin
        #to be used later for re-calibration during identification of mass peaks.
        
    return None
            
            
class MassCalibrationData():
    def __init__(self, filename: Path)-> None:
        """Read data from mass calibration file
        
        INPUT:
        - ``filename`` -- Path object, can be opened with .open('r')
        """
        mass_cal_file = filename.open('r')
        nb_mc_params = 9 #number of parameters saved in the mass cal file
        self.mass_cal_u = []
        self.mass_cal_u_exact_mass = []
        self.mass_cal_tof_values = []
        
        # line 1:
        line = mass_cal_file.readline(); raw=line.split()
        mass_cal_len = int(raw[0])
        # line 2
        line = mass_cal_file.readline(); raw=line.split()
        mass_cal_u_len = int(raw[0])
        # line 3
        line = mass_cal_file.readline(); raw=line.split()
        self.mass_cal_mode = int(raw[0])
        # line 4
        line = mass_cal_file.readline(); raw=line.split()
        self.sigma_tofidx_slope = float(raw[0])
        # line 5
        line = mass_cal_file.readline(); raw=line.split()
        self.sigma_tofidx_b = float(raw[0])
        # line 6
        line = mass_cal_file.readline(); raw=line.split()
        self.sigma_m_slope = float(raw[0])
        # line 7
        line = mass_cal_file.readline(); raw=line.split()
        self.sigma_m_b = float(raw[0])
        # line 8
        line = mass_cal_file.readline(); raw=line.split()
        self.alpha_pseudo_voigt = float(raw[0])
        # line 9
        line = mass_cal_file.readline(); raw=line.split()
        self.s_per_bin = float(raw[0])
        
        #load mass calibration parameters: rt, param1, param2, param3 (using mode 2)
        self.mass_cal_parameters = [None]*mass_cal_len
        for i in range(mass_cal_len):
            line = mass_cal_file.readline(); raw=line.split()
            self.mass_cal_parameters[i] = [float(raw[0]) * self.s_per_bin, float(raw[1]), float(raw[2]), float(raw[3])]
        #extract mass calibration uncertainty, varying over the mass domain:
        if mass_cal_u_len > 0:
            self.mass_cal_u_exact_mass = [None]*mass_cal_u_len
            self.mass_cal_u = [None]*mass_cal_u_len
            for i in range(mass_cal_u_len):
                line = mass_cal_file.readline(); raw=line.split()
                self.mass_cal_u_exact_mass[i] = float(raw[0])
                self.mass_cal_u[i] = float(raw[1])
        else:
            self.mass_cal_u_exact_mass = [0, 400]
            self.mass_cal_u = [50, 50]
            
        #load rt_idx, tof, exact_mass, u_calib (ppm)
        self.mass_cal_tof_values = []
        line = mass_cal_file.readline();
        while not len(line) == 0:
            raw=line.split()
            self.mass_cal_tof_values.append([float(raw[0])* self.s_per_bin, float(raw[1]), float(raw[2]), float(raw[3])])
            line = mass_cal_file.readline();
        mass_cal_file.close()
