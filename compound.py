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


Created on Mon Aug 17 11:18:06 2020

"""
import numpy as np
import math
import operator

from statistics import median

import io_tools
from io_tools import MassCalibrationData
import periodic_table as chem_data
from utils_data_extraction import interpol_mass_cal, f_pseudoVoigt, f_Gauss, fit_linear_np

run_type_tofwerk_tof = 0
run_type_NIST = 1
run_type_no_mass_cal = 2
run_type_unit_mass = 3

class Compound():
    """
    Compound is a group of fragments belonging to the same RT index.
    They co-elute and therefore are likely to belong to the same compound.
    Note: they could also belong to several co-eluting compounds.
    """
    mass_u_k_loose = float(3.0) # arbitrary choice. Should it be larger than dict_k['m_u_k']?

    def __init__(self, one_batch_fragment, compound_number, LOD):
        """
        
        INPUT:
        - ``one_batch_fragment``: a list of lists, output of io_tools:get_data_from_frag_file
        - ``compound_number``: a number ( bin number)
        - "str_ionisation": string that is a code to inform what type of ionisation was used

        """
        self.fragment_index = compound_number
        self.data = one_batch_fragment
        #self.data = [one_batch_fragment[i][:] for i in range(len(one_batch_fragment)) if one_batch_fragment[i][io_tools.idx_ionisation] == str_ionisation]
        #data sorted by increasing mass:
        #needed to speed up looking for isotopologues
        self.data.sort(key = lambda x: x[io_tools.idx_m])
        #self.data.sort(key= operator.itemgetter(1)) #order by increasing masses
        #print(self.data)
        
        self.meas_len = len(self.data)
        print("==== len of measured data  " + str(self.meas_len))
        

        # re-extract data (split by columns)
        #retention time of centre of peak, list
        self.meas_rt= [self.data[idx][io_tools.idx_rt] for idx in range(self.meas_len)]
        
        #Ionisation mode
        #In the future this may influence the minimum allowed DBE value
        self.ionisation = [self.data[idx][io_tools.idx_ionisation] for idx in range(self.meas_len)]

        
        #Measured mass
        self.meas_mass = [self.data[idx][io_tools.idx_m] for idx in range(self.meas_len)]
        
        #uncertainty due to measurement noise, type A
        self.meas_mass_u = [self.data[idx][io_tools.idx_m_u] for idx in range(self.meas_len)]
        
                
        #measured intensity of peak 
        #for Empa APRECON-GC-TOFMS: 
        #each peak integrated by pseudo-voigt fit in mass domain 
        #and Gaussian with tailing in rt domain
        self.meas_I = [self.data[idx][io_tools.idx_area] for idx in range(self.meas_len)]
        #self.meas_I_u = [self.data[idx][io_tools.idx_area_u] for idx in range(self.meas_len)]

        #LOD at the measured mass (depends on background noise, is mass dependent)
        #Not used at present, may be used in the future
        #to define a LOD specific for each measured mass
        self.meas_LOD = [self.data[idx][io_tools.idx_LOD] for idx in range(self.meas_len)]

        self.average_rt = sum([self.meas_rt[i]*self.meas_I[i] for i in range(self.meas_len)])/sum(self.meas_I)
        self.sum_I = sum(self.meas_I)
        self.meas_I_max = max(self.meas_I)
        self.meas_I_min = min(self.meas_I)
        
        #contains total assigned signal after optimisation, per measured mass
        self.assigned_I = [float(0)]*self.meas_len    
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        """
        Trying to define a more meaningful and flexible LOD
        """
        #signal: limit of detection
        self.LOD=float(min(LOD, self.meas_I_min))

        
        self.atom_list_likelihood = [None] 
        self.atom_list_likelihood_idx_sorted = [None]
        
        self.unidentified_mass_idx_list = []
        self.identified_mass_idx_list = []

    def initialise_mass_profile(self):
        self.meas_mass_profile_x = []
        self.delta_mass_u = float(0.0)
        self.meas_mass_profile_idx = []
        self.meas_mass_profile = []

    def do_meas_mass_profile_x(self, m_start: float, m_stop: float, mass_res: float) -> None:
        """Initialize self.meas_mass_profile_x
        
        Write np.array with x abcisse values, which is the mass axis,
        where the intensities for each theoretical isotopologue
        and for the measured signal 
        will be computed.
        
        INPUT:
        - ``m_start``: minimum mass abcissa
        - ``m_stop``: maximum mass abcissa
        - ``mass_res``: interval in between each mass point
        
        OUTPUT: None
        """
        self.meas_mass_profile_x = np.linspace(round(m_start) - self.mass_axis_offset, round(m_stop) + self.mass_axis_offset, round(mass_res*(round(m_stop) + self.mass_axis_offset - (round(m_start) - self.mass_axis_offset) + 1)))              


        
    def comp_list_identified_unidentified_mass(self, G):
        """
        Parameters
        ----------
        G : TYPE: networkx graph. Each node contains an instance of Isotopologues.
            DESCRIPTION.

        Returns
        -------
        List of indexes of identified masses.
        List of indexes of non-identified masses.
        Calculates total assigned signal for each measured mass. Total may be > measured.

        """
        #used 20200730
        self.identified_mass_idx_list = []
        self.unidentified_mass_idx_list = []
        self.assigned_I = [float(0)]*self.meas_len
        for i_node in G.nodes():
            for idx_m in range(len(G.nodes[i_node]['g_iso'].meas_mass_idx)):
                if G.nodes[i_node]['g_iso'].meas_mass_idx[idx_m] is not None:
                    self.assigned_I[idx_m] += min(G.nodes[i_node]['g_iso'].iso_list_I_rel[idx_m] * G.nodes[i_node]['g_iso'].k_guess, self.meas_I[idx_m])
                    if G.nodes[i_node]['g_iso'].meas_mass_idx[idx_m] not in self.identified_mass_idx_list:
                        self.identified_mass_idx_list.append(G.nodes[i_node]['g_iso'].meas_mass_idx[idx_m])
                    
        for idx_m in range(self.meas_len):
            if idx_m not in self.identified_mass_idx_list:
                self.unidentified_mass_idx_list.append(idx_m)

class CompoundTofwerkTOF(Compound):
    """A class dedicated to TofwerkTOF data that inherits from Compound
    """
    def __init__(self, one_batch_fragment: list, compound_number: int, mass_u_k: float, LOD: float): #, mass_calibration_data: MassCalibrationData
        # call superclass constructor
        super().__init__(one_batch_fragment, compound_number, LOD)
        self.run_type = run_type_tofwerk_tof
        u_threshold = 10.0 #10.0
        #u_threshold = [max(u_threshold_k*2, self.meas_mass_u_cal[i]) for i in self.meas_mass_u_cal if self.meas_mass[i] < 68.0]
        #mass_calibration_data.mass_cal_u = [max(u_threshold, mass_calibration_data.mass_cal_u[i]) for i in range(len(mass_calibration_data.mass_cal_u))]
        #uncertainty of mass calibration, list
        self.meas_mass_u_cal = [max(u_threshold*self.data[idx][io_tools.idx_m]/1000000.0, self.data[idx][io_tools.idx_m_cal_u]) for idx in range(self.meas_len)] #uncertainty due to calibration, type B
        #self.meas_mass_u_cal = np.interp(np.array(self.meas_mass), np.array(mass_calibration_data.mass_cal_u_exact_mass), np.array(mass_calibration_data.mass_cal_u)*np.array(mass_calibration_data.mass_cal_u_exact_mass)/10**6)
        #self.meas_mass_u_cal = [max(u_threshold, val) for val in self.meas_mass_u_cal]
        #mygu 20201210

        #we always have problem with small masses beeing outside the exact mass range
        #because mass cal below 60 is badly constrained
        #so we try to systematically increase u_cal for masses < 60
        #this should not be done here but beforehand, during mass cal process
        """
        for i in range(len(self.meas_mass_u_cal)):
            if self.meas_mass[i] < float(60.0):
                self.meas_mass_u_cal[i] = max ( self.meas_mass_u_cal[i], float(40.0)/float(1000000) * self.meas_mass[i])
        """
        self.meas_mass_u_loose = [self.mass_u_k_loose * math.sqrt(self.meas_mass_u[mass_idx]**2 + self.meas_mass_u_cal[mass_idx]**2) for mass_idx in range(self.meas_len)]

        #local mass calibration parameters around average rt:        
        #self.mass_cal_params_comp = interpol_mass_cal(self.average_rt, mass_calibration_data.mass_cal_parameters, mass_calibration_data.mass_cal_mode)
        self.mass_u_k = mass_u_k #2.5
        self.mass_axis_offset = 0.5

        self.meas_mass_U_combined = [self.mass_u_k * math.sqrt(self.meas_mass_u[idx_m]**2 + self.meas_mass_u_cal[idx_m]**2) for idx_m in range(self.meas_len)]
        self.meas_mass_min = [self.meas_mass[idx_m] - self.meas_mass_U_combined[idx_m] for idx_m in range(self.meas_len)]
        self.meas_mass_max = [self.meas_mass[idx_m] + self.meas_mass_U_combined[idx_m] for idx_m in range(self.meas_len)]
        
        self.meas_mass_u_median = median(self.meas_mass_U_combined)
        self.meas_mass_u_median_ppm = median([self.meas_mass_U_combined[idx_m]/self.meas_mass[idx_m]*1000000 for idx_m in range(self.meas_len)])
        print("meas_mass_u_median_ppm  " + str(self.meas_mass_u_median_ppm))

        #self.sigma_m_slope = mass_calibration_data.sigma_m_slope
        #self.sigma_m_b = mass_calibration_data.sigma_m_b
        self.sigma_m_slope, self.sigma_m_b, linear_fit, std_dev_linear_fit = fit_linear_np(np.array(self.meas_mass), np.array([self.data[idx][io_tools.idx_peak_width] for idx in range(self.meas_len)]))
        #self.alpha_pseudo_voigt = mass_calibration_data.alpha_pseudo_voigt
        self.alpha_pseudo_voigt = sum([self.data[idx][io_tools.idx_peak_alpha] for idx in range(self.meas_len)])/self.meas_len
        #print("sigma_m_slope, sigma_m_b, alpha ")
        #print([self.sigma_m_slope, self.sigma_m_b, self.alpha_pseudo_voigt])
        #print([self.data[idx][io_tools.idx_peak_width] for idx in range(self.meas_len)])
    def get_mass_range_looser_uncertainty(self, mass_idx: int, other_mass: float) -> float:
        #compute the sigma (standard deviation)
        #that represents the measured width of the peak in mass domain
        #(not the sigma of the mean position other_mass, which is much smaller)
        mass_sigma_peak = self.sigma_m_b + self.sigma_m_slope*other_mass
        mass_uncertainty_range = self.mass_u_k_loose * mass_sigma_peak + self.meas_mass_u_loose[mass_idx]
        return mass_uncertainty_range

    def do_meas_mass_profile_x(self, m_start: float, m_stop: float, ppm_mass) -> None:
        """Initialize self.meas_mass_profile_x
        
        TODO Myriam what does this function do?
        
        INPUT:
        - ``m_start``:
        - ``m_stop``:
        - ``ppm_mass``: TODO Myriam int or float?
        
        OUTPUT: None
        """
        mass_res = (float(1000000.0) / (ppm_mass * float(m_start+m_stop)/2.0))
        super().do_meas_mass_profile_x(m_start, m_stop, mass_res)

    def do_meas_mass_profile(self, delta_mass):
        """Do measured mass profile

        INPUT:
        - ``delta_mass``: potential mass axis offset, at initial step it is zero.
          Used to shift position of measured masses by a constant offset.
          delta_mass is optimised by the lmfit algo.
        """
        #indexes of measured masses within window range:        
        self.meas_mass_profile_idx = [idx for idx in range(self.meas_len) if self.meas_mass[idx] >=  self.meas_mass_profile_x[0] and self.meas_mass[idx] <= self.meas_mass_profile_x[-1]]
        #print("self.meas_mass_profile_idx")
        #print(self.meas_mass_profile_idx)
        #print([self.meas_mass_u_cal[i] for i in self.meas_mass_profile_idx])
        self.delta_mass_u = min([self.meas_mass_u_cal[i] for i in self.meas_mass_profile_idx])
        self.meas_mass_profile = np.zeros(len(self.meas_mass_profile_x))
        #delta_mass = ppm_mass
        for idx_m in range(len(self.meas_mass_profile_idx)): # and self.meas_mass[self.meas_mass_profile_idx[idx_m]] <= max(self.meas_mass_profile_x):
            #f_pseudoVoigt(x, A, mu, sigma, alpha)
            sigma = self.meas_mass[self.meas_mass_profile_idx[idx_m]] * self.sigma_m_slope + self.sigma_m_b
            self.meas_mass_profile += f_pseudoVoigt(self.meas_mass_profile_x, self.meas_I[self.meas_mass_profile_idx[idx_m]], self.meas_mass[self.meas_mass_profile_idx[idx_m]]+ delta_mass, sigma, self.alpha_pseudo_voigt)

    def do_iso_profiles(self, iso_list_idx: list, isotopologues_list: list):
        """
        Generate a list, containing for each row a numpy array
        corresponding to the sum profile of one isotopologue series.
        All isotopologue series are in the same mass window.
        
        INPUT:
        -``iso_list_idx``:
        -``isotopologues_list``:
        """
        #per isotopologue we need:
        #centre of mass (exact mass)
        #sigma (cf. mass cal)
        #Arel, determined knowing relative intensity
        #alpha, mass cal
        #we set a factor k to link Arel to measured intensity
        #at beginning k is calculated using Arel of abundant = 1 and I of matched frag abund.
        no_iso = len(iso_list_idx)
        #self.k_guesses_max = [None] * no_iso
        self.iso_profiles = [None] * no_iso
        for i_no_iso in range(no_iso):
            iso_profile = np.zeros(len(self.meas_mass_profile_x))
            #for i_iso in range(len(iso_list_idx[i_no_iso])):
            #we take one isotopologue series
            #and sum up the profile
            #at the required mass values
            #print('i_no_iso ' + str(i_no_iso))
            isotopologue = isotopologues_list[iso_list_idx[i_no_iso]]
            for i_iso_series in range(isotopologue.iso_list_len):
                #f_pseudoVoigt(x, A, mu, sigma, alpha)
                iso_mass = chem_data.get_mass(isotopologue.iso_list[i_iso_series])
                sigma = iso_mass * self.sigma_m_slope + self.sigma_m_b
                iso_profile += f_pseudoVoigt(self.meas_mass_profile_x, 
                                             isotopologue.iso_list_I_rel[i_iso_series],
                                             iso_mass,
                                             sigma,
                                             self.alpha_pseudo_voigt)
            self.iso_profiles[i_no_iso] = iso_profile.copy()


class CompoundNoMassCalibration(Compound):
    """A class dedicated to data without mass calibration that inherits from Compound
    """
    def __init__(self, one_batch_fragment: list, compound_number: int, LOD: float):
        super().__init__(one_batch_fragment, compound_number)
        self.run_type = run_type_no_mass_cal
        self.meas_mass_u_cal = np.zeros(len(self.meas_mass)) # numpy array zero vector
        self.meas_mass_u_loose = [self.mass_u_k_loose * self.meas_mass_u[mass_idx] for mass_idx in range(self.meas_len)]
        self.mass_u_k = 2.5
        self.mass_axis_offset = 0.5
        self.meas_mass_U_combined = [self.mass_u_k * math.sqrt(self.meas_mass_u[idx_m]**2 + self.meas_mass_u_cal[idx_m]**2) for idx_m in range(self.meas_len)]
        self.meas_mass_min = [self.meas_mass[idx_m] - self.meas_mass_U_combined[idx_m] for idx_m in range(self.meas_len)]
        self.meas_mass_max = [self.meas_mass[idx_m] + self.meas_mass_U_combined[idx_m] for idx_m in range(self.meas_len)]

        self.meas_mass_u_median = 0.0
        self.meas_mass_u_median_ppm = 0.0
        
    def get_mass_range_looser_uncertainty(self, mass_idx: int, other_mass: float=None) -> float:
        return self.meas_mass_u_loose[mass_idx]

    def do_meas_mass_profile_x(self, m_start: float, m_stop: float, ppm_mass) -> None:
        """Initialize self.meas_mass_profile_x
        
        TODO Myriam what does this function do?
        
        INPUT:
        - ``m_start``:
        - ``m_stop``:
        - ``ppm_mass``: TODO Myriam int or float?
        
        OUTPUT: None
        """
        mass_res = (float(1000000.0) / (ppm_mass * float(m_start+m_stop)/2.0))
        super().do_meas_mass_profile_x(m_start, m_stop, mass_res)

    def do_meas_mass_profile(self, delta_mass):
        """Do measured mass profile

        INPUT:
        - ``delta_mass``: potential mass axis offset, at initial step it is zero.
          Used to shift position of measured masses by a constant offset.
          delta_mass is optimised by the lmfit algo.
        """
        #indexes of measured masses within window range:        
        self.meas_mass_profile_idx = [idx for idx in range(self.meas_len) if self.meas_mass[idx] >=  self.meas_mass_profile_x[0] and  self.meas_mass[idx] <= self.meas_mass_profile_x[-1]]
        self.delta_mass_u = float(0.0)
        self.meas_mass_profile = np.zeros(len(self.meas_mass_profile_x))
        
        for idx_m in range(len(self.meas_mass_profile_idx)): # and self.meas_mass[self.meas_mass_profile_idx[idx_m]] <= max(self.meas_mass_profile_x):
            #we assume Gaussian shape of the peak.
            #try to define a sigma based on uncertainty of measured mass
            sigma = self.meas_mass_u[self.meas_mass_profile_idx[idx_m]]
            self.meas_mass_profile += f_Gauss(self.meas_mass_profile_x, self.meas_I[self.meas_mass_profile_idx[idx_m]], self.meas_mass[self.meas_mass_profile_idx[idx_m]]+ delta_mass, sigma)

    def do_iso_profiles(self, iso_list_idx: list, isotopologues_list: list):
        """
        Generate a list, containing for each row a numpy array
        corresponding to the sum profile of one isotopologue series.
        All isotopologue series are in the same mass window.
        
        INPUT:
        -``iso_list_idx``:
        -``isotopologues_list``:
        """
        no_iso = len(iso_list_idx)
        self.iso_profiles = [None] * no_iso
        for i_no_iso in range(no_iso):
            iso_profile = np.zeros(len(self.meas_mass_profile_x))
            isotopologue = isotopologues_list[iso_list_idx[i_no_iso]]
            for i_iso_series in range(isotopologue.iso_list_len):
                iso_mass = chem_data.get_mass(isotopologue.iso_list[i_iso_series])
                sigma = np.interp(iso_mass, self.meas_mass, self.meas_mass_u)
                iso_profile += f_Gauss(self.meas_mass_profile_x, 
                                       isotopologue.iso_list_I_rel[i_iso_series],
                                       iso_mass,
                                       sigma)
            self.iso_profiles[i_no_iso] = iso_profile.copy()


class CompoundNIST(Compound):
    """A class dedicated to NIST data that inherits from Compound
    """
    def __init__(self, one_batch_fragment: list, compound_number: int, LOD: float):
        # call superclass constructor
        super().__init__(one_batch_fragment, compound_number, LOD)
        self.run_type = run_type_NIST
        self.meas_mass_u_cal = np.zeros(len(self.meas_mass))
        self.meas_mass_u_loose = [self.mass_u_k_loose * self.meas_mass_u[mass_idx] for mass_idx in range(self.meas_len)]
        self.mass_u_k = 2.0
        self.mass_axis_offset = 0.0
        self.meas_mass_U_combined = [0.4995 for idx_m in range(self.meas_len)]
        self.meas_mass_min = [self.meas_mass[idx_m] - 0.4995 for idx_m in range(self.meas_len)]
        self.meas_mass_max = [self.meas_mass[idx_m] + 0.4995 for idx_m in range(self.meas_len)]

        self.meas_mass_u_median = 0.0
        self.meas_mass_u_median_ppm = 0.0
        


    def get_mass_range_looser_uncertainty(self, mass_idx: int, other_mass: float=None) -> float:
        return self.meas_mass_u_loose[mass_idx]

    def do_meas_mass_profile_x(self, m_start: float, m_stop: float, ppm_mass) -> None:
        """Initialize self.meas_mass_profile_x
        
        TODO Myriam what does this function do?
        
        INPUT:
        - ``m_start``:
        - ``m_stop``:
        - ``ppm_mass``: TODO Myriam int or float?
        
        OUTPUT: None
        """
        mass_res = float(1.0)
        super().do_meas_mass_profile_x(m_start, m_stop, mass_res)

    def do_meas_mass_profile(self, delta_mass):
        """Do measured mass profile

        INPUT:
        - ``delta_mass``: potential mass axis offset, at initial step it is zero.
          Used to shift position of measured masses by a constant offset.
          delta_mass is optimised by the lmfit algo.
          For unit-mass resolution data delta_mass is zero and therefore unused.
        """
        
        #indexes of measured masses within window range:        
        self.meas_mass_profile_idx = [idx for idx in range(self.meas_len) if self.meas_mass[idx] >=  self.meas_mass_profile_x[0] and  self.meas_mass[idx] <= self.meas_mass_profile_x[-1]]
        self.delta_mass_u = float(0.0)
        self.meas_mass_profile = np.zeros(len(self.meas_mass_profile_x))

        #for idx_m in range(len(self.meas_mass_profile_idx)):
        #    #here could do a while loop
        #    for idx_x in range(len(self.meas_mass_profile_x)):
        #        if self.meas_mass_profile_x[idx_x] == self.meas_mass[self.meas_mass_profile_idx[idx_m]]:
        #            self.meas_mass_profile[idx_x] += self.meas_I[self.meas_mass_profile_idx[idx_m]]
        
        idx_start = 0
        for idx_profile in range(len(self.meas_mass_profile_x)):
            idx = 0 #idx_start
            found = False
            while idx < len(self.meas_mass) and found == False:
                if round(self.meas_mass_profile_x[idx_profile]) == round(self.meas_mass[idx]):
                    self.meas_mass_profile[idx_profile] = self.meas_I[idx]
                    found = True
                    idx_start = idx
                else:
                    idx += 1
            
        #print(str(self.meas_mass_profile_x)  + "  .meas_mass_profile  " + str(self.meas_mass_profile))

    def do_iso_profiles(self, iso_list_idx: list, isotopologues_list: list):
        """
        Generate a list, containing for each row a numpy array
        corresponding to the sum profile of one isotopologue series.
        All isotopologue series are in the same mass window.
        Uni mass spectrum: signal of all same round(mass) added up together.
        
        INPUT:
        -``iso_list_idx``:
        -``isotopologues_list``:
        """
        no_iso = len(iso_list_idx)
        self.iso_profiles = [None] * no_iso
        for i_no_iso in range(no_iso):
            iso_profile = np.zeros(len(self.meas_mass_profile_x))
            isotopologue = isotopologues_list[iso_list_idx[i_no_iso]]
            for i_iso_series in range(isotopologue.iso_list_len):
                iso_mass = chem_data.get_mass(isotopologue.iso_list[i_iso_series])
                #print("iso_mass:  " + str(iso_mass))
                test_m = False
                idx_m = 0
                while idx_m < len(self.meas_mass_profile_x) and test_m == False:
                    if round(iso_mass) == self.meas_mass_profile_x[idx_m]:
                        iso_profile[idx_m] += isotopologue.iso_list_I_rel[i_iso_series]
                        test_m = True
                    idx_m += 1
            self.iso_profiles[i_no_iso] = iso_profile.copy()
            
        #print(".iso_profiles  " + str(self.iso_profiles))
