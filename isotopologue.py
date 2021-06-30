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


Created on Mon Aug 17 11:23:25 2020
"""
import numpy as np
import math
import operator

from compound import Compound, CompoundTofwerkTOF, CompoundNoMassCalibration, CompoundNIST
from runtype import RunType
import periodic_table as chem_data

from utils_identification import generate_all_isotopolog_fragments
from utils_identification import generate_all_isotopolog_fragments_with_relative_proba
from utils_identification import p_iso_rare_compute
from utils_identification import number_rare_isotopes, nb_abundant_isotopes

from utils_data_extraction import f_tofidx_to_mass , interpol_mass_cal, f_pseudoVoigt, f_Gauss

run_type_tofwerk_tof = 0
run_type_NIST = 1
run_type_no_mass_cal = 2
run_type_unit_mass = 3

class Isotopologues():
    """
    This is a class to gather, for each generated abundant fragment,
    all the info needed to then associate each family of fragment
    with the best estimated intensity.
    The aim is to eliminate all solutions that are likely not present,
    because no isotopes have been detected.
    In this step, uncertainties will be associated to each intensity value.
    """
    cl_comp = None
    

    def __init__(self, frag_abund: list, meas_mass_idx: list, cl_comp: Compound):

        """Initialise instance with minimal data

        INPUT:
        - `frag_abund`: encodes a fragment (as a list of integers), contains abundant atoms only.
        - `cl_comp`: instance of class Compound (already initialised)
        - `meas_mass_idx`: index of corresponding measured mass


        iso_rare: list of fragments, all fragments containing rare isotopes,
        being isotopologues of frag_abundant.
        
        I: intensity of frag_abundant.
        At start: I is the measured intensity for the measured mass
        to which frag_abundant is associated.
        u_I: uncertainty of I, 
        at start it is I_u (the standard deviation of data compared to fit).
        """


        #copy the address of the parameters:
        #(so that we dont need to pass all these parameters again later in the methods)
        #this is only an address! if somewherelse the variables are updated that's fine.
        self.set_cl_comp(cl_comp)

        #information about the abundant fragment, here meaning the fragment made of abundant atoms only
        #note: this 'abundant fragment' may not be the one with maximum intensity,
        #depending on the isotopic profile.
        self.frag_abund = frag_abund
        self.frag_abund_idx = meas_mass_idx # list of index of first reference fragment made of only abundant atoms
        self.frag_abund_mass = chem_data.get_mass(self.frag_abund)
        self.frag_abund_formula = chem_data.get_string_formula(self.frag_abund)
        self.DBE = chem_data.get_DBE_value(self.frag_abund)


    @classmethod
    def set_cl_comp(self, cl_comp: Compound):
        self.cl_comp = cl_comp

    def create_isotopologue_suite(self, test_new_variant:bool = True, verbose=False):
        """Generate the isotopologues.
        
        20200416 new: not part of init anymore.
        This part is expensive and done only if the node 
        is not a singleton after initialisation of the graph.
        
        INPUT: None
        OUTPUT: None
        """

        if verbose:
            print("generating isotopocules for {}".format(chem_data.get_string_formula(self.frag_abund)))
        #list of isotopocules
        self.iso_list = []
        #list of exact, theoretical intensities of each isotopocule relative to abundant fragment
        self.iso_list_I_rel = []
        
        #generate all rare isotopologues in case they exist
        frag_iso_list = []
        if number_rare_isotopes(self.frag_abund) == 0 \
        and nb_abundant_isotopes(self.frag_abund) > 0:
            #generate formula of isotopologues
            # improvement: generate isotopic patterns for each abundant atom,
            # compute respective probability, then mix all: combine the patterns,
            # multiply the probabilities
            if self.frag_abund_idx != None:
                LOD_threshold_ratio = float(self.cl_comp.LOD/sum([self.cl_comp.meas_I[self.frag_abund_idx[i]] for i in range(len(self.frag_abund_idx))]))
            else:
                #no idx of mass: this is a candidate molecular ion
                LOD_threshold_ratio = float(1)
            
            #print("LOD_threshold_ratio: " + str(LOD_threshold_ratio))
            arbitrary_choice_of_linear_factor = float(1.1)
            if test_new_variant:
                frag_iso_list, list_I_rel = generate_all_isotopolog_fragments_with_relative_proba(self.frag_abund, min_relative_proba = float(LOD_threshold_ratio/arbitrary_choice_of_linear_factor))
                self.iso_list = frag_iso_list
                self.iso_list_I_rel = list_I_rel
                max_I_rel = max(self.iso_list_I_rel)# should be at index 0
            else:
                frag_iso_list, list_I_rel = generate_all_isotopolog_fragments_with_relative_proba(self.frag_abund)
                #Eliminate fragments expected to be below LOD
                #the best would be to not generate them at all (pruning Loos et al. ?)
                for i_iso in range(len(frag_iso_list)):
                    #print('isp p test: ' + '{:.5f}'.format(self.iso_list_I_rel) + '\t' + '{:.5f}'.format(cl_comp.meas_I[self.frag_abund_idx]))
                    if list_I_rel[i_iso]*arbitrary_choice_of_linear_factor >= LOD_threshold_ratio:
                        self.iso_list.append(frag_iso_list[i_iso])
                        self.iso_list_I_rel.append(list_I_rel[i_iso])
                max_I_rel = max(self.iso_list_I_rel)
            # remove too low isotopologues (so start to compare at the end)
            idx = len(self.iso_list_I_rel)
            while idx > 0 and self.iso_list_I_rel[idx-1] < LOD_threshold_ratio:
                idx -= 1
            if idx < len(self.iso_list_I_rel): # delete isotopologues of relative intensity below LOD_threshold_ratio
                self.iso_list_I_rel = self.iso_list_I_rel[:idx]
                self.iso_list = self.iso_list[:idx]
            # normalise so that max_I_rel == 1
            if max_I_rel > 1.0:
                self.iso_list_I_rel = [val/max_I_rel for val in self.iso_list_I_rel]
        if verbose:
            print("result {} {}".format(", ".join([chem_data.get_string_formula(si) for si in self.iso_list]), ", ".join(["{:.4f}".format(fi) for fi in self.iso_list_I_rel]) ))
        #print('**********')
        #print(self.frag_abund_formula)
        #print('iso_list_I_rel: ' + str(self.iso_list_I_rel))
            
        #length of list of isotopologue presumably above LOD
        self.iso_list_len = len(self.iso_list)
        
        #for each mass of isotopologue, intensity of matching closest measured mass
        #this is the upper bound for assigned intensity
        self.iso_list_I_meas_max = [0.0] * self.iso_list_len


        #effective maximum I: min between measured matching mass and optimised signal based on isotopologue profile optimisation
        self.iso_list_I_opt = [0.0] * self.iso_list_len # PB computed in utils_graph not in this class

        
        #sum of maximum optimised signal for this fragment and all its isotopocules
        self.iso_list_sum_signal = float(0.0) # PB computed in utils_graph not in this class
        
        #calculate the boundaries of the mass domain covered by the isotopologue set:
        #this is used later to define sub-domain of the entire mass domain
        #to run the fit to the measured data, along a sub-domain of the mass axis.
        self.mass_iso_list = [chem_data.get_mass(vec) for vec in self.iso_list]
        self.iso_list_m_min = min(self.mass_iso_list)
        self.iso_list_m_max = max(self.mass_iso_list)

        #**************************************************************
        #now we look for mass index where the mass of each isotopologue
        #can be found, and we associate each iso with this index
        #we do not use LOD so far
        #note 20200416: weakness: allow only one measured mass to be associated with one isotopologue.
        #if mass uncertainties are large, there could be several.
        #20200723: this is not relevant as we are interested in the closest mass peak.
        self.meas_mass_idx = [None]*self.iso_list_len # it will contain a mass index idx so that the mass cl_comp.meas_mass[idx] matches the isotopologue mass
        #self.meas_mass_idx[0] = self.frag_abund_idx
        #self.iso_list_I_meas_max[0] = self.cl_comp.meas_I[self.meas_mass_idx[0]]
        
        #print('ini max signal frag abund: ' + self.frag_abund_formula + '  ' + str(self.cl_comp.meas_I[self.frag_abund_idx]) +';  ' + str(self.frag_abund_idx))

        #self.meas_mass_u = [None]*self.iso_list_len
        #self.meas_mass_u[0] = self.cl_comp.meas_mass_u[self.meas_mass_idx[0]]

        #if self.iso_list_len > 1:
        if self.frag_abund_idx != None:
            #this is not an added candidate molecular ion
            self.find_idx_meas_mass_from_exact_mass()
        
            
            

        #factor to multiply to theoretical isotopocule profile to fit measured profile
        self.k_guess_easy = float(0.0)
        self.k_guess = float(0.0)
        self.k_guesses_min = float(0.0)
        self.k_guesses_max = float(0.0)
        #k_guess but after optimisation
        self.k_opt = float(0.0)
        
        #compute k_guess, factor to multiply to relative iso profile to match measured profile
        #if this solution is taken individually

        self.k_guess_easy = self.compute_k_guess()

    

    def find_idx_meas_mass_from_exact_mass(self):
        """
        Fill in self.iso_list_I_meas_max. List of len of iso_len.
        !!! each line is also a list.
        Find the closest measured mass(es), matching uncertainty criteria, for all isotopocules.
        The idea is that since the intensities have been optimised for one set of isotopocules together,
        the mass uncertainty criteria should be loose. If a true mass is on the side of a measured mass, 
        its optimised intensity is low anyway.
        
        The information is the index of the closest measured mass, saved in self.meas_mass_idx
        Only a unique index per iso mass is returned, the index of the closest mass.
 
        !!! One exact mass may be the true fragment of several measured masses,
        in case the mass list is aggregated from different runs.
        We need to account for that.
        ----------

        """
        #print("=== Assign measured masses")
        #print("=== iso_list_len  " + str(self.iso_list_len))
        for i_iso in range(0, self.iso_list_len):
            #exact mass of isotopologue:
            iso_mass = self.mass_iso_list[i_iso]
            #print("=== iso_mass  " + str(iso_mass))
            #self.frag_abund_mass is used as starting point because we know already its index in the list of measured masses.
        
            #exact mass difference between abundant and rare isotopologues:
            iso_mass_diff = self.frag_abund_mass - iso_mass
            #print('ref mass: ' + str(self.frag_abund_mass) + 'iso_mass_diff  ' + str(iso_mass_diff))
            if iso_mass_diff == 0.0: # this is actually the abundant fragment (the fragment made of abundant atoms)
            
                self.meas_mass_idx[i_iso] = self.frag_abund_idx #this is a list
                #self.iso_list_I_meas_max[i_iso] = sum([self.cl_comp.meas_I[self.meas_mass_idx[i_iso][i]] for i in range(len(self.meas_mass_idx[i_iso]))])
                #print("abundant iso, assigned meas mass " + str(self.cl_comp.meas_mass[self.meas_mass_idx[i_iso]]) + "  " + str(self.meas_mass_idx[i_iso]) + "  " + str(round(self.iso_list_I_meas_max[i_iso])))
                #print()
            else: #iso_mass_diff != 0.0:
                #now loop over all measured masses to find the closest measured mass(es)
                #Does a measured mass match the mass of the candidate isotopologue?
                # start the seach at the index of the fragment of reference
                #the measured mass are already ordered by increasing mass
                
                m_i_list = [i_m for i_m in range(self.cl_comp.meas_len) if abs(self.cl_comp.meas_mass[i_m] - iso_mass) <= self.cl_comp.get_mass_range_looser_uncertainty(i_m, iso_mass)]
                if len(m_i_list) >0:
                    self.meas_mass_idx[i_iso] = m_i_list
                else:
                    self.meas_mass_idx[i_iso] = [None]
                    self.iso_list_I_meas_max[i_iso] = float(0.0)
                    
                """       
                ii = min(self.frag_abund_idx)  
                initial_mass_diff = abs(self.cl_comp.meas_mass[ii] - iso_mass) + 1 #artificially increase to force enter the while loop        
                possible_idx = ii
                while ii < self.cl_comp.meas_len and ii>=0 \
                and abs(self.cl_comp.meas_mass[ii] - iso_mass) <= initial_mass_diff:
                    initial_mass_diff = abs(self.cl_comp.meas_mass[ii] - iso_mass)
                    possible_idx = ii
                    #print("candidate meas mass " + str(self.cl_comp.meas_mass[ii]))
                    if iso_mass_diff < 0.0:
                        ii += 1
                    else:
                        ii -= 1
        
                #now check if the possible match measured mass is within uncertainties:
                if abs(self.cl_comp.meas_mass[possible_idx] - iso_mass) > self.cl_comp.get_mass_range_looser_uncertainty(possible_idx, iso_mass):
                    self.meas_mass_idx[i_iso] = [None]
                    self.iso_list_I_meas_max[i_iso] = float(0.0)
                    #print("No measured mass match")
                else:
                    self.meas_mass_idx[i_iso] = [possible_idx]
                    #self.iso_list_I_meas_max[i_iso] = self.cl_comp.meas_I[self.meas_mass_idx[i_iso][0]] #*1.05
                    
                #if self.meas_mass_idx[i_iso] != None:
                    #print("Assigned meas mass index  " + str(self.meas_mass_idx[i_iso]))
                """
        if self.meas_mass_idx != None:
            self.iso_list_I_meas_max = [sum([self.cl_comp.meas_I[self.meas_mass_idx[i_iso][i]] for i in range(len(self.meas_mass_idx[i_iso])) if self.meas_mass_idx[i_iso] != [None]]) for i_iso in range(self.iso_list_len)]
        else:
            #this is a candidate molecular ion
            self.iso_list_I_meas_max = [float(0) for i in range(self.iso_list_len)]
        #print('max signal frag abund: ' + self.frag_abund_formula + '  ' + str(self.iso_list_I_meas_max) )
        #print('idx frag abund: ' + str(self.meas_mass_idx[0]))
    
    def compute_k_guess(self) -> float:
        #compute max intensity for each isotopocule
        #if no meas mass is assigned, take LOD value as max
        #then compute Imeas/ Irel for each iso
        if self.iso_list_len != len(self.iso_list_I_meas_max) or self.iso_list_len != len(self.iso_list_I_rel):
            print("self.iso_list_len = {} len(self.iso_list_I_meas_max) = {} len(self.iso_list_I_rel) = {}".format(self.iso_list_len, len(self.iso_list_I_meas_max), len(self.iso_list_I_rel)))
        return min([max(self.cl_comp.LOD, self.iso_list_I_meas_max[i_iso])/self.iso_list_I_rel[i_iso] for i_iso in range(self.iso_list_len)])

    def init_k_guess(self, x_profile):
        """
        Compute a first rough estimate of k_guess, taking each individual set of isotopologues
        and fitting the profile to the measured profile, 
        on the mass domain covered by the candidate isotopologue.
        

        Parameters
        ----------
        x_profile : np array
            DESCRIPTION.
        measured_profile : np array
            DESCRIPTION.

        Returns
        -------
        None.

        """
        #x_profile_mass_interval
        
        #for each x abcisse, compute the theoretical profile of the set of isotopologue        
        iso_profile = self.compute_individual_iso_profile(x_profile)
        meas_profile = self.compute_meas_profile(x_profile)
        if self.cl_comp.run_type == run_type_NIST or self.cl_comp.run_type == run_type_unit_mass:
            iso_profile = np.where(iso_profile < 1.0, 1.0, iso_profile)
        self.k_guess = np.mean(meas_profile/iso_profile)
        #print("iso_profile" + str(iso_profile))
        #print("meas_profile" + str(meas_profile))
        #print("init k_guess l. 290  " + str(self.k_guess))



    def compute_individual_iso_profile(self, x_profile):
        
        #min and max mass provided by min and max masses of enumerated isotopologues
        iso_profile = np.zeros(len(x_profile))
        
        if self.cl_comp.run_type == run_type_NIST or self.cl_comp.run_type == run_type_unit_mass:
            for i_iso in range(self.iso_list_len):           
                test_m = False
                idx_m = 0
                while idx_m < len(x_profile) and test_m == False:
                    if round(self.mass_iso_list[i_iso]) == x_profile[idx_m]:
                        iso_profile[idx_m] += self.iso_list_I_rel[i_iso]
                        test_m = True
                    idx_m += 1

        elif self.cl_comp.run_type == run_type_tofwerk_tof:
            for i_iso in range(self.iso_list_len):
                #f_pseudoVoigt(x, A, mu, sigma, alpha)
                sigma = self.mass_iso_list[i_iso] * self.cl_comp.sigma_m_slope + self.cl_comp.sigma_m_b
                iso_profile += f_pseudoVoigt(x_profile,
                                             self.iso_list_I_rel[i_iso],
                                             self.mass_iso_list[i_iso],
                                             sigma,
                                             self.cl_comp.alpha_pseudo_voigt)
        return iso_profile
        
    def compute_meas_profile(self, x_profile):
        """
        This is a special case, it generates the measured profile using peaks matching the candidate,
        not using all measured peaks. This is to compute k_guess without overestimation:
        a measured peak not explain by a candidate is likely explained by a different candidate.       

        Parameters
        ----------
        x_profile : np array.
            DESCRIPTION. x values (mass) of abcisse where to compute the signal.

        Returns
        -------
        meas_profile : np array
            DESCRIPTION. Computed signal value at each abcisse x of the x_profile.

        """

        
        #min and max mass provided by min and max masses of enumerated isotopologues
        meas_profile = np.zeros(len(x_profile))
        for i_iso in range(self.iso_list_len):
            #f_pseudoVoigt(x, A, mu, sigma, alpha)
            for i_m in range(len(self.meas_mass_idx[i_iso])):
                if self.meas_mass_idx[i_iso][i_m] != None:
                                  
                    if self.cl_comp.run_type == run_type_tofwerk_tof:
                        sigma = self.cl_comp.meas_mass[self.meas_mass_idx[i_iso][i_m]] * self.cl_comp.sigma_m_slope + self.cl_comp.sigma_m_b
                        meas_profile += f_pseudoVoigt(x_profile,
                                                     self.cl_comp.meas_I[self.meas_mass_idx[i_iso][i_m]],
                                                     self.cl_comp.meas_mass[self.meas_mass_idx[i_iso][i_m]],
                                                     sigma,
                                                     self.cl_comp.alpha_pseudo_voigt)
    
                    elif self.cl_comp.run_type == run_type_NIST or self.cl_comp.run_type == run_type_unit_mass:
                        test_m = False
                        idx_m = 0
                        while idx_m < len(x_profile) and test_m == False:
                            if round(self.mass_iso_list[i_iso]) == x_profile[idx_m]:
                                meas_profile[idx_m] += self.cl_comp.meas_I[self.meas_mass_idx[i_iso][i_m]]
                                test_m = True
                            idx_m += 1
                    
        return meas_profile

    def update_k_opt(self, k_opt: float) -> None:
        """Set k_opt and update iso_list_I_opt and iso_list_sum_signal """
        self.k_opt = k_opt
        self.iso_list_I_opt = [min(self.iso_list_I_rel[i_series] * self.k_opt, self.iso_list_I_meas_max[i_series]) for i_series in range(self.iso_list_len)]
        self.iso_list_sum_signal = sum(self.iso_list_I_opt)

    def initialise_k_guess(self, dict_k_ppm_mass: float, cl_comp: Compound, run_type: RunType) -> None:
        """moved here from utils_graph"""



        #create x axis
        w_start = round(self.iso_list_m_min) #- cl_comp.mass_axis_offset
        w_stop  = round(self.iso_list_m_max) #+ cl_comp.mass_axis_offset
        
        cl_comp.initialise_mass_profile()
        #meas_mass_profile_x = []
        #delta_mass_u = float(0.0)
        #meas_mass_profile_idx = []
        #meas_mass_profile = []

        cl_comp.do_meas_mass_profile_x(w_start, w_stop, dict_k_ppm_mass)
        #meas_mass_profile_x = cl_comp.get_meas_mass_profile_x(w_start, w_stop, dict_k_ppm_mass)
        #compute measured mass profile at each x value on the x axis (mass axis)

        delta_mass = 0.0
        cl_comp.do_meas_mass_profile(delta_mass)
        #calculate measured sum signal on this portion of the mass axis
        #sum_signal_meas = sum([cl_comp.meas_I[i_m] for i_m in cl_comp.meas_mass_profile_idx])
        #if sum_signal_meas == 0:
        #    print(str(i_node) + '\t' + 'sum signal zero')



        #self.k_guess = 0
        if len(cl_comp.meas_mass_profile_idx) > 0 and len(cl_comp.meas_mass_profile_x) >= 2:
            #print("Compute first k_guess")
            #at least one measured mass within mass axis window

            self.init_k_guess(cl_comp.meas_mass_profile_x)
            self.update_k_opt(self.k_guess)
            #print()
            
            #print(" first k_guess:  " + str(cl_comp.k_guesses) + "  " + str(cl_comp.iso_profiles))
            #k_opt, delta_mass_opt = fit_iso_profile(cl_comp, delta_mass, run_type, [0], cl_comp.k_guesses, cl_comp.iso_profiles, fit_mode, graph = plot_mode)

            #print("=== k_guess_easy, k_guess, k_opt: " + "{:.6f}".format(self.k_guess_easy) + "\t" + "{:.6f}".format(self.k_guess) + "\t" +"{:.6f}".format(k_opt[0]))
            #self.k_guess = k_opt[0] #there is only one self to fit
            #self.k_opt = k_opt[0]
            
            #print("iso_list_sum_signal after k_guess:  " + str(self.iso_list_sum_signal))
        elif len(cl_comp.meas_mass_profile_idx) == 1 and len(cl_comp.meas_mass_profile_x) == 1:
            if run_type.is_NIST() or run_type.is_unit_mass():
                self.k_guess = sum(self.iso_list_I_meas_max) / sum(self.iso_list_I_rel)
                self.update_k_opt(self.k_guess)
            else:
                self.k_guess = 0
                self.k_opt = 0
                self.iso_list_sum_signal = 0
        else:
            #print("Not enough data point, k_guess = 0")
            self.k_guess = 0
            self.k_opt = 0
            self.iso_list_sum_signal = 0
        #print('-.-.-.-.-.-.-.-.-.-.-.-.-.-')
        #print(self.k_guess)
        #print(self.iso_list_I_opt)
        #print([self.iso_list_I_rel[i] * self.k_guess for i in range(self.iso_list_len)])
        #print(self.iso_list_sum_signal)

    def __repr__(self):
        s = "["
        for frag in self.iso_list:
            s += chem_data.get_string_formula(frag) + ", "
        s += "]"
        return s
