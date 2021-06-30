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


Created on Wed Aug 14 10:33:08 2019
"""

#All needed constants for identification

#Default list of potentially present atoms:
default_chem_elements = 'HCNOSFClBrI'

#Coverage factor for the mass uncertainty, 
#mutiply this to the mass uncertainty.
#unit: no unit
#See equation (1) in Guillevic et al.
m_u_k = float(2.5)

#Smallest detectable mass with Empa Aprecon-GC-ToF. 
#This is because we use a mass filter to prevent water from entering our detector.
#Change this value depending on your instrumentation.
#Unit: [m/z]
#This value will influence the computation of the likelihood estimator
#at Step 5
min_detect_mass = float(23.0)

#Value of the smallest detectable peak,
#this is the instrumental limit of detection
#Sets of isotopologues whose *total* signal is below this value will be eliminated
#at the end of Step 7
#(this is a conservative approach)
LOD = float(80.0)

#Percentage of measured signal to fit
#Ex.: 95 means that the sotfare will aim at reproducing 95% of the signal,
#and will stop the optimisation once this value is reached.
#Choose this value according to the measurement noise of your instrument,
#e.g. choose 90.0 (100% - 10%) if your noise level is around 10%.
#Or, choose a lower value to speed up the identification.
#If you wish to fit all measured mass peak, set this to e.g. 99%,
#but better get a correct list of present chemical elements first.
#threshold_id_sum_signal = 95.0
max_opt_signal = float(95.0)

#Fraction (or ratio) of measured masses to not fit.
#This parameter works in a similar way as "threshold_id_sum_signal".
#Choose a value according to how confident you are that each mass peak is "real",
#and not caused by e.g. electronic noise.
#If you suspect that some peaks are due to ringing,
#choose a value above zero, e.g. 0.1 means 10% of peaks may not be real.
#if you are sure that all peaks are real, set it to 0.0001.
#This will force the software to fit all peaks,
#but bette do this once you have identified the chemical elements present in the molecule.
max_ratio_not_opt_masses = float(0.100)


#This is a "safety" parameter to prevent the optimisation
#from staying too long in the loop, Steps XX to XX.
#Once this number of nodes have been considered for optimisation,
#the optimisation is stopped.
#You can increase this number of you have many measured peaks
#for the same molecule.
max_added_nodes = int(300) 

#This is the mass interval in between each "x" point on the x-axis.
#This value is used to compute theoretical isotopic profiles
#along the "x" mass axis, Step 4 and Step 7.
#We use 8.00 for a mass resolution of approx. 3500.
#Decrease this value if you have a better mass resolution, 
#computing time may increase.
#unit: ppm of mass.
ppm_mass_res = float(8.0)


#We have experienced that the software does not produce reliable results
#if very few masses are measured. If less masses than this number are measured,
#solutions are splited into groups of solutions belonging to
#the same maximal fragment, and the optimisation is run separately for
#each group. Multiple results may be produced.
#In such a case, more than 100% of the signal may be reconstructed.
min_no_masses = int(6)
    

