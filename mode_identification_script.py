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


Created on Mon Aug 24 09:30:00 2020

"""
import matplotlib.pyplot as plt
from pathlib import Path

#from MyClassRunParameters import PathParameters
from io_tools import make_dict_of_paths
from runtype import RunType
import periodic_table as chem_data
from mode_identification import make_identification
from plot_figures import plot_method_behaviour
plt.close('all')

#cl_path = PathParameters()
dict_path = make_dict_of_paths()
#**************************************
#LOAD EMPA TOF RESULTS

#200714.2200.cal6L.7.frag.1551
#200714.2200.cal6L.7.frag.1509_Oin_Thorin_test2
#toffile_list = list(dict_path['fragments'].glob('201101.0024.labair.11.frag.scan2700_unit_mass.txt'))

#toffile_list = list(dict_path['fragments'].glob('[1,2]*.frag.val_*.txt'))
#toffile_list1 = list(dict_path['fragments'].glob('[1,2]*.frag.val_*.txt'))
#toffile_list = toffile_list0 + toffile_list1

toffile_list = list(dict_path['fragments'].glob('[1,2]*.frag.train_*.txt'))
#for file in toffile_list:
#    print(file.name)
#toffile_list = list(dict_path['fragments'].glob('210623.2029.tank.11.frag.search_flurane_1386s.txt'))

#toffile_list = list(dict_path['fragments'].glob('200714.2200.cal6L.7.frag.train_NF3_.txt'))


#toffile_list = list(dict_path['fragments'].glob('200128.1735.tank.5.frag.1490s_MrMiralago.txt'))
#toffile_list = list(dict_path['fragments'].glob('200224.0907.tank.5.frag.val_PFPHP.txt'))

#COS: 463-58-1-Mass.jdx
#toffile_list = list(dict_path['NIST_EI'].glob('107-05-1-Mass.jdx'))

if len(toffile_list) == 0:
    raise ValueError('Error mode_identification_script.py line 28.\nNo file found. Check that path is correct and that filename is correct.')

#write results of no nodes and runtime in results file
method_performance_results = []
method_performance_file = dict_path['fragments'] / 'method_performance.txt'

list_likelihood_true_frag = []
list_likelihood_false_frag = []
list_likelihood_true_att = []
list_likelihood_false_att = []
list_ranking_true_frag = []
list_ranking_false_frag = []
list_ranking_true_att = []
list_ranking_false_att = []

for path_file in toffile_list:
    #print("Processing file {}".format(str(path_file)))
    filename = path_file.name
    run_type = RunType()
    fragfile_name = None
    if "Mass" in filename:
        run_type.set_NIST()
    elif "no_mass_cal" in filename:
        run_type.set_no_mass_cal()
    elif "unit_mass" in filename:
        run_type.set_unit_mass()
    else:        
        run_type.set_tofwerk_tof()
        toffile_name = path_file.stem # the filename without extension (.txt)
        print("toffile_name = {}".format(toffile_name))
        if '.frag.' in toffile_name:
            fragfile_name = toffile_name.split('.')[5]
            print(fragfile_name)
        #if '.val_' in toffile_name:
        #    val_name = toffile_name.split('.')[5]
        #    print(val_name)
    if 'train_' in filename or 'val_' in filename or 'known_' in filename:
        do_val = True
    else:
        do_val = False
        
    
    # go to one parent directory to find the directory "mass_spectra" where to save figures:
    fig_path = path_file.parents[1]/"mass_spectra"
    formulas_path = path_file.parents[1]/"formulas"

    chosen_list_frag_idx = [0] # index of chosen compound, zero by default
    min_frag_idx = None #15
    target_elements = None #list of chemical elements potentially present
    target_formula = None #"CCl2F2" #formula of molecular ion, if present
    # choose one of:
    #make_identification(run_type, path_file, fig_path, formulas_path, fragment_indices=chosen_list_frag_idx, fragfile_name=fragfile_name)
    # or
    
        
    
    method_performance_results_i, fragment_results, validation_dict_results = \
    make_identification(run_type,
                        path_file,
                        fig_path,
                        formulas_path,
                        fragment_indices =chosen_list_frag_idx,
                        min_fragment_index=min_frag_idx,
                        target_elements = target_elements,
                        target_formula = target_formula,
                        fragfile_name=fragfile_name,
                        show_plt_figures = False)
    

    if do_val:
        #method_performance_results.append(method_performance_results_i)
        method_performance_results.append(method_performance_results_i +
                                          [validation_dict_results["frac_correct_frag"]] +
                                          [validation_dict_results["frac_true_signal_to_total_signal"]] +
                                          [validation_dict_results["frac_true_frag_top10"]] +
                                          [validation_dict_results["frac_true_signal_top10"]])
    
    
        #print(method_performance_results)
        list_likelihood_true_frag = list_likelihood_true_frag + validation_dict_results["list_likelihood_true_frag"]
        list_likelihood_false_frag = list_likelihood_false_frag + validation_dict_results["list_likelihood_false_frag"]
        list_likelihood_true_att = list_likelihood_true_att + validation_dict_results["list_likelihood_true_att"]
        list_likelihood_false_att = list_likelihood_false_att + validation_dict_results["list_likelihood_false_att"]
        list_ranking_true_frag = list_ranking_true_frag + validation_dict_results["list_ranking_true_frag"]
        list_ranking_false_frag = list_ranking_false_frag + validation_dict_results["list_ranking_false_frag"]
        list_ranking_true_att = list_ranking_true_att + validation_dict_results["list_ranking_true_att"]
        list_ranking_false_att = list_ranking_false_att + validation_dict_results["list_ranking_false_att"]
        

        """
        #Extract results for article
        frag_result_file_name = 'frag_res.' + path_file.stem + '.txt'
        frag_result_file = dict_path['fragments'] / frag_result_file_name
        with frag_result_file.open('w') as file:
            file.write('formula'+ '\t' +
                       'correct'+ '\t' +
                       'mass'+ '\t' +
                       'iso contrib. to meas. signal'+ '\t' +
                       'att contrib. to meas. signal'+ '\t' +
                       'likelihood'+ '\t' +
                       'rank'+ '\t' +
                       'max. el.'+ '\n')
            for i in range(len(fragment_results)):
                file.write(str(fragment_results[i][0]) + '\t' +
                           str(fragment_results[i][1]) + '\t' +
                           '{:.5f}'.format(fragment_results[i][2]) + '\t' +
                           '{:.1f}'.format(fragment_results[i][3]) + '\t' +
                           '{:.1f}'.format(fragment_results[i][4]) + '\t' +
                           '{:.1f}'.format(fragment_results[i][5]) + '\t' +
                           str(fragment_results[i][6]) + '\t' +
                           str(fragment_results[i][7]) + '\n')
        """

    
    

#plots for article   
plot_method_behaviour(method_performance_results, list_likelihood_true_frag, list_likelihood_false_frag, list_likelihood_true_att, list_likelihood_false_att, list_ranking_true_frag, list_ranking_false_frag, list_ranking_true_att, list_ranking_false_att, fig_path)


i_met_res_filename = 0
i_met_res_no_meas_mass = 1
i_met_res_sum_I = 2                        #total measured signal, Volt
i_met_res_meas_mass_u_median = 3       
i_met_res_meas_mass_u_median_ppm = 4              
i_met_res_sum_fitted_signal = 5
i_met_res_no_nodes_knapsack = 6
i_met_res_percent_nodes_LOD = 7
i_met_res_percent_nodes_singletons = 8
i_met_res_percent_nodes_not_treated = 9
i_met_res_percent_nodes_validated = 10
i_met_res_runtime_step1_knapsack = 11
i_met_res_runtime_step2_graph = 12
i_met_res_runtime_step3_enum_iso = 13
i_met_res_runtime_step7_opt_loop = 14
i_met_res_runtime_total = 15
i_met_res_mol_ion = 16



#export_met_res_filename = [method_performance_results[i][i_met_res_filename] for i in range(len(method_performance_results))]
#export_met_res_runtime_step1_knapsack = [method_performance_results[i][i_met_res_runtime_step1_knapsack] for i in range(len(method_performance_results))]
#export_met_res_runtime_step2_graph = [method_performance_results[i][i_met_res_runtime_step2_graph] for i in range(len(method_performance_results))]
#export_met_res_runtime_step3_enum_iso = [method_performance_results[i][i_met_res_runtime_step3_enum_iso] for i in range(len(method_performance_results))]

filenameout = dict_path['formulas'] /"runtimes.txt"
fout = filenameout.open('w')

fout.write("Compound" + "\t" + "&" + "\t" +
           "No. knapsack formulae" + "\t" + "&" + "\t" +
           "Step 1: Knapsack" + "\t" + "&" + "\t" +
           "Step 2: Graph" + "\t" + "&" + "\t" +
           "Step 3: Iso. Enum" + "\t" + "&" + "\t" +
           "Step 7: Optimisation" + "\t" + "&" + "\t" +
           "Total " + "\\" + "\\" + "\n")
for i in range(len(method_performance_results)):
    fout.write(method_performance_results[i][i_met_res_filename] + "\t" + "&" + "\t" +
        str(method_performance_results[i][i_met_res_no_nodes_knapsack]) + "\t" + "&" + "\t" +               
        "{:.5f}".format(method_performance_results[i][i_met_res_runtime_step1_knapsack]) + "\t" + "&" + "\t" +
               "{:.5f}".format(method_performance_results[i][i_met_res_runtime_step2_graph]) + "\t" + "&" + "\t" +
               "{:.5f}".format(method_performance_results[i][i_met_res_runtime_step3_enum_iso]) + "\t" + "&" + "\t" +
               "{:.5f}".format(method_performance_results[i][i_met_res_runtime_step7_opt_loop]) + "\t" + "&" + "\t" +
               "{:.5f}".format(method_performance_results[i][i_met_res_runtime_total]) + "\t" + "&" + "\t" +
               str(method_performance_results[i][i_met_res_mol_ion]) + "  " + "\\" + "\\" + "\n")


fout.close()    

"""

i_met_res_frac_true_frag = 16
i_met_res_frac_true_signal_to_total = 17
i_met_res_frac_true_frag_top10 = 18
i_met_res_frac_true_signal_top10 = 19



with method_performance_file.open('w') as file:
    file.write('file ' + '\t' +
               'no_meas_mass'+ '\t' +
               'sum_fitted_signal'+ '\t' +
               'no_nodes_knapsack'+ '\t' +
               'percent_nodes_LOD'+ '\t' +
               'percent_nodes_singletons'+ '\t' +
               'percent_nodes_not_treated'+ '\t' +
               'percent_nodes_validated'+ '\t' +
               'frac_true_frag'+ '\t' +
               'frac_true_signal_to_total'+ '\t' +
               'top10_frac_true_frag_to_total' + '\t' + 
               'top10_frac_true_signal_to_total' + '\n')
    for i in range(len(method_performance_results)):
        #print(method_performance_results[i])
        file.write(method_performance_results[i][i_met_res_filename] + '\t' +
                   str(method_performance_results[i][i_met_res_no_meas_mass]) + '\t' +
                   '{:.2f}'.format(method_performance_results[i][i_met_res_sum_fitted_signal]) + '\t' +
                   str(method_performance_results[i][i_met_res_no_nodes_knapsack]) + '\t' +
                   '{:.2f}'.format(method_performance_results[i][i_met_res_percent_nodes_LOD]) + '\t' +
                   '{:.2f}'.format(method_performance_results[i][i_met_res_percent_nodes_singletons]) + '\t' +
                   '{:.2f}'.format(method_performance_results[i][i_met_res_percent_nodes_not_treated]) + '\t' +
                   '{:.2f}'.format(method_performance_results[i][i_met_res_percent_nodes_validated]) + '\t' +
                   '{:.2f}'.format(method_performance_results[i][i_met_res_frac_true_frag]) + '\t' +
                   '{:.2f}'.format(method_performance_results[i][i_met_res_frac_true_signal_to_total]) + '\t' +
                   '{:.2f}'.format(method_performance_results[i][i_met_res_frac_true_frag_top10]) + '\t' +
                   '{:.2f}'.format(method_performance_results[i][i_met_res_frac_true_signal_top10]) + 
                   '\n')

"""
