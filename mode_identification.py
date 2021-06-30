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


Created on Mon Feb 25 11:15:46 2019

Non-target screening
Use knapsack algorithm
Based on list of exact masses as measured by our TOF,
to produce a list of potential fragment formula,
within mass uncertainty.

Strategy for uncertainty: unit is same as the quantity it is related to (no percent, no ppm).
The uncertainty is also always 1 sigma (or k=1), because to sum up, values with k=1 should be used.
"""
import sys
import math
import numpy as np
import time
import operator
import matplotlib.pyplot as plt
#package needed to represent fragments as graphs:
import networkx as nx
from pathlib import Path
#https://docs.python.org/fr/3/library/pathlib.html#module-pathlib

import periodic_table as chem_data

from utils_graph import graph_optimise_isotopologues, graph_order_by_att_likelihood, graph_percent_sum_signal
from utils_graph import add_nodes_and_edges_to_graph_from_list_fragments, remove_singletons, set_str_formula_and_attribute_for_pydot
from utils_graph import compute_all_subfragments_of_nodes, compute_all_subfragments_of_nodes_with_lists
from utils_graph import initialise_isotopologue_set, initialise_isotopologue_k_guess
from utils_graph import get_list_identified_unidentified_mass
from utils_graph import remove_node_update_edges
from utils_graph import add_maximal_node_and_edges_to_graph_from_max_fragment
from utils_identification import define_knapsack_target, knapsack_double, binary_search_nearest_low_sorted_list_of_tuples, knapsack_double_with_lists

import const_identification as const_id

#=================================
#***Load Metadata***
from runtype import RunType
from io_tools import get_data_from_frag_file, get_data_from_NIST_file, get_data_from_frag_unit_mass_file

from compound import CompoundTofwerkTOF, CompoundNoMassCalibration, CompoundNIST

from plot_figures import plot_mass_defect_of_fragments, plot_graph_of_frag, plot_mass_spectrum
from plot_figures import plot_forest_of_frag_with_pydot

#--------------------------------------------
#for validation step
# import validation sub-module
import validation
from validation.validation_utils import validation, validate_frag_list, compute_validation_results, find_validation_idx_from_database
from validation.mol_validation import mol_validation
from validation.mol_validation import dic_labels as validation_labels
#--------------------------------------------

measured_masses_enough = 0
measured_masses_few = 1
measured_masses_many = 2

fit_mode_fixed_x_axis = 0
fit_mode_variable_x_axis = 1

def make_identification(run_type: RunType,
                        path_file: Path,
                        fig_path: Path,
                        formulas_path: Path,
                        fragment_indices: list=None,
                        min_fragment_index: int=None,
                        target_elements: str=None,
                        target_formula: str=None,
                        fragfile_name: str=None,
                        show_plt_figures: bool=False,
                        verbose_knapsack: bool=False)-> None:
    """The main function of non-target identification of fragments

    INPUT:
    - ``run_type``: run type is tofwerk_tof, NIST or no_mass_cal
    - ``path_file``: a Path filename to a data file of masses (see
      io_tools.py for format)
    - ``fig_path``: a Path filename where to save the figure of
      identified and non-identified fragments
    - ``formulas_path``: a Path filename where to save the text output
    - ``fragment_indices``: a subset of idx_comp_bin to consider; if a
      number is not valid, it is ignored
    - ``min_fragment_index``: minimum fragment index where to start
    - ``target_elements``: a string of target elements, possibly present,
      e.g. "CHOCl"
    - ``target_formula``:  the string of the target formula, e.g. "C6Br6"
      for NIST data
    - ``fragfile_name``: the part of filename with the fragment name,
      such as val_CCl4, toluene
    - ``show_plt_figures``: to show the figures. If False, the figures
      will be saved in fig_path
    - ``verbose_knapsack``: debug of the knapsack, print data
    """

    t0_total = time.time()

    print("processing file {}".format(str(path_file)))
    if fragfile_name is not None:
        print('Fragment: ' + fragfile_name)

    if run_type.is_tofwerk_tof() or run_type.is_no_mass_cal():
        batch_fragments = get_data_from_frag_file(path_file)
        nist_metadata = None
    elif run_type.is_NIST():
        batch_fragments, nist_metadata = get_data_from_NIST_file(path_file)
    elif run_type.is_unit_mass():
        batch_fragments = get_data_from_frag_unit_mass_file(path_file)
        nist_metadata = None


    # get the keys of the dictionary that are the fragment indices:
    list_frag_idx = [i for i in batch_fragments if i != -1]

    toffile_name = path_file.stem
    # set fig_file to some valid Path filename to save the figure

    """
    if run_type.is_tofwerk_tof() and show_plt_figures:
        try:
            plot_mass_defect_of_fragments(batch_fragments, toffile_name, fig_file=None, show=True)
        except:
            print("No figure of mass defect.")
    """
    const_id_max_opt_signal = const_id.max_opt_signal
    const_id_max_ratio_not_opt_masses = const_id.max_ratio_not_opt_masses
    #if fragfile_name is not None and 'known_' in fragfile_name:
    if fragfile_name is not None and 'known_' in fragfile_name:
        #use formula from validation database
        idx_val = find_validation_idx_from_database(fragfile_name, mol_validation, validation_labels)
        target_elements = None
        if idx_val is not None:
            target_formula = mol_validation[idx_val][validation_labels["Sum formula"]]
        else:
            target_formula = None
        const_id_max_opt_signal = min(const_id.max_opt_signal, 99.0)
        const_id_max_ratio_not_opt_masses = min(const_id.max_ratio_not_opt_masses, 0.00)



    #for test purpose: select one or several index of co-eluting compounds
    if (run_type.is_tofwerk_tof() or run_type.is_no_mass_cal() or run_type.is_unit_mass()) and fragment_indices is not None:
        print("using fragment_indices subset = {}".format(fragment_indices))
        list_frag_idx = [i for i in fragment_indices if i in list_frag_idx]
    elif (run_type.is_tofwerk_tof() or run_type.is_no_mass_cal() or run_type.is_unit_mass()) and min_fragment_index is not None:
        print("using minimal fragment index {}".format(min_fragment_index))
        list_frag_idx = [i for i in range(min_fragment_index, max(list_frag_idx)+1) if i in list_frag_idx]

    if len(list_frag_idx) == 0:
        raise ValueError('Chosen compound index does not exist.')



    method_performance_results = []
    for compound_number in list_frag_idx:
        if run_type.is_tofwerk_tof():
            print("==== Initialise Compound")
            cl_comp = CompoundTofwerkTOF(batch_fragments[compound_number], compound_number, const_id.m_u_k, const_id.LOD)
            const_id_min_detect_mass = min(min(cl_comp.meas_mass), const_id.min_detect_mass)

        elif run_type.is_NIST():
            cl_comp = CompoundNIST(batch_fragments[compound_number], compound_number, const_id.LOD)
            print('***LOD***: ' + str(cl_comp.LOD))
            const_id_min_detect_mass = min(10.0, min(cl_comp.meas_mass)) #24 if ToF APRECON
            const_id_max_ratio_not_opt_masses = 0.0 #force to iterate through all measured masses

        elif run_type.is_unit_mass():
            cl_comp = CompoundNIST(batch_fragments[compound_number], compound_number, const_id.LOD)
            print('***LOD***: ' + str(cl_comp.LOD))
            const_id_min_detect_mass = min(min(cl_comp.meas_mass), const_id.min_detect_mass)
        elif run_type.is_no_mass_cal():
            cl_comp = CompoundNoMassCalibration(batch_fragments[compound_number], compound_number, const_id.LOD)
            const_id_min_detect_mass = min(min(cl_comp.meas_mass), const_id.min_detect_mass)
        else:
            raise ValueError("run_type = {} not valid".format(run_type))

        if run_type.is_NIST():
            target_formula = str(nist_metadata.NIST_target_formula)
            target_elements = None
        else:
            print("target_elements " + str(target_elements))
            print("target_formula "+ str(target_formula))
            if target_elements is None and target_formula is None:
                target_elements = const_id.default_chem_elements
                print("target_elements " + str(target_elements))



        idx_list_kn, max_no_each_element = define_knapsack_target(target_elements, target_formula)
        #by construction, idx_list_kn contains elements ordered by increasing mass.
        print("idx_list_kn = {} = {}".format(idx_list_kn, [chem_data.element[i] for i in idx_list_kn]))

        #20200131: here, before the knpasack,
        #it would be good to index as 'not used'
        #masses that are impossible
        #and result from ringing of the machine.
        #this could be done at least for relatively small fragments,
        #where we know a mass domain where nothing is possible.


        #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        # ALPINAC: STEP 1: KNAPSACK (cf. Fig. 1 in Guillevic et al. xx)      +
        # Generate all fragment formulae matching the measured masses.       +
        #Use most abundant isotopes only.                                    +
        #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        print('&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&')
        print('*** Fragment index: '+ str(cl_comp.fragment_index) + '***')

        if verbose_knapsack:
            print("LOD = {}".format(cl_comp.LOD))
            print("masses at input of knapsack with uncertainty range and intensity:")
            print("mass           massmin        massmax        Intensity")
            for id_m in range(cl_comp.meas_len):
                print("{:12.8f}   {:12.8f}   {:12.8f}   {:8.2f}".format(cl_comp.meas_mass[id_m], cl_comp.meas_mass_min[id_m], cl_comp.meas_mass_max[id_m], cl_comp.meas_I[id_m]))
            for id_m in range(cl_comp.meas_len):
                print("meas_mass = {}\nmeas_mass_max = {}\nmeas_mass_min = {}".format(cl_comp.meas_mass[id_m], cl_comp.meas_mass_max[id_m], cl_comp.meas_mass_min[id_m]))
                print("idx_list_kn = {}".format(idx_list_kn))
                print("number_solutions{0}, valid_DBE_solutions{0}, list_solutions{0} = knapsack_double({1}, {2}, {3}, {4}, None, True)".format(id_m, cl_comp.meas_mass[id_m], cl_comp.meas_mass_max[id_m], cl_comp.meas_mass_min[id_m], idx_list_kn))
        print('***Knapsack: Generating fragments***')
        t0 = time.time()
        t0_knapsack = time.time()
        #find formula for abundant isotopologues for the measured masses:
        #this is the knapsack part
        all_solutions = [None]*cl_comp.meas_len
        idx_m = cl_comp.meas_len-1
        total_solutions, total_DBE_solutions, all_solutions[idx_m], list_multi, list_mono = knapsack_double(cl_comp.meas_mass[idx_m], cl_comp.meas_mass_max[idx_m], cl_comp.meas_mass_min[idx_m], idx_list_kn, max_no_each_element, return_lists=True, verbose=verbose_knapsack)
        for idx_m in range(cl_comp.meas_len-2, -1, -1):
            # filter the two lists
            if len(list_multi) > 0:
                idx_multi = binary_search_nearest_low_sorted_list_of_tuples(list_multi, cl_comp.meas_mass_max[idx_m], 0)
                list_multi = list_multi[:idx_multi+1]
            if len(list_mono) > 0:
                idx_mono = binary_search_nearest_low_sorted_list_of_tuples(list_mono, cl_comp.meas_mass_max[idx_m], 0)
                list_mono = list_mono[:idx_mono+1]

            number_solutions, valid_DBE_solutions, all_solutions[idx_m] = knapsack_double_with_lists(cl_comp.meas_mass_max[idx_m], cl_comp.meas_mass_min[idx_m], idx_list_kn, list_multi, list_mono, verbose=verbose_knapsack, double_check_DBE=False)
            total_solutions += number_solutions
            total_DBE_solutions += valid_DBE_solutions
        t1_knapsack = time.time()
        runtime_step1_knapsack = t1_knapsack - t0_knapsack
        print(str(total_DBE_solutions) + ' solutions for ' + str(cl_comp.meas_len) +' measured masses in ' + '{:.2f}'.format(runtime_step1_knapsack) + 's: ')
        #+ str([len(v) for v in all_solutions]))
        # print mass and number of solutions when it is not zero
        #print("\n".join(["{}: {} {}".format(cl_comp.meas_mass[i], len(all_solutions[i]), [chem_data.get_string_formula(frag) for frag in all_solutions[i]]) for i in range(len(all_solutions)) if len(all_solutions[i]) != 0]))
        #raise ValueError("test")

        non_identified_mass = [cl_comp.meas_mass[i] for i in range(len(all_solutions)) if len(all_solutions[i]) == 0]
        print("non-identified masses, abundant-isotope formulae only: {}".format(len(non_identified_mass)))
        #Note: a mass not identified at end of STEP 1 may be cause by a rare-isotope fragment, will be generated at STEP 3.

        if total_DBE_solutions == 0:
            #print("There were no knapsack solutions. Consider adding other elements.")
            raise ValueError("There were no knapsack solutions. Consider adding other elements or increasing your mass uncertainty.")

        #print("non-identified masses: {} {}".format(len(non_identified_mass), non_identified_mass))
        #print('*** Knapsack, candidate elements: ***')
        #print(str([chem_data.element[i] for i in idx_list_kn]))

        #print('*** Knapsack, used elements: ***')
        #cl_comp.update_comp_variable(all_solutions)
        #print(str([chem_data.element[i_element] for i_element in cl_comp.knapsack_list_elements_idx]))

        if total_DBE_solutions < 15:
            print("list_frag_str = {}".format([chem_data.get_string_formula(fi) for ri in all_solutions for fi in ri]))
            print("all_solutions = {}".format([[chem_data.get_string_formula(fi) for fi in ri] for ri in all_solutions]))

        if run_type.is_NIST():
            #we know the target molecule
            #so we force the code to fit the maximum of data signal and measured peaks.
            const_id.max_added_nodes = total_DBE_solutions + 1
            const_id_max_opt_signal = 99.990
            const_id_max_ratio_not_opt_masses = 0.001
        
        #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        # ALPINAC: STEP 2: Initialise directed graph                         +
        # Connect a fragment to closest, larger fragments.                   +
        # Eliminate singletons.                                              +
        #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        """---------------------------------------------------------------------
        Now we have all abundant-isotope formulae matching the mass criteria
        and the double-bound-equivalent criteria.

        We want to eliminate all wrong formulae.
        at step2, we build the graph with all nodes.
        On the directed graph with all formulas, the formula should not be a singleton
        otherwise it is eliminated.

        Now we prepare the directed graph. To minimise the comparisons when building the graph,
        the fragments are expected to be given by decreasing mass. It is very important,
        otherwise the algorithm of graph initialisation is wrong:
        add_nodes_and_edges_to_graph_from_list_fragments
        ---------------------------------------------------------------------"""


        #------------------------
        #ici ajout ligne calcul peak width for match isotopocules rares

        t0_graph = time.time()

        #create empty directed graph.
        # maximal_nodes is a list of maximal elements, "roots" in the graph.
        # add all elements in decreasing order of mass, so that the subfragments of any fragment are added after it.
        # first, add the heaviest nodes, there are all maximal elements. (there can be other maximal elements, of smaller mass).

        G_reversed_edges = nx.DiGraph(maximal_nodes = [],
                                      candidates_per_measured_mass = {},
                                      no_nodes_knapsack = total_DBE_solutions,
                                      no_nodes_singletons = 0,
                                      no_nodes_below_LOD = 0,
                                      no_nodes_validated = 0,
                                      no_nodes_not_optimised = 0,
                                      added_nodes = [],
                                      nodes_ordered_by_likelihood = [],
                                      nodes_validated = [],
                                      nodes_molecular_ions = [])
        # candidates_per_measured_mass: save the number of candidates per mass for remove_singletons
        for i in range(len(all_solutions)-1, -1, -1):
            list_frag = all_solutions[i]
            if len(list_frag) == 0:
                continue
            add_nodes_and_edges_to_graph_from_list_fragments(G_reversed_edges, list_frag, idx_measured_mass=i, verbose=False)
        t1_graph = time.time()
        runtime_step2_graph = t1_graph-t0_graph
        print("build graph and edges of all {} fragments in {:.6f}s".format(total_DBE_solutions, t1_graph-t0_graph))
        max_nodes = [G_reversed_edges.nodes[n]['node_lite'] for n in G_reversed_edges.graph['maximal_nodes']]
        singletons = [n for n in max_nodes if len(G_reversed_edges.succ[n.unique_id]) == 0]
        if len(max_nodes) <= 30:
            print("{} maximal element(s) {}".format(len(max_nodes), " ".join([str(ni) for ni in max_nodes]) ))
        else:
            print("{} maximal elements".format(len(max_nodes)))
        print("{} singleton(s) {}".format(len(singletons), " ".join([str(ni) for ni in singletons])))

        if show_plt_figures and total_DBE_solutions <= 200:
            try:
                set_str_formula_and_attribute_for_pydot(G_reversed_edges)
                plot_forest_of_frag_with_pydot(G_reversed_edges, graphviz_layout_prog='dot')
            except:
                print("an error occured while drawing graph, maybe pydot is not installed?")

        removed_singletons = remove_singletons(G_reversed_edges)
        print("{} removed singleton(s): {}".format(len(removed_singletons), " ".join([str(ni)+" ({})".format(ni.idx_measured_mass[0]) for ni in removed_singletons])))
        print("recompute maximal elements")
        max_nodes = [G_reversed_edges.nodes[n]['node_lite'] for n in G_reversed_edges.graph['maximal_nodes']]
        if len(max_nodes) <= 30:
            print("{} maximal element(s) {}".format(len(max_nodes), " ".join([str(ni) for ni in max_nodes]) ))
        else:
            print("{} maximal elements".format(len(max_nodes)))
        # now compute likelihood: k_guess, number of successors in the graph, number of subfragments given by knapsack
        # the singletons where already removed

        #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        # ALPINAC: STEP 3: Initialise isotopocule sets                       +
        # Generate rare-isotope formulae.                                    +
        # compute isotopocule profiles.                                      +
        #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        t0_step3_enum_iso = time.time()
        initialise_isotopologue_set(G_reversed_edges, cl_comp)
        print("initialise_isotopologue_set done.")
        runtime_step3_enum_iso = time.time() - t0_step3_enum_iso
        #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        # ALPINAC: STEP 4: Compute max. contribution of each set             +
        # for each isotopocule set, individually.                            +
        # Prepare for computation of likelihood estimator.                   +
        #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        #theoretical values needed to compute likelihood estimator,
        #these values will be constant until the end:
        print("compute the theoretical number of possible subfragments for each node")
        knapsack_subfragments_with_lists = True
        if knapsack_subfragments_with_lists:
            t0_graph = time.time()
            compute_all_subfragments_of_nodes_with_lists(G_reversed_edges, min_measured_mass=const_id_min_detect_mass) # , verbose=True
            t1_graph = time.time()
            print("compute all sub-fragments with lists: {:.6f}".format(t1_graph-t0_graph))
        else:
            t1_graph = time.time()
            compute_all_subfragments_of_nodes(G_reversed_edges, min_measured_mass=const_id_min_detect_mass)
            t2_graph = time.time()
            print("compute all sub-fragments:            {:.6f}".format(t1_graph-t2_graph))
        # now each node has a field G_reversed_edges.nodes[n]['node_lite'].subgraph_size and
        # G_reversed_edges.nodes[n]['node_lite'].no_all_subfragments
        # initialise isotopologue class for nodes and compute isotopologue profile


        initialise_isotopologue_k_guess(G_reversed_edges, const_id.ppm_mass_res, cl_comp, run_type)


        ########################################################################
        t0_graph = time.time()
        # Now, change G_reversed_edges to G_of_frag so that the edges are now pointing to the larger fragments

        G_of_frag = G_reversed_edges

        G_of_frag.graph['added_nodes'] = list(G_reversed_edges.nodes)
        G_of_frag.graph['list_of_fitted_nodes'] = []


        if cl_comp.meas_len < const_id.min_no_masses:
            #less than 6 measured masses. This is an under-constrained system.
            mode_identification = measured_masses_few
            print("***There are less than 6 measured masses. Treating each maximal fragment separately.***")
            #do not fit all solutions together, the risk that a bad one excludes the correct one is too high.
            #optimise each maximal fragment separately.


            nodes_maximal_fragments = G_of_frag.graph['maximal_nodes'] # this gives a list of node numbers, for maximal nodes
            print("maximal nodes are: {}".format(" ".join([G_of_frag.nodes[node_i]['node_lite'].iso.frag_abund_formula for node_i in nodes_maximal_fragments])))

            #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            # ALPINAC: MODE FEW MASSES, STEP 5 TO 8                          +
            # OPTIMISE EACH MAXIMAL FRAGMENT WITH ALL ITS CHILDREN           +
            #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            runtime_step7_opt_loop = 0.0
            for node_max in nodes_maximal_fragments:
                G_of_frag.graph['nodes_validated'] = []
                if node_max in G_of_frag and node_max is not None:
                    print('Maximal fragment: ' + G_of_frag.nodes[node_max]['node_lite'].iso.frag_abund_formula)
                    G_of_frag.graph['nodes_validated'].append(node_max)
                    set_subfragments = nx.descendants(G_of_frag, node_max)
                    for val in set_subfragments:
                        if val not in G_of_frag.graph['nodes_validated'] and val is not None:
                            G_of_frag.graph['nodes_validated'].append(val)
                    print('all nodes of maximal fragment, incl. max.: ' + str(len(G_of_frag.graph['nodes_validated'])))

                    iso_fit_mode = fit_mode_fixed_x_axis
                    if not (run_type.is_NIST() or run_type.is_unit_mass()):
                        iso_fit_mode = fit_mode_variable_x_axis

                    t0_step7_opt_loop = time.time()
                    graph_optimise_isotopologues(G_of_frag, run_type, cl_comp, const_id.ppm_mass_res, iso_fit_mode, plot_mode = False)
                    runtime_step7_opt_loop += time.time() - t0_step7_opt_loop

            nodes_maximal_fragments = G_of_frag.graph['maximal_nodes'] # this gives a list of node numbers, for maximal nodes
            print("+++maximal fragments are: {}".format(" ".join([G_of_frag.nodes[node_i]['node_lite'].iso.frag_abund_formula for node_i in nodes_maximal_fragments])))

            #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            # ALPINAC: MODE FEW MASSES, STEP 6: Select most likely max fragment  +
            # and all its children, mark all others as "not validated"           +
            # Add max frag and children until signal reaches thereshold.         +
            #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            graph_order_by_att_likelihood(G_of_frag, cl_comp)

            G_of_frag.graph['nodes_validated'] = []
            sum_signal = float(0.0)
            percent_sum_signal = float(0.0)
            idx_node = 0

            while idx_node < len(G_of_frag.nodes_ordered_by_likelihood) and percent_sum_signal <= const_id_max_opt_signal:
                #take the most likely maximal node
                node_max = G_of_frag.nodes_ordered_by_likelihood[idx_node]
                print("node_max  " + str(node_max))
                if node_max in nodes_maximal_fragments and node_max not in G_of_frag.graph['nodes_validated']:
                    G_of_frag.graph['nodes_validated'].append(node_max)
                    #print("nodes_validated   " + str(G_of_frag.graph['nodes_validated']))

                    set_subfragments = nx.descendants(G_of_frag, node_max)
                    for val in set_subfragments:
                        if val not in G_of_frag.graph['nodes_validated'] and val is not None:
                            G_of_frag.graph['nodes_validated'].append(val)
                    

                    sum_signal, percent_sum_signal, no_optimised_masses, no_masses_to_optimise = graph_percent_sum_signal(G_of_frag, cl_comp.sum_I)
                    print("Mode few masses, percent_sum_signal  " + str(percent_sum_signal))
                idx_node += 1
            print("+++Mode few masses, nodes_validated:   " + str(G_of_frag.graph['nodes_validated']))

        else:
            #^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            #for cases with at least 6 measured masses
            #do not work with maximal fragments, start with most likely fragments.
            mode_identification =  measured_masses_enough

            #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            # ALPINAC: STEP 5: Rank fragments                                    +
            # Per fragment including all its children:                           +
            # rank according to likelihood estimator                             +
            #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

            graph_order_by_att_likelihood(G_of_frag, cl_comp) # TODO Aurore rename as graph_order_nodes_by_likelihood

            #all existing nodes are now ordered by decreasing likelihood in G.nodes_ordered_by_likelihood


            #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            # ALPINAC: NOW ENTERING LOOP,                                    +
            # STEPS 5 TO 8                                                   +
            #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            runtime_step7_opt_loop = 0

            #now fit several fragment formulas together,
            #each fragment added to the list of node one after the other,
            #taking all existing children per fragment

            #INITIALISATION OF LOOP
            list_optimised_masses = []
            no_iteration = 0
            no_added_nodes = 0
            sum_signal = 0.0
            percent_sum_signal = 0.0
            ratio_not_opt_masses = sum([1 for v in all_solutions if len(v)>0])
            not_optimised_masses = True


            #the loop goes on as long as:
                #the number of "validated" nodes is less than the number of nodes in G
                #as safety: the number of added nodes is less than a fixed value (to improve: threshold set as function of no. meas masses?)
                #the minimum fraction of fitted signal is not reached
                #the minimum fraction of not optimised measured masses

            while len(G_of_frag.graph['nodes_validated']) < G_of_frag.number_of_nodes() \
            and not_optimised_masses \
            and no_added_nodes < const_id.max_added_nodes\
            and(percent_sum_signal <= const_id_max_opt_signal\
            or ratio_not_opt_masses > const_id.max_ratio_not_opt_masses):


                #^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
                #mygu 20201005
                #if explained signal has reached set threshold,
                #but fraction of explained masses has not reached set threshold,
                #we modify likelihood list to take first nodes of not explained masses.
                #Re-write likelihood list, keep order but put first nodes corresponding to a non-explained mass,
                #and after all other nodes (still in their own likelihood order).
                if percent_sum_signal > const_id_max_opt_signal \
                and not_optimised_masses:
                    print("===Sum signal > " + str(const_id_max_opt_signal) + " , now adding masses of non-fitted measured masses")
                    nodes_of_not_opt_masses = []

                    for key_node in G_of_frag.nodes_ordered_by_likelihood:
                        test_add = True
                        idx_m_to_add = 0

                        iso = G_of_frag.nodes[key_node]['node_lite'].iso
                        while idx_m_to_add < len(iso.meas_mass_idx) and test_add:
                            for i_m in range(len(iso.meas_mass_idx[idx_m_to_add])):
                                if iso.meas_mass_idx[idx_m_to_add][i_m]  is not None and iso.meas_mass_idx[idx_m_to_add] in list_optimised_masses:
                                    test_add = False
                            idx_m_to_add += 1
                        if test_add:
                            nodes_of_not_opt_masses.append(key_node)


                    #print("Added nodes of non-opt masses: " + str(len(nodes_ordered_by_likelihood_new)))
                    #now add next all remaining nodes:
                    if len(nodes_of_not_opt_masses) > 0:
                        print('{:.2f}'.format(ratio_not_opt_masses)  + " frac. not opt. masses; " + str(len(nodes_of_not_opt_masses)) + " node(s) of non-optimised masses to be added.")
                        for i_add_node in range(G_of_frag.number_of_nodes()):
                            if G_of_frag.nodes_ordered_by_likelihood[i_add_node] not in nodes_of_not_opt_masses:
                                nodes_of_not_opt_masses.append(G_of_frag.nodes_ordered_by_likelihood[i_add_node])
                    else:
                        print("WARNING! There are no nodes of non-optimise masses to add.")


                        G_of_frag.nodes_ordered_by_likelihood = nodes_of_not_opt_masses.copy()
                #^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

                #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                # ALPINAC STEP 6: Select most likely fragments                       +
                # and all their children                                             +
                #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

                #the list of nodes ordered by likelihood potentially changes at each iteration
                #we want to add the most likely node not yet added
                #we need to check the entire list of added nodes
                added_node = 0
                i_add = 0
                while i_add < G_of_frag.number_of_nodes() \
                and added_node == 0:
                    node_added = G_of_frag.nodes_ordered_by_likelihood[i_add]

                    if node_added in G_of_frag and node_added not in G_of_frag.graph['nodes_validated']  and node_added is not None:
                        added_node += 1
                        no_iteration += 1
                    else:
                        i_add +=1


                if added_node > 0:

                    G_of_frag.graph['nodes_validated'].append(node_added)
                    #print('Nodes valid. at l. 478: ' + str(G_of_frag.graph['nodes_validated']))
                    #print('list_optimised_masses: ' + str(list_optimised_masses))

                    #update list of validated masses
                    #not sure this is the right place, what if the node is deleted after?
                    for key_node in G_of_frag.graph['nodes_validated']:
                        iso = G_of_frag.nodes[key_node]['node_lite'].iso
                        for i_iso in range(iso.iso_list_len):
                            if iso.meas_mass_idx[i_iso] is not None and iso.meas_mass_idx[i_iso] not in list_optimised_masses:
                                list_optimised_masses.append(iso.meas_mass_idx[i_iso])

                    set_subfragments = nx.descendants(G_of_frag, node_added)
                    for val in set_subfragments:
                        if val not in G_of_frag.graph['nodes_validated'] and val is not None:
                            G_of_frag.graph['nodes_validated'].append(val)
                            added_node += 1

                    #print(str([G_of_frag.nodes[val]['node_lite'].iso.frag_abund_formula for val in set_subfragments]))
                    no_added_nodes += added_node
                    print('****************************************')
                    print('Added node: ' + str(G_of_frag.nodes[node_added]['node_lite'].iso.frag_abund_formula) + " and " + str(added_node-1) + " childrens. Total: " + str(no_added_nodes))


                    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                    # ALPINAC STEP 7: Optimise multiple isotopocule sets                 +
                    # Optimise contribution of isotopocule sets to fit measured profile  +
                    # Use Python package lmfit                                           +
                    # Eliminate sets < LOD.                                              +
                    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


                    if not (run_type.is_NIST() or run_type.is_unit_mass()):
                        iso_fit_mode = fit_mode_variable_x_axis
                    else:
                        iso_fit_mode = fit_mode_fixed_x_axis

                    t0_opt_loop = time.time()
                    graph_optimise_isotopologues(G_of_frag, run_type, cl_comp, const_id.ppm_mass_res, iso_fit_mode, plot_mode=False)
                    runtime_step7_opt_loop += (time.time() - t0_opt_loop)

                    #update ranking of nodes
                    graph_order_by_att_likelihood(G_of_frag, cl_comp)


                    #***************************************
                    #Compute explained sum signal
                    sum_signal, percent_sum_signal, no_optimised_masses, no_masses_to_optimize = graph_percent_sum_signal(G_of_frag, cl_comp.sum_I)
                    ratio_not_opt_masses = float(no_masses_to_optimize-no_optimised_masses)/float(no_masses_to_optimize)
                    if no_optimised_masses >= no_masses_to_optimize:
                        #if no_masses_to_optimize == 0:
                        not_optimised_masses = False

                    #**************************************
                    #compute likelihood of being present for each chemical element
                    #chem_elem_likelihood= [float(0)]*len(G_of_frag.nodes[max_node]['node_lite'].iso.frag_abund
                    

                    print('***Fraction of not optimised masses: ' + '{:.4f}'.format(ratio_not_opt_masses))
                    print('***% Sum signal after ' + str(no_iteration) + ' iterations: ' + '{:.2f}'.format(percent_sum_signal))
                    print('masses: [optimised; solution(s) exist; measured] ' + str([no_optimised_masses, no_masses_to_optimize, cl_comp.meas_len]))



        #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        # ALPINAC STEP 9: Eliminate not optimised sets                       +
        #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        #remove all nodes that have not been optimised
        #This is done for all cases, with more or less than 6 measured masses.
        nodes_to_be_removed = [i_remain_node for i_remain_node in G_of_frag.nodes if i_remain_node not in G_of_frag.graph['nodes_validated']]
        for i_remain_node in nodes_to_be_removed:
            remove_node_update_edges(G_of_frag, G_of_frag.nodes[i_remain_node]['node_lite'], cl_comp)
        G_of_frag.graph['no_nodes_not_optimised'] += len(nodes_to_be_removed)
        print("graph max nodes " + str(G_of_frag.graph['maximal_nodes']))

        #check again for singletons!
        old_percent_sum_signal = percent_sum_signal
        #print("old_percent_sum_signal " + str(old_percent_sum_signal))
        list_nodes_singletons = [node.unique_id for node in remove_singletons(G_of_frag)]

        sum_signal, percent_sum_signal, no_optimised_masses, no_masses_to_optimize = graph_percent_sum_signal(G_of_frag, cl_comp.sum_I)
        #print("percent_sum_signal " + str(percent_sum_signal))
        if ((percent_sum_signal-old_percent_sum_signal) < float(-5.0)) or percent_sum_signal < const_id_max_opt_signal:
            print("***More singleton removed, optimise signal one more time***")
            #at least 5% of signal now lost because validated nodes have now been removed as singleton
            #because all remaining nodes have been removed.
            #update signal one last time

            graph_optimise_isotopologues(G_of_frag, run_type, cl_comp, const_id.ppm_mass_res, iso_fit_mode, plot_mode=False)


            graph_order_by_att_likelihood(G_of_frag, cl_comp)

            sum_signal, percent_sum_signal, no_optimised_masses, no_masses_to_optimize = graph_percent_sum_signal(G_of_frag, cl_comp.sum_I)



        #**********************************************
        #optimisation of signal for all selected nodes is done.

        #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        # ALPINAC STEP 10: create candidate molecular ion(s)                 +
        #                                                                    +
        #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        #Guessing the molecular ion is not a piece of cake.
        #We do that only now, using validated maximal elements
        #We create extra molecular formula and add to the graph,
        #without searching to match measured masses.

        #A node may be a max element in the set of validated but not in the total set.
        #So we need to remove all non-validated nodes first,
        #and update the list of maximal nodes.
        
        #take one max element with max likelihood.
        #If DBE is an even number,
        #this may be a molecular ion, do not add more atoms.       
        #If DBE is an odd number,
        #it needs at least one more atom to be a molecular ion.
        
        
        old_max_nodes = G_of_frag.graph['maximal_nodes'].copy()
        mol_ion_nodes = []
        frag_candidate_mol_ion = []
        max_nodes_odd_DBE = []
        for max_node in G_of_frag.graph['maximal_nodes']:
            #apply Kind and Fiehn 2007, rule 2
            #test SENIOR rule: "i) The sum of valences or the total number of atoms having
            #odd valences is even"
            if (G_of_frag.nodes[max_node]['node_lite'].iso.DBE)-int(G_of_frag.nodes[max_node]['node_lite'].iso.DBE) ==0:
                #even DBE number, this may be a molecular ion                                    
                #Check: The sum of valences is greater than or equal to twice the maximum valence
                #cf. Kind and Fiend 2007, p. 8, 1st column, rule (ii)
                #with this test, e.g. CFCl cannot be a molecular ion.
                if chem_data.test_SENIOR_rule_ii(G_of_frag.nodes[max_node]['node_lite'].iso.frag_abund):
                    mol_ion_nodes.append(max_node)
                    print("Candidate mol. ion even DBE: " + str(G_of_frag.nodes[max_node]['node_lite'].iso.frag_abund_formula))
                    #is already on the graph, does not need to be added to the graph
            else:
                #odd DBE number, at least one atom missing
                max_nodes_odd_DBE.append(max_node)
                print("Max node odd DBE:  " + str(G_of_frag.nodes[max_node]['node_lite'].iso.frag_abund_formula))
                #print("Added candidate mol. ion odd DBE:  " + str(chem_data.get_string_formula(G_of_frag.nodes[max_node]['node_lite'].iso.frag_abund)))
                #print(G_of_frag.nodes[max_node]['node_lite'].iso.frag_abund)

        print(str(len(max_nodes_odd_DBE)) + " maximal nodes with odd DBE found")            
        #now in the list of odd DBE max nodes,
        #we want to enumerate all solution of adding one chemical element
        #to all possible max nodes
        max_nodes_odd_and_even = mol_ion_nodes + max_nodes_odd_DBE
        
        if len(max_nodes_odd_DBE) > 0:
            #list max number for all validated atoms
            #and their maximum number in the entire set of validated molecules
            frag_max = [max([G_of_frag.nodes[max_node_odd]['node_lite'].iso.frag_abund[i_a] for max_node_odd in max_nodes_odd_and_even]) for i_a in range(chem_data.len_frag)]
            for i_a in range(chem_data.len_frag):
                if frag_max[i_a]>0:
                    frag_max[i_a] +=1            
            #print("constructed max frag: " + chem_data.get_string_formula(frag_max))

            #now for each odd DBE max node, add one more monovalent atom at a time
            for max_node in max_nodes_odd_DBE:
                for i_a in chem_data.idx_mono_valent:
                    candidate_max_node = list(G_of_frag.nodes[max_node]['node_lite'].iso.frag_abund) # deep copy!
                    if candidate_max_node[i_a] < frag_max[i_a]:
                        candidate_max_node[i_a] += 1
                        if candidate_max_node not in frag_candidate_mol_ion\
                        and chem_data.test_SENIOR_rule_ii(candidate_max_node):
                            frag_candidate_mol_ion.append(candidate_max_node)
                            print("Candidate mol. ion: " + str(chem_data.get_string_formula(candidate_max_node)))
                            
                    

        #add candidate molecular ions to graph
        if len(frag_candidate_mol_ion)>0:
            print("===Add molecular ion(s) to graph===")
            add_maximal_node_and_edges_to_graph_from_max_fragment(G_of_frag, 
                                                                  list_max_fragments = frag_candidate_mol_ion, 
                                                                  molecular_ion = True, 
                                                                  verbose = True)
            
            #indexes of added nodes (candidate molecular ions)
            idx_new_max_nodes = [i_node for i_node in G_of_frag.graph['maximal_nodes'] if i_node not in old_max_nodes]
            mol_ion_nodes = mol_ion_nodes + idx_new_max_nodes

            initialise_isotopologue_set(G_of_frag, cl_comp, add_mol_ion=True, mol_ion_nodes = idx_new_max_nodes)
            
            compute_all_subfragments_of_nodes_with_lists(G_of_frag, min_measured_mass=const_id_min_detect_mass)
            
            graph_order_by_att_likelihood(G_of_frag, cl_comp)

        #for key_node in G_of_frag.graph['maximal_nodes']:           
        #    print("iso.iso_list_sum_signal: " + str(G_of_frag.nodes[key_node]['node_lite'].iso.iso_list_sum_signal))
        #    print("subgraph_sum_I: " + str(G_of_frag.nodes[key_node]['node_lite'].subgraph_sum_I))
        #    print("subgraph_percent_I: " + str(G_of_frag.nodes[key_node]['node_lite'].subgraph_percent_I))
        #    print("no_all_subfragments: " + str(G_of_frag.nodes[key_node]['node_lite'].no_all_subfragments))
            
        


        G_of_frag.graph['no_nodes_validated'] = G_of_frag.number_of_nodes()
        print('*********************************')
        if mode_identification == measured_masses_few:
            print('Identification done separately for each maximal fragment, total % ' + '{:.2f}'.format(percent_sum_signal))
        else:
            print('***%Sum signal after ' + str(no_iteration) + ' iterations: ' + '{:.2f}'.format(percent_sum_signal))
        print('knapsack solutions:        ' + str(G_of_frag.graph['no_nodes_knapsack']) + '\t'      + str(100) + ' %')
        print('nodes removed, < LOD:      ' + str(G_of_frag.graph['no_nodes_below_LOD']) + '\t'     + str(G_of_frag.graph['no_nodes_below_LOD']     /G_of_frag.graph['no_nodes_knapsack']*100) + ' %')
        print('nodes removed, singletons: ' + str(G_of_frag.graph['no_nodes_singletons']) + '\t'    + str(G_of_frag.graph['no_nodes_singletons']    /G_of_frag.graph['no_nodes_knapsack']*100) + ' %')
        print('nodes not optimised:       ' + str(G_of_frag.graph['no_nodes_not_optimised']) + '\t' + str(G_of_frag.graph['no_nodes_not_optimised'] /G_of_frag.graph['no_nodes_knapsack']*100) + ' %')
        print('nodes validated:           ' + str(G_of_frag.graph['no_nodes_validated']) + '\t'     + str(G_of_frag.graph['no_nodes_validated']     /G_of_frag.graph['no_nodes_knapsack']*100) + ' %')
        print('*********************************')


        print('***Most likely fragments***')

        graph_order_by_att_likelihood(G_of_frag, cl_comp)
        runtime_total = time.time() - t0_total

        
        candidate_mol_ion_formulae = [G_of_frag.nodes[key_node]['node_lite'].iso.frag_abund_formula for key_node in mol_ion_nodes]
        candidate_mol_ion_formulae_str = ""

        if run_type.is_tofwerk_tof():
            outfilename = fragfile_name
        else:
            outfilename = toffile_name

        export_name = outfilename

        fragment_results = []
        validation_dict_results = {}
        validation_dict_results["frac_correct_frag"]  = 0
        validation_dict_results["frac_true_signal_to_total_signal"] = 0
        validation_dict_results["frac_true_frag_top10"] = 0
        validation_dict_results["frac_true_signal_top10"] = 0
        if G_of_frag.number_of_nodes() == 0:
            print("There are no solutions. Consider adding other elements or increasing the mass uncertainty.")

        else:

            #print the remaining maximal nodes
            nodes_maximal_fragments = G_of_frag.graph['maximal_nodes']

            print('****************************************')
            print(str(len(nodes_maximal_fragments)) + ' most likely maximal fragments(s)')
            print(nodes_maximal_fragments)

            print('formula' + '\t' +
                  'DBE' + '\t' +
                  'likelihood'+ '\t' +
                  'ranking')

            for idx_node in range(G_of_frag.number_of_nodes()):
                if G_of_frag.nodes_ordered_by_likelihood[idx_node] in nodes_maximal_fragments:

                    node = G_of_frag.nodes[G_of_frag.nodes_ordered_by_likelihood[idx_node]]['node_lite']
                    print("{}\t{}\t{:.1f}\t{}".format(node.iso.frag_abund_formula, node.iso.DBE, node.subgraph_likelihood, idx_node+1))

            
            print('****************************************')
            print(str(len(mol_ion_nodes)) + " Most likely molecular ion(s)")
            for idx_node in range(G_of_frag.number_of_nodes()):
                if G_of_frag.nodes_ordered_by_likelihood[idx_node] in mol_ion_nodes:

                    node = G_of_frag.nodes[G_of_frag.nodes_ordered_by_likelihood[idx_node]]['node_lite']
                    print("{}\t{}\t{:.1f}\t{}".format(node.iso.frag_abund_formula, node.iso.DBE, node.subgraph_likelihood, idx_node+1))
                    candidate_mol_ion_formulae_str += str(chem_data.get_string_formula_sub(node.iso.frag_abund)) + "; "
                
                
                
            if mode_identification != measured_masses_few and show_plt_figures:
                plot_graph_of_frag(G_of_frag)


            #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            # EXTRA STEP: VALIDATION (if filename part of validation set)    +
            # The reconstructed fragments are compared                       +
            # to theoretically possible fragments, with re-arrangements.     +
            #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

            list_identified_fragments = [G_of_frag.nodes[key_node]['node_lite'].iso.frag_abund  for key_node in G_of_frag.nodes_ordered_by_likelihood]
            #print(list_identified_fragments)

            list_checked = [None]*len(list_identified_fragments)
            do_val = False
            val_name = "---"

            if fragfile_name is not None:
                if 'train_' in fragfile_name or 'val_' in fragfile_name or 'known_' in fragfile_name:
                    do_val = True 
            if run_type.is_tofwerk_tof() and do_val:
                #find correct index in mol_validation, index 6 should be same as text string after val_ in filename


                idx_val = find_validation_idx_from_database(fragfile_name, mol_validation, validation_labels)
                if idx_val is not None:

                    #for mol_tuple in mol_validation[idx_val]:
                    #for mol_tuple in mol_validation:
                    #load SMILES code of corresponding true molecule
                    smiles = mol_validation[idx_val][validation_labels["SMILES code"]]
                    mol_name = mol_validation[idx_val][validation_labels["Sum formula"]]
    
                    str_sum_formula = chem_data.get_string_formula(chem_data.get_fragment_from_string_formula(mol_name))
    
                    list_checked = validate_frag_list(smiles, mol_name, list_identified_fragments, draw_graph= show_plt_figures)
                    #for i in range(len(list_identified_fragments)):
                    #    print(str(chem_data.get_string_formula(list_identified_fragments[i])) + "\t" + str(list_checked[i]))
    
                    validation_dict_results, fragment_results = compute_validation_results(list_checked, G_of_frag, sum_signal, cl_comp.sum_I)
                    validation_dict_results["val_name"] = val_name
                    export_name = mol_validation[idx_val][validation_labels["Compound"]]

            
            
            method_performance_results = [export_name,
                                          cl_comp.meas_len,
                                          cl_comp.sum_I,
                                          cl_comp.meas_mass_u_median,
                                          cl_comp.meas_mass_u_median_ppm,
                                          percent_sum_signal,
                                          G_of_frag.graph['no_nodes_knapsack'],
                                          G_of_frag.graph['no_nodes_below_LOD']/G_of_frag.graph['no_nodes_knapsack']*100,
                                          G_of_frag.graph['no_nodes_singletons']/G_of_frag.graph['no_nodes_knapsack']*100,
                                          G_of_frag.graph['no_nodes_not_optimised']/G_of_frag.graph['no_nodes_knapsack']*100,
                                          G_of_frag.graph['no_nodes_validated']/G_of_frag.graph['no_nodes_knapsack']*100,
                                          runtime_step1_knapsack,
                                          runtime_step2_graph,
                                          runtime_step3_enum_iso,
                                          runtime_step7_opt_loop,
                                          runtime_total,
                                          candidate_mol_ion_formulae_str
                                          ]

            #chemSpider
            #https://matt-swain.com/blog/2012-02-02-chemspipy-python-wrapper-chemspider-api
            #https://chemspipy.readthedocs.io/en/stable/

            #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            # WRITE RESULTS IN TEXT FILE                                     +
            #                                                                +
            #                                                                +
            #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

            ires_idx_m = 0
            ires_rt = 1
            ires_I = 2
            ires_I_norm = 3
            ires_m_meas = 4
            ires_m_exact = 5
            ires_m_meas_u = 6
            ires_m_diff = 7
            ires_formula = 8
            ires_DBE = 9
            ires_I_frag = 10
            ires_I_frag_norm = 11
            ires_ionisation = 12

            #group all generated formulas by index of measured mass:
            results_found_list = []#*len(meas_mass)
            results_not_found_list = []
            #print('max mass index ' + str(cl_comp.meas_len))

            for i_node in G_of_frag.nodes():
                #if i_list in iso_set_idx_archive:

                iso = G_of_frag.nodes[i_node]['node_lite'].iso
                #print(iso.meas_mass_idx)
                for i_iso in range(len(iso.iso_list)):

                    #compute sum of given isotopocule in case it matches several measured masses
                    #sum_meas_I = sum([cl_comp.meas_I[idx_m] for idx_m in iso.meas_mass_idx[i_iso] if idx_m is not None])
                    #iso.iso_list_I_meas_max
                    I_meas_max = iso.iso_list_I_meas_max[i_iso]
                    if I_meas_max == 0:
                        I_meas_max = float(1)
                    #for i_m in range(len(iso.meas_mass_idx[i_iso])):
                    #    idx_m = iso.meas_mass_idx[i_iso][i_m]

                    if iso.meas_mass_idx[i_iso] is not None:
                    #    sum_meas_I = sum([cl_comp.meas_I[idx_m] for idx_m in iso.meas_mass_idx[i_iso]])
                        for i_m in range(len(iso.meas_mass_idx[i_iso])):
                            idx_m = iso.meas_mass_idx[i_iso][i_m]


                            candidate_mass = chem_data.get_mass(iso.iso_list[i_iso])

                            #print('idx_m ' + str(idx_m))
                            #print(chem_data.get_string_formula(isotopologues_list[i_list].iso_list[i_iso]))
                            if idx_m is not None: # and isotopologues_list[i_list].iso_list_I_opt[i_iso] >= cl_comp.LOD/10 :
                                #print(chem_data.get_string_formula(isotopologues_list[i_list].iso_list[i_iso]))
                                results_found_list.append([
                                    idx_m,
                                    cl_comp.meas_rt[idx_m],
                                    cl_comp.meas_I[idx_m],
                                    cl_comp.meas_I[idx_m]*100.00 / cl_comp.meas_I_max,
                                    cl_comp.meas_mass[idx_m],
                                    candidate_mass,
                                    cl_comp.meas_mass_U_combined[idx_m]/cl_comp.meas_mass[idx_m] * 10**6,
                                    (candidate_mass - cl_comp.meas_mass[idx_m])/ cl_comp.meas_mass[idx_m] * 10**6,
                                    chem_data.get_string_formula(iso.iso_list[i_iso]),
                                    iso.DBE,
                                    iso.iso_list_I_opt[i_iso]*cl_comp.meas_I[idx_m]/I_meas_max,
                                    iso.iso_list_I_opt[i_iso] / I_meas_max,
                                    cl_comp.ionisation[idx_m]
                                    ])
                                #how to still keep results were no solution was found?
                            else:
                                results_not_found_list.append([
                                    -1,
                                    0.0,
                                    iso.iso_list_I_opt[i_iso],
                                    iso.iso_list_I_opt[i_iso]*100.00 / cl_comp.meas_I_max,
                                    0.0,
                                    candidate_mass,
                                    0.0,
                                    0.0,
                                    chem_data.get_string_formula(iso.iso_list[i_iso]) + ", not found",
                                    iso.DBE,
                                    iso.iso_list_I_opt[i_iso],
                                    0.0,
                                    '--',
                                    ])

            #results_list.sort(key=operator.itemgetter(4, 8))

            #average retentio time:
            rt_average = sum([results_found_list[i_res][ires_rt] * results_found_list[i_res][ires_I] for i_res in range(len(results_found_list))])/ \
            sum([results_found_list[i_res][ires_I] for i_res in range(len(results_found_list))])



            results_found_list.sort(key=operator.itemgetter(ires_I_frag, ires_I_frag_norm), reverse=True)


            #add unidientified masses to the list of fragments:

            identified_mass_idx_list, unidentified_mass_idx_list = get_list_identified_unidentified_mass(G_of_frag, cl_comp)


            for idx_m in unidentified_mass_idx_list:
                results_not_found_list.append([
                        idx_m,
                        cl_comp.meas_rt[idx_m],
                        cl_comp.meas_I[idx_m],
                        cl_comp.meas_I[idx_m]*100.00 / cl_comp.meas_I_max,
                        cl_comp.meas_mass[idx_m],
                        0.0,
                        cl_comp.meas_mass_u[idx_m]/cl_comp.meas_mass[idx_m] * 10**6,
                        0.0,
                        'unidentified',
                        0.0,
                        cl_comp.meas_I[idx_m],
                        1.0,
                        cl_comp.ionisation[idx_m]
                        ])

            results_list = results_found_list + results_not_found_list
            results_list.sort(key=operator.itemgetter(ires_I_frag, ires_I_frag_norm), reverse=True)

            #****************************************************************
            if run_type.is_NIST():
                suptitle = "NIST mass spectrum: {}, CAS: {}".format(nist_metadata.NIST_target_formula, nist_metadata.NIST_CAS) # nist_metadata.NIST_chem_name
                fig_name = str(nist_metadata.NIST_CAS) + "_NIST_mass_spectrum_high_res.png"
            else:
                suptitle = toffile_name + ": mass spectrum at {:.2f} s".format(rt_average)
                fig_name = toffile_name + "_mass_spectrum_{:.2f}s_{}_valid.png".format(rt_average, compound_number)

            fig_file = str(fig_path/fig_name)

            if show_plt_figures:
                plot_mass_spectrum(run_type, rt_average, results_found_list, results_not_found_list, chem_data.electron_mass, cl_comp.meas_I_max, suptitle, fig_file)
            #else:
            #    print('***There were less than 6 measured masses. No spectrum given. Check most likely fragments.')
            #****************************************************************

            #^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

            #^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            #WRITE OUTPUT FILE
            if run_type.is_NIST():
                out_formula_file_name = toffile_name + '.highres'  + '.txt'

            elif run_type.is_tofwerk_tof() or run_type.is_no_mass_cal() or run_type.is_unit_mass():

                #out_formula_file_name = toffile_name + '.comp.'  + str(int(round(cl_run.rt_min))) + '.' + str(int(round(cl_run.rt_max))) + '_' + str(list_frag_idx[0]) + '.txt'
                out_formula_file_name = toffile_name + '.comp.' + '{:.2f}'.format(rt_average) + '_' + str(compound_number) + '.txt'

            out_file = formulas_path / out_formula_file_name
            with out_file.open('w') as out_formula_file:
                out_formula_file.write('Run file: ' + toffile_name + '\n')

            #out_formula_file_name = 'formula_' + str(fragment_index) + '_' + str(filename)
            #out_formula_file_name = toffile_name + '.comp.' + str(cl_comp.fragment_index) + '.' + '{:.2f}'.format(cl_comp.average_rt) + '.txt'
            #out_formula_path = cl_path.dict_path['formulas']

            #file = out_formula_path / out_formula_file_name
            #with out_file.open('w') as out_formula_file:
            with out_file.open('a') as out_formula_file:
                out_formula_file.write('#############################################\n')
                if run_type.is_NIST():
                    #out_formula_file.write('NIST EI library' +  '\n')
                    out_formula_file.write('NIST EI library, compound formula: ' + str(nist_metadata.NIST_target_formula) +  '\n')
                elif run_type.is_tofwerk_tof():
                    out_formula_file.write('Compound index: ' + str(cl_comp.fragment_index) +  '\n')
                    out_formula_file.write('Average retention time, seconds: ' + '{:.2f}'.format(rt_average) +  '\n')
                out_formula_file.write('Largest detected fragment series:'  +  '\n')
                out_formula_file.write('Formula' + '\t' + 'Likelihood indicator (%)' + '\t' + 'DBE' + '\t' +  'Molecular ion?' + '\n')


                
                if G_of_frag.number_of_nodes() > 0:

                    #write suggested molecular ions
                    for idx_node in range(G_of_frag.number_of_nodes()):
                        if G_of_frag.nodes_ordered_by_likelihood[idx_node] in nodes_maximal_fragments:
        
                            node = G_of_frag.nodes[G_of_frag.nodes_ordered_by_likelihood[idx_node]]['node_lite']
                            is_mol_ion = node.iso.frag_abund_formula in candidate_mol_ion_formulae
                            #print("{}\t{}\t{:.1f}\t{}".format(node.iso.frag_abund_formula, node.iso.DBE, node.subgraph_likelihood, idx_node+1))
                            out_formula_file.write("{}+: \t{:.2f}\t{}\t{}\n".format(
                                node.iso.frag_abund_formula,
                                node.subgraph_likelihood,
                                node.iso.DBE,
                                is_mol_ion))


                out_formula_file.write('% norm I, exact' + '\t' +
                                       'RT [s]' + '\t' +
                                       'Total Intensity' + '\t' +
                                       '% of Total Intensity' + '\t' +
                                       'Meas. mass [m/z]' + '\t' +
                                       'Exact mass [m/z]' + '\t' +
                                       'Exact mass, u [ppm]' + '\t' +
                                       'Mass diff [ppm]' + '\t' +
                                       'Formula' + '\t' +
                                       'DBE'+ '\t' +
                                       'Intensity' + '\t' +
                                       'Intensity, fraction'  + '\t' +
                                       'Ionisation mode' +
                                       '\n')

                if run_type.is_tofwerk_tof() or run_type.is_no_mass_cal():
                    for line in results_list:

                        """
                        ires_idx_m = 0
                        ires_rt = 1
                        ires_I = 2
                        ires_I_norm = 3
                        ires_m_meas = 4
                        ires_m_exact = 5
                        ires_m_meas_u = 6
                        ires_m_diff = 7
                        ires_formula = 8
                        ires_DBE = 9
                        ires_I_frag = 10
                        ires_I_frag_norm = 11
                        """
                        if line[ires_idx_m] not in unidentified_mass_idx_list and line[ires_idx_m] != -1 and line[ires_I] >= cl_comp.LOD:
                            #found and identified fragment
                            out_formula_file.write('{:.4f}'.format(line[ires_I]/cl_comp.sum_I*line[ires_I_frag_norm]*100) + '\t' #idx measured mass
                                                   + '{:.2f}'.format(line[ires_rt]) +'\t' #rt
                                                   + '{:.2f}'.format(line[ires_I]) + '\t' #I total of peak
                                                   + '{:.1f}'.format(line[ires_I_norm]) + '\t' #I total of peak, %
                                                   + '{:.6f}'.format(line[ires_m_meas]-chem_data.electron_mass) + '\t' #measured mass
                                                   + '{:.6f}'.format(line[ires_m_exact]-chem_data.electron_mass) + '\t' #exact mass
                                                   + '{:.1f}'.format(line[ires_m_meas_u]) + '\t'
                                                   + '{:.1f}'.format(line[ires_m_diff]) + '\t' #mass diff, ppm
                                                   + str(line[ires_formula]) + '+'  + '\t' #formula
                                                   + str(line[ires_DBE]) + '\t' #DBE
                                                   + '{:.2f}'.format(line[ires_I_frag]) + '\t' #I frag
                                                   + '{:.3f}'.format(line[ires_I_frag_norm]) + '\t'
                                                   + str(line[ires_ionisation]) +'\n') #I frag, fraction

                        elif line[ires_idx_m] in unidentified_mass_idx_list and line[ires_idx_m] != -1 and line[ires_I] >= cl_comp.LOD:
                            #found fragment, not identified
                            out_formula_file.write('{:.2f}'.format(line[ires_I]/cl_comp.sum_I*line[ires_I_frag_norm]*100) + '\t' #idx measured mass
                                                   + '{:.2f}'.format(line[ires_rt]) +'\t' #rt
                                                   + '{:.2f}'.format(line[ires_I]) + '\t' #I total
                                                   + '{:.1f}'.format(line[ires_I_norm]) + '\t' #I total, %
                                                   + '{:.6f}'.format(line[ires_m_meas]-chem_data.electron_mass) + '\t' #measured mass
                                                   + str(line[ires_m_exact]) + '\t' #exact mass
                                                   + '{:.1f}'.format(line[ires_m_meas_u]) + '\t'
                                                   + str(line[ires_m_diff]) + '\t' #mass diff, ppm
                                                   + str(line[ires_formula]) + '\t' #formula
                                                   + str(line[ires_DBE]) + '\t' #DBE
                                                   + '{:.2f}'.format(line[ires_I_frag]) + '\t' #I frag
                                                   + '{:.3f}'.format(line[ires_I_frag_norm]) + '\t'
                                                   + str(line[ires_ionisation]) +'\n') #I frag, fraction
                        elif line[ires_I] >= cl_comp.LOD:
                            #expected fragment, not found
                            out_formula_file.write('{:.2f}'.format(line[ires_I]/cl_comp.sum_I*line[ires_I_frag_norm]*100) + '\t' #idx measured mass
                                                   + '{:.2f}'.format(rt_average) +'\t' #rt
                                                   + '{:.2f}'.format(line[ires_I]) + '\t' #I total
                                                   + '{:.1f}'.format(line[ires_I_norm]) + '\t' #I total, %
                                                   + str(line[ires_m_meas]) + '\t' #measured mass
                                                   + '{:.6f}'.format(line[ires_m_exact]-chem_data.electron_mass) +'\t' #exact mass
                                                   + '{:.1f}'.format(line[ires_m_meas_u]) + '\t'
                                                   + str(line[ires_m_diff]) + '\t' #mass diff, ppm
                                                   + str(line[ires_formula]) + '\t' #formula
                                                   + str(line[ires_DBE]) + '\t' #DBE
                                                   + '{:.2f}'.format(line[ires_I_frag]) + '\t' #I frag
                                                   + '{:.3f}'.format(line[ires_I_frag_norm]) + '\t'
                                                   + str(line[ires_ionisation]) +'\n') #I frag, fraction

                    #out_formula_file.write('\n')

                elif run_type.is_NIST():
                    for line in results_list:

                        if line[ires_idx_m] not in unidentified_mass_idx_list and line[ires_idx_m] != -1 and line[ires_I] >= cl_comp.LOD:
                            #found and identified fragment
                            out_formula_file.write('{:.2f}'.format(line[ires_I]/cl_comp.sum_I*line[ires_I_frag_norm]*100) + '\t' #idx measured mass
                                                   + '{:.2f}'.format(line[ires_rt]) +'\t' #rt
                                                   + '{:.1f}'.format(line[ires_I]) + '\t' #I total
                                                   + '{:.1f}'.format(line[ires_I_norm]) + '\t' #I total, %
                                                   + '{:.1f}'.format(line[ires_m_meas]) + '\t' #measured mass
                                                   + '{:.6f}'.format(line[ires_m_exact]-chem_data.electron_mass) + '\t' #exact mass
                                                   + '{:.1f}'.format(line[ires_m_meas_u]) + '\t'
                                                   + '{:.1f}'.format(line[ires_m_diff]) + '\t' #mass diff, ppm
                                                   + str(line[ires_formula]) + '+'  + '\t' #formula
                                                   + str(line[ires_DBE]) + '\t' #DBE
                                                   + '{:.2f}'.format(line[ires_I_frag]) + '\t' #I frag
                                                   + '{:.3f}'.format(line[ires_I_frag_norm]) +'\n') #I frag, fraction

                        elif line[ires_idx_m] in unidentified_mass_idx_list and line[ires_idx_m] != -1 and line[ires_I] >= cl_comp.LOD:
                            #found fragment, not identified
                            out_formula_file.write('{:.2f}'.format(line[ires_I]/cl_comp.sum_I*line[ires_I_frag_norm]*100) + '\t' #idx measured mass
                                                   + '{:.2f}'.format(line[ires_rt]) +'\t' #rt
                                                   + '{:.1f}'.format(line[ires_I]) + '\t' #I total
                                                   + '{:.1f}'.format(line[ires_I_norm]) + '\t' #I total, %
                                                   + '{:.1f}'.format(line[ires_m_meas]) + '\t' #measured mass
                                                   + str(line[ires_m_exact]) + '\t' #exact mass
                                                   + '{:.1f}'.format(line[ires_m_meas_u]) + '\t'
                                                   + str(line[ires_m_diff]) + '\t' #mass diff, ppm
                                                   + str(line[ires_formula]) + '\t' #formula
                                                   + str(line[ires_DBE]) + '\t' #DBE
                                                   + '{:.2f}'.format(line[ires_I_frag]) + '\t' #I frag
                                                   + '{:.3f}'.format(line[ires_I_frag_norm]) +'\n') #I frag, fraction
                        elif line[ires_I] >= cl_comp.LOD:
                            #expected fragment, not found
                            out_formula_file.write('{:.2f}'.format(line[ires_I]/cl_comp.sum_I*line[ires_I_frag_norm]*100) + '\t' #idx measured mass
                                                   + '{:.2f}'.format(rt_average) +'\t' #rt
                                                   + '{:.1f}'.format(line[ires_I]) + '\t' #I total
                                                   + '{:.1f}'.format(line[ires_I_norm]) + '\t' #I total, %
                                                   + str(line[ires_m_meas]) + '\t' #measured mass
                                                   + '{:.6f}'.format(line[ires_m_exact]-chem_data.electron_mass) + '\t' #exact mass
                                                   + '{:.1f}'.format(line[ires_m_meas_u]) + '\t'
                                                   + str(line[ires_m_diff]) + '\t' #mass diff, ppm
                                                   + str(line[ires_formula]) + '\t' #formula
                                                   + str(line[ires_DBE]) + '\t' #DBE
                                                   + '{:.2f}'.format(line[ires_I_frag]) + '\t' #I frag
                                                   + '{:.3f}'.format(line[ires_I_frag_norm]) +'\n') #I frag, fraction
                    #out_formula_file.write('\n')


    #t_stop = time.time()
    #print('Duration: ' + '{:.2f}'.format(t_stop - t_start))

    return method_performance_results, fragment_results, validation_dict_results



if __name__ == "__main__":
    plt.close('all')

    # path_file is one filename
    args=sys.argv
    print("usage: python3 {} <filename> [--target-elements <elements> --target-formula <formula> --verbose-knapsack]".format(args[0]))
    if len(args) <= 1:
        raise ValueError("please provide on command-line input a string as filename")

    filename = args[1]
    current_working_dir = Path.cwd()
    path_file = current_working_dir/Path(filename)
    print("file {}".format(filename))

    i=2
    target_elements=None
    target_formula=None
    verbose_knapsack=False
    cmd_options = ["--target-elements", "--target-formula", "--verbose-knapsack"]

    while i < len(args) and target_elements is None and target_formula is None:
        if args[i] == "--target-elements" and i<len(args)-1:
            target_elements = args[i+1]
            i = i+2
            continue
        if args[i] == "--target-formula" and i<len(args)-1:
            target_formula = args[i+1]
            i = i+2
            continue
        if args[i] == "--verbose-knapsack":
            verbose_knapsack=True
        i = i+1

    if target_elements is not None:
        print("target_elements={}".format(target_elements))
    if target_formula is not None:
        print("target_formula={}".format(target_formula))
    if verbose_knapsack:
        print("verbose_knapsack is True")

    run_type = RunType()
    fragfile_name = None
    if "Mass" in filename:
        run_type.set_NIST()
        # how to set the fragfile_name?
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

    # go to one parent directory to find the directory "mass_spectra" where to save figures:
    fig_path = path_file.parents[1]/"mass_spectra"
    formulas_path = path_file.parents[1]/"formulas"

    print("run_type = {}".format(run_type))

    show_plt_figures = False

    method_performance_results = \
    make_identification(run_type,
                        path_file,
                        fig_path,
                        formulas_path,
                        target_elements=target_elements,
                        target_formula=target_formula,
                        fragfile_name=fragfile_name,
                        show_plt_figures= show_plt_figures,
                        verbose_knapsack=verbose_knapsack)

# def make_identification(run_type: RunType,
#                         path_file: Path,
#                         fig_path: Path,
#                         formulas_path: Path,
#                         fragment_indices: list=None,
#                         min_fragment_index: int=None,
#                         target_elements: str=None,
#                         target_formula: str=None,
#                         fragfile_name: str=None,
#                         show_plt_figures: bool=True,
#                         verbose_knapsack: bool=False)-> None:
