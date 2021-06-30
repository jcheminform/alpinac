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


Validation of the data: does the knapsack give chimically possible answers?
"""

import math

import sys
from sys import version_info
if sys.version_info[0] < 3:
    from exceptions import ValueError

import pathlib
from pathlib import Path, PurePath

import networkx as nx
import matplotlib
import matplotlib.pyplot as plt


# https://github.com/pckroon/pysmiles
# pip3 install pysmiles
import pysmiles
from pysmiles.read_smiles import read_smiles

import itertools

import periodic_table as chem_data
from utils_identification import generate_all_isotopolog_fragments, p_iso_rare_compute
from io_tools import get_data_from_frag_file, MassCalibrationData

def graph_to_sum_formula(g):
    """Given a graph whose nodes are chemical elements represented as pairs (integer, string code),
    returns the list of elements as a list [n0,n1,...,ni] where ni = 0 if element i is not in the graph, or ni = the number of elements i in the graph (with multiplicity).
    """
    a = [0] * chem_data.len_frag
    nodes = g.nodes(data='element')
    for (ni,label) in nodes:
        a[chem_data.dict_element_idx[label]] += 1
    return a

def get_subgraphs(g, final_check_connectivity=False,verbose=False):
    # 0. make list of single atoms and pairs of atoms
    all_connected_subgraphs = [g.subgraph(gi) for gi in list(g.nodes()) + list(g.edges())]
    # 1. consider subgraphs made of multi-valent atoms
    multi_val = []
    mono_val = []
    dict_is_mono_val = {}
    #for ni in list(g.nodes()):
    for (ni,label) in list(g.nodes(data='element')):
        if verbose:
            if chem_data.valence[chem_data.dict_element_idx[label]] > 1 and len(list(g.adj[ni])) <= 1:
                print("-----------------------------------------------")
                print("node {} is chimically multi-valent but has {} edge".format((ni,label), len(list(g.adj[ni]))) )
        if len(list(g.adj[ni])) > 1:
            multi_val.append(ni)
            dict_is_mono_val[ni] = False
        else:
            mono_val.append(ni)
            dict_is_mono_val[ni] = True
    if verbose:
        print("single-edge nodes: {}\nmulti-edge nodes: {}".format(mono_val, multi_val))
    
    # we do not want to consider the entire graph itself, only proper subgraphs, strictly smaller    
    for nb_nodes in range(1, 1+min(len(multi_val), g.number_of_nodes()-1)):
        for selected_nodes in itertools.combinations(multi_val, nb_nodes):
            sg = g.subgraph(selected_nodes)
            if nx.is_connected(sg):
                if nb_nodes > 2:
                    all_connected_subgraphs.append(sg)
                # now try to add mono-valent atoms
                # 1. list all adjacent nodes that are monovalent atoms
                mono_val_sg = []
                for ni in list(selected_nodes):
                    for neighbor in list(g.adj[ni]):
                        if dict_is_mono_val[neighbor]: # if monovalent, no need to double-check uniqueness and (neighbor not in mono_val_sg):
                            mono_val_sg.append(neighbor)
                #if verbose:
                #    print("multi_val_nodes = {}, mono_val_neighbor_nodes = {}".format(selected_nodes, mono_val_sg))
                if nb_nodes == 1:
                    # all sub-graphs of 1 and 2 vertices were already considered,
                    # start at nb_nodes + min_mono_val = 3
                    min_mono_val = 2
                else:
                    min_mono_val = 1
                max_mono_val = min(g.number_of_nodes()-1-nb_nodes, len(mono_val_sg))
                # 2. add mono-valent atoms
                for no_mono_val in range(min_mono_val, max_mono_val+1):
                    for nodes_mono in itertools.combinations(mono_val_sg, no_mono_val):
                        sg2 = g.subgraph(list(selected_nodes) + list(nodes_mono))
                        if final_check_connectivity:
                            if not nx.is_connected(sg2):
                                raise ValueError("Error in algorithm")
                        all_connected_subgraphs.append(sg2)
            # if it is not connected, mono-valent atoms will not help on this. so do nothing.
    return all_connected_subgraphs


def get_mass_sum_formula_abundant_fragments(smiles, verbose=False):
    """
    smiles: a string that is a SMILES code of a molecule
    verbose: print details
    
    Returns a list of pairs (mass, sum_formula) for each possible fragment of the molecule.
    List is sorted by increasing mass.
    Duplicates are removed (same sum-formula obtained from a different piece of the molecule)
    Returns molecular_ion, the sum_formula (vector) of the entire molecular ion
    """
    mol = read_smiles(smiles, explicit_hydrogen=True, zero_order_bonds=True, reinterpret_aromatic=False)
    molecular_ion = graph_to_sum_formula(mol)
    str_formula = chem_data.get_string_formula(molecular_ion)
    print("\nmolecule {}".format(str_formula))
    is_forest = nx.is_forest(mol)
    is_tree = nx.is_tree(mol)
    is_connected = nx.is_connected(mol)
    if verbose:
        print("connected: {} tree: {} forest: {}".format(is_connected, is_tree, is_forest))
    #print("Nodes: {} Edges: {}".format(mol.nodes(data='element'),mol.edges()))
    
    # compute sub-graphs
    all_connected_subgraphs = get_subgraphs(mol, verbose=True)
    all_connected_subgraphs.append(mol)
    # compute sum formula
    all_sum_formulas = [graph_to_sum_formula(gi) for gi in all_connected_subgraphs]
    all_sum_formulas.sort()
    distinct_sum_formulas = [all_sum_formulas[0]]
    for i in range(1,len(all_sum_formulas)):
        if all_sum_formulas[i-1] != all_sum_formulas[i]:
            distinct_sum_formulas.append(all_sum_formulas[i])
    
    sum_formulas = [(chem_data.get_mass(a),a) for a in distinct_sum_formulas]
    sum_formulas.sort()
    mass_ambiguity = ""
    (si,ai) = sum_formulas[0]
    ci = chem_data.get_string_formula(ai)
    s = "{}:{:.2f}".format(ci,si)
    rj = round(si)
    cj = ci
    for i in range(1,len(sum_formulas)):
        (si,ai) = sum_formulas[i]
        ci = chem_data.get_string_formula(ai)
        s += ", {}:{:.2f}".format(ci,si)
        ri = round(si)
        if ri == rj:
            mass_ambiguity += "({},{},{}) ".format(cj,ci,ri)
        rj = ri
        cj = ci
    if verbose:
        print(s)
        #if len(mass_ambiguity) > 0:
        #    print("+++++++++++++++++++++++++++++++++++++++++++++++")
        #    print("integer_mass_ambiguity: "+mass_ambiguity)
        print("")
    return sum_formulas, molecular_ion

def get_rare_isotopologues_all_fragments(fragments_by_mass):
    """
    fragments_by_mass: a list of pairs (mass (float), sum_formula (as vector)) that are monoisotopic
    comment myriam 20201006: monoisotopic means that there are no exisitng rare isotopes (e.g., F, I have no rare isotopes)
    Do you mean the sum formula only have abundant isotopes? (e.g. C and no [13C])

    Returns a list of lists of all isotopologues, each one is a triple (mass (float), sum_formula (as vector), relative intensity (float))
    Returns the maximum mass of everything (float)
    """
    # generate isotopologues of fragments as sets of isotopologues, with their relative intensity
    mass_max = 0
    all_fragments_and_iso = [None] * len(fragments_by_mass)
    i = 0
    for (mass, frag_abundant) in fragments_by_mass:
        iso_set = generate_all_isotopolog_fragments(frag_abundant)
        # the abundant one is not in the set
        # compute relative intensity for each fragment, wrt the one made of abundant atoms only
        mass_iso = [(mass,frag_abundant, float(1.0))] + [(chem_data.get_mass(fi),fi, p_iso_rare_compute(frag_abundant, fi)) for fi in iso_set]
        mass_iso.sort() # sort by increasing mass
        mass_max = max(mass_max, mass_iso[-1][0])
        all_fragments_and_iso[i] = mass_iso
        i += 1
    return all_fragments_and_iso, mass_max

def draw_graph_theoric_fragments(all_fragments_and_iso_by_mass, mass_max, str_sum_formula, mol_name, block=None):
    """
    draw a graph with maplotlib
    x-axis: masses
    y-axis: relative intensity
    Each set of isotopologues of a fragment is in a distinct color
    The relative intensity of the fragment made of abundant atoms only is set to 1.0 and the others have positive relative intensity (can be > 1.0)
    """
    # matplotlib.colors.cnames
    #colornames = list(matplotlib.colors.cnames.keys())
    colornames = ['aqua', 'blue', 'cadetblue', 'cornflowerblue', 'cyan', 'darkblue', 'darkcyan', 'darkslateblue', 'darkturquoise', 'dodgerblue', 'indigo', 'mediumblue', 'mediumslateblue', 'midnightblue', 'navy', 'royalblue', 'slateblue', 'steelblue']
    if mol_name == str_sum_formula:
        label_fig = mol_name
    else:
        label_fig = mol_name + " ("+str_sum_formula+")"
    fig, ax1 = plt.subplots(figsize=(8,4))
    plt.xlabel('Exact mass [m/z] of theoretical fragments of '+label_fig)
    plt.ylabel('Relative intensity')
    ax1.ticklabel_format(axis = 'both', style = 'plain', useMathText=False)
    ax1.grid(True, color='xkcd:grey', linestyle='--', linewidth=0.5)
    i = 0
    for mass_iso in all_fragments_and_iso_by_mass:
        # draw graph
        frag_abundant = mass_iso[0][1]
        mass = mass_iso[0][0]
        col = colornames[i % len(colornames)]
        plt.bar([x[0] for x in mass_iso], [x[2] for x in mass_iso], 0.5, color=col)
        t = ax1.text(mass, 1.0 + 0.05*(i % 2), chem_data.get_string_formula(frag_abundant))
        i += 1
    plt.xlim(0.0,mass_max+5)
    plt.show(block=block)
    #plt.close()
    return label_fig

def draw_graph_measured_masses(measured_masses, measured_intensity, mass_max, label_fig,block=None):
    """
    Draw a graph, x-axis: masses from 0 to mass_max, y-axis: intensity
    """
    fig2, ax2 = plt.subplots(figsize=(8,4))
    plt.xlim(0.0,mass_max+5)
    plt.xlabel('Measured mass [m/z] of TOF fragments of '+label_fig)
    plt.ylabel('Intensity')
    ax2.ticklabel_format(axis = 'both', style = 'plain', useMathText=False)
    ax2.grid(True, color='xkcd:grey', linestyle='--', linewidth=0.5)
    plt.bar(measured_masses, measured_intensity, 0.5, color='indigo')
    plt.show(block=block)

def get_filename(pattern: str, mol_name: str, str_sum_formula: str) -> Path:
    """
    given a pattern string, search for a file with that pattern in the name, in
    data/nontarget_screening/fragments
    Returns filename or None
    """
    current_working_dir = Path.cwd()
    fdir = current_working_dir/'data'/'nontarget_screening'/'fragments'
    list_files = []
    if pattern != None:
        list_files = list(fdir.glob('*'+pattern+'.txt'))
    
    if pattern == None or len(list_files) == 0:
        print("No valid filename pattern provided for {} ({}) (filename pattern {} given), trying sum formula".format(mol_name, str_sum_formula, pattern))
        list_files = list(fdir.glob('*'+mol_name+'.txt'))
        if len(list_files) == 0:
            print("No valid filename matching {} found".format(mol_name))
            if mol_name != str_sum_formula:
                list_files = list(fdir.glob('*'+str_sum_formula+'.txt'))
            if len(list_files) == 0:
                #raise FileNotFoundError("file for molecule {} ({}) was not found".format(mol_name, str_sum_formula))
                print("file for molecule {} ({}) was not found".format(mol_name, str_sum_formula))
    if len(list_files) >=1 :
        filename = list_files[0]
    else:
        filename = None
    return filename

def match_masses_with_fragments(list_sum_formulas, all_iso_fragments):
    """
    all_iso_fragments a list, each item is a set of isotopologues with mass, sum formula, relative intensity.
    list_sum_formulas: a list of sum_formulas (as vectors) to match
    1. do a list with all iso of all fragments, ordered by increasing mass (use merge sort)
    2. order by increasing mass the list of sum formulas
    3. read onward a list, downward the other and match.
    """

    return

def match_fragfromgraph_with_fragfromidentifcation(list_sum_formulas, all_iso_fragments):
    """
    This function Myriam 20201006
    Trying to compare fragments reconstituted by identification routine with list of fragments from graph.
    Fragments resulting from re-arrangement are counted as false (they are not on the graph).
    This function is probably a draft of match_masses_with_fragments

    Returns
    -------
    None.

    """
    return

def match_fragfromsumformula_with_fragfromidentifcation():
    """
    This function Myriam 20201006    
    Trying to compare fragments reconstituted by identification routine with list of fragments from correct sum formula.
    This allows to count re-arrangements as true (as long as they are a subset of the sumformula).
    Returns
    -------
    None.

    """
    return

def find_validation_idx_from_database(fragfile_name, mol_validation, validation_labels):

    val_name = fragfile_name.split('_')[1]
    #print(val_name)  
    found_val = False
    idx_val = -1

    while idx_val < len(mol_validation)-1 and found_val == False:
        idx_val += 1
        
        i_aliases = 0
        while i_aliases < len(mol_validation[idx_val][validation_labels["Aliases"]]) and found_val == False:
            found_val = val_name == mol_validation[idx_val][validation_labels["Aliases"]][i_aliases]
            i_aliases += 1

    
    if found_val == True:
        print("Compound: " + str(val_name) + "; Database: " + str(mol_validation[idx_val][6]))
        print(mol_validation[idx_val])
    else:
        print("Compound " + str(val_name) + " not found in database, validation not possible. Please check your database.")
        idx_val = None
    
    return idx_val
    
    
def validation(smiles, mol_name, path_file: Path, draw_graph=True):
    """
    smiles: a string that is a valid SMILES code
    mol_name: a string that is close to a sum formula
    path_file: a valid Path filename of experimental data where to read data, or None
    """
    fragments_by_mass, molecular_ion = get_mass_sum_formula_abundant_fragments(smiles)
    all_fragments_and_iso, mass_max = get_rare_isotopologues_all_fragments(fragments_by_mass)

    #str_sum_formula = chem_data.get_string_formula(molecular_ion)
    mol_ion = chem_data.get_fragment_from_string_formula(mol_name)
    if mol_ion != molecular_ion:
        print("Error mol_name = {} mol_ion = {} molecular_ion = {}".format(mol_name, mol_ion, molecular_ion))
        raise ValueError("Error string: mol_name = {} mol_ion = {} molecular_ion = {}".format(mol_name, mol_ion, molecular_ion))
    str_sum_formula = chem_data.get_string_formula(mol_ion)
    if draw_graph:
        label_fig = draw_graph_theoric_fragments(all_fragments_and_iso, mass_max, str_sum_formula, mol_name, block=False)
    # now compare to experimental mass spectrum
    # 1. open file of experimental mass spectrum
    if path_file == None:
        return
    print("path_file = {}".format(str(path_file)))
    batches_fragments = get_data_from_frag_file(path_file)
    # this is a dictionary indexed by comp_bin and each entry is sorted by increasing mass
    frag_data = batches_fragments[0]
    
    # how to match measured masses with masses of possible fragments?
    # for each abundant fragment, consider the most abundant isotopologue, and try to match
    
    #comment myriam 20201006: I am not sure matching measured masses and reconstituted masses from identified fragments 
    #is the best strategy. Indeed, one measured mass can be made up of several fragments, and the apparent measured mass
    #may be slightly different.
    #Maybe it is better to match reconstituted fragments with list of theoretically possible fragments.
    
    measured_masses = [fi[3] for fi in frag_data]
    measured_intensity = [fi[5] for fi in frag_data]
    if draw_graph:
        draw_graph_measured_masses(measured_masses, measured_intensity, mass_max, label_fig,block=False)
    #plt.close('all')

    # measured mass uncertainty
    measured_mass_u = [fi[4] for fi in frag_data]
    # a mass is within the range [m-u,m+u] for m in measured_masses and u in measured_mass_u
    for i in range(len(measured_masses)):
        m = measured_masses[i]
        #U = u_k * math.sqrt(measured_mass_u[i]**2 + measured_mass_u_cal[i]**2)
        #m_max= m + U
        #m_min= m - U
        # try matching the mass
    
    print("")
    # TODO Aurore: matching algorithm
    
    
def validate_frag_list(smiles, mol_name, list_identified_fragments, draw_graph=True):
    """
    Check if a knapsack gnerated fragment formula, still part of the graph at the end,
    is indeed a plausible fragment.
    Here we assume that re-organisations are possible,
    so any fragment that is a sub-fragment of the molecular ion formula is assumed correct.
    This is a simplification.
    smiles: a string that is a valid SMILES code
    mol_name: a string that is close to a sum formula
    path_file: a valid Path filename of experimental data where to read data, or None
    """
    list_checked = [False]*len(list_identified_fragments)
    mol_ion = chem_data.get_fragment_from_string_formula(mol_name)
    str_sum_formula = chem_data.get_string_formula(mol_ion)


    #myriam 20201006 matching identified fragments with theoretically possible fragments
    
    #list of sum formulas of identified fragments
    
# =============================================================================
#     #=============================================================
#     #this would be the test to test if frag is a subgraph
#     #list of theoretically generated fragments from graph
#     #fragments_by_mass, molecular_ion = get_mass_sum_formula_abundant_fragments(smiles)
#     #all_fragments_and_iso, mass_max = get_rare_isotopologues_all_fragments(fragments_by_mass)
#     #str_sum_formula = chem_data.get_string_formula(molecular_ion)
#     #if mol_ion != molecular_ion:
#     #    print("Error mol_name = {} mol_ion = {} molecular_ion = {}".format(mol_name, mol_ion, molecular_ion))
#     #    raise ValueError("Error string: mol_name = {} mol_ion = {} molecular_ion = {}".format(mol_name, mol_ion, molecular_ion))
#     #if draw_graph:
#     #    label_fig = draw_graph_theoric_fragments(all_fragments_and_iso, mass_max, str_sum_formula, mol_name, block=False)
#
#     #fragments_by_mass: a list of pairs (mass (float), sum_formula (as vector)) that are monoisotopic
#     #print("len(list_identified_fragments) " + str(len(list_identified_fragments)))
#     for i_frag in range(len(list_identified_fragments)):
#         identified = False
#         idx_test = -1
#         while idx_test < len(fragments_by_mass)-1 and identified == False:
#             idx_test += 1
#             #print(idx_test)
#             #print(list_identified_fragments[i_frag])
#             #print(fragments_by_mass[idx_test])
#             
#             identified = list_identified_fragments[i_frag] == fragments_by_mass[idx_test][1]
#         
#         list_checked[i_frag] = identified
#         #print(list_checked)
# =============================================================================

    for i_frag in range(len(list_identified_fragments)):
        identified = False
        less_than_mol_ion = True
        i_atom = 0
        while i_atom < len(mol_ion) and less_than_mol_ion == True:
            less_than_mol_ion = list_identified_fragments[i_frag][i_atom] <= mol_ion[i_atom]
            i_atom += 1
            #print(idx_test)
            #print(list_identified_fragments[i_frag])
            #print(fragments_by_mass[idx_test])
            
            #identified = list_identified_fragments[i_frag] == fragments_by_mass[idx_test][1]
        
        list_checked[i_frag] = less_than_mol_ion
        #print(list_checked)
    
    return list_checked

def compute_validation_results(list_checked, G, sum_signal, measured_sum_signal):
    """

    Parameters
    ----------
    list_checked : TYPE
        DESCRIPTION.
    
    measured_sum_signal: this is cl_comp.sum_I, type float.

    Returns
    -------
    None.

    """
    nodes_maximal_fragments = G.graph['maximal_nodes']
    
    fragment_results = []
    validation_dict_results = {}
    validation_dict_results["frac_correct_frag"]  = 0
    validation_dict_results["frac_true_signal_to_total_signal"] = 0
    validation_dict_results["frac_true_frag_top10"] = 0
    validation_dict_results["frac_true_signal_top10"] = 0

    list_likelihood_true_frag = []
    list_likelihood_false_frag = []
    list_likelihood_true_att = []
    list_likelihood_false_att = []
    list_ranking_true_frag = []
    list_ranking_false_frag = []
    list_ranking_true_att = []
    list_ranking_false_att = []

    #for article: create list with
    #frag formula, Exact mass, Assigned signal, Likelihood, ranking & Max. element
    rank = 0
    for key_node in G.nodes_ordered_by_likelihood:
        rank += 1
        #if key_node in nodes_maximal_fragments
        #is_maximal_frag = 'False'
        if key_node in nodes_maximal_fragments:
            is_maximal_frag = 'True'
        else:
            is_maximal_frag = 'False'

        node = G.nodes[key_node]['node_lite']
        iso = node.iso
        fragment_results.append([iso.frag_abund_formula,
                                 list_checked[rank-1],
                                 iso.frag_abund_mass,
                                 iso.iso_list_sum_signal/measured_sum_signal*100.0,
                                 node.subgraph_sum_I/measured_sum_signal*100.0,
                                 node.subgraph_likelihood,
                                 rank,
                                 is_maximal_frag])

                     
    frac_correct_frag = sum([1 for i in range(len(list_checked)) if list_checked[i] == True ])/len(list_checked)
    print("frac_correct_frag: " + "{:.2f}".format(frac_correct_frag))



    #fraction of correct signal compared to total signal                
    frac_true_signal_to_total_signal = 0
    
    #now iterating over nodes starting by most likely, so this is rank 1
    i_rank = 1
    for idx_node in range(len(G.nodes_ordered_by_likelihood)):

        node = G.nodes[G.nodes_ordered_by_likelihood[idx_node]]['node_lite']
        if list_checked[idx_node] == True:
            list_likelihood_true_frag.append(node.subgraph_likelihood)
            list_ranking_true_frag.append(i_rank)
            frac_true_signal_to_total_signal += node.iso.iso_list_sum_signal

            if node.unique_id in nodes_maximal_fragments:
                list_likelihood_true_att.append(node.subgraph_likelihood)

                list_ranking_true_att.append(i_rank)
            #for i_iso in range(iso.iso_list_len):
            #    frac_true_signal_to_total_signal += iso.iso_list_I_opt[i_iso]

        else:

            list_likelihood_false_frag.append(node.subgraph_likelihood)
            list_ranking_false_frag.append(i_rank)
            if node.unique_id in nodes_maximal_fragments:
                list_likelihood_false_att.append(node.subgraph_likelihood)

                list_ranking_false_att.append(i_rank)
                
        i_rank += 1

    #print("no. att.: " + str(len(nodes_maximal_fragments)) + "\t" + "no. true att. " + str(len(list_ranking_true_att)) +"\t" + "no. false att. " + str(len(list_ranking_false_att)))
    frac_true_signal_to_total_signal /= sum_signal #cl_comp.sum_I
    print("Frac true signal: " + "{:.2f}".format(frac_true_signal_to_total_signal))

    #For top-10 likelihood data:
    frac_true_frag_top10 = 0
    frac_true_signal_top10 = 0
    sum_signal_top10 = 0
            
    for idx_node in range(min(10, len(G.nodes_ordered_by_likelihood))):

        iso = G.nodes[G.nodes_ordered_by_likelihood[idx_node]]['node_lite'].iso

        sum_signal_top10 += iso.iso_list_sum_signal

        if list_checked[idx_node] == True:
            frac_true_frag_top10 += 1
            frac_true_signal_top10 += iso.iso_list_sum_signal
    frac_true_frag_top10 /= float(min(10, len(G.nodes_ordered_by_likelihood)))
    frac_true_signal_top10 /= sum_signal_top10
                
    validation_dict_results["frac_correct_frag"]  = frac_correct_frag
    validation_dict_results["frac_true_signal_to_total_signal"] = frac_true_signal_to_total_signal
    validation_dict_results["frac_true_frag_top10"] = frac_true_frag_top10
    validation_dict_results["frac_true_signal_top10"] = frac_true_signal_top10            


    validation_dict_results["list_likelihood_true_frag"] = list_likelihood_true_frag
    validation_dict_results["list_likelihood_false_frag"] = list_likelihood_false_frag
    validation_dict_results["list_likelihood_true_att"] = list_likelihood_true_att
    validation_dict_results["list_likelihood_false_att"] = list_likelihood_false_att
    validation_dict_results["list_ranking_true_frag"] = list_ranking_true_frag
    validation_dict_results["list_ranking_false_frag"] = list_ranking_false_frag
    validation_dict_results["list_ranking_true_att"] = list_ranking_true_att
    validation_dict_results["list_ranking_false_att"] = list_ranking_false_att

                
    return validation_dict_results, fragment_results




