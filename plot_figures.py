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


Created on Tue Aug 18 17:01:47 2020
"""
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
try:
    import pydot
    from networkx.drawing.nx_pydot import graphviz_layout
    has_pydot = True
except:
    has_pydot = False

from pathlib import Path
import io_tools
from runtype import RunType

from periodic_table import nice_formula_from_frag_formula
#import re

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
                'xkcd:mauve', 'xkcd:dark red', 'xkcd:taupe', 'xkcd:cerulean', 'xkcd:eggplant']
color2_names = ['xkcd:sky blue', 'xkcd:teal', 'xkcd:peach', 'xkcd:mauve', 'xkcd:dark red', 'xkcd:taupe', 'xkcd:cerulean', 'xkcd:eggplant']

#grey scale
#color3_names = ['xkcd:grey', 'xkcd:slate', 'xkcd:black', 'xkcd:grey', 'xkcd:slate', 'xkcd:black', 'xkcd:grey', 'xkcd:slate', 'xkcd:black']

color3_names = ['xkcd:dusty pink', 'xkcd:grey', 'xkcd:wine red', 'xkcd:grey purple', 'xkcd:dark rose', 
                'xkcd:dark salmon', 'xkcd:slate', 'xkcd:black', 'xkcd:eggplant', 'xkcd:taupe', 'xkcd:mauve']

def plot_mass_defect_of_fragments(frag_data: dict, toffile_name: str, fig_file: Path=None, show: bool=True)->None:
    """plot mass defect vs RT, to pre-group fragments before identification

    INPUT:
    - ``frag_data``: a dictionary indexed by bin_number (see io_tools.py), this is to label the graph
    - ``toffile_name``: a string filename without extension (i.e. path_file.stem with Path)
    - ``fig_file``: a Path filename where to save the graph (as .png), if None, do not save the graph
    - ``show``: a boolean to control the figure
    """
    
    
    #plot mass defect vs RT:
    fig, ax1 = plt.subplots(figsize=(8,4))
    fig.suptitle(toffile_name + ': Mass defect vs RT', fontsize=12)
    plt.xlabel('RT [s]')
    plt.ylabel('Mass defect [m/z]')
    ax1.ticklabel_format(axis = 'both', style = 'plain', useMathText=False)
    ax1.grid(True, color='xkcd:grey', linestyle='--', linewidth=0.5)

    for j in list(frag_data):
        if j == -1:
            continue
        frag_all_mass_defect = [fji[io_tools.idx_m] - int(round(fji[io_tools.idx_m])) for fji in frag_data[j]]
        #mass defect vs mass:
        #ax1.scatter([frag_data[j][i][io_tools.idx_m] for i in range(len(frag_data[j]))], frag_all_mass_defect, c = str(color1_names[j])) #, s = 2* isotopologues_list[i_list].iso_list_I_opt[i_iso]*100.00 / cl_comp.meas_I_max
        #mass defect/mass vs mass:
        #ax1.scatter([frag_data[j][i][io_tools.idx_m] for i in range(len(frag_data[j]))], [frag_all_mass_defect[i]/frag_data[j][i][io_tools.idx_m] for i in range(len(frag_data[j]))], c = str(color1_names[j])) #, s = 2* isotopologues_list[i_list].iso_list_I_opt[i_iso]*100.00 / cl_comp.meas_I_max
        #mass defect vs RT:
        #ax1.scatter([frag_data[j][i][io_tools.idx_rt] for i in range(len(frag_data[j]))], frag_all_mass_defect, c = str(color1_names[j])) #, s = 2* isotopologues_list[i_list].iso_list_I_opt[i_iso]*100.00 / cl_comp.meas_I_max
        #mass defect/mass vs RT:
        ax1.scatter([fji[io_tools.idx_rt] for fji in frag_data[j]], [frag_all_mass_defect[i]/frag_data[j][i][io_tools.idx_m] for i in range(len(frag_data[j]))], c = str(color1_names[j])) #, s = 2* isotopologues_list[i_list].iso_list_I_opt[i_iso]*100.00 / cl_comp.meas_I_max

    ax1.set_ylim(-0.002, 0.002)
    if fig_file is not None:
        plt.savefig(fig_file, bbox_inches='tight', transparent = True, dpi = 600)
    if show == True:
        plt.show()
    plt.close('all')

def plot_graph_of_frag(G, show_plt_figures = True):

    """
    ['bipartite_layout',
     'circular_layout',
     'fruchterman_reingold_layout',
     'kamada_kawai_layout',
     'multipartite_layout',
     'planar_layout',
     'random_layout',
     'rescale_layout',
     'shell_layout',
     'spectral_layout',
     'spiral_layout',
     'spring_layout']
    """
    fig = plt.subplots(figsize=(8,4))                
    label_formula = {}
    node_color_values = [0]*G.number_of_nodes()
    idx= 0
    for key_node in G.nodes():
        node = G.nodes[key_node]['node_lite']
        label_formula[key_node] = node.iso.frag_abund_formula
        node_color_values[idx] = node.iso.frag_abund_mass
        idx+=1
    pos=nx.spectral_layout(G) # positions for all nodes
    nx.draw_networkx_labels(G, pos, labels=label_formula)
    #node_size = [0.0005 * nx.get_node_attributes(G, 'population')[v] for v in G]
    nx.draw_networkx_nodes(G, pos, cmap=plt.get_cmap('jet'), alpha=0.3,
                            node_size = 500, node_color = node_color_values) #,
    #nx.draw(G_of_frag, pos, node_size=50, vmin=0.0, vmax=5.0, labels=label_formula, with_labels = True, font_size = 8) #, font_weight='bold'
    nx.draw_networkx_edges(G, pos, alpha=0.5, arrows=True)
    if show_plt_figures:
        plt.show()
    plt.close('all')
    
    return None

def plot_forest_of_frag_with_pydot(G, show_plt_figures = True, graphviz_layout_prog="dot"):
    #fig = plt.subplots(figsize=(24,4))
    if not has_pydot:
        print("Package pydot not installed, please install it first. Continuing...")
        return
    fig = plt.figure(figsize=(24,4))
    label_formula = {}
    node_color_values = [0]*G.number_of_nodes()
    idx= 0
    for key_node in G.nodes():
        label_formula[key_node] = G.nodes[key_node]['str_formula']
        node_color_values[idx] = G.nodes[key_node]['attribute']
        idx+=1
    label_edges = {}
    for key_edge in G.edges():
        label_edges[key_edge] = G.edges[key_edge]['label']
    # 'dot', 'twopi', 'fdp', 'sfdp', 'circo'
    pos = graphviz_layout(G, prog=graphviz_layout_prog)
    nx.draw_networkx_labels(G, pos, labels=label_formula, font_size=6)
    nx.draw_networkx_nodes(G, pos, cmap=plt.get_cmap('jet'), alpha=0.3,
                            node_size = 500, node_color = node_color_values)# , font_size=6)
    nx.draw_networkx_edge_labels(G, pos, edge_labels=label_edges, label_pos=0.5, font_size=6)
    nx.draw_networkx_edges(G, pos, alpha=0.5, arrows=True)
    if show_plt_figures:
        fig.tight_layout()
        plt.show()
    # for saving, consider savefig(filename, bbox_inches=0) or bbox_inches='tight'

def plot_graph_with_pydot(G, label_nodes, label_edges, node_colors, show_plt_figures = True, graphviz_layout_prog="dot"):
    if not has_pydot:
        print("Package pydot not installed, please install it first. Continuing...")
        return
    fig = plt.figure(figsize=(24,4))
    # 'dot', 'twopi', 'fdp', 'sfdp', 'circo'
    pos = graphviz_layout(G, prog=graphviz_layout_prog)
    nx.draw_networkx_labels(G, pos, labels=label_nodes, font_size=6)
    nx.draw_networkx_nodes(G, pos, cmap=plt.get_cmap('jet'), alpha=0.3,
                            node_size = 500, node_color = node_colors)# , font_size=6)
    nx.draw_networkx_edge_labels(G, pos, edge_labels=label_edges, label_pos=0.5, font_size=6)
    nx.draw_networkx_edges(G, pos, alpha=0.5, arrows=True)
    if show_plt_figures:
        fig.tight_layout()
        plt.show()
    # for saving, consider savefig(filename, bbox_inches=0) or bbox_inches='tight'

def plot_mass_spectrum(run_type:RunType, rt_average:float, results_found_list:list, results_not_found_list:list, electron_mass:float, meas_I_max:float, suptitle:str, fig_file:str, show: bool=True):

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
 

    #MASS SPECTRUM
    #max_volume = meas_I_max
    #print('Retention time: ' + str(rt_centre_weighted))
    #rt_centre_weighted_std_dev = 2*np.sqrt(sum([(mass_spectrum_list[i][3]*(mass_spectrum_list[i][0] - rt_centre_weighted))**2 for i in range(len(mass_spectrum_list))]) / sum([(mass_spectrum_list[i][3])**2 for i in range(len(mass_spectrum_list))]))
    
    #Plot of extracted mass spectrum
    #to plot relative intensities:
    #y_positions_rel = chromato_res_np[indexes, ic_area]/meas_I_max*100 + 0.7
    #to plot absolute intensities:
    #y_positions = [line[3] + 0.7 for line in mass_spectrum_list]
    #if idx_rt_group == 3:
    #    y_positions[5] +=3
    #y_positions[13] -=3
    #y_positions[12] +=2
    #y_positions[17] +=3
    #y_positions[15] +=5
    #y_positions[2] +=4
    #y_positions[11] +=5

    text_size = 12
    text_threshold = 10
    fig, ax1 = plt.subplots(figsize=(8.5,4))
    ax1.tick_params(axis='x', labelsize = text_size)
    ax1.tick_params(axis='y', labelsize = text_size)

    fig.suptitle(suptitle) # , fontsize=12
    #fig.suptitle("Reconstructed mass spectrum for carbon tetrachloride")
    if run_type.is_NIST():
        ms_color = 'xkcd:scarlet'
    else:
        ms_color = 'xkcd:light orange'
        ms_color_EI = 'xkcd:light orange'
        ms_color_CI = 'xkcd:bright blue'
        
        
    plt.ylabel('Intensity, normalised', fontsize=text_size) #plt.ylabel('Volume of peak, normalised')
    plt.xlabel('Measured mass [m/z]', fontsize=text_size)
    ax1.grid(True, color='xkcd:grey', linestyle='--', linewidth=0.5)
    #ax1.set_xlim(0, 110)
    ax1.set_ylim(0, 110)
    width = 0.6
    #width = 0.003 * (max([results_list[i][ires_m_exact] for i in range(len(results_list))]) - min([results_list[i][ires_m_exact] for i in range(len(results_list))]))

    
    #to plot normalised values:

    for i in range(len(results_found_list)):
        if results_found_list[i][ires_ionisation] == 'EI':
            ms_color = ms_color_EI
        elif results_found_list[i][ires_ionisation] == 'CI':
            ms_color = ms_color_CI
        plt.bar(results_found_list[i][ires_m_exact]-electron_mass, results_found_list[i][ires_I_frag]/meas_I_max*100, width, color = ms_color, alpha=0.5) #xkcd:crimson, yerr=results_found_list[i]*2/meas_I_max*100)
        if results_found_list[i][ires_I_frag]/meas_I_max*100 > text_threshold:
            #t = ax1.text(results_found_list[i][ires_m_exact], results_found_list[i][ires_I_frag]/meas_I_max*100 + 0.7, '{:.3f}'.format(results_found_list[i][ires_m_exact]), fontsize=text_size)
            #t = ax1.text(results_found_list[i][ires_m_exact], results_found_list[i][ires_I_frag]/meas_I_max*100 + 0.7, str(int(results_found_list[i][ires_m_exact])), fontsize=text_size)
            #write integer masses:
            #t = ax1.text(results_found_list[i][ires_m_exact]-electron_mass, results_found_list[i][ires_I_frag]/meas_I_max*100 + 0.7,  str(round(results_found_list[i][ires_m_exact]-electron_mass)), fontsize=text_size)
            #write mass only:
            #t = ax1.text(results_found_list[i][ires_m_exact]-electron_mass, results_found_list[i][ires_I_frag]/meas_I_max*100 + 0.7,  '{:.4f}'.format(results_found_list[i][ires_m_exact]-electron_mass), fontsize=text_size)
            #write mass + formula:       
            text_formula = nice_formula_from_frag_formula(str(results_found_list[i][ires_formula]))
            #write formula and exact mass
            t = ax1.text(results_found_list[i][ires_m_exact]-electron_mass, results_found_list[i][ires_I_frag]/meas_I_max*100 + 0.7, text_formula + '\n' + '{:.3f}'.format(results_found_list[i][ires_m_exact]-electron_mass), fontsize=text_size)
            #write formula:
            #t = ax1.text(results_found_list[i][ires_m_exact]-electron_mass, results_found_list[i][ires_I_frag]/meas_I_max*100 + 0.7, text_formula, fontsize=text_size)

    for i in range(len(results_not_found_list)):
        if results_not_found_list[i][ires_idx_m] == -1:
            #expected isotopologue, not found, exact mass is known.
            plt.bar(results_not_found_list[i][ires_m_exact]-electron_mass, results_not_found_list[i][ires_I_frag]/meas_I_max*100, width, color = 'xkcd:taupe', alpha=0.5) #, yerr=results_found_list[i]*2/meas_I_max*100)
            #if results_not_found_list[i][ires_I_frag]/meas_I_max*100 > text_threshold:
                #t = ax1.text(results_not_found_list[i][ires_m_exact], results_not_found_list[i][ires_I_frag]/meas_I_max*100 + 0.7, '{:.3f}'.format(results_not_found_list[i][ires_m_exact]), fontsize=text_size)

        else:
            #measured mass but unidentified, measured mass is known.
            plt.bar(results_not_found_list[i][ires_m_meas]-electron_mass, results_not_found_list[i][ires_I_frag]/meas_I_max*100, width, color = 'xkcd:terracotta', alpha=0.5) #, yerr=results_found_list[i]*2/meas_I_max*100)                 
            #if results_not_found_list[i][ires_I_frag]/meas_I_max*100 > text_threshold:
                #t = ax1.text(results_not_found_list[i][ires_m_meas], results_not_found_list[i][ires_I_frag]/meas_I_max*100 + 0.7, '{:.3f}'.format(results_not_found_list[i][ires_m_meas]), fontsize=text_size)
    
    #ax1.set_ylim(0, 10)
    #ax1.set_yscale('log')
    
    #ax1.legend(['Calculated, type A', 'Calculated, type A + B', 'Measured after identification'])
    plt.savefig(fig_file, bbox_inches='tight', transparent = True, dpi = 600)
    if show == True:
        plt.show()
    plt.close('all')
    
    return None

    """            
    #plot mass vs RT for all found:
    if run_type.is_tofwerk_tof():
        fig, ax1 = plt.subplots(figsize=(8,4))
        fig.suptitle(suptitle, fontsize=12)
        plt.xlabel('RT [s]')
        plt.ylabel('Mass [m/z]')
        ax1.ticklabel_format(axis = 'both', style = 'plain', useMathText=False)
        ax1.grid(True, color='xkcd:grey', linestyle='--', linewidth=0.5)
        idx_color = 0
        for i_list in range(len(isotopologues_list)):
            for i_iso in range(len(isotopologues_list[i_list].iso_list)):
                idx_m = isotopologues_list[i_list].meas_mass_idx[i_iso]
                candidate_mass = chem_data.get_mass(isotopologues_list[i_list].iso_list[i_iso])
                if idx_m is not None:
                    #ax1.plot(results_found_list[i][ires_rt], results_found_list[i][ires_m_exact], '.', color = 'xkcd:blue')
                    ax1.scatter(cl_comp.meas_rt[idx_m], candidate_mass, c = str(color1_names[idx_color]), s = 2* isotopologues_list[i_list].iso_list_I_opt[i_iso]*100.00 / cl_comp.meas_I_max)
            idx_color += 1
            if idx_color > len(color1_names)-1: idx_color = 0
        #fig.tight_layout()
        #plt.savefig(fig_file, bbox_inches='tight', transparent = True, dpi = 300)
        #plt.show()
    """
    
    """
    #plot calculated mass uncertainty and measured mass diff vs area for identified compounds
    fig, ax1 = plt.subplots(figsize=(8,4))
    fig.suptitle('Uncertainty: reconstructed and measured, vs area', fontsize=12)
    plt.xlabel('Area [V * m/z * s]')
    plt.ylabel('Mass uncertainty [ppm]')
    ax1.ticklabel_format(axis = 'both', style = 'plain', useMathText=False)
    ax1.grid(True, color='xkcd:grey', linestyle='--', linewidth=0.5)
    for i in range(len(results_found_list)):
        
        u_mass_cal = np.interp(results_found_list[i][ires_m_meas], np.array(cl_run.mass_cal_u_exact_mass), np.array(cl_run.mass_cal_u))
        
        ax1.plot(results_found_list[i][ires_I], results_found_list[i][ires_m_meas_u]*2, '.', color = 'xkcd:blue')
        ax1.plot(results_found_list[i][ires_I], results_found_list[i][ires_m_meas_u]*2 + u_mass_cal * 2, '.', color = 'xkcd:green')
        ax1.plot(results_found_list[i][ires_I], abs(results_found_list[i][ires_m_diff]), '.', color = 'xkcd:orange')
        
    ax1.set_xscale('log')
    ax1.legend(['Calculated, type A', 'Calculated, type A + B', 'Measured after identification'])
    fig.tight_layout()
    plt.show()
    """
    
    """
    #plot calculated mass uncertainty vs measured mass diff for identified compounds
    fig, ax1 = plt.subplots(figsize=(6,6))
    fig.suptitle('Uncertainty: reconstructed and measured, vs area', fontsize=12)
    plt.xlabel('Measured mass difference [ppm]')
    plt.ylabel('Calculated mass uncertainty, 2 sigma [ppm]')
    ax1.ticklabel_format(axis = 'both', style = 'plain', useMathText=False)
    ax1.grid(True, color='xkcd:grey', linestyle='--', linewidth=0.5)
    #datax_plot = [d for d in range(0, max(abs(results_list[:][ires_m_diff])), 100)]
    ax1.plot([abs(results_list[i][ires_m_diff]) for i in range(len(results_list))], [abs(results_list[i][ires_m_diff]) for i in range(len(results_list))], linestyle='--', color = 'xkcd:slate' ) 
    
    ax1.plot(color = 'xkcd:slate')
    for i in range(len(results_list)):
        ax1.plot(abs(results_list[i][ires_m_diff]), results_list[i][ires_m_meas_u]*2, '.', color = 'xkcd:blue')
    fig.tight_layout()
    plt.show()
    """

    #**************************************************
    #PLOT MASS SPECTRUM per set of co-eluting fragments
    #**************************************************
    """
    idx_frag = 0
    idx_comp_bin_old = -1
    while idx_frag < len(chromato_res_np):
        idx_comp_bin = chromato_res_np[idx_frag, ic_comp_bin]

        
        if idx_comp_bin > idx_comp_bin_old:
            
            indexes = np.flatnonzero(chromato_res_np[:, ic_comp_bin] == idx_comp_bin)
            
            if len(indexes)>2:   
                #plot mass spectrum for this compound
                max_volume = max(chromato_res_np[indexes, ic_area])
                rt_centre_weighted = sum(chromato_res_np[indexes, ic_rt] * chromato_res_np[indexes, ic_area]) / sum(chromato_res_np[indexes, ic_area])
                #print('Retention time: ' + str(rt_centre_weighted))
                #rt_centre_weighted_std_dev = 2*np.sqrt(sum([(mass_spectrum_list[i][3]*(mass_spectrum_list[i][0] - rt_centre_weighted))**2 for i in range(len(mass_spectrum_list))]) / sum([(mass_spectrum_list[i][3])**2 for i in range(len(mass_spectrum_list))]))
                
                #Plot of extracted mass spectrum
                #to plot relative intensities:
                y_positions_rel = chromato_res_np[indexes, ic_area]/max_volume*100 + 0.7
                #to plot absolute intensities:
                #y_positions = [line[3] + 0.7 for line in mass_spectrum_list]
                #if idx_rt_group == 3:
                #    y_positions[5] +=3
                #y_positions[13] -=3
                #y_positions[12] +=2
                #y_positions[17] +=3
                #y_positions[15] +=5
                #y_positions[2] +=4
                #y_positions[11] +=5
                
            
                fig, ax1 = plt.subplots(figsize=(8.5,4))
                fig.suptitle(cl_run.toffile_name + ': mass spectrum at ' + '{:.2f}'.format(rt_centre_weighted) + ' s') #, fontsize=12
                plt.ylabel('Intensity, normalised') #plt.ylabel('Volume of peak, normalised')
                plt.xlabel('Measured mass [m/z]')
                ax1.grid(True, color='xkcd:grey', linestyle='--', linewidth=0.5)
                width = 0.0125 * (max(chromato_res_np[indexes, ic_m]) - min(chromato_res_np[indexes, ic_m]))
                i = 0
                #to plot normalised values:
                plt.bar(chromato_res_np[indexes, ic_m], chromato_res_np[indexes, ic_area]/max_volume*100, width, color = 'xkcd:light orange', yerr=chromato_res_np[indexes, ic_area_u]*2/max_volume*100)
                #to plot absolute values:
                #plt.bar(chromato_res_np[indexes, ic_m], chromato_res_np[indexes, ic_area], width, color = 'xkcd:light orange', yerr=chromato_res_np[indexes, ic_area_u]*2)
                
                for index in indexes:
                    if chromato_res_np[index, ic_area]/max_volume*100 > 5:
                        t = ax1.text(chromato_res_np[index, ic_m], chromato_res_np[index, ic_area]/max_volume*100 + 0.7, '{:.3f}'.format(chromato_res_np[index, ic_m]))
                    i+=1 
                #fig.tight_layout() #h_pad=0
                #fig_name = cl_run.toffile_name + '_mass_spectrum_'+'{:.2f}'.format(rt_centre_weighted) + 's' +   '.png'
                #fig_path = dict_path['mass_spectra']
                #fig_file = fig_path/fig_name
                #plt.savefig(fig_file, bbox_inches='tight', transparent = True, dpi = 300)
                plt.show()
    
            idx_frag += max(indexes) - min(indexes) +1
            idx_comp_bin_old = idx_comp_bin
            
        else: idx_frag += 1
    """ 

def plot_method_behaviour(method_performance_results, 
                          list_likelihood_true_frag, 
                          list_likelihood_false_frag, 
                          list_likelihood_true_att, 
                          list_likelihood_false_att,
                          list_ranking_true_frag,
                          list_ranking_false_frag,
                          list_ranking_true_att, 
                          list_ranking_false_att,
                          fig_path):


    
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
    i_met_res_frac_true_frag = 17
    i_met_res_frac_true_signal_to_total = 18
    i_met_res_frac_true_frag_top10 = 19
    i_met_res_frac_true_signal_top10 = 20


                                      
    """
    #plot model behaviour vs number of knapsack solutions
    text_size = 10
    fig, ax1 = plt.subplots(figsize=(8.5,4))
    ax1.tick_params(axis='x', labelsize = text_size)
    ax1.tick_params(axis='y', labelsize = text_size)
    fig.suptitle('Method behaviour') #, fontsize=12   
    plt.ylabel('Percentage of chemical formulae', fontsize=text_size) #plt.ylabel('Volume of peak, normalised')
    plt.xlabel('Number of knapsack-generated chemical formulae', fontsize=text_size)
    ax1.grid(True, color='xkcd:grey', linestyle='--', linewidth=0.5)
    
    x_abcisse = [method_performance_results[i][i_met_res_no_nodes_knapsack] for i in range(len(method_performance_results))]
    
    ax1.plot(x_abcisse, [method_performance_results[i][i_met_res_percent_nodes_LOD] for i in range(len(method_performance_results))], '*', color = 'xkcd:mauve')
    ax1.plot(x_abcisse, [method_performance_results[i][i_met_res_percent_nodes_singletons] for i in range(len(method_performance_results))], '*', color = 'xkcd:grey')
    ax1.plot(x_abcisse, [method_performance_results[i][i_met_res_percent_nodes_not_treated] for i in range(len(method_performance_results))], '*', color = 'xkcd:brick red')
    ax1.plot(x_abcisse, [method_performance_results[i][i_met_res_percent_nodes_validated] for i in range(len(method_performance_results))], '+', color = 'xkcd:blue')
    ax1.set_xscale('log') 
    ax1.set_yscale('log')
    ax1.legend(['< LOD', 'Singletons', 'Not opt.', 'Validated'])
    plt.show()
    """
    
    
    text_size = 10
    fig, ax1 = plt.subplots(figsize=(5,5))
    ax1.tick_params(axis='x', labelsize = text_size)
    ax1.tick_params(axis='y', labelsize = text_size)
    fig.suptitle('Method behaviour') #, fontsize=12   
    plt.ylabel('Number of chemical formulae', fontsize=text_size) #plt.ylabel('Volume of peak, normalised')
    plt.xlabel('Number of knapsack-generated chemical formulae', fontsize=text_size)
    ax1.grid(True, color='xkcd:grey', linestyle='--', linewidth=0.5)
    
    x_abcisse = [method_performance_results[i][i_met_res_no_nodes_knapsack] for i in range(len(method_performance_results))]

    ax1.plot(x_abcisse, [method_performance_results[i][i_met_res_percent_nodes_not_treated] * method_performance_results[i][i_met_res_no_nodes_knapsack] /100 for i in range(len(method_performance_results))], 'o', color = 'xkcd:brick red')
    ax1.plot(x_abcisse, [method_performance_results[i][i_met_res_percent_nodes_LOD] * method_performance_results[i][i_met_res_no_nodes_knapsack] /100 for i in range(len(method_performance_results))], '*', color = 'xkcd:mauve')
    ax1.plot(x_abcisse, [method_performance_results[i][i_met_res_percent_nodes_singletons] * method_performance_results[i][i_met_res_no_nodes_knapsack] /100 for i in range(len(method_performance_results))], 'x', color = 'xkcd:grey')
    ax1.plot(x_abcisse, [method_performance_results[i][i_met_res_percent_nodes_validated] * method_performance_results[i][i_met_res_no_nodes_knapsack] /100 for i in range(len(method_performance_results))], '+', color = 'xkcd:blue')
    ax1.set_xscale('symlog')
    ax1.set_yscale('symlog') 
    #ax1.set_xlim(0.0)
    #ax1.set_ylim(0)
    ax1.legend(['Rejected: Not optimised', 'Rejected: < LOD', 'Rejected: Singletons', 'Validated'])
    
    
    fig_name = 'method_behaviour.png'
    #fig_path = cl_path.dict_path['fragments']
    fig_file = fig_path/fig_name
    plt.savefig(fig_file, bbox_inches='tight', transparent = True, dpi = 600)
    plt.show()
    

    #Plot runtimes

    """    
    i_met_res_runtime_step1_knapsack = 11
    i_met_res_runtime_step2_graph = 12
    i_met_res_runtime_step3_enum_iso = 13
    i_met_res_runtime_step7_opt_loop = 14
    i_met_res_runtime_total = 15
    """
    text_size = 10
    fig, ax1 = plt.subplots(figsize=(5,5))
    ax1.tick_params(axis='x', labelsize = text_size)
    ax1.tick_params(axis='y', labelsize = text_size)
    fig.suptitle('Runtimes') #, fontsize=12   
    plt.ylabel('Runtime / s', fontsize=text_size) #plt.ylabel('Volume of peak, normalised')
    plt.xlabel('Number of knapsack-generated chemical formulae', fontsize=text_size)
    ax1.grid(True, color='xkcd:grey', linestyle='--', linewidth=0.5)
    
    x_abcisse = [method_performance_results[i][i_met_res_no_nodes_knapsack] for i in range(len(method_performance_results))]

    ax1.plot(x_abcisse, [method_performance_results[i][i_met_res_runtime_step1_knapsack] for i in range(len(method_performance_results))], 'o', color = 'xkcd:brick red')
    ax1.plot(x_abcisse, [method_performance_results[i][i_met_res_runtime_step2_graph] for i in range(len(method_performance_results))], '*', color = 'xkcd:mauve')
    ax1.plot(x_abcisse, [method_performance_results[i][i_met_res_runtime_step3_enum_iso] for i in range(len(method_performance_results))], 'x', color = 'xkcd:grey')
    ax1.plot(x_abcisse, [method_performance_results[i][i_met_res_runtime_step7_opt_loop] for i in range(len(method_performance_results))], '+', color = 'xkcd:blue')
    ax1.plot(x_abcisse, [method_performance_results[i][i_met_res_runtime_total] for i in range(len(method_performance_results))], '+', color = 'xkcd:black')

    ax1.set_xscale('symlog')
    ax1.set_yscale('symlog') 
    #ax1.set_xlim(0.0)
    #ax1.set_ylim(0)
    ax1.legend(['Step1: knapsack', 'Step2: graph', 'Step3: enum. iso.', 'Step7: opt.', 'Total'])
    
    
    fig_name = 'method_runtimes.png'
    #fig_path = cl_path.dict_path['fragments']
    fig_file = fig_path/fig_name
    plt.savefig(fig_file, bbox_inches='tight', transparent = True, dpi = 600)
    plt.show()
    
    

    
    """
    #**********************************    
    #2D-histogram of fractions of true signal vs fraction of true fragments
    n_results=len(method_performance_results)

    n_bins = 20
    fig, ax = plt.subplots(figsize=(4,4)) #, tight_layout=True
    hist = ax.hist2d([method_performance_results[i][i_met_res_frac_true_frag] for i in range(len(method_performance_results))], [method_performance_results[i][i_met_res_frac_true_signal_to_total] for i in range(len(method_performance_results))], bins = n_bins)

    plt.xlabel('Fraction of true fragments', fontsize=text_size) #plt.ylabel('Volume of peak, normalised')
    plt.ylabel('Fraction of true signal', fontsize=text_size)
    
    #cbar = plt.colorbar()
    #cbar.set_label('Fraction of results')
    
    plt.show()
    """

    """
    #**********************************
    #plot histograms of correct signal and correct fragments
    fig, ax1 = plt.subplots(figsize=(5,4)) #, tight_layout=True
    fig.suptitle('Method accuracy') #, scenario "sum formula"
    #axs.grid(True, color='xkcd:grey', linestyle='--', linewidth=0.5)    
    # We can set the number of bins with the `bins` kwarg
    #weights_met_res_frac_true_frag = 
    bins = np.linspace(0, 1, num = 20)
    ax1.hist([method_performance_results[i][i_met_res_frac_true_frag] for i in range(len(method_performance_results))], bins, color = 'xkcd:grey', alpha =0.5) #density = True
    ax1.hist([method_performance_results[i][i_met_res_frac_true_signal_to_total] for i in range(len(method_performance_results))], bins, color = 'xkcd:salmon', alpha =0.5)
    ax1.set_xlabel('Fraction of correct fragments or signal', fontsize=text_size) #plt.ylabel('Volume of peak, normalised')
    #axs[1].set_xlabel('Fraction of correct signal', fontsize=text_size)
    ax1.set_ylabel('Number of compounds', fontsize=text_size)
    ax1.grid(True, color='xkcd:grey', linestyle='--', linewidth=0.5)
    ax1.legend(['fragments', 'signal'])
    #axs[1].grid(True, color='xkcd:grey', linestyle='--', linewidth=0.5)

    fig_name = 'method_accuracy_scenario_sum_formula.png'
    #fig_path = cl_path.dict_path['fragments']
    fig_file = fig_path/fig_name
    plt.savefig(fig_file, bbox_inches='tight', transparent = True, dpi = 600)

    plt.show()
    """

    #**********************************
    #plot histograms of correct signal and correct fragments, for total and top 10
    fig, axs = plt.subplots(figsize=(6,3.5), nrows=1, ncols=2, sharex=True, sharey = True)
    fig.suptitle('Method performance: fraction of correct results') #, scenario "sum formula"
    #axs.grid(True, color='xkcd:grey', linestyle='--', linewidth=0.5)    
    # We can set the number of bins with the `bins` kwarg
    #weights_met_res_frac_true_frag = 
    bins = np.linspace(-0.01, 1.005, num = 20)

    ax1 = axs[0]

    ax1.hist([method_performance_results[i][i_met_res_frac_true_frag] for i in range(len(method_performance_results))], bins, color = 'xkcd:grey', alpha =0.5) #density = True
    ax1.hist([method_performance_results[i][i_met_res_frac_true_signal_to_total] for i in range(len(method_performance_results))], bins, color = 'xkcd:salmon', alpha =0.5)
    ax1.set_xlabel('Correct fraction, all results', fontsize=text_size) #plt.ylabel('Volume of peak, normalised')
    #axs[1].set_xlabel('Fraction of correct signal', fontsize=text_size)
    ax1.set_ylabel('Number of compounds', fontsize=text_size)
    ax1.grid(True, color='xkcd:grey', linestyle='--', linewidth=0.5)
    ax1.legend(['fragments', 'signal'])
    
    ax1 = axs[1]

    ax1.hist([method_performance_results[i][i_met_res_frac_true_frag_top10] for i in range(len(method_performance_results))], bins, color = 'xkcd:grey', alpha =0.5) #density = True
    ax1.hist([method_performance_results[i][i_met_res_frac_true_signal_top10] for i in range(len(method_performance_results))], bins, color = 'xkcd:salmon', alpha =0.5)
    ax1.set_xlabel('Correct fraction, top 10', fontsize=text_size) #plt.ylabel('Volume of peak, normalised')
    #axs[1].set_xlabel('Fraction of correct signal', fontsize=text_size)
    #ax1.set_ylabel('Number of compounds', fontsize=text_size)
    ax1.grid(True, color='xkcd:grey', linestyle='--', linewidth=0.5)
    ax1.legend(['fragments, top 10', 'signal, top 10'])
    #axs[1].grid(True, color='xkcd:grey', linestyle='--', linewidth=0.5)

    fig_name = 'method_accuracy_scenario_sum_formula_top10.png'
    #fig_path = cl_path.dict_path['fragments']
    fig_file = fig_path/fig_name
    plt.savefig(fig_file, bbox_inches='tight', transparent = True, dpi = 600)

    plt.show()


    


    #plot histograms of likelihood values for true and false fragments/maximal elements
    fig, axs = plt.subplots(figsize=(6,3.5), nrows=1, ncols=2, sharex=True)
    fig.suptitle('Distribution of likelihood values for fragments and max. elements')
    
    bins = np.linspace(0.0, 100.0, num = 40)
    
    ax1 = axs[0]
    ax1.hist(list_likelihood_true_frag, bins, color = 'xkcd:blue', alpha =0.5)
    ax1.hist(list_likelihood_false_frag, bins, color = 'xkcd:brick red', alpha =0.5)
    ax1.set_xlabel('Likelihood value of fragments', fontsize=text_size)
    ax1.set_ylabel('Number of fragments or max. elements', fontsize=text_size)
    ax1.legend(['correct fragments', 'incorrect fragments'])

    ax1 = axs[1]
    ax1.hist(list_likelihood_true_att, bins, color = 'xkcd:blue', alpha =0.5)
    ax1.hist(list_likelihood_false_att, bins, color = 'xkcd:brick red', alpha =0.5)
    ax1.set_xlabel('Likelihood value of max. elements', fontsize=text_size)
    #ax1.set_ylabel('Number of results', fontsize=text_size)
    ax1.legend(['correct max. el.', 'incorrect max. el.'])


    fig_name = 'hist_likelihood_frag_att.png'
    #fig_path = cl_path.dict_path['fragments']
    fig_file = fig_path/fig_name
    plt.savefig(fig_file, bbox_inches='tight', transparent = True, dpi = 600)

    plt.show()    


    #plot histograms of ranking values for true and false fragments/maximal elements
    fig, axs = plt.subplots(figsize=(6,3.5), nrows=1, ncols=2, sharex=True)
    fig.suptitle('Distribution of ranking values for fragments and max. elements')
    if len(list_ranking_true_frag) > 0 and len(list_ranking_false_frag) > 0:
        max_ranking = max(max(list_ranking_true_frag), max(list_ranking_false_frag))
    elif len(list_ranking_true_frag) > 0 and len(list_ranking_false_frag) == 0:
        max_ranking = max(list_ranking_true_frag)
    elif len(list_ranking_true_frag) == 0 and len(list_ranking_false_frag) > 0:
        max_ranking = max(list_ranking_false_frag)
    
    bins = np.linspace(0.0, max_ranking, num = max_ranking+1) + 0.5
    #print("bins   " + str((bins)))
    
    ax1 = axs[0]
    ax1.hist(list_ranking_true_frag, bins, color = 'xkcd:blue', alpha =0.5)
    ax1.hist(list_ranking_false_frag, bins, color = 'xkcd:brick red', alpha =0.5)
    ax1.set_xlabel('Ranking value of fragments', fontsize=text_size)
    ax1.set_ylabel('Number of fragments or max. elements', fontsize=text_size)
    ax1.legend(['correct fragments', 'incorrect fragments'])

    ax1 = axs[1]
    ax1.hist(list_ranking_true_att, bins, color = 'xkcd:blue', alpha =0.5)
    ax1.hist(list_ranking_false_att, bins, color = 'xkcd:brick red', alpha =0.5)
    ax1.set_xlabel('Ranking value of max. elements', fontsize=text_size)
    ax1.legend(['correct max. el.', 'incorrect max. el.'])


    fig_name = 'hist_ranking_frag_att.png'
    #fig_path = cl_path.dict_path['fragments']
    fig_file = fig_path/fig_name
    plt.savefig(fig_file, bbox_inches='tight', transparent = True, dpi = 600)

    plt.show()    


    """
    
    #colormap scatter plot of frac of correct formulas vs frac of correct % signal with number of measured masses as colormap
    fig, ax1 = plt.subplots(figsize=(5,5))
    ax1.tick_params(axis='x', labelsize = text_size)
    ax1.tick_params(axis='y', labelsize = text_size)
    fig.suptitle('Accuracy of results') #, fontsize=12   
    plt.xlabel('Fraction of true fragments', fontsize=text_size) #plt.ylabel('Volume of peak, normalised')
    plt.ylabel('Fraction of true signal', fontsize=text_size)
    plt.scatter( [method_performance_results[i][i_met_res_frac_true_frag] for i in range(len(method_performance_results))], 
                [method_performance_results[i][i_met_res_frac_true_signal_to_total] for i in range(len(method_performance_results))], 
                c = [method_performance_results[i][i_met_res_no_meas_mass] for i in range(len(method_performance_results))],
                cmap = 'rainbow', alpha =0.5)    
    cbar = plt.colorbar()
    cbar.set_label('Number of measured masses')

    fig_name = 'method_accuracy_vs_no_masses.png'
    #fig_path = cl_path.dict_path['fragments']
    fig_file = fig_path/fig_name
    plt.savefig(fig_file, bbox_inches='tight', transparent = True, dpi = 600)
    plt.show()


    #*********************************************
    #colormap scatter plot of frac of correct formulas vs frac of correct % signal with total measured I as colormap
    fig, ax1 = plt.subplots(figsize=(5,5))
    ax1.tick_params(axis='x', labelsize = text_size)
    ax1.tick_params(axis='y', labelsize = text_size)
    fig.suptitle('Accuracy of results') #, fontsize=12   
    plt.xlabel('Fraction of true fragments', fontsize=text_size) #plt.ylabel('Volume of peak, normalised')
    plt.ylabel('Fraction of true signal', fontsize=text_size)
    plt.scatter( [method_performance_results[i][i_met_res_frac_true_frag] for i in range(len(method_performance_results))], 
                [method_performance_results[i][i_met_res_frac_true_signal_to_total] for i in range(len(method_performance_results))], 
                c = [method_performance_results[i][i_met_res_sum_I] for i in range(len(method_performance_results))],
                cmap = 'rainbow', alpha =0.5)    
    cbar = plt.colorbar()
    cbar.set_label('Total measured signal per substance, Volt')

    fig_name = 'method_accuracy_vs_sum_I.png'
    #fig_path = cl_path.dict_path['fragments']
    fig_file = fig_path/fig_name
    plt.savefig(fig_file, bbox_inches='tight', transparent = True, dpi = 600)
    plt.show()

    #*********************************************
    #colormap scatter plot of frac of correct formulas vs frac of correct % signal with median mass ppm as colormap
    fig, ax1 = plt.subplots(figsize=(5,5))
    ax1.tick_params(axis='x', labelsize = text_size)
    ax1.tick_params(axis='y', labelsize = text_size)
    fig.suptitle('Accuracy of results') #, fontsize=12   
    plt.xlabel('Fraction of true fragments', fontsize=text_size) #plt.ylabel('Volume of peak, normalised')
    plt.ylabel('Fraction of true signal', fontsize=text_size)
    plt.scatter( [method_performance_results[i][i_met_res_frac_true_frag] for i in range(len(method_performance_results))], 
                [method_performance_results[i][i_met_res_frac_true_signal_to_total] for i in range(len(method_performance_results))], 
                c = [method_performance_results[i][i_met_res_meas_mass_u_median_ppm] for i in range(len(method_performance_results))],
                cmap = 'rainbow', alpha =0.5)    
    cbar = plt.colorbar()
    cbar.set_label('Median mass uncertainty, ppm')

    fig_name = 'method_accuracy_vs_mass_median_ppm.png'
    #fig_path = cl_path.dict_path['fragments']
    fig_file = fig_path/fig_name
    plt.savefig(fig_file, bbox_inches='tight', transparent = True, dpi = 600)
    plt.show()

    """

    """
    fig, axs = plt.subplots(1, 2, figsize=(5,4), sharey=True, tight_layout=True)
    nx, xbins, ptchs = axs[0].hist([method_performance_results[i][i_met_res_frac_true_frag] for i in range(len(method_performance_results))], bins=n_bins)
    plt.clf() # Get rid of this histogram since not the one we want.
    
    nx_frac = nx/float(len(nx)) # Each bin divided by total number of objects.
    width = xbins[1] - xbins[0] # Width of each bin.
    x = np.ravel(zip(xbins[:-1], xbins[:-1]+width))
    y = np.ravel(zip(nx_frac,nx_frac))
    
    axs[0].plot(x, y, linestyle="dashed", label="MyLabel")
    
    plt.show()
    """
    
    
    return None
