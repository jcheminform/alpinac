# -*- coding: utf-8 -*-
"""
Initial Software, Myriam Guillevic and Aurore Guillevic,
Copyright Empa and Inria, 2019 - 2021.

Created on Wed Apr 15 11:27:12 2020
"""

import numpy as np
import math
import time

import networkx as nx
import operator


from utils_data_extraction import fit_iso_profile
from utils_identification import enumerate_all_valid_subfragments
from utils_identification import enumerate_all_valid_subfragments_from_lists
import periodic_table as chem_data
from runtype import RunType
from compound import Compound, CompoundTofwerkTOF, CompoundNoMassCalibration, CompoundNIST


# https://pubchem.ncbi.nlm.nih.gov/
from isotopologue import Isotopologues

fit_mode_fixed_x_axis = 0
fit_mode_variable_x_axis = 1

class NodeLite():
    """A class for the minimal data to encode a node. """
    total_number_nodes_lite = 0 # a class parameter
    total_number_duplicates = 0 # class parameter
    def __init__(self, fragment: list, unique_id: int=None, visited: bool=False, idx_measured_mass: int=None) -> None:
        """Constructor of NodeLite, initialize a node with fragment

        INPUT:
        - `fragment`: a list of non-negative integers encoding a sum-formula
        - `unique_id`: the number-id of the node, should be unique for each node
        - `visited`: flag for efficient graph traversal
        - `idx_measured_mass`: the index of the mass for which fragment is a candidate sum-formula
        
        the children will be stored as adjacent nodes of the graph
        """
        self.fragment = fragment
        #self.children = [] # a list of other NodeFragment that are children (when it will be safe to do it with adjacent nodes of the graph, remove it: use G.successors(n) or G.neighbors(n))
        self.subgraph_size = 1 # the number of nodes of distinct fragment in the dag (can be a tree) rooted at self (#self=1)
        if unique_id is not None:
            self.unique_id = unique_id
        else:
            self.unique_id = self.get_total_number_nodes()
        self.increment_node_counter() # do it with a class method otherwise it creates an instance field
        self.visited = visited
        self.added_somewhere_as_subfragment = False # needed when building the graph
        if idx_measured_mass is not None:
            self.idx_measured_mass = [idx_measured_mass] # index of mass
        else:
            self.idx_measured_mass = None
        # Predicted number of formula solutions using all subsets of abundant fragment formula:
        self.no_all_subfragments = None # will be filled later in compute_all_subfragments_of_nodes
        self.no_all_subfragments_check = None # for test
        self.iso = None # for the class Isotopologue

        #sum signal of this fragment and all its subfragment nodes in the graph (successors)
        # this is the sum of self.iso.iso_list_sum_signal and all node_lite.iso.iso_list_sum_signal for all node_lite that are subfragments of self
        self.subgraph_sum_I = float(0.0)

        #percentage of signal explained by this fragment and all its subfragment nodes in the graph (successors)
        self.subgraph_percent_I = float(0.0)

        #likelohood value for this fragment (captures data from all the subfragment nodes in the graph)
        self.subgraph_likelihood = float(0.0)
        
        #is the node validated, i.e. has the k factor of the isotopologue profile 
        #at least been once optimised using lmfit?
        self.optimised = False


    def __repr__(self):
        return chem_data.get_string_formula(self.fragment) + " id {}".format(self.unique_id) + " visited: {}".format(self.visited)

    def __str__(self):
        return chem_data.get_string_formula(self.fragment)

    @classmethod
    def increment_node_counter(self):
        self.total_number_nodes_lite += 1
    @classmethod
    def reset_node_counter(self):
        self.total_number_nodes_lite = 0
        self.total_number_duplicates = 0
    @classmethod
    def get_total_number_nodes(self):
        return self.total_number_nodes_lite
    @classmethod
    def increment_number_duplicates(self):
        self.total_number_duplicates += 1
    @classmethod
    def decrement_number_duplicates(self):
        self.total_number_duplicates -= 1
    @classmethod
    def get_total_number_duplicates(self):
        return self.total_number_duplicates

    def is_subfragment(self, other: list) -> bool:
        i=0
        while i < chem_data.len_frag and self.fragment[i] <= other[i]:
            i += 1
        return i >= chem_data.len_frag
    def is_same_fragment(self, other: list) -> bool:
        i=0
        while i < chem_data.len_frag and self.fragment[i] == other[i]:
            i += 1
        return i >= chem_data.len_frag
    def is_supfragment(self, other: list) -> bool:
        i=0
        while i < chem_data.len_frag and self.fragment[i] >= other[i]:
            i += 1
        return i >= chem_data.len_frag
    def is_strict_supfragment(self, other: list) -> bool:
        return sum(self.fragment) > sum(other) and self.is_supfragment(other)
    def is_strict_subfragment(self, other: list) -> bool:
        return sum(self.fragment) < sum(other) and self.is_subfragment(other)

    def lost_fragment(self, fragment):
        if self.is_supfragment(fragment):
            lost = [self.fragment[i] - fragment[i] for i in range(chem_data.len_frag)]
            return lost
        if self.is_subfragment(fragment):
            lost = [fragment[i] - self.fragment[i] for i in range(chem_data.len_frag)]
            return lost
        return [0]*self.len_frag

# these are not methods of class NodeLite anymore, but functions
def add_edges_to_node(G: nx.DiGraph, start_node: NodeLite, node: NodeLite) -> (bool, int):
    """Find all the minimal ancestors of `node` and add edges
    
    search for edges to node, starting search at start_node (search for
    ancestors of node). Do not add edges from node (assume there is no
    smaller element, i.e. no successors). Assume that
    G.add_node(unique_id, node) was done before.

    INPUT:
    - `G`: an initialized graph G, containing start_node and node
    - `start_node`: a node in `G` from where starting to look for
      ancestors of `node`
    - `node`: a new node, in G, without ancestors

    OUTPUT:
    Returns True if the node was added somewhere (or a twin was found),
    and the number of occurences, that is the number of new edges to it.
    If a twin was found (a node with the same fragment but another mass
    index), the new mass index is added to the existing note. It returns
    True, 0.
    
    FUTURE potential development:
    One edge could be added only if the neutral loss 
    (the fragment difference between parent and child) has a positive
    DBE value. This way, e.g. H3 could not be a valid neutral loss.
    """
    if start_node.visited:
        return start_node.added_somewhere_as_subfragment, 0 # 0 because if it is already present, do not count it twice
    if start_node.is_supfragment(node.fragment): # yes, it is an ancestor
        if start_node.is_same_fragment(node.fragment): # it should happen rarely...
            # how to distinguish between a same fragment being a possibility for two distinct masses,
            # and in a DAG variant, another occurence of the node? -> with the flag 'visited'
            # if start_node.visited: there was a return statement above, so here this is a twin node, not a duplicate
            start_node.visited = True
            start_node.added_somewhere_as_subfragment = True
            if start_node.idx_measured_mass is not None and node.idx_measured_mass is not None:
                start_node.idx_measured_mass.append(node.idx_measured_mass[0])
                #print("+++There is a duplicate: " + chem_data.get_string_formula(start_node.fragment))
            start_node.increment_number_duplicates()
            return True, 0
        # ok, now check each child.
        # if none of the children is visited and none can have it as subfragment, add it as child of start_node
        # it works even if start_node has no successor: it will add node as sucessor
        occ = 0
        added_somewhere = False
        for child_id in G.successors(start_node.unique_id):
            child = G.nodes[child_id]['node_lite']
            added_as_subfragment, occurences = add_edges_to_node(G, child, node) # recursive call
            # if the child is already visited, it returns child.added_somewhere_as_subfragment, 0
            occ += occurences
            added_somewhere = added_somewhere or added_as_subfragment
            if added_as_subfragment:
                start_node.added_somewhere_as_subfragment = True
        if not added_somewhere:
            # add the fragment as child here: add edge to it
            G.add_edge(start_node.unique_id, node.unique_id)
            start_node.subgraph_size += 1 # assume node as no successors, otherwise add node.subgraph_size instead of 1
            occ += 1
        start_node.visited = True
        start_node.added_somewhere_as_subfragment = True
        return True, occ
    start_node.visited = True
    start_node.added_somewhere_as_subfragment = False
    return False, 0

def reset_unvisited(G: nx.DiGraph, start_node: NodeLite=None) -> None:
    """Reset recursively flag NodeLite.visited to False """
    if start_node is None:
        for max_node_id in G.graph['maximal_nodes']:
            max_node = G.nodes[max_node_id]['node_lite']
            reset_unvisited(G, max_node)
    else:
        if start_node.visited:
            # 1. reset sub-graph
            for child_id in G.successors(start_node.unique_id):
                child = G.nodes[child_id]['node_lite']
                reset_unvisited(G, child)
            # 2. reset itself
            start_node.visited = False
        # if start_node is unvisited, do nothing.

def add_maximal_node_and_edges_to_graph_from_max_fragment(G: nx.DiGraph, list_max_fragments:list, idx_measured_mass: int=None, molecular_ion: bool=False, verbose: bool=False) -> None:
    """Add one maximal node and its edges from it to other former maximal nodes

    Assume that the fragment is a sub-fragment of a node already in the graph, or it is incomparable (it is a singleton)
    It can be used to add a possible molecular ion, in that case `molecular_ion` is set to True

    The new fragment is added to G.graph['maximal_nodes']
    and the former nodes in G.graph['maximal_nodes'] that are
    now its children are removed from G.graph['maximal_nodes'].

    INPUT:
    - `G`: a directed acyclic graph
    - `list_max_fragments`: a list of incomparable fragments that are all incomparable or >= each max fragment of G
    - `idx_measured_mass`: the index of the mass for which the fragments are a candidate sum-formula
    - `molecular_ion`: the fragment is a possible molecular ion, and does not correspond to a measured mass
    OUTPUT: None
    """
    if verbose:
        print("add_maximal_node_and_edges_to_graph_from_max_fragment(G, list_max_fragments={}, idx_measured_mass={}, molecular_ion={}, verbose={})".format(len(list_max_fragments), idx_measured_mass, molecular_ion, verbose))
    if idx_measured_mass is not None:
        if idx_measured_mass in G.graph['candidates_per_measured_mass']:
            # just in case there is already another max element of same mass w.r.t. uncertainty
            G.graph['candidates_per_measured_mass'][idx_measured_mass] += 1
        else: # new fragment with  idx_measured_mass
            G.graph['candidates_per_measured_mass'][idx_measured_mass] = 1
    flag_keep_max_nodes = [True]*len(G.graph['maximal_nodes']) # a flag to say if the node stay maximal after processing all the list of new fragments
    new_max_nodes = []
    for i in range(len(list_max_fragments)):
        max_fragment = list_max_fragments[i]
        if verbose:
            print("adding #{} {}".format(i, chem_data.get_string_formula(max_fragment)))
    
        if len(G.graph['maximal_nodes']) == 0: # G has no root, it should mean it is empty
            node = NodeLite(max_fragment, visited = False, idx_measured_mass=idx_measured_mass)
            # set visited to False do avoid the step reset_unvisited
            G.add_node(node.unique_id, node_lite = node)
            new_max_nodes.append(node.unique_id)
        else:
            children_of_new_max_node = []
            has_twin = False
            # now check the former maximal nodes if they can be possible children
            for j in range(len(G.graph['maximal_nodes'])):# counter of max node number
                former_max_node_id = G.graph['maximal_nodes'][j]
                former_max_node = G.nodes[former_max_node_id]['node_lite']
                # max_fragment should be larger than former_max_node.fragment or incomparable
                if former_max_node.is_subfragment(max_fragment): # sub or equal
                    if former_max_node.is_same_fragment(max_fragment):
                        # duplicate case
                        # do not add a new node, this is a twin (it should not happen though)
                        if former_max_node.idx_measured_mass is not None and idx_measured_mass is not None:
                            former_max_node.idx_measured_mass.append(idx_measured_mass)
                        former_max_node.increment_number_duplicates()
                        has_twin = True
                    else: # former_max_node.is_strict_subfragment(max_fragment):
                        children_of_new_max_node.append(former_max_node)
                        flag_keep_max_nodes[j] = False
            # now update the list of maximal nodes
            if not has_twin:
                new_max_node = NodeLite(max_fragment, visited = False, idx_measured_mass=idx_measured_mass)
                G.add_node(new_max_node.unique_id, node_lite = new_max_node)
                G.graph['nodes_validated'].append(new_max_node.unique_id)
                # now add edges from the the node to its children given in children_of_new_max_node
                for child in children_of_new_max_node:
                    G.add_edge(new_max_node.unique_id, child.unique_id)
                    new_max_node.subgraph_size += child.subgraph_size
                new_max_nodes.append(new_max_node.unique_id)
    
    G.graph['maximal_nodes'] = [G.graph['maximal_nodes'][j] for j in range(len(flag_keep_max_nodes)) if flag_keep_max_nodes[j]] + new_max_nodes



def add_nodes_and_edges_to_graph_from_list_fragments(G: nx.DiGraph, list_fragments:list, idx_measured_mass: int=None, add_mol_ion: bool=False, verbose: bool=False) -> None:
    """Add nodes for each of the fragments, and the edges to them
    
    Assume that the fragments are incomparable or given in decreasing
    order of mass (very important)
    if there is no edge to a fragment, its unique_id is added to
    G.graph['maximal_nodes'] (it is actually a singleton)
    if the node is a twin (there is already a same fragment encoding for
    a different mass of overlapping uncertainty range), it is not added
    to the Graph, but its mass index is added to the twin node.

    INPUT:
    - `G`: a directed acyclic graph
    - `list_fragments`: a list of fragments incomparable or ordered >=
    - `idx_measured_mass`: the index of the mass for which all the fragments are candidate sum-formulas
    OUTPUT: None
    """
    if verbose:
        print("add_nodes_and_edges_to_graph_from_fagments, len(list_fragments) = {}".format(len(list_fragments)))
    if idx_measured_mass is not None:
        G.graph['candidates_per_measured_mass'][idx_measured_mass] = len(list_fragments)
        # save the number of candidates per mass for remove_singletons
    for i in range(len(list_fragments)):
        fragment = list_fragments[i]
        if verbose:
            print("adding #{} {}".format(i, chem_data.get_string_formula(fragment)))
        if len(G.graph['maximal_nodes']) == 0: # G has no root, it should mean it is empty
            node = NodeLite(fragment, visited = False, idx_measured_mass=idx_measured_mass)
            # set visited to False do avoid the step reset_unvisited
            G.add_node(node.unique_id, node_lite = node)
            G.graph['maximal_nodes'].append(node.unique_id)
        else:
            node = NodeLite(fragment, visited = True, idx_measured_mass=idx_measured_mass)
            G.add_node(node.unique_id, node_lite = node)
            added_somewhere = False
            total_occ = 0
            
            if add_mol_ion == False:
                for max_node_id in G.graph['maximal_nodes']:
                    max_node = G.nodes[max_node_id]['node_lite']
                    
                    #node must be smaller than max_node
                    added, occ = add_edges_to_node(G, max_node, node)
                    
                    added_somewhere = added_somewhere or added
                    total_occ += occ
                if added_somewhere and total_occ == 0: # this is a twin, remove the node
                    print("node is a twin, removed: {}".format(node))
                    G.remove_node(node.unique_id)
                elif not added_somewhere:
                    # add it to the list of max nodes
                    G.graph['maximal_nodes'].append(node.unique_id)
                    


            else:
                #we deal with adding a candidate molecular ion
                
                
                print("listed max nodes")
                print([nod for nod in G.graph['maximal_nodes']])
                print([G.nodes[nod]['node_lite'].iso.frag_abund_formula for nod in G.graph['maximal_nodes']])
                #check if previously existing maximal ions should be removed
                #add edges with previously existing maximal nodes, if applicable
                old_max_nodes_list = G.graph['maximal_nodes']
                for old_max_node_id in old_max_nodes_list:
                    old_max_node = G.nodes[old_max_node_id]['node_lite']
                    #print(G.nodes[old_max_node_id]['node_lite'].iso.frag_abund_formula)
                    if node.is_supfragment(old_max_node.iso.frag_abund): # yes, it is an ancestor
                        print(str(node) + " is supfrag from " + str(old_max_node.iso.frag_abund_formula))
                        # add the fragment as child here: add edge to it
                        #Note: here old_max_node is a smaller fragment than node
                        G.add_edge(node.unique_id, old_max_node.unique_id)
                print("nx.descendants "+ str(nx.descendants(G, node.unique_id)))
                print("nx.ancestors " + str(nx.ancestors(G, node.unique_id)))
                print("Edges: " + str(nx.edges(G, node.unique_id)))
                
                G.graph['maximal_nodes'].append(node.unique_id)
                G.graph['nodes_validated'].append(node.unique_id)
                print("max. nodes: " + str([nod for nod in G.graph['maximal_nodes']]))
                                    
            reset_unvisited(G)
            
                

def update_graph_data_after_removing_nodes(G: nx.DiGraph, removed_nodes: list) -> None:
    """Remove the node numbers from nodes_validated and added_nodes """
        # update data on nodes
    if 'nodes_validated' in G.graph:
        for s in removed_nodes:
            if s in G.graph['nodes_validated']:
                G.graph['nodes_validated'].remove(s)
    if 'added_nodes' in G.graph:
        for s in removed_nodes:
            if s in G.graph['added_nodes']:
                G.graph['added_nodes'].remove(s)

def remove_all_singletons(G: nx.DiGraph):
    """remove all maximal nodes that have no child
    """
    max_nodes = [G.nodes[n]['node_lite'] for n in G.graph['maximal_nodes']]
    not_singletons = [node for node in max_nodes if len(G.succ[node.unique_id]) > 0]
    singletons = [node for node in max_nodes if len(G.succ[node.unique_id]) == 0]
    G.graph['maximal_nodes'] = [node_i.unique_id for node_i in not_singletons]
    if len(G.graph['candidates_per_measured_mass']) > 0:
        for s in singletons:
            if s.idx_measured_mass is not None:
                for i_m in s.idx_measured_mass:
                    G.graph['candidates_per_measured_mass'][i_m] -= 1
    singleton_numbers = [s.unique_id for s in singletons]
    G.remove_nodes_from(singleton_numbers)
    if 'no_nodes_singletons' in G.graph:
        G.graph['no_nodes_singletons'] += len(singleton_numbers)
    update_graph_data_after_removing_nodes(G, singleton_numbers)
    return singletons

def remove_singletons(G: nx.DiGraph):
    """remove maximal nodes that have no child

    Morever, check that all the maximal nodes of one mass index are not
    all the nodes of that mass index. If a mass index has only singleton
    nodes, do not remove them (otherwise there is no candidate anymore
    for that mass index).

    """
    max_nodes = [G.nodes[n]['node_lite'] for n in G.graph['maximal_nodes']]
    not_singletons = [node for node in max_nodes if len(G.succ[node.unique_id]) > 0]
    singletons = [node for node in max_nodes if len(G.succ[node.unique_id]) == 0]
    singleton_numbers = [s.unique_id for s in singletons]
    # if idx_measured_mass is not provided, remove all singletons
    if len(G.graph['candidates_per_measured_mass']) == 0:
        G.remove_nodes_from(singleton_numbers)
        G.graph['maximal_nodes'] = [node_i.unique_id for node_i in not_singletons]
        update_graph_data_after_removing_nodes(G, singleton_numbers)
        return singletons
    # now check that a set of singletons is not all for a mass
    # that is, do not remove singletons of a mass which has only singletons
    # make batches of singletons according to their idx_measured_mass
    singletons_dict = {}
    removed_singletons = []
    G.graph['maximal_nodes'] = [node_i.unique_id for node_i in not_singletons] # and add the singletons that are not removed
    for s in singletons:
        if s.idx_measured_mass is None:
            G.remove_node(s.unique_id)
            removed_singletons.append(s)
        else:
            for i_m in s.idx_measured_mass:
                if i_m in singletons_dict:
                    singletons_dict[i_m].append(s)
                else:
                    singletons_dict[i_m] = [s]
    # check that removing the singletons will not remove all the candidates for a mass
    for i_m in singletons_dict:
        if len(singletons_dict[i_m]) < G.graph['candidates_per_measured_mass'][i_m]:
            #G.remove_nodes_from([s.unique_id for s in singletons_dict[i_m]])
            #G.graph['candidates_per_measured_mass'][i_m] -= len(singletons_dict[i_m])
            #removed_singletons += singletons_dict[i_m]
            # in case there are 'twin' nodes:
            for s in singletons_dict[i_m]:
                if len(s.idx_measured_mass) == 1: # no twin node
                    G.remove_node(s.unique_id)
                    removed_singletons.append(s)
                else:
                    s.idx_measured_mass.remove(i_m)
                    s.decrement_number_duplicates()
                    G.graph['maximal_nodes'].append(s.unique_id) # actually the node itself is still there
                G.graph['candidates_per_measured_mass'][i_m] -= 1
        else: # re-add the nodes of singletons_dict[i_m] to the list of maximal elements
            for s in singletons_dict[i_m]:
                G.graph['maximal_nodes'].append(s.unique_id)
    singleton_numbers = [s.unique_id for s in removed_singletons]
    if 'no_nodes_singletons' in G.graph:
        G.graph['no_nodes_singletons'] += len(singleton_numbers)
    update_graph_data_after_removing_nodes(G, singleton_numbers)
    return removed_singletons

def update_parents_after_removed_node(G: nx.DiGraph, parent_numbers: list, removed_iso_list_sum_signal: float, cl_comp: Compound) -> None:
    """Ipdate data of parents recursively
    """
    for parent_number in parent_numbers:
        parent = G.nodes[parent_number]['node_lite']
        if parent.visited:
            continue
        parent.visited = True
        # update size of subgraph
        parent.subgraph_size -= 1
        parent.subgraph_sum_I -= removed_iso_list_sum_signal
        parent.subgraph_percent_I = float(100.0)*parent.subgraph_sum_I/cl_comp.sum_I
        parent.subgraph_likelihood = (float(100.0) - abs(float(100.0)-parent.subgraph_percent_I))* parent.subgraph_size / parent.no_all_subfragments
        parents = G.predecessors(parent_number)
        update_parents_after_removed_node(G, parents, removed_iso_list_sum_signal, cl_comp)

def reset_parents_unvisited(G: nx.DiGraph, parent_numbers: list) -> None:
    """set NodeLite.visited to False recursively upwards"""
    for parent_number in parent_numbers:
        parent = G.nodes[parent_number]['node_lite']
        if parent.visited:
            parent.visited = False
            parents = G.predecessors(parent_number)
            reset_parents_unvisited(G, parents)

def remove_node_update_edges(G: nx.DiGraph, node: NodeLite, cl_comp: Compound, update_data: bool=False) -> None:
    """Remove one node and update graph data

    INPUT:
    - `G`: directed acyclic graph
    - `node`: the node to be removed
    - `cl_comp`: Compound data (to update the likelihood of the parent nodes)

    OUTPUT: None (the graph itself is updated in place)

    This function updates the edges before removing a node in the graph.
    A node encodes a fragment, it has parent nodes encoding sup-fragments, and child nodes encoding sub-fragments.
    When removing a node, its childen (the sub-fragments) should be reconnected to the parents only in certain cases.

    1. remove all edges (p, n) and (n, c) where p are parents and c are children
    2. for each child, its parents cannot be grand-parents of n (the removed node) by construction of G
       set edges between parents of n and children of n in appropriate cases
    """
    # 1. lists the parents of node
    parents = list(G.predecessors(node.unique_id))
    if len(parents) == 0:
        G.graph['maximal_nodes'].remove(node.unique_id)
    children = list(G.neighbors(node.unique_id))
    for parent_number in parents:
        G.remove_edge(parent_number, node.unique_id)
    if update_data:
        update_parents_after_removed_node(G, parents, node.iso.iso_list_sum_signal, cl_comp)
        reset_parents_unvisited(G, parents)
    #else:
    #    for parent_number in parents:
    #        G.nodes[parent_number]['node_lite'].subgraph_size -= 1

    for child_number in children:
        G.remove_edge(node.unique_id, child_number)
        # compare the child's parents to the node's parents (some of the child's grand-parents)
        # if one of the **other** parents of the child is a sub-fragment of one of the node's parents,
        # nothing to so.
        # otherwise, add the child as a child of the considered node's parent
        childs_parents_number = list(G.predecessors(child_number))
        # compare childs_parent to all nodes'parents in parents
        for parent_number in parents:
            is_subfragment = False
            for childs_parent_number in childs_parents_number:
                is_subfragment = is_subfragment or G.nodes[childs_parent_number]['node_lite'].is_subfragment(G.nodes[parent_number]['node_lite'].fragment)
            if not is_subfragment:
                G.add_edge(parent_number, child_number)

        if len(list(G.predecessors(child_number))) == 0:
            G.graph['maximal_nodes'].append(child_number)

    G.remove_node(node.unique_id)

def print_maximal_nodes(G: nx.DiGraph) -> None:
    print("maximal fragments are {}".format(" ".join([G.nodes[node_i]['node_lite'].iso.frag_abund_formula for node_i in G.graph['maximal_nodes']])))

def check_maximal_nodes(G: nx.DiGraph) -> bool:
    """Check that G.graph['maximal_nodes'] is made of exactly all the maximal nodes"""
    effective_maximal_nodes = [node_number for node_number in G.nodes if len(list(G.predecessors(node_number))) == 0]
    missing_max_node = [node_i for node_i in effective_maximal_nodes if node_i not in G.graph['maximal_nodes']]
    false_max_node = [node_i for node_i in G.graph['maximal_nodes'] if node_i not in effective_maximal_nodes]
    if len(false_max_node) > 0:
        print("Error maximal nodes: {} ({}) in G.graph['maximal_nodes'] are not in effective_maximal_nodes".format(false_max_node, " ".join([G.nodes[node_i]['node_lite'].iso.frag_abund_formula for node_i in false_max_node])))
    if len(missing_max_node) > 0:
        print("Error maximal nodes: {} ({}) in effective_maximal_nodes are not in G.graph['maximal_nodes']".format(missing_max_node, " ".join([G.nodes[node_i]['node_lite'].iso.frag_abund_formula for node_i in missing_max_node])))
    #return len(missing_max_node) == 0 and len(false_max_node) == 0
    if not (len(missing_max_node) == 0 and len(false_max_node) == 0):
        raise ValueError("Error maximal nodes: missing_max_node = {} ({}) false_max_node = {} ({})".format(missing_max_node, " ".join([G.nodes[node_i]['node_lite'].iso.frag_abund_formula for node_i in missing_max_node]), false_max_node, " ".join([G.nodes[node_i]['node_lite'].iso.frag_abund_formula for node_i in false_max_node])))

def compute_all_subfragments_of_nodes(G: nx.DiGraph, min_measured_mass:float, start_node: NodeLite=None, check: bool=False) -> None:
    """Computes the number of theoretical sub-fragments of each node
    iteratively with the knapsack algorithm

    INPUT:
    - `G` a directed acyclic graph whose maximal nodes (numbers) are
    given in G.graph['maximal_nodes']
    - `min_measured_mass`: the minimal mass measured by the machine
    (such as 23.0 m/z)
    - `start_node`: start at this node and recursively inspect sub-nodes
    - `check`: if True, set field mynodelite.no_all_subfragments_check,
    otherwise, set mynodelite.no_all_subfragments (should be the default)

    OUTPUT: None, the information is stored in
    G.nodes[node_number]['node_lite'].no_all_subfragments<_check>
    """
    if start_node is None:
        for max_node_id in G.graph['maximal_nodes']:
            max_node = G.nodes[max_node_id]['node_lite']
            compute_all_subfragments_of_nodes(G, min_measured_mass, max_node, check=check)
        reset_unvisited(G)
        return
    if start_node.visited:
        return
    no_all_subfragments, list_subfragments = enumerate_all_valid_subfragments(start_node.fragment, min_mass=min_measured_mass)
    #print(str([chem_data.get_string_formula(list_subfragments[i]) for i in range(len(list_subfragments))]))
    if check:
        start_node.no_all_subfragments_check = no_all_subfragments
    else:
        start_node.no_all_subfragments = no_all_subfragments
    # recursive call
    start_node.visited = True
    for child_id in G.successors(start_node.unique_id):
        child = G.nodes[child_id]['node_lite']
        compute_all_subfragments_of_nodes(G, min_measured_mass, child, check=check)

def compute_all_subfragments_of_nodes_with_lists(G: nx.DiGraph, min_measured_mass:float, start_node: NodeLite=None, parent_node: NodeLite=None, multi_valent: list=None, mono_valent: dict=None, verbose: bool=False) -> None:
    """Computes the number of theoretical sub-fragments of each node
    recursively with the knapsack algorithm

    INPUT:
    - `G` a directed acyclic graph whose maximal nodes (numbers) are
    given in G.graph['maximal_nodes']
    - `min_measured_mass`: the minimal mass measured by the machine
    (such as 23.0 m/z)
    - `start_node`: start at this node and recursively inspect sub-nodes
    - `parent_node`: the parent of `start_node`
    - `multi_valent`: list of partial fragments of the parent node made
    of multi-valent atoms only, None if start_node is None or maximal
    - `mono_valent`: list of partial fragments of the parent node made
    of mono-valent atoms only, None if start_node is None or maximal

    OUTPUT: None, the information is stored in
    G.nodes[node_number]['node_lite'].no_all_subfragments

    Recursive call, with a Depth-first-search algorithm.
    If start_node is None, start with all the maximal elements.
    Set to visited the nodes that have been processed.
    To improve efficiency, try to share some parts of knapsack
    enumeration. With inputs max_mass = fragment.mass and mass_min = 0.0
    the knapsack outputs all the subfragments. It computed internally a
    list of subfragments made of multi-valent atoms, and a list of sets
    of at most v mono-valent atoms, where v is the highest valence of
    the fragments in the multi-valent list.
    Sort by increasing order of mass both lists, and cut each list at
    the mass of the fragment made of all the multi-valent atoms of
    start_node, then refine the list by discarding all fragments that
    are not subfagments of multi-valent(start_node.fragment).
    Do the same for the list of mono-valents.
    In this way, the re-run of the knapsack is replaced by list
    filtering. It could be faster. Need to be checked.
    """
    if start_node is None:
        for max_node_id in G.graph['maximal_nodes']:
            max_node = G.nodes[max_node_id]['node_lite']
            if verbose:
                print("processing maximal node {}".format(str(max_node)))
            compute_all_subfragments_of_nodes_with_lists(G, min_measured_mass, max_node, verbose=verbose)
        reset_unvisited(G)
        return
    if start_node.visited:
        return
    if verbose:
        print("processing node {}".format(str(start_node)))
    if multi_valent is None or mono_valent is None:
        # we are at a maximal element. Call knapsack.
        no_all_subfragments, list_subfragments, multi_valent, mono_valent = enumerate_all_valid_subfragments(start_node.fragment, min_mass=min_measured_mass, return_partial_lists=True)
    else: # filter the partial lists and compute the number of subfragments of start_node
        no_all_subfragments, list_subfragments, multi_valent, mono_valent = enumerate_all_valid_subfragments_from_lists(start_node.fragment, multi_valent, mono_valent, parent_node.fragment, min_mass=min_measured_mass)
    start_node.no_all_subfragments = no_all_subfragments
    if verbose:
        print("processing node {}, {} subfragments".format(str(start_node), no_all_subfragments))
    # recursive call
    start_node.visited = True
    for child_id in G.successors(start_node.unique_id):
        child = G.nodes[child_id]['node_lite']
        if not child.visited:
            compute_all_subfragments_of_nodes_with_lists(G, min_measured_mass, child, start_node, multi_valent, mono_valent, verbose=verbose)

def node_edge_labels_colors_for_pydot(G: nx.DiGraph):
    node_labels = {}
    node_colors = [0]*G.number_of_nodes()
    idx= 0
    for key_node in G.nodes():
        node_labels[key_node] = chem_data.get_string_formula(G.nodes[key_node]['node_lite'].fragment)
        node_colors[idx] = chem_data.get_mass(G.nodes[key_node]['node_lite'].fragment)
        idx+=1
    edge_labels = {}
    for key_edge in G.edges():
        (u,v) = key_edge
        edge_labels[key_edge] = chem_data.get_string_formula(G.nodes[u]['node_lite'].lost_fragment(G.nodes[v]['node_lite'].fragment))
    return node_labels, node_colors, edge_labels

def set_str_formula_and_attribute_for_pydot(G: nx.DiGraph, start_node: NodeLite=None) -> None:
    if start_node is None:
        for max_node_id in G.graph['maximal_nodes']:
            set_str_formula_and_attribute_for_pydot(G, G.nodes[max_node_id]['node_lite'])
        reset_unvisited(G)
    else:
        if not start_node.visited:
            G.nodes[start_node.unique_id]['str_formula'] = chem_data.get_string_formula(start_node.fragment)
            G.nodes[start_node.unique_id]['attribute'] = chem_data.get_mass(start_node.fragment)
            start_node.visited = True
            for child_id in G.successors(start_node.unique_id):
                child = G.nodes[child_id]['node_lite']
                G.edges[(start_node.unique_id, child.unique_id)]['label'] = chem_data.get_string_formula(start_node.lost_fragment(child.fragment))
                set_str_formula_and_attribute_for_pydot(G, child)

def initialise_isotopologue_set(G: nx.DiGraph, cl_comp, add_mol_ion: bool=False, mol_ion_nodes: list = None) -> None:
    """Initialise the Isotopologue class of each node.

    Do it with a depth-first-search graph traversal so that it will be
    possible in an improved version to avoid non-promising nodes.
    TODO, Work in progress
    """
    
    list_nodes_to_initialise = []
    if add_mol_ion == False:
        list_nodes_to_initialise = G.nodes()
    else:
        list_nodes_to_initialise = mol_ion_nodes

    for node_id in list_nodes_to_initialise:
        node = G.nodes[node_id]['node_lite']
        node.iso = Isotopologues(node.fragment, node.idx_measured_mass, cl_comp)
        node.iso.create_isotopologue_suite()

