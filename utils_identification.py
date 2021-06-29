# -*- coding: utf-8 -*-
"""
Initial Software, Myriam Guillevic and Aurore Guillevic,
Copyright Empa and Inria, 2019 - 2021.

created 23 04 2019

"""
import math
#import time
#from pathlib import Path, PurePath

#import numpy as np
from periodic_table import formula_to_e_list, formula_to_frag
import periodic_table as chem_data

class Binomial():
    tab_binomial = [
        1,
        [1, 2],
        [1, 3],
        [1, 4, 6],
        [1, 5, 10],
        [1, 6, 15, 20],
        [1, 7, 21, 35],
        [1, 8, 28, 56, 70],
        [1, 9, 36, 84, 126],
        [1, 10, 45, 120, 210, 252],
        [1, 11, 55, 165, 330, 462],
        [1, 12, 66, 220, 495, 792, 924],
        [1, 13, 78, 286, 715, 1287, 1716],
        [1, 14, 91, 364, 1001, 2002, 3003, 3432],
        [1, 15, 105, 455, 1365, 3003, 5005, 6435],
        [1, 16, 120, 560, 1820, 4368, 8008, 11440, 12870],
        [1, 17, 136, 680, 2380, 6188, 12376, 19448, 24310],
        [1, 18, 153, 816, 3060, 8568, 18564, 31824, 43758, 48620],
        [1, 19, 171, 969, 3876, 11628, 27132, 50388, 75582, 92378],
        [1, 20, 190, 1140, 4845, 15504, 38760, 77520, 125970, 167960, 184756],
        [1, 21, 210, 1330, 5985, 20349, 54264, 116280, 203490, 293930, 352716],
        [1, 22, 231, 1540, 7315, 26334, 74613, 170544, 319770, 497420, 646646, 705432],
        [1, 23, 253, 1771, 8855, 33649, 100947, 245157, 490314, 817190, 1144066, 1352078],
        [1, 24, 276, 2024, 10626, 42504, 134596, 346104, 735471, 1307504, 1961256, 2496144, 2704156],
        [1, 25, 300, 2300, 12650, 53130, 177100, 480700, 1081575, 2042975, 3268760, 4457400, 5200300],
        [1, 26, 325, 2600, 14950, 65780, 230230, 657800, 1562275, 3124550, 5311735, 7726160, 9657700, 10400600],
        [1, 27, 351, 2925, 17550, 80730, 296010, 888030, 2220075, 4686825, 8436285, 13037895, 17383860, 20058300],
        [1, 28, 378, 3276, 20475, 98280, 376740, 1184040, 3108105, 6906900, 13123110, 21474180, 30421755, 37442160, 40116600],
        [1, 29, 406, 3654, 23751, 118755, 475020, 1560780, 4292145, 10015005, 20030010, 34597290, 51895935, 67863915, 77558760],
        [1, 30, 435, 4060, 27405, 142506, 593775, 2035800, 5852925, 14307150, 30045015, 54627300, 86493225, 119759850, 145422675, 155117520],
        [1, 31, 465, 4495, 31465, 169911, 736281, 2629575, 7888725, 20160075, 44352165, 84672315, 141120525, 206253075, 265182525, 300540195],
        [1, 32, 496, 4960, 35960, 201376, 906192, 3365856, 10518300, 28048800, 64512240, 129024480, 225792840, 347373600, 471435600, 565722720, 601080390],
        [1, 33, 528, 5456, 40920, 237336, 1107568, 4272048, 13884156, 38567100, 92561040, 193536720, 354817320, 573166440, 818809200, 1037158320, 1166803110],
        [1, 34, 561, 5984, 46376, 278256, 1344904, 5379616, 18156204, 52451256, 131128140, 286097760, 548354040, 927983760, 1391975640, 1855967520, 2203961430, 2333606220],
        [1, 35, 595, 6545, 52360, 324632, 1623160, 6724520, 23535820, 70607460, 183579396, 417225900, 834451800, 1476337800, 2319959400, 3247943160, 4059928950, 4537567650],
        [1, 36, 630, 7140, 58905, 376992, 1947792, 8347680, 30260340, 94143280, 254186856, 600805296, 1251677700, 2310789600, 3796297200, 5567902560, 7307872110, 8597496600, 9075135300],
        [1, 37, 666, 7770, 66045, 435897, 2324784, 10295472, 38608020, 124403620, 348330136, 854992152, 1852482996, 3562467300, 6107086800, 9364199760, 12875774670, 15905368710, 17672631900],
        [1, 38, 703, 8436, 73815, 501942, 2760681, 12620256, 48903492, 163011640, 472733756, 1203322288, 2707475148, 5414950296, 9669554100, 15471286560, 22239974430, 28781143380, 33578000610, 35345263800],
        [1, 39, 741, 9139, 82251, 575757, 3262623, 15380937, 61523748, 211915132, 635745396, 1676056044, 3910797436, 8122425444, 15084504396, 25140840660, 37711260990, 51021117810, 62359143990, 68923264410],
        [1, 40, 780, 9880, 91390, 658008, 3838380, 18643560, 76904685, 273438880, 847660528, 2311801440, 5586853480, 12033222880, 23206929840, 40225345056, 62852101650, 88732378800, 113380261800, 131282408400, 137846528820],
    ]

    @classmethod
    def build_tab_binom(self,n = 30):
        tab_binomial = [1]
        for n in range(2,n+1):
            tab_binomial.append([1,n] + [None]*(n//2 - 1))
            for i in range(2, n//2+1):
                if i-1 > (n-1)//2:
                    r = tab_binomial[n-2][(n-1) - (i-1)]
                else:
                    r = tab_binomial[n-2][i-1]
                if i > (n-1)//2:
                    s = tab_binomial[n-2][(n-1)-i]
                else:
                    s = tab_binomial[n-2][i]
                tab_binomial[n-1][i] =  r + s
        return tab_binomial

    @classmethod
    def print_tab_binom(self):
        print("tab_binomial = [")
        for n in range(len(self.tab_binomial)-1):
            print("    {},".format(self.tab_binomial[n]))
        print("    {}\n]".format(self.tab_binomial[-1]))

    @classmethod
    def print_tab_binom(self, tab_binom):
        print("tab_binomial = [")
        for n in range(len(tab_binom)-1):
            print("    {},".format(tab_binom[n]))
        print("    {}\n]".format(tab_binom[-1]))


    @classmethod
    def binom(self, i, n):
        if i > n or i < 0:# this case should not happen
            return 0
        if i > n//2:
            i = n-i
        if i == 0:
            return 1
        if i == 1:
            return n
        if n <= len(self.tab_binomial) and self.tab_binomial[n-1][i] is not None:
            return self.tab_binomial[n-1][i]
        res = self.binom(i-1,n-1) + self.binom(i,n-1)
        for j in range(len(self.tab_binomial), n):
            self.tab_binomial.append([1,j+1] + [None]*((j-1)//2))
        self.tab_binomial[n-1][i] = res
        return res

def p_iso_rare_compute_wrt_pure_abundant(abundant_frag: list, fragment: list)-> float:
    """Assume abundant_frag is made of abundant atoms only
    """
    p_iso=float(1.0)
    # consider only rare atoms
    for i in [_ for _ in chem_data.idx_abun_isotopes if abundant_frag[_] > 0 and len(chem_data.dict_abundant_rare[_]) > 0]:
        n_e = abundant_frag[i]
        for idx_rare in [_ for _ in chem_data.dict_abundant_rare[i] if fragment[_] > 0]:
            n_ei = fragment[idx_rare]
            bi = float(Binomial.binom(n_ei,n_e))
            #if bi == 0:
            #    print("Error Binomial({},{}) -> 0, initial n_e = {}, n_ei = {}".format(n_ei, n_e, abundant_frag[i], [fragment[_] for _ in chem_data.dict_abundant_rare[i]]))
            if n_ei > 1:
                p_iso *= (chem_data.abundance[idx_rare]/chem_data.abundance[i])**n_ei * bi
            else:
                p_iso *= (chem_data.abundance[idx_rare]/chem_data.abundance[i]) * bi
            n_e -= n_ei
    return p_iso

def _get_patterns_aux(sum_: int, vec_a: list, idx: int, res: list) -> None:
    """Recursive enumeration of patterns (multinomial)

    INPUT:
    - ``sum_``: the number of items left to attribute
    - ``vec_a``: the current fragment as a vector
    - ``idx``: the current index in vec_a
    - ``res``: a list of results that grows along the recursive calls

    OUTPUT: None, the result is stored in ``res``
    """
    if idx >= len(vec_a):
        return
    if sum_ == 0 and idx < len(vec_a):
        for i in range(idx, len(vec_a)):
            vec_a[i] = 0
        res.append(vec_a.copy())
        return
    if idx == len(vec_a)-1:
        vec_a[idx] = sum_
        res.append(vec_a.copy())
        return
    for ai in range(sum_, -1, -1):
        vec_a[idx] = ai
        _get_patterns_aux(sum_-ai, vec_a, idx+1, res)


def get_all_isotopic_patterns_one_element(atom_idx: int, n: int) -> list:
    """ given one atom index and an integer n, generate all possible isotopic combinations
    input: atom, an index number

    INPUT:
    - ``atom_idx``: an abundant atom index between 0 and chem_data.len_frag
    - ``
    """
    list_all_isotopologues = []
    indices_iso = [atom_idx] + chem_data.dict_abundant_rare[atom_idx]
    iso = len(indices_iso)
    # enumerating the possible choices is like enumerating the possible ways to sum up to n
    # do a knapsack-like algo
    a = [0] * iso
    res = []
    _get_patterns_aux(n, a, 0, res)
    for c in res:
        vec = [0 for i in range(chem_data.len_frag)]
        for i_ in range(iso):
            vec[indices_iso[i_]] = c[i_]
        list_all_isotopologues.append(vec)
    return list_all_isotopologues

def generate_all_isotopolog_fragments(frag):
    """Given a fragment, generate all possible isotopologues.
    """
    list_all_fragments = []
    for idx in [ _ for _ in chem_data.idx_abun_isotopes if frag[_] > 0]:
        iso_atoms = get_all_isotopic_patterns_one_element(idx, frag[idx])
        if len(list_all_fragments) == 0: # first step
            list_all_fragments = iso_atoms
        else:
            list_all_fragments = [[frg[i_] + vec[i_] for i_ in range(chem_data.len_frag)] for vec in iso_atoms for frg in list_all_fragments]
    list_all_fragments.remove(list_all_fragments[0])
    return list_all_fragments

def generate_all_isotopolog_fragments_with_relative_proba(frag:list, min_relative_proba:float=None, verbose:bool=False):
    """Generate all possible isotopologues and relative probabilities

    INPUT:
    - ``frag``: a fragment (list of coefficients) made of abundant atoms only
    - ``min_relative_proba``: do not list fragments of lower relative probability
    
    OUTPUT: a list of pairs [(fragment, relative_proba)] where the
    relative probability is the ratio proba(frag_rare)/proba(frag_abun)
    """
    # the list of abundant atoms in the fragment
    idx_abundant = [ _ for _ in chem_data.idx_abun_isotopes if frag[_] > 0]
    # the patterns of abundant and rare elements for each atom (independently)
    list_iso_atoms = [get_all_isotopic_patterns_one_element(idx, frag[idx]) for idx in idx_abundant]

    #list_proba = [float(1.0)] + [p_iso_rare_compute(frg_abun, iso_atoms[i_]) for i_ in range(1, len(iso_atoms))]
    list_fragments_with_proba_all_elements = [[(float(1.0), iso_atoms[0])]
                                              + [(p_iso_rare_compute_wrt_pure_abundant(iso_atoms[0], iso_atoms[i_]), iso_atoms[i_]) for i_ in range(1, len(iso_atoms))]
                                              for iso_atoms in list_iso_atoms]
    # now find the maximum abundant fragment for each element
    if min_relative_proba is not None:
        max_rel_abundance = float(1.0)
        list_max_rel_abun = [1.0]*len(list_fragments_with_proba_all_elements)
        j_ = 0
        for li in list_fragments_with_proba_all_elements:
            max_i = max([pi for (pi,_) in li])
            list_max_rel_abun[j_] = max_i
            max_rel_abundance *= max_i # this is >= 1 and yes it can be > 1 for example Br[81]Br wrt Br2
            j_ += 1
    list_all_fragments_with_proba = []
    j_ = 0
    for idx in idx_abundant:
        iso_atoms = list_iso_atoms[j_]
        list_fragments_with_proba = list_fragments_with_proba_all_elements[j_]
        frg_abun = iso_atoms[0]
        if min_relative_proba is not None:
            min_rel_pr = min_relative_proba/max_rel_abundance*list_max_rel_abun[j_]
            # 1. partition list_fragments_with_proba in two subsets according to the threshold min_relative_proba and cut
            l=0
            r=len(list_fragments_with_proba)
            while l < r:
                if list_fragments_with_proba[l][0] >= min_rel_pr:
                    l += 1
                else:# swap
                    list_fragments_with_proba[l],list_fragments_with_proba[r-1] = list_fragments_with_proba[r-1],list_fragments_with_proba[l]
                    r -= 1
            #print("")
            cut = r-1 # all items at indices > r-1 are < min_rel_pr
            if cut+1 < len(list_fragments_with_proba) and list_fragments_with_proba[cut+1][0] >= min_rel_pr:
                cut += 1
            if cut < len(list_fragments_with_proba)-1:
                list_fragments_with_proba = list_fragments_with_proba[:cut+1]
        list_fragments_with_proba.sort(reverse=True)
        if len(list_fragments_with_proba) == 0:
            raise ValueError("Error list_fragments_with_proba is empty, abundant frag is {}, threshold is {}".format(chem_data.get_string_formula(frg_abun), min_rel_pr))
        
        if len(list_all_fragments_with_proba) == 0:
            list_all_fragments_with_proba = list_fragments_with_proba
        else:
            if min_relative_proba is not None:
                list_all_fragments_with_proba_tmp = []
                for (pa,frg) in list_all_fragments_with_proba:
                    j=0
                    (pi,vec) = list_fragments_with_proba[j]
                    pr = pi*pa
                    pr_ok = pr >= min_rel_pr
                    while j+1 < len(list_fragments_with_proba) and pr_ok:
                        list_all_fragments_with_proba_tmp.append((pr, [frg[i_] + vec[i_] for i_ in range(chem_data.len_frag)]))
                        j += 1
                        (pi,vec) = list_fragments_with_proba[j]
                        pr = pi*pa
                        pr_ok = pr >= min_rel_pr
                    if j == len(list_fragments_with_proba)-1 and pr_ok:
                        list_all_fragments_with_proba_tmp.append((pr, [frg[i_] + vec[i_] for i_ in range(chem_data.len_frag)]))
                list_all_fragments_with_proba = list_all_fragments_with_proba_tmp
                list_all_fragments_with_proba.sort(reverse=True)
            else:
                list_all_fragments_with_proba = [(pa*pi, [frg[i_] + vec[i_] for i_ in range(chem_data.len_frag)]) for (pi,vec) in list_fragments_with_proba for (pa,frg) in list_all_fragments_with_proba]
        j_ += 1
    
    list_all_fragments = [frg for (_,frg) in list_all_fragments_with_proba]
    list_all_probas = [pr for (pr,_) in list_all_fragments_with_proba]
    return list_all_fragments, list_all_probas

def generate_each_isotopolog_pattern(frag):
    """Given a fragment, generate the list of isotopologues per element
    """
    idx_abundant = [ _ for _ in chem_data.idx_abun_isotopes if frag[_] > 0]
    list_patterns = [get_all_isotopic_patterns_one_element(idx, frag[idx]) for idx in idx_abundant]
    return idx_abundant, list_patterns

def unbounded_knapsack_rec_aux(mass_measured: float, mass_max: float, mass_min: float, vec_fragment: list, id_kn: int, idx_list_kn: list, list_results: list, max_no_each_atom: list=None, max_sum_no_remain_atoms: int=None, check_DBE: bool=False, ctr: int=0, ctr_DBE: int=0):
    """auxiliary recursive function for knapsack

    INPUT:
    - `mass_measured`: float (the target mass)
    - `mass_max`: upper bound on target mass
    - `mass_min`: lower bound on target mass
    - `vec_fragment`: current possibility (starts at [0,...,0] in the header function)
    - `id_kn`: index in idx_list_kn, and vec_fragment[idx_list_kn[idx_kn]] is >= 0, one step updates the value vec_fragment[idx_list_kn[id_kn]], min 0, max len(idx_list_kn)-1
    - `idx_list_kn`: list of indices to use in the list of masses of atoms (which of the periodic table to use)
       idx_list_kn[0], ..., idx_list_kn[id_kn], ..., idx_list_kn[len(idx_list_kn)-1] correspond to index in periodic table
    - `list_results`: a list of all possible vec_fragment
    - `max_no_each_atom`: array of positive int, max number for each atom (in case the initial molecule is known for example), array is None if no constraint
       max_no_each_atom[i] corresponds to an element in the fragment at index idx_list_kn[i]
    - `max_sum_no_remain_atoms`: positive integer or None. Maximum sum of numbers of atoms to set (because of valence constraints for example)
      so 0 <= a_i <= max_sum_no_remain_atoms and  0 <= sum(vec_fragment) <= max_sum_no_atoms
    - `check_DBE`: do the DBE test
    - `ctr`: a counter for enumerated possibilities so far
    - `ctr_DBE`: a counter for enumerated possibilities of valid DBE so far

    the masses are in increasing order in the lists, the index id_kn starts at the end (with the heaviest atom)
    """
    #id_kn is the index that is decreased by step of -1 (needed for recursive function)
    
    #change 20190827
    #idx_list_kn is a list of indices to use to pick masses for the knapsack algo.
    #id_kn should successively take all values of idx_list_kn (and no other values)
    # the values of mass_max and mass_min decrease (or increase) at each update of vec_fragment
    if id_kn < 0 or mass_max < 0 or (mass_min > 0 and mass_max < chem_data.mass[0]):
        # in this 'if', this is a dead-end: no solution found. That is why we
        # really need to check that mass_min > 0 AND mass_max < m0: it ensure
        # that the current mass is not yet in the valid interval (mass_min > 0)
        # but there is no more room to add any atom (mass_max < m0)
        return ctr, ctr_DBE
    # enumerate the possible a_i at step i
    id_m = idx_list_kn[id_kn]
    mass_i = chem_data.mass[id_m]

    if id_kn == 0: # this is the last and lightest atom to consider

        # check for each possible value of a_0 if the solution is valid
        max_ai = math.floor(mass_max/mass_i)
        min_ai = max(0,math.ceil(mass_min/mass_i))
        if min_ai < 0:
            print("max_ai = {}, min_ai = {}".format(max_ai, min_ai))
        if max_no_each_atom is not None and max_no_each_atom[id_kn] is not None:
            max_ai = min(max_ai, max_no_each_atom[id_kn])
            #if min_ai > max_ai:
            #    the minimal number of atoms required to fill the remaining mass
            #    is higher than the allowed number of atoms:
            #    the target mass cannot be reached
            #    return ctr,ctr_DBE
        if max_sum_no_remain_atoms is not None:
            max_ai = min(max_ai, max_sum_no_remain_atoms)
        #print("max_ai = {} min_ai = {}".format(max_ai, min_ai))
        for a_i in range(max_ai, min_ai-1, -1):
            vec_fragment[id_m] = a_i
            # now vec_fragment might be a solution, in particular, a_i = 0 is allowed

            if float(a_i)*mass_i <= mass_max and float(a_i)*mass_i >= mass_min:
                # this is a valid solution
                #mass = chem_data.get_mass(vec_fragment)
                vec = [ai for ai in vec_fragment]
                #print(vec)
                # copy of the row-solution otherwise it does not work within a recursive function
                ctr += 1
                if (not check_DBE) or (chem_data.get_DBE_value(vec) >= 0.0): 
                    ctr_DBE += 1
                    list_results.append(vec)
                    #if (chem_data.get_DBE_value(vec) < 0.0):
                        #print("l.442 {} added (current vec_fragment = {})".format(chem_data.get_string_formula(vec), chem_data.get_string_formula(vec_fragment)))
                    
                    #mass_diff=(mass-mass_measured)/mass_measured*10**6
                    #str_formula = chem_data.get_string_formula(vec_fragment)
                    #print(str_formula+'\t'+str(mass)+'\t'+str(DBE))

        return ctr,ctr_DBE

    if id_kn > 0:

        # choice 0 means we do not take this atom
        max_ai = math.floor(mass_max/mass_i)
        if max_no_each_atom is not None and max_no_each_atom[id_kn] is not None:
            max_ai = min(max_ai, max_no_each_atom[id_kn])
        if max_sum_no_remain_atoms is not None:
            if max_sum_no_remain_atoms <= 0: # there is no more room for more atoms (valence)
                # should we fill with zeros the remaining a_i and return the vector?
                if mass_min <= 0.0 and mass_max >= 0.0:# indeed somethimes we enter this branch
                    vec = [ai for ai in vec_fragment]
                    for ii in range(id_kn, -1, -1):
                        vec[idx_list_kn[ii]] = 0
                    ctr += 1
                    # I am not sure it is needed.
                    if (not check_DBE) or (chem_data.get_DBE_value(vec_fragment) >= 0.0): 
                        ctr_DBE += 1
                        list_results.append(vec)
                        #if (chem_data.get_DBE_value(vec) < 0.0):
                            #print("l.469 {} added (current vec_fragment = {})".format(chem_data.get_string_formula(vec), chem_data.get_string_formula(vec_fragment)))

                return ctr,ctr_DBE
            max_ai = min(max_ai, max_sum_no_remain_atoms)
            for a_i in range(max_ai, -1, -1):
                vec_fragment[id_m] = a_i
                ctr, ctr_DBE = unbounded_knapsack_rec_aux(mass_measured, mass_max-float(a_i)*mass_i, mass_min-float(a_i)*mass_i, vec_fragment, id_kn-1, idx_list_kn, list_results, max_no_each_atom=max_no_each_atom, max_sum_no_remain_atoms=max_sum_no_remain_atoms-a_i, check_DBE=check_DBE, ctr=ctr, ctr_DBE=ctr_DBE)
        else:
            for a_i in range(max_ai, -1, -1):
                vec_fragment[id_m] = a_i
                ctr, ctr_DBE = unbounded_knapsack_rec_aux(mass_measured, mass_max-float(a_i)*mass_i, mass_min-float(a_i)*mass_i, vec_fragment, id_kn-1, idx_list_kn, list_results, max_no_each_atom=max_no_each_atom, check_DBE=check_DBE, ctr=ctr, ctr_DBE=ctr_DBE)

        return ctr,ctr_DBE

def unbounded_knapsack_rec(mass_measured: float, mass_max: float, mass_min: float, idx_list_kn: list, max_no_each_atom: list=None, max_sum_no_atoms: int=None, check_DBE: bool=True):
    """Recursive knapsack function
    
    INPUT:
    - `mass_measured`: float (the target mass)
    - `mass_max`: upper bound on target mass
    - `mass_min`: lower bound on target mass
    - `idx_list_kn`: list of indices to use in the list of masses of atoms (which of the periodic table to use)
    - `max_no_each_atom`: array of positive int, max number for each atom (in case the initial molecule is known for example), array is None if no constraint
    - `max_sum_no_atoms`: positive integer or None. Maximum sum of numbers of atoms to set (because of valence constraints for example)
      0 <= sum(vec_fragment) <= max_sum_no_atoms
    - `check_DBE`: do the DBE test
    """
    vec_frag = [0 for _ in range(chem_data.len_frag)] #changed 20190827, the final lengh right at the beginning of the algo
    list_results = []
    #we start to loop at idx_list_kn[len(idx_list_kn)-1], index of heaviest atom to use for the knapsack algo
    #we need to give the list idx_list_kn as info as well so that it can find the next index to look at
    number_solutions, number_valid_DBE = unbounded_knapsack_rec_aux(mass_measured, mass_max, mass_min, vec_frag, len(idx_list_kn)-1, idx_list_kn, list_results, max_no_each_atom=max_no_each_atom, max_sum_no_remain_atoms=max_sum_no_atoms, check_DBE=check_DBE, ctr=0, ctr_DBE=0)
    
    return number_solutions, number_valid_DBE, list_results

def binary_search_sorted_list_of_tuples(A, m, t0):
    """Return the nearest index i such that A[i-1][t0] < m <= A[i][t0] for a sorted list of tuples A
    """
    i = 0
    j = len(A)-1
    if j == 0 or m > A[j][t0]:
        return j
    while i < j:
        k = (i+j)//2
        if A[k][t0] < m:
            i = k+1
        else:
            j = k
    return i

def binary_search_nearest_low_sorted_list_of_tuples(A, m, t0):
    """Return the nearest index i such that A[i][t0] <= m < A[i+1][t0] for a sorted list of tuples A
    """
    i = 0
    j = len(A)-1
    if j == 0 or m < A[0][t0]:
        return 0
    while i < j:
        k = (i+j)//2
        if m <= A[k][t0]:
            j = k
        else:
            i = k+1
    return i

def knapsack_double_with_lists(mass_max: float, mass_min: float, idx_list_kn: list, list_multi: list, list_mono: list, verbose: bool=False, double_check_DBE: bool=False) -> (int, int, list):
    """Knapsack optimized (twos knapsacks on multi-valent, then mono-valent atoms)

    INPUT:
    - `mass_max`: upper bound on target mass
    - `mass_min`: lower bound on target mass
    - `idx_list_kn`: list of indices to use in the list of masses of atoms (which of the periodic table to use)
    - `list_multi`: multi-valent list of tuples (mass, 2*DBE, fragment) sorted in increasing order of mass and covering the entire range of masses from 1.0 to mass_max
    - `list_mono`: mono-valent list of tuples (mass, number of atoms, fragment) sorted in increasing order of mass and covering the entire range of masses from 1.0 to mass_max
    - `verbose`: print intermediate data
    - `double_check_DBE`: only >= DBE values are expected, but still re-compute DBE and double-check >= 0

    OUTPUT:
    - `ctr`: counter of fragments of valid mass
    - `ctr_DBE`: counter of fragments of valid mass and >= DBE (non-negative BDE)
    - `list_results`: the list of fragments of valid mass and non-negative DBE
    """
    # now read the two lists and find pairs. Note: the zero-vector is not in the respective lists.
    # the two lists are sorted in increasing order of mass.
    idx_list_kn_multi_val = [ei for ei in idx_list_kn if chem_data.valence[ei] > 1]
    idx_list_kn_mono_val = [ei for ei in idx_list_kn if chem_data.valence[ei] == 1]
    mass_atoms_frag_mono_val = list_mono
    mass_DBE2_frag_multi_val = list_multi
    list_results = []
    ctr = 0
    ctr_DBE = 0
    i_multi = 0
    i_mono = len(mass_atoms_frag_mono_val)-1
    while i_multi < len(mass_DBE2_frag_multi_val) and i_mono >= 0:
        # decrease i_mono until the mass is lower than mass_max
        while i_mono >= 0 and mass_DBE2_frag_multi_val[i_multi][0] + mass_atoms_frag_mono_val[i_mono][0] > mass_max:
            i_mono -= 1
        # here, mi <= mass_max: loop over i_mono (decreasing), or i_mono < 0 --> end of the loop
        while i_mono >= 0 and mass_DBE2_frag_multi_val[i_multi][0] + mass_atoms_frag_mono_val[i_mono][0] >= mass_min:
            # now there is a match, check DBE compatibility
            if mass_atoms_frag_mono_val[i_mono][1] <= mass_DBE2_frag_multi_val[i_multi][1]:
                vec_fragment = [0]*chem_data.len_frag
                for j in idx_list_kn_mono_val:
                    vec_fragment[j] = mass_atoms_frag_mono_val[i_mono][2][j]
                for j in idx_list_kn_multi_val:
                    vec_fragment[j] = mass_DBE2_frag_multi_val[i_multi][2][j]
                #vec_fragment = [mass_DBE2_frag_multi_val[i_multi][2][j] + mass_atoms_frag_mono_val[i_mono][2][j] for j in range(chem_data.len_frag)]
                ctr += 1
                DBE2 = 0
                if double_check_DBE:
                    DBE2 = chem_data.get_DBE_value2(vec_fragment)
                if DBE2 >= 0:
                    ctr_DBE += 1
                    list_results.append(vec_fragment)
            i_mono -= 1
        # now, increase i_multi: it will increase mi (total paired mass).
        # Store i_mono, then increase i_mono as long as the total mass stays below the max mass, find all matches,
        # then reset i_mono to the saved value, and restart the whole process (go back to the beginning of the top while loop)
        i_multi += 1
        if i_multi >= len(mass_DBE2_frag_multi_val): # end of the while loop
            continue
        while i_mono < len(mass_atoms_frag_mono_val) and (mass_DBE2_frag_multi_val[i_multi][0] + mass_atoms_frag_mono_val[i_mono][0]) < mass_min:
            i_mono += 1
        i_mono_tmp = i_mono
        while i_mono < len(mass_atoms_frag_mono_val) and (mass_DBE2_frag_multi_val[i_multi][0] + mass_atoms_frag_mono_val[i_mono][0]) <= mass_max:
            # now there is a match, check DBE compatibility
            if mass_atoms_frag_mono_val[i_mono][1] <= mass_DBE2_frag_multi_val[i_multi][1]:
                vec_fragment = [0]*chem_data.len_frag
                for j in idx_list_kn_mono_val:
                    vec_fragment[j] = mass_atoms_frag_mono_val[i_mono][2][j]
                for j in idx_list_kn_multi_val:
                    vec_fragment[j] = mass_DBE2_frag_multi_val[i_multi][2][j]
                #vec_fragment = [mass_DBE2_frag_multi_val[i_multi][2][j] + mass_atoms_frag_mono_val[i_mono][2][j] for j in range(chem_data.len_frag)]
                ctr += 1
                DBE2 = 0
                if double_check_DBE:
                    DBE2 = chem_data.get_DBE_value2(vec_fragment)
                if DBE2 >= 0:
                    ctr_DBE += 1
                    list_results.append(vec_fragment)
            i_mono += 1
        i_mono = i_mono_tmp - 1
        # now explore the matches while decreasing i_mono: it is done by starting again the while loop
    # now consider the multi-valents but without mono-valent: like if i_mono = -1
    i_multi = len(mass_DBE2_frag_multi_val)-1
    while i_multi >= 0 and mass_DBE2_frag_multi_val[i_multi][0] > mass_max:
        i_multi -= 1
    while i_multi >= 0 and mass_DBE2_frag_multi_val[i_multi][0] >= mass_min: # valid fragment
        vec_fragment = list(mass_DBE2_frag_multi_val[i_multi][2]) # copy
        ctr += 1
        DBE2 = 0
        if double_check_DBE:
            DBE2 = chem_data.get_DBE_value2(vec_fragment)
        if DBE2 >= 0:
            ctr_DBE += 1
            list_results.append(vec_fragment)
        i_multi -= 1
    i_mono = len(mass_atoms_frag_mono_val)-1
    while i_mono >= 0 and mass_atoms_frag_mono_val[i_mono][0] > mass_max:
        i_mono -= 1
    while i_mono >= 0 and mass_atoms_frag_mono_val[i_mono][0] >= mass_min: # valid fragment if DBE >= 0
        if mass_atoms_frag_mono_val[i_mono][1] <= 2:
            vec_fragment = list(mass_atoms_frag_mono_val[i_mono][2]) # copy
            ctr += 1
            DBE2 = 0
            if double_check_DBE:
                DBE2 = chem_data.get_DBE_value2(vec_fragment)
            if DBE2 >= 0:
                ctr_DBE += 1
                list_results.append(vec_fragment)
        i_mono -= 1

    if verbose:
        print("no_solutions = {}".format(len(list_results)))
        print("list_results = {}".format(", ".join([chem_data.get_string_formula(vec) for vec in list_results])))
    return ctr, ctr_DBE, list_results

def knapsack_double(mass_measured: float, mass_max: float, mass_min: float, idx_list_kn: list, max_no_each_atom: bool=None, return_lists:bool=False, verbose: bool=False, double_check_DBE: bool=False):
    """Knapsack optimized (twos knapsacks on multi-valent, then mono-valent atoms)
    
    INPUT:
    - `mass_measured`: float (the target mass)
    - `mass_max`: upper bound on target mass
    - `mass_min`: lower bound on target mass
    - `idx_list_kn`: list of indices to use in the list of masses of atoms (which of the periodic table to use)
    - `max_no_each_atom`: array of positive int, max number for each atom (in case the initial molecule is known for example), array is None if no constraint
    - `return_lists`: returns the list of multi-valent and the list of mono-valent atoms
    - `verbose`: print intermediate data
    - `double_check_DBE`: only >= DBE values are expected, but still re-compute DBE and double-check >= 0
    """
    #0. identify the multi-valent atoms
    #1. call knapsack on multi-valent atoms
    #2. call knapsack on mono-valent atoms
    idx_list_kn_multi_val = [ei for ei in idx_list_kn if chem_data.valence[ei] > 1]
    idx_list_kn_mono_val = [ei for ei in idx_list_kn if chem_data.valence[ei] == 1]
    if max_no_each_atom is None:
        max_no_each_atom_multi_val = None
        max_no_each_atom_mono_val = None
    else:
        max_no_each_atom_multi_val = [max_no_each_atom[i] for i in range(len(max_no_each_atom)) if chem_data.valence[idx_list_kn[i]] > 1]
        max_no_each_atom_mono_val = [max_no_each_atom[i] for i in range(len(max_no_each_atom)) if chem_data.valence[idx_list_kn[i]] == 1]
    
    # the first run is called with minumum mass 1.0 instead of 0.0 to avoid the zero-vector in the output list
    no_solutions_multi_val, no_valid_DBE_multi_val, list_results_multi_val = unbounded_knapsack_rec(mass_measured, mass_max, float(1.0), idx_list_kn_multi_val, max_no_each_atom=max_no_each_atom_multi_val, max_sum_no_atoms=None, check_DBE=False)
    if verbose:
        print("mass_measured: {} mass_max: {} mass_min: 0.0".format(mass_measured, mass_max))
        print("idx_list_kn_multi_val       = {} = {}".format(idx_list_kn_multi_val, ", ".join([chem_data.element[vi] for vi in idx_list_kn_multi_val])))
        print("idx_list_kn_mono_val        = {} = {}".format(idx_list_kn_mono_val, ", ".join([chem_data.element[vi] for vi in idx_list_kn_mono_val])))
        print("no_solutions_multi_val      = {}".format(no_solutions_multi_val))
        print("no_valid_DBE_multi_val      = {}".format(no_valid_DBE_multi_val))
        #print("len(list_results_multi_val) = {}".format(len(list_results_multi_val)))
        print("list_results_multi_val = {}".format(", ".join([chem_data.get_string_formula(vec) for vec in list_results_multi_val])))
    
    if verbose:
        list_DBEm = [chem_data.get_DBE_value2_at_indices(vec_fragment_multi_val, idx_list_kn_multi_val) for vec_fragment_multi_val in list_results_multi_val]
        list_DBEm.sort()
        list_DBEm_unique = []
        i = 0
        while i < len(list_DBEm):
            if list_DBEm[i] not in list_DBEm_unique:
                list_DBEm_unique.append(list_DBEm[i])
            i += 1
        print("list_DBEm = {}".format(list_DBEm_unique))
    # compute DBE and mass of each fragment made of multi-valent atoms, then sort according to mass
    mass_DBE2_frag_multi_val = [(chem_data.get_mass_at_indices(vec_fragment_multi_val, idx_list_kn_multi_val), chem_data.get_DBE_value2_at_indices(vec_fragment_multi_val, idx_list_kn_multi_val), vec_fragment_multi_val) for vec_fragment_multi_val in list_results_multi_val]
    mass_DBE2_frag_multi_val.sort()
    # now compute maximum DBE
    # future optimization: cut in 3 or more sub-sets with maximal mass and maximal DBE according to the first list.
    # That would need to find pivot masses at which the maximum DBE (up to that mass) changes.
    # then call one time the computation of the second list for each sub-interval and concatenate the results.
    if len(mass_DBE2_frag_multi_val) < 50:
        # one list
        dbe1 = max(2, max([mass_DBE2_frag_multi_val[i][1] for i in range(len(mass_DBE2_frag_multi_val))]))
        no_solutions_mono_val, _, list_results_mono_val = unbounded_knapsack_rec(mass_measured, mass_max, float(1.0), idx_list_kn_mono_val, max_no_each_atom=max_no_each_atom_mono_val, max_sum_no_atoms=dbe1, check_DBE=False)
        mass_atoms_frag_mono_val = [(chem_data.get_mass_at_indices(vec_fragment_mono_val, idx_list_kn_mono_val), sum([vec_fragment_mono_val[i] for i in idx_list_kn_mono_val]), vec_fragment_mono_val) for vec_fragment_mono_val in list_results_mono_val]
        mass_atoms_frag_mono_val.sort()
        if verbose:
            print("no_solutions_mono_val = {}".format(no_solutions_mono_val))
            print("list_results_mono_val = {}".format(", ".join([chem_data.get_string_formula(vec) for vec in list_results_mono_val])))
    else:
        delta_m = (mass_max-mass_DBE2_frag_multi_val[0][0])/float(3.0)
        m1 = mass_DBE2_frag_multi_val[0][0] + delta_m
        m2 = mass_DBE2_frag_multi_val[0][0] + 2*delta_m
        # find where the list mass_DBE2_frag_multi_val splits in 3 ranges of masses
        nearest_index1 = binary_search_sorted_list_of_tuples(mass_DBE2_frag_multi_val, m1, 0) # with tuple_index=0
        nearest_index2 = binary_search_sorted_list_of_tuples(mass_DBE2_frag_multi_val, m2, 0) # with tuple_index=0
        dbe1 = max(2, max([mass_DBE2_frag_multi_val[i][1] for i in range(nearest_index1+1)]))
        while (nearest_index1+1) < len(mass_DBE2_frag_multi_val) and mass_DBE2_frag_multi_val[nearest_index1+1][1] <= dbe1:
            nearest_index1 += 1
        m1 = mass_DBE2_frag_multi_val[nearest_index1][0]
        dbe2 = max(2, max([mass_DBE2_frag_multi_val[i][1] for i in range(nearest_index1+1, nearest_index2+1)]))
        while (nearest_index2+1) < len(mass_DBE2_frag_multi_val) and mass_DBE2_frag_multi_val[nearest_index2+1][1] <= dbe2:
            nearest_index2 += 1
        m2 = mass_DBE2_frag_multi_val[nearest_index2][0]
        dbe3 = max(2, max([mass_DBE2_frag_multi_val[i][1] for i in range(nearest_index2+1, len(mass_DBE2_frag_multi_val))]))
        if verbose:
            print("binary_search_sorted_list_of_tuples m1 = {} dbe1 = {} m2 = {} dbe2 = {} max = {} dbemax = {}".format(m1, dbe1, m2, dbe2, mass_DBE2_frag_multi_val[-1][0], dbe3))
        # now we have 3 ranges:
        # 1.0 -- mass_DBE2_frag_multi_val[nearest_index1][0], of DBE max = dbe1
        # mass_DBE2_frag_multi_val[nearest_index1+1][0] -- mass_DBE2_frag_multi_val[nearest_index2][0], of DBE max = dbe2
        # mass_DBE2_frag_multi_val[nearest_index2+1][0] -- mass_DBE2_frag_multi_val[-1][0], of DBE max = max(all)
        # now get the complement mass and call the routine to list the mono-valents
        no_solutions_mono_val0, _, list_results_mono_val0 = unbounded_knapsack_rec(mass_measured, mass_max, mass_max-m1, idx_list_kn_mono_val, max_no_each_atom=max_no_each_atom_mono_val, max_sum_no_atoms=dbe1, check_DBE=False)
        no_solutions_mono_val1, _, list_results_mono_val1 = unbounded_knapsack_rec(mass_measured-m1, mass_max-m1, mass_max-m2, idx_list_kn_mono_val, max_no_each_atom=max_no_each_atom_mono_val, max_sum_no_atoms=dbe2, check_DBE=False)
        no_solutions_mono_val2, _, list_results_mono_val2 = unbounded_knapsack_rec(mass_measured-m2, mass_max-m2, float(1.0), idx_list_kn_mono_val, max_no_each_atom=max_no_each_atom_mono_val, max_sum_no_atoms=dbe3, check_DBE=False)
        no_solutions_mono_val = no_solutions_mono_val0 + no_solutions_mono_val1 + no_solutions_mono_val2
        # now compute the mass, sort, compute also the number of atoms per fragment: it should be <= DBE2(fragment multi val)
        mass_atoms_frag_mono_val0 = [(chem_data.get_mass_at_indices(vec_fragment_mono_val, idx_list_kn_mono_val), sum([vec_fragment_mono_val[i] for i in idx_list_kn_mono_val]), vec_fragment_mono_val) for vec_fragment_mono_val in list_results_mono_val0]
        mass_atoms_frag_mono_val0.sort()
        mass_atoms_frag_mono_val1 = [(chem_data.get_mass_at_indices(vec_fragment_mono_val, idx_list_kn_mono_val), sum([vec_fragment_mono_val[i] for i in idx_list_kn_mono_val]), vec_fragment_mono_val) for vec_fragment_mono_val in list_results_mono_val1]
        mass_atoms_frag_mono_val1.sort()
        mass_atoms_frag_mono_val2 = [(chem_data.get_mass_at_indices(vec_fragment_mono_val, idx_list_kn_mono_val), sum([vec_fragment_mono_val[i] for i in idx_list_kn_mono_val]), vec_fragment_mono_val) for vec_fragment_mono_val in list_results_mono_val2]
        mass_atoms_frag_mono_val2.sort()
        mass_atoms_frag_mono_val = mass_atoms_frag_mono_val2 + mass_atoms_frag_mono_val1 + mass_atoms_frag_mono_val0

        if verbose:
            print("no_solutions_mono_val      = {}, {}, {}".format(no_solutions_mono_val0, no_solutions_mono_val1, no_solutions_mono_val2))
            print("list_results_mono_val = {}\n    {}\n    {}".format(", ".join([chem_data.get_string_formula(vec) for vec in list_results_mono_val0]),
                                                                    ", ".join([chem_data.get_string_formula(vec) for vec in list_results_mono_val1]),
                                                                    ", ".join([chem_data.get_string_formula(vec) for vec in list_results_mono_val2])))
    ctr, ctr_DBE, list_results = knapsack_double_with_lists(mass_max, mass_min, idx_list_kn, mass_DBE2_frag_multi_val, mass_atoms_frag_mono_val, verbose=verbose, double_check_DBE=double_check_DBE)
    if return_lists:
        return ctr, ctr_DBE, list_results, mass_DBE2_frag_multi_val, mass_atoms_frag_mono_val
    else:
        return ctr, ctr_DBE, list_results

def bounded_integer_knapsack_rec_aux(vec_fragment: list, id_kn: int, idx_list_kn: list, list_results: list, max_no_each_atom: list, max_sum_no_remain_atoms: int=None, min_mass: float=None, ctr: int=0) -> int:
    """Recursive sub-function for bounded_integer_knapsack_rec

    Enumerate all valid integer values at index idx_list_kn[id_kn] in vectors of length chem_data.len_frag which have same indices as vec_fragment[idx_list_kn[id_kn]+1:]
    and value 0 <= v_s(i) <= max_no_each_atom[i] for indices s(i)=idx_list_kn[i]

    INPUT:
    - `vec_fragment`: common data at indices [idx_list_kn[id_kn]+1:]
    - `id_kn`: index where to enumerate possible values, the index decreases
    - `idx_list_kn`: list of indices of atoms to consider
    - `list_results`: append the new results to this list
    - `max_no_each_atom`: maximum nuber of each atom
    - `max_sum_no_atoms`: total maximal number of atoms
    - `min_mass`: minimum mass or None
    - `ctr`: a counter of the number of solutions

    OUTPUT: the number of solutions
    """
    if id_kn < 0 or (max_sum_no_remain_atoms is not None and max_sum_no_remain_atoms < 0):
        return ctr
    # enumerate the possible a_i at index i
    id_m = idx_list_kn[id_kn]
    max_ai = max_no_each_atom[id_kn]
    if max_sum_no_remain_atoms is not None:
        max_ai = min(max_ai, max_sum_no_remain_atoms)
    if id_kn == 0: # this is the last index to consider
        # check non-zero vector:
        a_i = 0
        vec_fragment[id_m] = 0
        if sum(vec_fragment) > 0 and (min_mass is None or chem_data.get_mass(vec_fragment) >= min_mass):
            ctr += 1
            vec = [ai for ai in vec_fragment]
            list_results.append(vec)
        min_ai = 1
        if min_mass is not None:
            mass_i = chem_data.get_mass(vec_fragment)
            if mass_i < min_mass:
                min_ai = int(math.ceil((min_mass - mass_i)/chem_data.mass[id_m])) # ceil, not floor
        for a_i in range(max_ai, min_ai-1, -1):
            vec_fragment[id_m] = a_i
            vec = [ai for ai in vec_fragment]
            ctr += 1
            list_results.append(vec)
        return ctr
    if id_kn > 0:
        if max_sum_no_remain_atoms is not None and max_sum_no_remain_atoms == 0:
            # finish here
            vec = [ai for ai in vec_fragment]
            for ii in range(id_kn, -1, -1):# fill with 0 the remaining indices, including the current one
                vec[idx_list_kn[ii]] = 0
            if sum(vec) > 0 and (min_mass is None or chem_data.get_mass(vec) >= min_mass):
                ctr += 1
                list_results.append(vec)
            return ctr
        for a_i in range(max_ai, -1, -1):
            vec_fragment[id_m] = a_i
            if max_sum_no_remain_atoms is not None:
                max_rem = max_sum_no_remain_atoms - a_i
            else:
                max_rem = None
            ctr = bounded_integer_knapsack_rec_aux(vec_fragment, id_kn-1, idx_list_kn, list_results, max_no_each_atom, max_sum_no_remain_atoms=max_rem, min_mass=min_mass, ctr=ctr)
        return ctr

def bounded_integer_knapsack_rec(idx_list_kn: list, max_no_each_atom: list, max_sum_no_atoms: int=None, min_mass: float=None) -> (list, int):
    """Enumerate all vectors of length chem_data.len_frag and value 0 <= v_s(i) <= max_no_each_atom[i] for indices s(i)=idx_list_kn[i]

    INPUT:
    - `idx_list_kn`: list of indices of atoms to consider
    - `max_no_each_atom`: maximum nuber of each atom
    - `max_sum_no_atoms`: total maximal number of atoms
    - `min_mass`: minimum mass or None

    OUTPUT: the number of solutions and the list of solutions without particular order
    """
    vec_frag = [0 for _ in range(chem_data.len_frag)] # we could restrict to len(idx_list_kn)
    list_results = []
    # no need of a recursive function anymore, an enumeration can do it
    # See knuth
    number_solutions = bounded_integer_knapsack_rec_aux(vec_frag, len(idx_list_kn)-1, idx_list_kn, list_results, max_no_each_atom, max_sum_no_remain_atoms=max_sum_no_atoms, min_mass=min_mass, ctr=0)
    return number_solutions, list_results

def enumerate_all_valid_subfragments(fragment: list, min_mass: float=None, return_partial_lists: bool=False, verbose: bool=False) -> (int, list):
    """Enumerate all valid (positive DBE) subfragments of mass >= `min_mass`

    Uses a knapsack-double-like algorithm but with integers (the number
    of occurences of each atom) instead of masses as float.

    INPUT:
    - `fragment`: a list of integers encoding a chemical fragment
    - `min_mass`: the minimum positive mass (e.g. 24.0), or None
    - `return_partial_lists`: return the partial lists of multi-valents
    and of mono-valents

    OUTPUT: the list of DBE-valid subfragments and their number, of
    mass >= min_mass
    """
    if sum(fragment) == 1:
        if not return_partial_lists:
            return 1, [fragment]
        else:
            return 1, [fragment], None, None
    idx_list_kn = [i for i in range(len(fragment)) if fragment[i] > 0]
    # 1. identify the multi-valent atoms
    idx_list_kn_multi_val = [ei for ei in idx_list_kn if chem_data.valence[ei] > 1]
    idx_list_kn_mono_val = [ei for ei in idx_list_kn if chem_data.valence[ei] == 1]
    max_no_each_atom_multi_val = [fragment[i] for i in idx_list_kn_multi_val]
    max_no_each_atom_mono_val = [fragment[i] for i in idx_list_kn_mono_val]
    # 2. enumerate all subfragments of multi-valent atoms
    # if there is no mono-valent atom in the fragment, consider mass >= min_mass
    if len(max_no_each_atom_mono_val) == 0 or sum(max_no_each_atom_mono_val) == 0:
        no_solutions_multi_val, list_results_multi_val = bounded_integer_knapsack_rec(idx_list_kn_multi_val, max_no_each_atom_multi_val, max_sum_no_atoms=None, min_mass=min_mass)
    else:
        no_solutions_multi_val, list_results_multi_val = bounded_integer_knapsack_rec(idx_list_kn_multi_val, max_no_each_atom_multi_val, max_sum_no_atoms=None, min_mass=None)
    zero_found = False
    for vec_fragment in list_results_multi_val:
        zero_found = zero_found or sum(vec_fragment) == 0
    if zero_found:
        raise ValueError("bounded_integer_knapsack_rec returned the zero-vector on input {}, max_number_each_atom = {}".format([chem_data.element[i] for i in idx_list_kn_multi_val], max_no_each_atom_multi_val))
    # 3. for each atom, compute DBE, deduce max number of mono-valent
    if len(max_no_each_atom_mono_val) == 0 or sum(max_no_each_atom_mono_val) == 0:
        # no mono-valent atom in the fragment, skip this step
        ctr_DBE = no_solutions_multi_val
        if not return_partial_lists:
            return ctr_DBE, list_results_multi_val
        else:
            # is the list sorted in any order? No. Sort it in decreasing order of mass
            partial_list_multi_val = [(chem_data.get_mass(fi), fi) for fi in list_results_multi_val]
            partial_list_multi_val.sort()
            return ctr_DBE, list_results_multi_val, partial_list_multi_val, {}
    list_results = []
    partial_list_multi_valent = [] # a list of tuples (mass, fragment)
    # for mono-valent, dict indexed by total number of atoms (that is min valence of associated multi-valent piece)
    ctr = 0
    ctr_DBE = 0
    if len(list_results_multi_val) > 0:
        #DBEf = chem_data.get_DBE_value(vec_fragment_multi_val)
        #DBEm = int(math.ceil(DBEf*2.0))
        DBEm_tab = [chem_data.get_DBE_value2_at_indices(vec_fragment_multi_val, idx_list_kn_multi_val) for vec_fragment_multi_val in list_results_multi_val]
        max_DBEm = max(DBEm_tab)
        #run knapsack on mono-valent atoms with two bounds: their total number and the max number of each of them
        no_solutions_mono_val, list_results_mono_val = bounded_integer_knapsack_rec(idx_list_kn_mono_val, max_no_each_atom_mono_val, max_sum_no_atoms=max_DBEm, min_mass=None)
        list_results_mono_val = [(sum(fi), fi) for fi in list_results_mono_val]
        list_results_mono_val.sort() # sort by increasing total number of atoms
        # make batches, one for each total_number_atoms
        dict_results_mono_val = {}
        i_ = 0
        total_number_atoms = list_results_mono_val[i_][0]
        j_ = i_
        while i_ < len(list_results_mono_val):
            j_ = i_
            while j_+1 < len(list_results_mono_val) and list_results_mono_val[j_+1][0] == total_number_atoms:
                j_ += 1
            dict_results_mono_val[total_number_atoms] = [(chem_data.get_mass(fi), fi) for (_, fi) in list_results_mono_val[i_:j_+1]]
            dict_results_mono_val[total_number_atoms].sort(reverse=True) # sort by decreasing mass
            j_ += 1
            i_ = j_
            if i_ < len(list_results_mono_val):
                total_number_atoms = list_results_mono_val[i_][0]
        # now dict_results_mono_val is a dictionary indexed by the number of atoms of fragments, and each batch is in decreasing order of mass
        for i_ in range(len(list_results_multi_val)):
            vec_fragment_multi_val = list_results_multi_val[i_]
            used_multi_val_fragment = False
            DBEm = DBEm_tab[i_]
            if min_mass is not None or return_partial_lists:
                mass_i = chem_data.get_mass_at_indices(vec_fragment_multi_val, idx_list_kn_multi_val)
            else:
                mass_i = 0.0
            if min_mass is None or mass_i >= min_mass:
                ctr += 1
                ctr_DBE += 1
                list_results.append(vec_fragment_multi_val) # without mono-valent atoms
                used_multi_val_fragment = True
            # matches partial multi-valent fragment with sets of mono-valent atoms
            #if DBEm <= 0.0: no assemblage possible, but it should not append with only multi-valent atoms
            if DBEm > 0:
                if min_mass is not None:
                    min_mass_i = max(0.0, min_mass - mass_i) # min_mass for the mono-valent part added later
                else:
                    min_mass_i = 0.0
                for DBEi in [_ for _ in range(1, DBEm+1) if _ in dict_results_mono_val]:
                    j_ = 0
                    while j_ < len(dict_results_mono_val[DBEi]) and dict_results_mono_val[DBEi][j_][0] >= min_mass_i: # the items are ordered by decreasing mass
                        # append the solution
                        used_multi_val_fragment = True
                        vec_fragment_mono_val = dict_results_mono_val[DBEi][j_][1]
                        vec_fragment = [vec_fragment_multi_val[i] + vec_fragment_mono_val[i] for i in range(chem_data.len_frag)]
                        ctr += 1
                        DBE = chem_data.get_DBE_value(vec_fragment) # double-check that DBE is >= 0
                        if DBE >= 0.0:
                            ctr_DBE += 1
                            list_results.append(vec_fragment)
                        j_ += 1
            if return_partial_lists and used_multi_val_fragment:
                partial_list_multi_valent.append((mass_i, vec_fragment_multi_val, DBEm))
        partial_list_multi_valent.sort()
        # 5. Now consider mono-valent only atoms: they can form a fragment of at most 2 elements
        # because bounded_integer_knapsack_rec is not supposed to return the zero-fragment
        if min_mass is None:
            min_mass_i = 0.0
        else:
            min_mass_i = min_mass
        for DBEi in [_ for _ in [1, 2] if _ in dict_results_mono_val]: # mono-valent atoms can form at most a pair
            j_ = 0
            while j_ < len(dict_results_mono_val[DBEi]) and dict_results_mono_val[DBEi][j_][0] >= min_mass_i:
                ctr += 1
                vec_fragment_mono_val = dict_results_mono_val[DBEi][j_][1]
                DBE = chem_data.get_DBE_value2_at_indices(vec_fragment_mono_val, idx_list_kn_mono_val) # double-check that DBE is >= 0
                if DBE >= 0:
                    ctr_DBE += 1
                    list_results.append(vec_fragment_mono_val)
                j_ += 1
    else: # only mono-valent atoms
        no_solutions_mono_val, list_results_mono_val = bounded_integer_knapsack_rec(idx_list_kn_mono_val, max_no_each_atom_mono_val, max_sum_no_atoms=2, min_mass=min_mass)
        for vec_fragment_mono_val in list_results_mono_val:
            ctr += 1
            DBE = chem_data.get_DBE_value2_at_indices(vec_fragment_mono_val, idx_list_kn_mono_val) # double-check that DBE is >= 0
            if DBE >= 0.0:
                ctr_DBE += 1
                list_results.append(vec_fragment_mono_val)
        if return_partial_lists: # make batches like above, one for each total_number_atoms
            list_results_mono_val = [(sum(fi), fi) for fi in list_results_mono_val]
            list_results_mono_val_1 = [(chem_data.get_mass_at_indices(fi, idx_list_kn_mono_val),fi) for (si, fi) in list_results_mono_val if si == 1]
            list_results_mono_val_2 = [(chem_data.get_mass_at_indices(fi, idx_list_kn_mono_val),fi) for (si, fi) in list_results_mono_val if si == 2]
            dict_results_mono_val = {}
            if len(list_results_mono_val_1) > 0:
                list_results_mono_val_1.sort(reverse=True) # reverse order-> decreasing mass
                dict_results_mono_val[1] = list_results_mono_val_1
            if len(list_results_mono_val_2) > 0:
                list_results_mono_val_2.sort(reverse=True) # reverse order-> decreasing mass
                dict_results_mono_val[2] = list_results_mono_val_2

    if not return_partial_lists:
        return ctr_DBE, list_results
    else:
        return ctr_DBE, list_results, partial_list_multi_valent, dict_results_mono_val

def enumerate_all_valid_subfragments_from_lists(fragment: list, multi_valent: list, mono_valent: dict, parent_fragment: list=None, min_mass: float=None, verbose: bool=False) -> (int, list, list, dict):
    """Computes from partial lists the subfragments of fragment
    
    Filter the list and dict of subfragments of multi-val and mono-val
    and obtain the list of sub-fragments of fragment

    INPUT:
    - `fragment`: a list of integers encoding a chemical fragment
    - `multi_valent`: a list of pairs (mass, fragment) sorted by increasing mass
    - `mono_valent`: a dictionary indexed by the total number of atoms is the set, then each batch sorted by decreasing mass
    - `min_mass`: the minimum positive mass (e.g. 24.0), or None

    OUTPUT: the list of DBE-valid subfragments and their number, of
    mass >= min_mass, the list of partial fragments made of multi-
    valents, and the dict of mono-valents
    """
    if verbose:
        print("fragment: {}".format(chem_data.get_string_formula(fragment)))
        print("multi_valent: {}".format([chem_data.get_string_formula(fi[1]) for fi in multi_valent]))
        print("mono_valent: {}".format({dbe: [chem_data.get_string_formula(fi[1]) for fi in mono_valent[dbe]] for dbe in mono_valent}))
    idx_list_kn = [i for i in range(len(fragment)) if fragment[i] > 0]
    # total numer of atoms in the fragment
    number_atoms = sum([fragment[i] for i in idx_list_kn])
    if number_atoms <= 1:
        return 1, [fragment], None, None
    # identify the multi-valent atoms
    idx_multi_val = [ei for ei in idx_list_kn if chem_data.valence[ei] > 1]
    idx_mono_val = [ei for ei in idx_list_kn if chem_data.valence[ei] == 1]

    number_atoms_multi_val = sum([fragment[i] for i in idx_multi_val])
    number_atoms_mono_val = sum([fragment[i] for i in idx_mono_val])
    
    if parent_fragment is not None:
        idx_list_kn_parent = [i for i in range(len(parent_fragment)) if parent_fragment[i] > 0]
        idx_multi_val_parent = [ei for ei in idx_list_kn_parent if chem_data.valence[ei] > 1]
        idx_mono_val_parent = [ei for ei in idx_list_kn_parent if chem_data.valence[ei] == 1]
    else:
        idx_multi_val_parent = chem_data.idx_multi_valent
        idx_mono_val_parent = chem_data.idx_mono_valent

    # 0. fragment multi-val, fragment mono-val
    frag_multi_val = [0]*chem_data.len_frag
    frag_mono_val = [0]*chem_data.len_frag
    for i in idx_multi_val:
        frag_multi_val[i] = fragment[i]
    for i in idx_mono_val:
        frag_mono_val[i] = fragment[i]
    mass_multi_val = chem_data.get_mass_at_indices(frag_multi_val, idx_multi_val)
    mass_mono_val = chem_data.get_mass_at_indices(frag_mono_val, idx_mono_val)
    # 1. filter the list of multi-valent
    if number_atoms_multi_val <= 1:
        partial_list_multi_val = [(chem_data.get_mass_at_indices(frag_multi_val, idx_multi_val), frag_multi_val, chem_data.get_DBE_value2_at_indices(frag_multi_val, idx_multi_val)) ]
    else:
        partial_list_multi_val = []
        j = 0
        while j < len(multi_valent) and multi_valent[j][0] <= mass_multi_val:
            # compare if subfragment
            k = 0
            # do to make the test with idx_multi_val of fragment because there are atoms in the list that are not in fragment
            while k < len(idx_multi_val_parent) and multi_valent[j][1][idx_multi_val_parent[k]] <= fragment[idx_multi_val_parent[k]]:
                k += 1
            if k >= len(idx_multi_val_parent):
                partial_list_multi_val.append(multi_valent[j])
            j += 1
    if verbose:
        print("filtered list of multi-valent: {}".format([chem_data.get_string_formula(fi[1]) for fi in partial_list_multi_val]))
    # if there is no mono-valent in the fragment, stops here.
    if len(idx_mono_val) == 0:
        result = [ri[1] for ri in partial_list_multi_val] # ri[0] is the mass
        return len(result), result, partial_list_multi_val, {}
    if verbose and len(idx_multi_val) > 0 and len(partial_list_multi_val) == 0:
        print("case len(idx_multi_val) = {} > 0 but len(partial_list_multi_val) == 0".format(len(idx_multi_val)))
        print("fragment = {}".format(chem_data.get_string_formula(fragment)))
        print("idx_multi_val <-> {}".format([chem_data.element[i] for i in idx_multi_val]))
        print("input multi_valent = {}".format([chem_data.get_string_formula(fi[1]) for fi in multi_valent]))
    if len(idx_multi_val) == 0 or len(partial_list_multi_val) == 0:
        max_DBE = 2
    else:
        # compute max DBE
        max_DBE = max(2, max([ri[2] for ri in partial_list_multi_val]))
    # 2. filter the list of mono-valent -> needs DBE. store DBE in the list.
    dict_mono_val = {}
    for dbe in [_ for _ in range(1, max_DBE+1) if _ in mono_valent]:
        j = len(mono_valent[dbe])-1 # start with lightest fragment
        mono_val_dbe = []
        while j >= 0 and mono_valent[dbe][j][0] <= mass_mono_val:
            if min_mass is not None and len(idx_multi_val) == 0 and mono_valent[dbe][j][0] < min_mass:
                j -= 1
                continue
            # compare if subfragment
            k = 0
            while k < len(idx_mono_val_parent) and mono_valent[dbe][j][1][idx_mono_val_parent[k]] <= fragment[idx_mono_val_parent[k]]:
                k += 1
            if k >= len(idx_mono_val_parent):
                mono_val_dbe.append(mono_valent[dbe][j])
            j -= 1
        if len(mono_val_dbe) > 0:
            mono_val_dbe.reverse()
            dict_mono_val[dbe] = mono_val_dbe
    # now multi-val and mono-val are filtered. Now combine both.
    if len(idx_multi_val) == 0: # only mono-vals, return dict
        results = [dict_mono_val[dbe][i][1] for dbe in dict_mono_val for i in range(len(dict_mono_val[dbe]))]
        return len(results), results, [], dict_mono_val

    list_results = []
    res_partial_list_multi_val = []
    for i_ in range(len(partial_list_multi_val)):
        used_multi_val_fragment = False
        mass_i, vec_fragment_multi_val, DBEm = partial_list_multi_val[i_]
        if min_mass is None or mass_i >= min_mass:
            list_results.append(vec_fragment_multi_val) # without mono-valent atoms
            used_multi_val_fragment = True
        if DBEm > 0:
            if min_mass is not None:
                min_mass_i = min_mass - mass_i # min_mass for the mono-valent part added later
            else:
                min_mass_i = 0.0
            for DBEi in [_ for _ in range(1, DBEm+1) if _ in dict_mono_val]:
                j_ = 0
                while j_ < len(dict_mono_val[DBEi]) and dict_mono_val[DBEi][j_][0] >= min_mass_i: # the items are ordered by decreasing mass
                    # append the solution
                    used_multi_val_fragment = True
                    vec_fragment_mono_val = dict_mono_val[DBEi][j_][1]
                    vec_fragment = [vec_fragment_multi_val[i] + vec_fragment_mono_val[i] for i in range(chem_data.len_frag)]
                    #DBE = chem_data.get_DBE_value(vec_fragment) # double-check that DBE is >= 0
                    #if DBE >= 0.0:
                    list_results.append(vec_fragment)
                    j_ += 1
        if used_multi_val_fragment:
            res_partial_list_multi_val.append((mass_i, vec_fragment_multi_val, DBEm))
    # add monovalent-only subfragments
    if min_mass is None:
        min_mass_i = 0.0
    else:
        min_mass_i = min_mass
    for DBEi in [_ for _ in [1, 2] if _ in dict_mono_val]: # mono-valent atoms can form at most a pair
        j_ = 0
        while j_ < len(dict_mono_val[DBEi]) and dict_mono_val[DBEi][j_][0] >= min_mass_i:
            vec_fragment_mono_val = dict_mono_val[DBEi][j_][1]
            #DBE = chem_data.get_DBE_value(vec_fragment_mono_val) # double-check that DBE is >= 0
            #if DBE >= 0.0:
            list_results.append(vec_fragment_mono_val)
            j_ += 1

    return len(list_results), list_results, res_partial_list_multi_val, dict_mono_val

