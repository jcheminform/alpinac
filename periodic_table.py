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


Periodic table of elements that can be in the air.
The exact mass and the isotopic relative abundances are from IUPAC, see ref.:
Meija, J., Coplen, T.B., Berglund, M., Brand, W.A., Bievre,
P.D., Greoning, M., Holden, N.E., Irrgeher, J., Loss,
R.D., Walczyk, T., Prohaska, T.:
Atomic weights of the elements 2013 (IUPAC technical report).
Pure andApplied Chemistry 88(3), 265–291 (2016).
doi:10.1515/pac-2015-0305
"""
labels = ("element","mass","unit_mass","valence_index","occurence","atomic_number")
dic_labels = {"element":0,"mass":1,"unit_mass":2,"valence_index":3,"occurence":4,"atomic_number":5}
periodic_table = [
    ("H",         1.0078250319,  1,  1, 0.99984426, 1),  #0
    ("[2H]",      2.0141017779,  2,  1, 0.00015574, 1),  #1
    ("[10B]",    10.0129371,    10,  3, 0.1982,     5),  #2
    ("B",        11.00930536,   11,  3, 0.8018,     5),  #3
    ("C",        12,            12,  4, 0.988922,   6),  #4
    ("[13C]",    13.00335484,   13,  4, 0.011078,   6),  #5
    ("N",        14.0030740074, 14,  3, 0.996337,   7),  #6
    ("[15N]",    15.000108973,  15,  3, 0.003663,   7),  #7
    ("O",        15.9949146223, 16,  2, 0.9976206,  8),  #8
    ("[17O]",    16.9991315,    17,  2, 0.000379,   8),
    ("[18O]",    17.9991604,    18,  2, 0.0020004,  8),
    ("F",        18.99840316,   19,  1, 1,          9),  #11
    ("Ne",       19.99244018,   20,  0, 0.9048,    10),
    ("[21Ne]",   20.99384669,   21,  0, 0.0027,    10),
    ("[22Ne]",   21.99138511,   22,  0, 0.0925,    10),
    ("Si",       27.97692649,   28,  4, 0.9222968, 14),  #15
    ("[29Si]",   28.97649468,   29,  4, 0.0468316, 14),
    ("[30Si]",   29.97377018,   30,  4, 0.0308716, 14),
    ("P",        30.973762,     31,  3, 1,         15),
    ("S",        31.97207073,   32,  6, 0.9504074, 16),  #19,  we need valence = 6
    ("[33S]",    32.97145854,   33,  6, 0.0074869, 16),
    ("[34S]",    33.96786687,   34,  6, 0.0419599, 16),  #21
    ("Cl",       34.96885271,   35,  1, 0.757647,  17),  #22
    ("[36S]",    35.96708088,   36,  6, 0.0001458, 16),
    ("[36Ar]",   35.9675451,    36,  0, 0.00334,   18),
    ("[37Cl]",   36.9659026,    37,  1, 0.242353,  17),
    ("[38Ar]",   37.96273211,   38,  0, 0.00062,   18),
    ("Ar",       39.96238312,   40,  0, 0.99604,   18),
    ("[78Kr]",   77.92036494,   78,  0, 0.0036,    36),
    ("Br",       78.9183376,    79,  1, 0.50686,   35),  #29
    ("[80Kr]",   79.91637808,   80,  0, 0.0229,    36),
    ("[81Br]",   80.9162897,    81,  1, 0.49314,   35),  #31
    ("[82Kr]",   81.91348273,   82, -1, 0.1159,    36),
    ("[83Kr]",   82.91412716,   83, -1, 0.115,     36),
    ("Kr",       83.91149773,   84, -1, 0.5699,    36),
    ("[86Kr]",   85.91061063,   86, -1, 0.1728,    36),
    ("[92Mo]",   91.90680796,   92, -1, 0.1465,    42),
    ("[94Mo]",   93.9050849,    94, -1, 0.0919,    42),
    ("[95Mo]",   94.90583877,   95, -1, 0.1587,    42),
    ("[96Mo]",   95.90467612,   96, -1, 0.1667,    42),
    ("[97Mo]",   96.90601812,   97, -1, 0.0958,    42),
    ("Mo",       97.90540482,   98, -1, 0.2429,    42),  #41
    ("[100Mo]",  99.9074718,   100, -1, 0.0974,    42),
    ("Ag",      106.9050916,   107, -1, 0.51839,   47),
    ("[109Ag]", 108.9047553,   109, -1, 0.48161,   47),
    ("[124Xe]", 123.905892,    124, -1, 0.00095,   54),
    ("[126Xe]", 125.9042983,   126, -1, 0.00089,   54), 
    ("I",       126.9044719,   127,  1, 1,         53),  #47
    ("[128Xe]", 127.903531,    128, -1, 0.0191,    54),
    ("[129Xe]", 128.90478086,  129, -1, 0.26401,   54),
    ("[130Xe]", 129.90350935,  130, -1, 0.04071,   54),
    ("[131Xe]", 130.90508406,  131, -1, 0.21232,   54),  #51
    ("Xe",      131.90415509,  132, -1, 0.26909,   54),  #52
    ("[134Xe]", 133.90539466,  134, -1, 0.10436,   54),  #53
    ("[136Xe]", 135.90721448,  136, -1, 0.08857,   54),  #54
    ("[185Re]", 184.9529545,   185, -1, 0.374,     75),  #55
    ("Re",      186.9557501,   187, -1, 0.626,     75)   #56
]

#define constants
#mass of a lost electron during ionisation process
electron_mass=float(0.00054858)
#mass of H+:
PRT_H = float(1.00782503) - electron_mass
len_frag = len(periodic_table) # 57
# abundance
abundance = [elem[dic_labels["occurence"]] for elem in periodic_table]
# mass
mass = [elem[dic_labels["mass"]] for elem in periodic_table]
# label (string code) of the element
element = [elem[dic_labels["element"]] for elem in periodic_table]
# valence
valence = [elem[dic_labels["valence_index"]] for elem in periodic_table]
atomic_number = [elem[dic_labels["atomic_number"]] for elem in periodic_table]
unit_mass = [elem[dic_labels["unit_mass"]] for elem in periodic_table]
# correspondance between string code and index
dict_element_idx = {e: i for (i,e) in enumerate(element)}
# correspondance between string code and mass
dict_element_mass = {e: mass[i] for (i,e) in enumerate(element)}

# correspondance of indices: given the abundant one, which ones are rare
dict_abundant_rare = {0: [1], 3: [2], 4: [5], 6: [7], 8: [10, 9], 11: [], 12: [14, 13], 15: [16, 17], 18: [], 19: [21, 20, 23], 22: [25], 27: [24, 26], 34: [35, 32, 33, 30, 28], 29: [31], 41: [39, 38, 36, 42, 40, 37], 43: [44], 47: [], 52: [49, 51, 53, 54, 50, 48, 45, 46], 56: [55]}
# correspondance of indices of rare elements to their abundant
dict_rare_abundant = {1: 0, 2: 3, 5: 4, 7: 6, 10: 8, 9: 8, 14: 12, 13: 12, 16: 15, 17: 15, 21: 19, 20: 19, 23: 19, 25: 22, 24: 27, 26: 27, 35: 34, 32: 34, 33: 34, 30: 34, 28: 34, 31: 29, 39: 41, 38: 41, 36: 41, 42: 41, 40: 41, 37: 41, 44: 43, 49: 52, 51: 52, 53: 52, 54: 52, 50: 52, 48: 52, 45: 52, 46: 52, 55: 56}

# the indices of abundant isotopes. Mono-isotopic atoms of abundance 1 are also given (F (11), P (18), I (47)).
idx_abun_isotopes = [0, 3, 4, 6, 8, 11, 12, 15, 18, 19, 22, 27, 34, 29, 41, 43, 47, 52, 56]
idx_rare_isotopes_i = [1, 2, 5, 7, 10, 9, 14, 13, 16, 17, 21, 20, 23, 25, 24, 26, 35, 32, 33, 30, 28, 31, 39, 38, 36, 42, 40, 37, 44, 49, 51, 53, 54, 50, 48, 45, 46, 55]
idx_rare_isotopes = [1, 2, 5, 7, 9, 10, 13, 14, 16, 17, 20, 21, 23, 24, 25, 26, 28, 30, 31, 32, 33, 35, 36, 37, 38, 39, 40, 42, 44, 45, 46, 48, 49, 50, 51, 53, 54, 55]
idx_mono_valent = [i for i in range(len(periodic_table)) if periodic_table[i][3] == 1]
idx_multi_valent = [i for i in range(len(periodic_table)) if periodic_table[i][3] > 1]
idx_zero_valent = [i for i in range(len(periodic_table)) if periodic_table[i][3] <= 0]

#order of atoms to display nice-looking formulae
ordered_elements_string_formula_sub = [4, 5, 6, 7, 8, 9, 10, 19, 20, 21, 23, 15, 16, 17, 18, 0, 1, 2, 3, 11, 22, 25, 29, 31, 47,
                                       12, 13, 14, 24, 26, 27, 28, 30, 32, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 
                                       48, 49, 50, 51, 52, 53, 54, 55, 56]

for i in range(len_frag):
    if i not in ordered_elements_string_formula_sub:
        ordered_elements_string_formula_sub.append(i)
    
def get_string_formula(vec_frag: list) -> str:
    """Standardized string encoding of chemical formula

    INPUT:
    - ``vec_frag``: a list of positive integers of length `len_frag`
    
    OUPTUT: a string encoding the chemical molecule
    """
    s = ""
    for i in range(len(vec_frag)):
        vi = vec_frag[i]
        if vi == 1:
            s+=element[i]
        elif vi > 1:
            s += element[i] + str(vi)
    return s

def get_string_formula_sub(vec_frag: list) -> str:
    """Standardized string encoding of chemical formula

    INPUT:
    - ``vec_frag``: a list of positive integers of length `len_frag`

    OUPTUT: a string encoding the chemical molecule with \Sub and \Sup
    
    The order of chemical elements is modified as well, to try to reflect order recommended by IUPAC.
    Multivalent atoms first, and then alphabetically (more or less).
    """
    s = ""
    #for i in range(len(vec_frag)):
    for i in ordered_elements_string_formula_sub:    
        vi = vec_frag[i]
        if vi == 1:
            s+=element[i]
        elif vi > 1:
            s += element[i] + "\\Sub{" + str(vi) + "}"
    return s

def get_mass(vec_frag: list) -> float:
    """Exact mass of the fragment

    INPUT:
    - ``vec_frag``: a list of positive integers of length `len_frag`
    
    OUPTUT: the exact mass of `vec_frag`
    """
    return float(sum([float(vec_frag[i])*mass[i] for i in range(len(vec_frag)) if vec_frag[i] > 0]))

def get_mass_at_indices(vec_frag: list, indices: list) -> float:
    """Exact mass of the fragment

    INPUT:
    - ``vec_frag``: a list of positive integers of length `len_frag`
    - ``indices``: a list of all the indices to consider for computing the mass
    
    OUPTUT: the exact mass of `vec_frag` considering only the given indices
    """
    return float(sum([float(vec_frag[i])*mass[i] for i in indices]))

def get_DBE_value(fragment: list) -> float:
    """Double-Bond-Equivalent (DBE) formula

    INPUT:
    - ``fragment``: a list of positive integers of length `len_frag`
    
    OUPTUT: the DBE value of the molecule `vec_frag`

    cf. p. 254 in Mass spectrometry, a textbook by Juergen H. Gross, 2007, eq. (6.9)
    'double bound equivalent': calculate number of rings or double bounds in a formula
    """
    free_bond = float(2.0)
    for i in range(len(fragment)):
        free_bond += fragment[i] * (valence[i] - 2)
    #to actually calculate and return DBE:
    free_bond /= float(2.0)
    return free_bond

def get_DBE_value_at_indices(fragment: list, indices: list) -> float:
    """Double-Bond-Equivalent (DBE) formula

    INPUT:
    - ``fragment``: a list of positive integers of length `len_frag`
    - ``indices``: a list of all the indices to consider for computing the DBE
    
    OUPTUT: the DBE value of the molecule `vec_frag` at indices

    cf. p. 254 in Mass spectrometry, a textbook by Juergen H. Gross, 2007, eq. (6.9)
    'double bound equivalent': calculate number of rings or double bounds in a formula
    """
    free_bond = float(2.0)
    for i in indices:
        free_bond += fragment[i] * (valence[i] - 2)
    #to actually calculate and return DBE:
    free_bond /= float(2.0)
    return free_bond

def get_DBE_value2(fragment: list) -> int:
    """Double-Bond-Equivalent (DBE) formula (times 2)

    INPUT:
    - ``fragment``: a list of positive integers of length `len_frag`
    
    OUPTUT: 2*(the DBE value of the molecule `vec_frag`)

    cf. p. 254 in Mass spectrometry, a textbook by Juergen H. Gross, 2007, eq. (6.9)
    'double bound equivalent': calculate number of rings or double bounds in a formula
    """
    return 2 + sum([fragment[i] * (valence[i] - 2) for i in range(len(fragment))])

def get_DBE_value2_at_indices(fragment: list, indices: list) -> int:
    """Double-Bond-Equivalent (DBE) formula (times 2)

    INPUT:
    - ``fragment``: a list of positive integers of length `len_frag`
    - ``indices``: a list of all the indices to consider for computing the DBE
    
    OUPTUT: 2*(the DBE value of the molecule `vec_frag`)

    cf. p. 254 in Mass spectrometry, a textbook by Juergen H. Gross, 2007, eq. (6.9)
    'double bound equivalent': calculate number of rings or double bounds in a formula
    """
    return 2 + sum([fragment[i] * (valence[i] - 2) for i in indices])

# an old comment:
#    This function does two calculations to estimate
#    if the candidate formula is chemically possible or not.
#    If the formula is chemically possible,
#    it calculates the DBE value (float number) and returns it.
#    If the formula is impossible, it returns (-1) which causes later on
#    the candidate formula to be eliminated.
#    First test: (i) use all multivalent (V>1) atoms to form a chain
#    and calculate number of remaining free slots to bind another atom.
#    (ii) Calculate number of monovalent atoms (V = 1) and then assign DBE = -1
#    if number of monovalent > free slots for atoms (i.e., all monovelents
#    cannot be bound to a unique chain).
#    Second test: DBE: cf. p. 254 in Mass spectrometry, a textbook by Juergen H. Gross, 2007, eq. (6.9)
#    'double bound equivalent': calculate number of rings or double bounds in a formula
#    Note 20190821: using DBE < 0 as criteria for elimination of the candidate is not enough
#    because e.g. in H2F2, DBE > 0 but chemically this is impossible.
#    
#    For compound identification (Gross p. 254):
#    the DBE algorithm produces integers for odd-electron ions and molecules
#    and non-integers for even-electron ions (that have to rounded to the next lower integer)
#    The molecular ion (if present) has to be an odd-electron ion, with integer DBE value.

def get_total_valence(fragment: list) -> int:
    """
    Compute the total valence of a given fragment.
    This is used to evaluate if a fragment can be a molecular ion or not.

    Parameters
    ----------
    fragment : list
        a list of positive integers of length `len_frag`.

    Returns
    -------
    int
        Value of the total valence of the fragment.

    """
    total_valence = sum([valence[i]*fragment[i] for i in range(len(fragment))])
    return total_valence

def get_max_valence(fragment: list) -> int:
    """
    Compute the maximum valence of a given fragment,
    or the valence of the chemical element having maximum valence.
    This is used to evaluate if a fragment can be a molecular ion or not.

    Parameters
    ----------
    fragment : list
        a list of positive integers of length `len_frag`.

    Returns
    -------
    int
        Value of the maximum valence of the fragment.

    """
    max_valence = max([valence[i] for i in range(len(fragment)) if fragment[i]>0])
    return max_valence

def test_SENIOR_rule_ii(fragment: list) -> bool:
    """test the SENIOR theorem, item (ii),
    as described in Kind and Fiehn 2007:
        The sum of valences is greater than or equal to twice the
        maximum valence.
    This prevents fragments such as CFCl to be considered as
    valid molecular ion."""
    
    test_SENIOR_rule = sum([valence[i]*fragment[i] for i in range(len(fragment))]) >= 2 * max([valence[i] for i in range(len(fragment)) if fragment[i]>0])
    return test_SENIOR_rule
    
    
    
    
def get_fragment_from_string_formula(string: str) -> list:
    """The standardised string formula of a fragment"""
    # a code is a upper case letter + maybe a lower case letter + a number
    i = 0
    elements = []
    vec = [0]*len_frag
    while i < len(string):
        c = string[i]
        i += 1
        if c == '[':
            # get number
            c += string[i]
            i += 1
            while string[i].isnumeric():
                c += string[i]
                i += 1
            # get upper case letter
            c += string[i]
            i += 1
            while string[i].islower():
                c += string[i]
                i += 1
            # should be ']'
            assert string[i] == ']'
            c += string[i]
            i += 1
        else:
            while i < len(string) and string[i].islower():
                c += string[i]
                i += 1
        n = ""
        while i < len(string) and string[i].isnumeric():
            n += string[i]
            i += 1
        if len(n) > 0:
            I = int(n)
        else:
            I = 1
        elements.append((c,I))
        vec[dict_element_idx[c]] += I
    #print(elements)
    return vec

import re #for formula_to_mass

def read_formula_pattern(formula: str) -> list:
    """
    This function: 
        handles parenthesis if any, for example (CF3)3 becomes CF3CF3CF3, using "pattern",
        then identify atom + number, as many times as needed, using "pattern2".
        It basically cuts the string into element pieces.

    INPUT:
    - ``formula``: a string formula with possible isotopologues (with brakets)

    OUTPUT: a list of atoms
    """
    pattern = r"(\()([\[\]\w]*)(\))(\d*)"
    all_matches = re.findall(pattern, formula)

    while all_matches:
        for match in all_matches:
            count = match[3]
            text =""
            if (count == ""):
                count = 1
            else:
                count = int(match[3])
            while (count >= 1):
                text = text + match[1]
                count -= 1
            formula = formula.replace('(' + match[1] + ')' + match[3], text)
    
    #at this point, we end up with no parenthesis, but still brackets from isotopologues.
    #we cannot remove the brackets right now otherwise we end up mixing up isotope with numbers.
    #pattern2 = "\[?(\d*[A-Z]{1}[a-z]?)\]?(\d*)"
    pattern2 = "(\[?\d*[A-Z]{1}[a-z]?\]?)(\d*)" #do not strip the brackets here!
    all_matches2 = re.findall(pattern2, formula)
    
    return all_matches2



def nice_formula_from_frag_formula(text_formula):
    """
    

    Parameters
    ----------
    text_formula : TYPE string.
        format: chemical formula written from a fragment. No parenthesis. Isotopes may be present.

    Returns
    -------
    nice_text: tring. Formula written in a Latex compatible format,
    with subscript and super script as should be.
    Use this as text string as nice label on graph with Mathplotlib.

    """

    pattern = "(\[?)(\d{0,3})([A-Z]{1}[a-z]{0,1})(\]?)(\d*)"
    text = re.findall(pattern, text_formula)
    nice_text = ""
    for match in text:
        if match[0] == "[":
            #this is an isotope notation as [37Cl]
            #for match in all_matches:
            nice_text += match[0] + '$^{' + match[1] +  '}$' + match[2] + match[3]
        else:
            nice_text += match[2] #+ 
        if match[4] != "":
            nice_text += '$_{' + match[4] + '}$'
            
    text_formula = nice_text 
   
    return nice_text
    
    
    
    
def formula_to_mass(chem_formula: str) -> float:
    """return the exact, ionized mass (minus one electron)

    INPUT:
    - ``chem_formula``: a chemical formula (string) corresponding to the
      formula of a fragment

    OUTPUT: the total mass of the chemical formula minus one electron
    """
    all_matches = read_formula_pattern(chem_formula)
    molecularFormula =""
    mass = 0
    
    for match in all_matches:
        #Search symbol
        symbol = match[0]
        #Search numbers
        number = match[1]
        #print(symbol,number)
        if (number == ""):
            number = 1
        else:
            number = int(match[1])
        mass += float(dict_element_mass[symbol])*float(number)
        while (number >= 1):
            molecularFormula = molecularFormula + symbol
            number -= 1
    #print(molecularFormula)
    mass -= electron_mass
    #print("formula: " + formula +  "+, m/z = " + str(mass))
    return mass

def formula_to_frag(chem_formula: str) -> list:
    """return the fragment (a list of positive integers)

    INPUT:
    - ``chem_formula``: a chemical formula (string) corresponding to the
      formula of a fragment
    
    OUTPUT: the fragment as a list of integers
    """
    all_matches = read_formula_pattern(chem_formula)
    frag = [0]*len_frag
    for match in all_matches:
        #Search symbol
        symbol = match[0]
        #print(symbol)
        #Search numbers
        number = match[1]
        #print(symbol,number)
        if (number == ""):
            number = 1
        else:
            number = int(match[1])
        frag[dict_element_idx[symbol]] += number
    return frag

def formula_to_e_list(chem_formula: str) -> list:
    """return the list of present elements (atoms)

    INPUT:
    - ``chem_formula``: a chemical formula (string) corresponding to the
      formula of a fragment

    OUTPUT:
    This function is based on the read_formula_pattern function
    """
    e_list = []
    all_matches = read_formula_pattern(chem_formula)
    for match in all_matches:
        #Search symbol
        symbol = match[0]
        #print(symbol)
        #if symbol not in e_list: ????? TOCHECK MYRIAM
        if symbol != [e_list[i] for i in range(len(e_list))]:
            e_list.append(symbol)
    return e_list

# =============================================================================
# def count_multivalent(frag):
#     """This function is used only in `find_frag_v_count` which is not used
#     TODO Myriam remove the function if it is useless
#     """
#     multivalent_count = 0
#     for i_atom in range(len(frag)):
#         if valence[i_atom] > 1:
#             multivalent_count += frag[i_atom]
#     return multivalent_count
# 
# def find_frag_v_count(list_structure, multivalent_count):
#     """This function is not used anywhere
#     TODO Myriam remove the function if it is useless
#     """
#     list_frag = []
#     for i_fragment in range(len(list_structure)):
#         fragment = list_structure[i_fragment]
#         multivalent_atom = count_multivalent(i_fragment)
#         if multivalent_atom == multivalent_count: # and get_DBE_value(fragment) == 0.5
#             list_frag.append(fragment)
#     return list_frag
# =============================================================================
