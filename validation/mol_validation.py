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

"""
labels = ("CAS", "Compound", "Chemical name", "Boiling point", "Sum formula",  "SMILES code", "Aliases")
dic_labels = {"CAS":0, "Compound":1, "Chemical name":2, "Boiling point":3, "Sum formula":4, "SMILES code":5, "Aliases":6}

# SMILES code from
# https://pubchem.ncbi.nlm.nih.gov/
# http://www.chemspider.com

#Molecules to use for validation



mol_validation = [
    #CAS number, Compound,  Chemical name,                                 boiling point, Sum formula, SMILES code,                    aliases
    #All chemical names, SMILES codes, and boiling points are from ChemSpider: http://www.chemspider.com/
    ("7783-54-2", "NF3",        "nitrogen trifluoride",                       -128.7, "NF3",       "N(F)(F)F" ,                    ["NF3"]),
    ("75-73-0", "CF4",        "tetrafluoromethane",                           -128.0, "CF4",       "C(F)(F)(F)F",                  ["CF4"]),
    ("74-84-0", "C2H6",       "ethane",                                      -88.6, "C2H6",      "CC",                           ["C2H6", "ethane"]),
    ("76-16-4", "PFC-116",    "perfluoroethane",                               -78.1, "C2F6",      "C(C(F)(F)F)(F)(F)F",           ["PFC-116", "PFC116"]),
    ("2551-62-4", "SF6",        "sulphur hexafluoride",                        -63.8, "SF6",       "FS(F)(F)(F)(F)F",              ["SF6"]),
    ("75-72-9", "CFC-13",     "chlorotrifluoromethane",                        -81.4, "CF3Cl",     "C(F)(F)(F)Cl",                 ["CFC-13", "CFC13"]),
    ("463-58-1", "COS",        "carbonyl sulphide",                            -50.0, "COS",       "C(=O)=S",                      ["COS"]),
    ("75-63-8", "Halon-1301", "bromo(trifluoro)methane",                       -57.8, "CF3Br",     "C(F)(F)(F)Br",                 ["Halon1301", "Halon-1301", "halon-1301", "H-1301"]),
    ("76-19-7", "PFC-218",    "perfluoropropane",                              -36.7, "C3F8",      "C(C(F)(F)F)(C(F)(F)F)(F)F",    ["PFC218", "PFC-218"]),
    ("74-98-6", "C3H8",       "propane",                                       -42.0, "C3H8",      "CCC",                          ["C3H8", "Propane", "propane"]),
    ("2699-79-8", "SO2F2",      "sulphuryl difluoride",                        -55.4, "SO2F2",     "O=S(=O)(F)F",                  ["SO2F2"]),
    ("76-15-3", "CFC-115",    "1-chloro-1,1,2,2,2-pentafluoroethane",          -39.1, "C2F5Cl",    "C(C(F)(F)Cl)(F)(F)F",          ["CFC115", "CFC-115"]),
    ("373-80-8", "SF5CF3",     "pentafluoro(trifluoromethyl)sulfur",           -20.4, "SF5CF3",    "C(F)(F)(F)S(F)(F)(F)(F)F",     ["SF5CF3"]),
    ("75-71-8", "CFC-12",     "dichlorodifluoromethane",                       -29.8, "CF2Cl2",    "C(F)(F)(Cl)Cl",                ["CFC12", "CFC-12"]),
    ("75-45-6", "HCFC-22",    "chlorodifluoromethane",                         -40.8, "HCF2Cl",    "C(F)(F)Cl",                    ["HCFC22", "HCFC-22"]),
    ("115-25-3", "PFC-c318",   "octafluorocyclobutane",                         -6.0, "C4F8",      "C1(C(C(C1(F)F)(F)F)(F)F)(F)F", ["PFCc318", "PFC-c318", "PFC-318"]),
    ("74-87-3", "CH3Cl",      "chloromethane",                                 -24.2, "CH3Cl",     "CCl",                          ["CH3Cl"]),
    ("353-59-3", "Halon-1211", "bromochlorodifluoromethane",                    -2.5, "CF2ClBr",   "C(F)(F)(Cl)Br",                ["Halon1211", "Halon-1211", "halon-1211", "H-1211"]),
    ("76-14-2", "CFC-114",    "1,2-dichloro-1,1,2,2-tetrafluoroethane",          3.8, "C2F4Cl2",   "C(C(F)(F)Cl)(F)(F)Cl",         ["CFC114", "CFC-114"]),
    ("2837-89-0", "HCFC-124",   "2-chloro-1,1,1,2-tetrafluoroethane",          -12.0, "HC2F4Cl",   "C(C(F)(F)F)(F)Cl",             ["HCFC124", "HCFC-124"]),
    ("75-68-3", "HCFC-142b",  "1-chloro-1,1-difluoroethane",                    -9.2, "H3C2F2Cl",  "CC(F)(F)Cl",                   ["HCFC142b", "HCFC-142b"]),
    ("74-83-9", "CH3Br",      "bromomethane",                                    4.0, "CH3Br",     "CBr",                          ["CH3Br"]),
    ("75-69-4", "CFC-11",     "trichlorofluoromethane",                         23.8, "CFCl3",     "C(F)(Cl)(Cl)Cl",               ["CFC11", "CFC-11"]),
    ("74-88-4", "CH3I",       "iodomethane",                                    42.5, "CH3I",      "CI",                           ["CH3I"]),
    ("355-42-0", "C6F14",      "perfluorohexane",                               59.5, "C6F14",     "C(C(C(C(F)(F)F)(F)F)(F)F)(C(C(F)(F)F)(F)F)(F)F", ["C6F14", "perfluorohexane"]),
    ("75-09-2", "CH2Cl2",     "dichloromethane",                                39.6, "CH2Cl2",    "C(Cl)Cl",                      ["CH2Cl2"]),
    ("1717-00-6", "HCFC-141b",  "1,1-Dichloro-1-fluoroethane",                  32.0, "H3C2FCl2",  "CC(F)(Cl)Cl",                  ["HCFC141b", "HCFC-141b"]),
    ("76-13-1", "CFC-113",    "1,1,2-trichloro-1,2,2-trifluoroethane",          48.0, "C2F3Cl3",   "C(C(F)(Cl)Cl)(F)(F)Cl",        ["CFC113", "CFC-113"]),
    ("354-58-5", "CFC-113a", "1,1,1-Trichlorotrifluoroethane",                  45.7, "C2Cl3F3",   "C(C(Cl)(Cl)Cl)(F)(F)F", ["CFC-113a", "CFC113a"]),
    ("124-73-2", "Halon-2402", "1,2-dibromo-1,1,2,2-tetrafluoroethane",         47.3, "C2F4Br2",   "C(C(F)(F)Br)(F)(F)Br",         ["Halon2402", "Halon-2402", "halon-2402", "H-2402"]),
    ("67-66-3", "CHCl3",      "chloroform",                                     61.3, "CHCl3",     "C(Cl)(Cl)Cl",                  ["CHCl3", "chloroform"]),
    ("56-23-5", "CCl4",       "tetrachloromethane",                             76.7, "CCl4",      "C(Cl)(Cl)(Cl)Cl",              ["CCl4"]),
    ("79-01-6", "C2HCl3",     "1,1,2-trichloroethene",                          86.7, "HC2Cl3",    "C(=C(Cl)Cl)Cl",                ["C2HCl3", "TCE"]),
    ("71-55-6", "CH3CCl3",    "1,1,1-trichloroethane",                          74.0, "H3C2Cl3",   "CC(Cl)(Cl)Cl",                 ["CH3CCl3"]),
    ("74-95-3", "CH2Br2",     "dibromomethane",                                 97.0, "CH2Br2",    "C(Br)Br",                      ["CH2Br2"]),
    ("127-18-4", "C2Cl4",      "1,1,2,2-tetrachloroethene",                    121.0, "C2Cl4",     "C(=C(Cl)Cl)(Cl)Cl",            ["C2Cl4", "PCE"]),
    ("71-43-2", "benzene",       "benzene",                                     80.0, "C6H6",      "C1=CC=CC=C1",                  ["C6H6", "benzene"]),
    ("75-25-2", "CHBr3",      "bromoform",                                     149.0, "CHBr3",     "C(Br)(Br)Br",                  ["CHBr3", "bromoform"]),
    ("108-88-3", "toluene",       "toluene",                                   111.0, "C7H8",      "CC1=CC=CC=C1",                 ["C7H8", "toluene"]),
    ("593-53-3", "HFC-41",       "fluoromethane",                                     -78.5, "CH3F",      "CF",                           ["HFC-41", "HFC41"]),
    ("75-46-7", "HFC-23",       "fluoroform",                                        -82.0, "CHF3",      "C(F)(F)F",                     ["HFC-23", "HFC23"]),
    ("75-10-5", "HFC-32",       "difluoromethane",                                   -51.6, "CH2F2",     "C(F)F",                        ["HFC-32", "HFC32"]),
    ("420-46-2", "HFC-143a",     "1,1,1-trifluoroethane",                             -47.2, "C2H3F3",    "CC(F)(F)F",                    ["HFC-143a", "HFC143a"]),
    ("354-33-6", "HFC-125",      "1,1,1,2,2-pentafluoroethane",                       -48.1, "C2HF5",     "C(C(F)(F)F)(F)F",              ["HFC-125", "HFC125"]),
    ("811-97-2", "HFC-134a",     "1,1,1,2-tetrafluoroethane norflurane",              -26.5, "C2H2F4",    "C(C(F)(F)F)F",                 ["HFC-134a", "HFC134a"]),
    ("75-37-6", "HFC-152a",     "1,1-difluoroethane",                                -24.7, "C2H4F2",    "CC(F)F",                       ["HFC-152a", "HFC152a"]),
    ("359-35-3", "HFC-134",      "1,1,2,2-tetrafluoroethane",                         -19.7, "C2H2F4",    "C(C(F)F)(F)F",                 ["HFC-134", "HFC134"]),
    ("430-66-0", "HFC-143",      "1,1,2-trifluoroethane",                               5.0, "C2H3F3",    "C(C(F)F)F",                    ["HFC-143", "HFC143"]),
    ("431-89-0", "HFC-227ea",    "1,1,1,2,3,3,3-heptafluoropropane",                  -16.3, "C3HF7",     "C(C(F)(F)F)(C(F)(F)F)F",       ["HFC-227ea", "HFC227ea"]),
    ("624-72-6", "HFC-152",      "1,2-difluoroethane",                                 30.7, "C2H4F2",    "C(CF)F",                       ["HFC-152", "HFC152"]),
    ("677-56-5", "HFC-236cb",    "1,1,1,2,2,3-hexafluoropropane",                      -1.4, "C3H2F6",    "C(C(C(F)(F)F)(F)F)F",          ["HFC-236cb", "HFC236cb"]),
    ("690-39-1", "HFC-236fa",    "1,1,1,3,3,3-hexafluoropropane",                      -1.1, "C3H2F6",    "C(C(F)(F)F)C(F)(F)F",          ["HFC-236fa", "HFC236fa"]),
    ("431-63-0", "HFC-236ea",    "1,1,1,2,3,3-hexafluoropropane",                       6.5, "C3H2F6",    "C(C(F)F)(C(F)(F)F)F",          ["HFC-236ea", "HFC236ea"]),
    ("460-73-1", "HFC-245fa",    "1,1,1,3,3-pentafluoropropane",                       15.0, "C3H3F5",    "C(C(F)F)C(F)(F)F",             ["HFC-245fa", "HFC245fa"]),
    ("679-86-7", "HFC-245ca",    "1,1,2,2,3-pentafluoropropane",                       25.5, "C3H3F5",    "C(C(C(F)F)(F)F)F",             ["HFC-245ca", "HFC245ca"]),
    ("406-58-6", "HFC-365mfc",   "1,1,1,3,3-pentafluorobutane",                        40.0, "C4H5F5",    "CC(CC(F)(F)F)(F)F",            ["HFC-365mfc", "HFC365mfc"]),
    ("138495-42-8", "HFC-43-10mee", "1,1,1,2,2,3,4,5,5,5-decafluoropentane",              54.0, "C5H2F10",   "C(C(C(F)(F)F)F)(C(C(F)(F)F)(F)F)F", ["HFC-43-10mee", "HFC-4310mee", "HFC4310mee"]),
    ("754-12-1", "HFO-1234yf",     "2,3,3,3-tetrafluoroprop-1-ene",                   -28.3, "H2C3F4",    "C=C(C(F)(F)F)F",               ["HFO-1234yf", "HFO1234yf"]),
    ("29118-24-9", "HFO-1234ze(E)",  "(E)-1,3,3,3-tetrafluoroprop-1-ene",               -19.0, "H2C3F4",    "C(=CF)C(F)(F)F",               ["HFO-1234zeE", "HFO1234zeE", "HFO-1234ze(E)"]),
    ("102687-65-0", "HCFO-1233zd(E)", "(E)-1-chloro-3,3,3-trifluoro prop-1-ene",          21.0, "H2C3F3Cl",  "C(=CCl)C(F)(F)F",              ["HCFO-1233zdE", "HCFO1233zdE", "HCFO-1233zd(E)"]),
    ("375-45-1", "TCHFB", "tetrachlorohexafluorobutane",                       134.2, "C4Cl4F6",   "C(C(C(F)(F)Cl)(F)Cl)(C(F)(F)Cl)(F)Cl", ["TCHFB", "tetrachlorohexafluorobutane", "tetrachlorohexafluoro-butane"]),
    ("87-68-3", "HCBD", "1,1,2,3,4,4-hexachlorobuta-1,3-diene",                215.0, "C4Cl6",     "C(=C(Cl)Cl)(C(=C(Cl)Cl)Cl)Cl", ["hexachlorobutadiene", "HCBD"]),
]




