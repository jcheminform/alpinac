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
    (122, "PFDMCH", "Perfluorodimethylcyclohexane",                            102.0, "C8F16",     "C1(C(C(C(C(C1(F)F)(F)F)(F)F)(C(F)(F)F)F)(F)F)(C(F)(F)F)F", ["PFDMCH", "Perfluorodimethylcyclohexane", "Perfluorodimethyl-cyclohexane"]),
    (123, "perfluorodecaline", "1,1,2,2,3,3,4,4,4a,5,5,6,6,7,7,8,8,8a-octadecafluoronaphthalene", 142.0, "C10F18", "C12(C(C(C(C(C1(F)F)(F)F)(F)F)(F)F)(C(C(C(C2(F)F)(F)F)(F)F)(F)F)F)F", ["Perfluorodecaline", "perfluorodecaline", "(cis/trans)-Perfluorodecaline", "PFD"]),
    ("375-45-1", "TCHFB", "tetrachlorohexafluorobutane",                       134.2, "C4Cl4F6",   "C(C(C(F)(F)Cl)(F)Cl)(C(F)(F)Cl)(F)Cl", ["TCHFB", "tetrachlorohexafluorobutane", "tetrachlorohexafluoro-butane"]),
    (125, "DCTFP", "dichloro-trifluoro-pyridine",                              159.5, "C5Cl2F3N",  "c1(c(c(nc(c1Cl)F)F)Cl)F", ["DCTFP", "dichloro-trifluoro-pyridine"]),
    (126, "PFTBA", "Tris(perfluorobutyl)-amine",                               179.0, "C12F27N",   "C(C(C(F)(F)F)(F)F)(C(N(C(C(C(C(F)(F)F)(F)F)(F)F)(F)F)C(C(C(C(F)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F", ["PFTBA", "Tris(perfluorobutyl)-amine"]),
    (127, "PFPHP", "1,1,2,2,3,3,4,4,4a,4b,5,5,6,6,7,7,8,8,8a,9,9,10,10,10a-tetracosafluorophenanthrene", 215.0, "C14F24",    "C12(C3(C(C(C(C1(C(C(C(C2(F)F)(F)F)(F)F)(F)F)F)(F)F)(F)F)(C(C(C(C3(F)F)(F)F)(F)F)(F)F)F)F)F", ["PFPHP"]),
    ("87-68-3", "HCBD", "1,1,2,3,4,4-hexachlorobuta-1,3-diene",                215.0, "C4Cl6",     "C(=C(Cl)Cl)(C(=C(Cl)Cl)Cl)Cl", ["hexachlorobutadiene", "HCBD"]),
    (200, "CFC-112",     "1,1,2,2-tetrachloro-1,2-difluoro-ethane",             92.8, "C2F2Cl4",   "C(C(F)(Cl)Cl)(F)(Cl)Cl",       ["CFC-112", "CFC112"]),
    (201, "HCFC-133a",     "2-chloro-1,1,1-trifluoroethane",                     6.9, "H2C2F3Cl",  "C(C(F)(F)F)Cl",                ["HCFC-133a", "HCFC133a"]),
    (202, "C4F10", "1,1,1,2,2,3,3,4,4,4-decafluorobutane",                      -2.1, "C4F10",     "C(C(C(F)(F)F)(F)F)(C(F)(F)F)(F)F", ["C4F10"]),
    (203, "HCFC-132b", "1,2-dichloro-1,1-difluoroethane",                       46.8, "H2C2F2Cl2", "C(C(F)(F)Cl)Cl",               ["HCFC-132b", "HCFC132b"]),
    (204, "HCFC-132a", "1,1-dichloro-2,2-difluoroethane",                       60.0, "H2C2F2Cl2", "C(C(Cl)Cl)(F)F",               ["HCFC-132a", "HCFC132a"]),
    (205, "HCFC-132c", "1,1-dichloro-1,2-difluoroethane",                       45.0, "H2C2F2Cl2", "C(C(F)(Cl)Cl)F",               ["HCFC-132c", "HCFC132c"]),
    (206, "HCFC-132", "1,2-dichloro-1,1-difluoroethane",                        58.5, "H2C2F2Cl2", "C(C(F)Cl)(F)Cl",               ["HCFC-132", "HCFC132"]),
    (207, "1,1-Dichloroethane", "1,1-dichloroethane",                           57.0, "C2H4Cl2",   "CC(Cl)Cl",                     ["1-1-Dichloroethane", "1,1-Dichloroethane", "1,1-dichloroethane", "1-1-dichloroethane"]),
    (208, "1,2-Dichloroethane", "1,2-dichloroethane",                           83.0, "C2H4Cl2",   "C(CCl)Cl",                     ["1-2-Dichloroethane", "1,2-Dichloroethane", "1,2-dichloroethane", "1-2-dichloroethane"]),
    (209, "HCFO-1224ydz", "Z-1-chloro-2,3,3,3-tetrafluoropropene",               1.9, "C3HClF4",   "C(=[CH](C(F)(F)F)F)Cl",        ["HCFO-1224ydz", "HFO-1224ydz", "HCFO1224ydz"]),
    (210, "HCFC-235fa", "1-chloro-1,1,3,3,3-pentafluoropropane",                28.4, "C3H2ClF5",  "C(C(F)(F)F)C(F)(F)Cl",         ["HCFC-235fa", "HCFC235fa"]),
    (211, "ethyne", "ethyne",                                                  -84.7, "C2H2",      "C#C",                          ["C2H2", "ethyne", "Ethyne", "acetylene"]),
    (212, "ethene", "ethene",                                                 -103.8, "C2H4",      "C=C",                          ["ethylene", "ethene"]),
    (213, "c-propane", "cyclopropane"  ,                                       -33.0, "C3H6",      "C1CC1",                        ["c-propane", "cpropane", "cyclopropane"]),
    (214, "desflurane",  "2-(difluoromethoxy)-1,1,1,2-tetrafluoroethane",       23.5, "C3H2F6O",   "C(C(F)(F)F)(OC(F)F)F",         ["desflurane"]),
    (215, "isoflurane",  "2-chloro-2-(difluoromethoxy)-1,1,1-trifluoroethane",  48.5, "C3H2ClF5O", "C(C(F)(F)F)(OC(F)F)Cl",        ["isoflurane"]),
    (216, "sevoflurane", "1,1,1,3,3,3-hexafluoro-2-(fluoromethoxy)propane",     58.1, "C4H3F7O",   "C(OC(C(F)(F)F)C(F)(F)F)F",     ["sevoflurane"]),
    (217, "C2F5I", "1,1,1,2,2-pentafluoro-2-iodoethane",                        12.5, "C2F5I",     "C(C(F)(F)I)(F)(F)F",           ["C2F5I"]),
    (218, "HCFC-225ca", "3,3-dichloro-1,1,1,2,2-pentafluoropropane",            51.1, "HC3Cl2F5",  "C(C(C(F)(F)F)(F)F)(Cl)Cl",     ["HCFC225ca", "HCFC-225ca"]),
    (219, "HCFC-225cb", "1,3-dichloro-1,1,2,2,3-pentafluoropropane",            54.0, "HC3Cl2F5",  "C(C(C(F)(F)Cl)(F)F)(F)Cl",     ["HCFC225cb", "HCFC-225cb"]),
    (220, "CH2BrCl", "bromochloromethane",                                      68.0, "CH2BrCl",   "ClCBr",                        ["CH2BrCl", "bromochloromethane"]),
    (221, "CHBrCl2", "bromodichloromethane",                                    89.0, "CHBrCl2",   "C(Cl)(Cl)Br",                  ["CHBrCl2", "Bromodichloromethane"]),
    (222, "HCFC-21", "dichlorofluoromethane",                                    8.9, "CHFCl2",    "C(F)(Cl)Cl",                   ["HCFC-21", "dichlorofluoromethane"]),
    (223, "PFC-216", "1,1,2,2,3,3-hexafluorocyclopropane",                     -31.0, "C3F6",      "C1(C(C1(F)F)(F)F)(F)F",        ["PFC-216"]),
    (224, "hexafluoropropene", "hexafluoropropene",                            -29.0, "C3F6",      "C(=C(F)F)(C(F)(F)F)F",         ["hexafluoropropene"]),
    (225, "tetrafluoroethene", "tetrafluoroethene",                            -76.3, "C2F4",      "C(=C(F)F)(F)F",                ["tetrafluoroethene", "C2F4"]),
    (226, "fluoroethene", "fluoroethene",                                      -72.0, "H3C2F",     "C=CF",                         ["fluoroethene", "vinyl fluoride", "C2H3F"]),
    (227, "C2F3Cl", "chlorotrifluoroethene",                                   -27.0, "C2F3Cl",    "C(=C(F)Cl)(F)F",               ["C2F3Cl", "C2ClF3", "CFC-1113"]),
    (228, "i-C4F10", "1,1,1,2,3,3,3-heptafluoro-2-(trifluoromethyl)propane",     0.0, "C4F10",     "C(C(F)(F)F)(C(F)(F)F)(C(F)(F)F)F", ["i-C4F10", "perfluoroisobutane"]),
    (229, "H-1201", "bromodifluoromethane",                                    -14.5, "CHBrF2",    "C(F)(F)Br",                    ["H-1201", "halon-1201", "h-1201"]),
    (230, "HCFC-31", "chlorofluoromethane",                                     -9.1, "H2CFCl",    "FCCl",                         ["HCFC-31", "HCFC31"]),
    (231, "CS2", "carbon disulfide",                                            46.0, "CS2",       "C(=S)=S",                      ["CS2"]),
    (232, "CF3I", "iodotrifluoromethane",                                      -22.5, "CF3I",      "C(F)(F)(F)I",                  ["CF3I"]),
    (233, "H-1202", "dibromodifluoromethane",                                   24.5, "CBr2F2",    "", ["H-1202", "h1202", "h-1202", "H1202"]),
    (234, "HFO-1225yeZ", "(Z)-1,1,1,2,3-pentafluoropropene",                   -18.0, "C3HF5",     "", ["HFO-1225yeZ", "HFO1225yeZ"]),
    (235, "HFO-1225yeE", "(E)-1,2,3,3,3-Pentafluoropropene",                   -10.0, "C3HF5",     "", ["HFO-1225yeE", "HFO1225yeE"]),
    (236, "HFO-1336mzzE", "(2E)-1,1,1,4,4,4-Hexafluoro-2-butene",                8.5, "H2C4F6",    "", ["HFO-1336mzzE", "HFO1336mzzE"]),
    (237, "HFO-1234zeZ", "(Z)-1,3,3,3-tetrafluoropropene",                      11.5, "C3H2F4",    "", ["HFO-1234zeZ", "HFO1234zeZ"]),
    (238, "HFO-1224ydZ", "Z-1-Chloro-2,3,3,3-tetrafluoropropene",                1.9, "C3HClF4",   "", ["HFO-1224ydZ", "HFO1224ydZ", "HCFO-1224ydZ", "HCFO1224ydZ"]),
    (239, "HFC-356mff", "1,1,1,4,4,4-Hexafluorobutane",                         24.6, "C4H4F6",    "", ["HFC-356mff", "HFC356mff"]),
    (240, "HCFC-123", "2,2-Dichloro-1,1,1-trifluoroethane",                     27.5, "C2HCl2F3",  "", ["HCFC-123", "HCFC123"]),
    (241, "HBFO-1233xf", "2-Bromo-3,3,3-trifluoropropene",                      33.5, "C3H2BrF3",  "", ["HBFO-1233xf", "HBFO1233xf", "HBFO-1233xfB"]),
    (242, "HFO-1336mzzZ", "(Z)-1,1,1,4,4,4-hexafluoro-2-butene",                33.0, "H2C4F6",    "", ["HFO-1336mzzZ", "HFO1336mzzZ"]),
    (243, "H-2311", "2-bromo-2-chloro-1,1,1-trifluoroethane",                   50.0, "C2HBrClF3", "", ["H-2311", "H2311", "h-2311", "h2311", "halothane"]),
    (244, "HCFC-122", "1,2,2-trichloro-1,1-difluoro-ethane",                    72.0, "C2HCl3F2",  "", ["HCFC-122", "HCFC122"]),
    (245, "HCFC-122a", "1,1,2-trichloro-1,2-difluoroethane",                    72.0, "C2HCl3F2",  "", ["HCFC-122a", "HCFC122a"]),
    (246, "CFC-112a", "1,1-Difluoro-1,2,2,2-tetrachloroethane",                 91.5, "C2F2Cl4",   "", ["CFC-112a", "CFC112a"]),
    (247, "HCFC-121", "1,1,2,2-tetrachloro-1-fluoro-ethane",                   116.6, "HC2FCl4",   "", ["HCFC-121", "HCFC121"]),
    (248, "HCFC-131", "1,1,2-Trichloro-2-fluoroethane",                        102.5, "C2H2Cl3F",  "", ["HCFC-131", "HCFC131"]),
    (249, "PCBTF", "1-chloro-4-(trifluoromethyl)benzene",                      128.5, "C7H4ClF3",  "", ["PCBTF"]),
    (250, "CHBr2Cl", "Dibromochloromethane",                                   119.0, "CHBr2Cl",   "", ["CHBr2Cl"]),
    (251, "propene","1-propene",                                               -47.7, "C3H6",      "CC=C", ["propene","1-propene", "propylene"]),
    (252, "i-butane","2-methylpropane",                                        -11.7, "C4H10",     "CC(C)C", ["i-butane", "isobutane", "2-methylpropane"]),
    (253, "n-butane","n-butane",                                                -0.5, "C4H10",     "CCCC", ["n-butane", "butane"]),
    (254, "1-butene","1-butene",                                                -6.3, "C4H8",      "", ["1-butene"]),
    (255, "neopentane","2,2-dimethyl-propane",                                   9.5, "C5H12",     "CC(C)(C)C", ["neopentane"]),
    (256, "propyne","propyne",                                                 -23.2, "C3H4",      "", ["propyne"]),
    (257, "2-buteneE","2-butene(E)",                                             1.0, "C4H8",      "", ["2-buteneE", "2buteneE"]),
    (258, "1,3-butadiene","1,3-butadiene",                                      -5.0, "C4H6",      "", ["1-3-butadiene", "1,3-butadiene"]),
    (259, "2-buteneZ","2-buteneZ",                                               3.7, "C4H8",      "", ["2-buteneZ", "2buteneZ"]),
    (260, "ethylbenzene","ethylbenzene",                                       136.0, "C8H10",     "", ["ethylbenz", "ethylbenzene"]),
    (261, "CFC-114a", "1,1,1,2-Tetrafluoro-2,2-dichloroethane",                  4.0, "CF3CCl2F",  "", ["CFC-114a", "CFC114a"]),
    (262, "HFC-161", "Fluoroethane",                                           -37.1, "C2H5F",     "", ["HFC-161", "HFC161"]),
    (263, "NOVEC-4710", "2,3,3,3-Tetrafluoro-2-(trifluoromethyl)propanenitrile",26.1, "C4F7N",     "C(#N)C(C(F)(F)F)(C(F)(F)F)F", ["NOVEC-4710", "NOVEC4710"]),
    (264, "C5F12", "dodecafluoropentane",                                       29.0, "C5F12",     "", ["C5F12", "dodecafluoropentane", "perfluoropentane"]),
    (265, "C-5-ketone", "1,1,1,3,4,4,4-Heptafluoro-3-(trifluoromethyl)-2-butanone", 24.0, "C5F10O", "", ["C-5-Ketone", "C5Ketone", "C-5-ketone", "C5ketone"]),
    (266, "PFNMM", "2,2,3,3,5,5,6,6-Octafluoro-4-(trifluoromethyl)morpholine",  50.0, "C5F11NO",   "", ["PFNMM"]),
    (267, "HCFC-151a", "1-Chloro-1-fluoroethane",                               16.1, "C2H4ClF",   "", ["HCFC-151a", "HCFC151a"]),
    (268, "i-C6F14pent", "1,1,1,2,2,3,3,4,5,5,5-Undecafluoro-4-(trifluoromethyl)pentane", 58.2, "C6F14", "", ["i-C6F14pent", "iC6F14pent"]),
    (269, "i-C6F14but", "1,1,1,2,3,4,4,4-Octafluoro-2,3-bis(trifluoromethyl)butane", 58.0, "C6F14", "", ["i-C6F14but", "iC6F14but"]),
    (271, "PFTEA", "1,1,2,2,2-Pentafluoro-N,N-bis(pentafluoroethyl)ethanamine", 68.5, "C6F15N",    "C(C(F)(F)F)(N(C(C(F)(F)F)(F)F)C(C(F)(F)F)(F)F)(F)F", ["PFTEA"]),
    (272, "hexafluorobenzene", "hexafluorobenzene",                             80.5, "C6F6",      "c1(c(c(c(c(c1F)F)F)F)F)F", ["hexafluorobenzene"]),
    (273, "perfluorodimethyl-cyclohexane", "1,1,2,2,3,3,4,5,5,6-decafluoro-4,6-bis(trifluoromethyl)cyclohexane", 102.0, "C8F16", "C1(C(C(C(C(C1(F)F)(F)F)(F)F)(C(F)(F)F)F)(F)F)(C(F)(F)F)F", ["perfluorodimethylcyclohexane", "perfluorodimethyl-cyclohexane"]),
    (274, "enflurane", "2-Chloro-1-(difluoromethoxy)-1,1,2-trifluoroethane",    56.8, "C3H2ClF5O", "C(C(OC(F)F)(F)F)(F)Cl", ["enflurane"]),
    (275, "perfluorodecaline", "1,1,2,2,3,3,4,4,4a,5,5,6,6,7,7,8,8,8a-octadecafluoronaphthalene", 142.0, "C10F18", "C12(C(C(C(C(C1(F)F)(F)F)(F)F)(F)F)(C(C(C(C2(F)F)(F)F)(F)F)(F)F)F)F", ["perfluorodecaline"]),
    (276, "chlorobenzene", "chlorobenzene",                                    131.7, "C6H5Cl",    "c1ccc(cc1)Cl", ["chlorobenzene"]),
    (277, "tetrachlorohexafluoro-butane", "1,2,3,4‐Tetrachlorohexafluorobutane", 134.2, "C4Cl4F6", "C(C(C(F)(F)Cl)(F)Cl)(C(F)(F)Cl)(F)Cl", ["TCHFB", "tetrachlorohexafluoro-butane", "tetrachlorohexafluorobutane"]),
    (278, "MrEfflamm", "1,1,2-trichloroethane",                               112.5, "C2H3Cl3",   "C(C(Cl)Cl)Cl", ["112Trichloroethane", "1,1,2-Trichloroethane", "1,1,2-trichloroethane", "MrEfflamm"]),
    ("suspect", "MsAlwena", "1,2-dichloropropane",                                   95.0, "H6C3Cl2",   "CC(CCl)Cl", ["12dichloropropane", "1,2-dichloropropane", "MsAlwena"]),
    (280, "MrBrieuc", "bromoethane",                                           38.5, "C2H5Br",    "CCBr", ["bromoethane", "MrBrieuc", "C2H5Br"]),
    (281, "MrBombur", "2-Bromo-1,1-difluoroethene",                             6.1, "HC2F2Br",   "C(=C(F)F)Br", ["MrBombur"]),
    (282, "MrBernina", "Vinyl bromide",                                        16.0, "C2H3Br",    "C=CBr", ["MrBernina"]),
    (283, "MrBofur", "Bromofluoromethane",                                     17.5, "H2CFBr",    "C(F)Br", ["MrBofur"]),
    (284, "MrGloin", "Mr_Gloin",                                                5.5, "C3F5Cl",    "", ["MrGloin"]),
    (285, "MsDominique", "1-chloro-1-fluoroethene",                           -24.0, "H2C2FCl",   "C=C(F)Cl", ["MsDominique"]),
    (286, "MrAofred", "(Z)1-chloro-2-fluoroethene",                            15.5, "H2C2FCl",   "C(=C\Cl)\F", ["MrAofred"]),
    ("suspect", "Mr_Poschiavo", "Pentafluoroethyliodide",                             12.5, "C2F5I",     "C(C(F)(F)I)(F)(F)F", ["Mr_Poschiavo", "C2F5I"]),
    (288, "MrBifur", "1-Bromo-1-chloro-2,2-difluoroethene",                    41.0, "C2F2ClBr",  "C(=C(Cl)Br)(F)F", ["MrBifur"]),
    (289, "Mr_Brusio", "1,1-Dibromo-1,2,2,2-tetrafluoroethane",                 47.0, "C2F4Br2",   "C(C(F)(Br)Br)(F)(F)F", ["Mr_Brusio"]),
    (290, "o-xylene", "o-xylene",                                              144.0, "C8H10",     "Cc1ccccc1C", ["o-xylene"]),
    (291, "Mr_Fili", "1,1-Difluoroethene",                                     -83.0, "H2C2F2",    "C=C(F)F", ["Mr_Fili", "1,1-Difluoroethene", "1,1-difluoroethene"]),
    (292, "Mr_Kili", "trifluoroethene",                                        -51.0, "HC2F3",     "C(=C(F)F)F", ["Mr_Kili", "trifluoroethene"]),
    (293, "2-methylbutane", "2-methylbutane",                                   30.0, "C5H12",     "CCC(C)C", ["2-methylbutane", "2methylbutane", "i-pentane"]),
    (294, "n-pentane", "n-pentane",                                             36.0, "C5H12",     "CCCCC", ["n-pentane", "pentane"]),
    (295, "c-pentane", "c-pentane",                                             50.0, "C5H10",     "C1CCCC1", ["c-pentane", "cyclopentane"]),
    (296, "3-methyl-1-butene", "3-methyl-1-butene",                             20.0, "C5H10",     "CC(C)C=C", ["3-methyl-1-butene"]),
    (297, "1-pentene", "1-pentene",                                             30.0, "C5H10",     "CCCC=C", ["1-pentene", "1pentene"]),
    (298, "2-penteneE", "2-penteneE",                                           37.0, "C5H10",     "CC/C=C/C", ["2-penteneE", "2-pentene(E)"]),
    (299, "acetone", "acetone",                                                 56.0, "C3H6O",     "CC(=O)C", ["acetone"]),
    (300, "2,2-dimethylbutane", "2,2-dimethylbutane",                           50.0, "C6H14",     "CCC(C)(C)C", ["2-2-dimethylbutane", "2,2-dimethylbutane"]),
    (301, "furan", "furan",                                                     32.0, "C4H4O",     "c1ccoc1", ["furan", "furane"]),
    (301, "alpha-methylstyrene", "alpha-methylstyrene",                        167.0, "C9H10",     "CC(=C)c1ccccc1", ["alpha-methylstyrene-c"]),
    (302, "2-3-dimethylbutane", "2-3-dimethylbutane",                           58.0, "C6H14",     "CC(C)C(C)C", ["2-3-dimethylbutane", "2,3-dimethylbutane"]),
    (303, "isoprene", "isoprene",                                               34.0, "C5H8",      "CC(=C)C=C", ["isoprene"]),
    (304, "2-methylpentane", "2-methylpentane",                                 62.0, "C6H14",     "CCCC(C)C", ["2-methylpentane", "2methylpentane"]),
    (305, "3-methylpentane", "3-methylpentane",                                 62.0, "C6H14",     "CCC(C)CC", ["3-methylpentane", "3methylpentane"]),
    (306, "n-hexane", "n-hexane",                                               69.0, "C6H14",     "CCCCCC", ["n-hexane", "hexane"]),
    (307, "c-hexane", "c-hexane",                                               80.7, "C6H12",     "C1CCCCC1", ["c-hexane", "cyclohexane"]),
    (308, "n-heptane", "n-heptane",                                             98.0, "C7H16",     "CCCCCCC", ["n-heptane", "heptane"]),
    (309, "n-PB", "1-bromopropane",                                             71.0, "C3H7Br",    "CCCBr", ["n-PB", "1-bromopropane"]),
    (310, "i-PB", "2-bromopropane",                                             59.0, "C3H7Br",    "CC(C)Br", ["i-PB", "2-bromopropane"]),
    (311, "2,2,4-trimethylpentane", "2,2,4-trimethylpentane",                   98.5, "C8H18",     "CC(C)CC(C)(C)C", ["224trimethylpentane", "2,2,4-trimethylpentane", "2-2-4-trimethylpentane"]),
    (312, "octane", "octane",                                                  126.0, "C8H18",     "CCCCCCCC", ["octane", "n-octane"]),
    (313, "dichloro-trifluoro-pyridine", "3,5-Dichloro-2,4,6-trifluoropyridine", 159.5, "C5Cl2F3N", "c1(c(c(nc(c1Cl)F)F)Cl)F", ["DCTFP", "dichloro-trifluoro-pyridine", "dichlorotrifluoropyridine"]),
    (314, "mp-xylene", "mp-xylene",                                            138.5, "C8H10",     "Cc1ccc(cc1)C", ["mp-xylene", "m+p-xylene"]),
    (315, "m-xylene", "m-xylene",                                              139.0, "C8H10",     "Cc1cccc(c1)C", ["m-xylene"]),
    (316, "p-xylene", "p-xylene",                                              138.0, "C8H10",     "Cc1ccc(cc1)C", ["p-xylene"]),
    (317, "perfluorotributylamine", "perfluorotributylamine",                  179.0, "C12F27N",   "C(C(C(F)(F)F)(F)F)(C(N(C(C(C(C(F)(F)F)(F)F)(F)F)(F)F)C(C(C(C(F)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F", ["PFTBA", "perfluorotributylamine"]),
    (318, "H-1201-suspect", "Bromodifluoromethane",                            -14.5, "CHBrF2",    "C(F)(F)Br", ["H-1201-suspect", "H-1201", "CHBrF2"]),
    (319, "Mr_Oin", "chloroethene",                                            -13.4, "C2H3Cl",    "C=CCl", ["Mr_Oin", "chloroethene"]),
    (320, "Ms_Lizenn", "chloroethane",                                          12.3, "H5C2Cl",    "CCCl", ["Ms_Lizenn", "chloroethane"]),
    (321, "c-C4F8O", "Perfluorotetrahydrofuran",                                 0.6, "C4F8O",     "C1(C(C(OC1(F)F)(F)F)(F)F)(F)F", ["c-C4F8O", "C4F8O"]),
    (322, "iodoethane", "iodoethane",                                           71.5, "C2H5I",     "CCI", ["iodoethane", "iodoethane_suspect", "C2H5I"]),
    (323, "Mr_Thorin", "1-Chloro-1,1,2,2,3,3,3-heptafluoropropane",             -2.0, "C3F7Cl",    "C(C(F)(F)F)(C(F)(F)Cl)(F)F", ["Mr_Thorin"]),
    (324, "H-2301-suspect", "2-Bromo-1,1,1-trifluoroethane",                    26.0, "CH2BrCF3",  "C(C(F)(F)F)Br", ["H-2301-suspect"]),
    (325, "Mr_Matthias", "Nonafluoropentanoyl fluoride",                        52.0, "C5F10O",    "C(=O)(C(C(C(C(F)(F)F)(F)F)(F)F)(F)F)F", ["Mr_Matthias"]),
    (236, "MsLuthien", "1,1,3,3,3-Pentafluoro-1-propene",                      -21.0, "C3HF5",     "C(=C(F)F)C(F)(F)F", ["MsLuthien", "HFO-1225zc", "HFO1225zc"]),
    (237, "MrDain", "1,1-dichloro-2,2-difluoroethene",                          19.0, "C2F2Cl2",   "C(=C(Cl)Cl)(F)F", ["MrDain", "Mr_Dain"]),
    (238, "MrDori", "1,2-Dichloro-1,2-difluoroethene",                          21.5, "C2F2Cl2",   "", ["MrDori", "Mr_Dori"]),
    (239, "MsMyriam", "1,1-dichloroethene",                                     31.0, "H2C2Cl2",   "C=C(Cl)Cl", ["MsMyriam", "Ms_Myriam", "1,1-Dichloroethene", "1,1-dichloroethene", "1-1-Dichloroethene"]),
    (340, "MsSofia", "(E)1,2-Dichloroethylene",                                 48.0, "H2C2Cl2",   "C(=C/Cl)\Cl", ["MsSofia", "Ms_Sofia", "trans-1,2-Dichloroethylene", "(E)1,2-Dichloroethylene", "(E)1,2-dichloroethylene"]),
    (241, "MrCleden", "(Z)1,2-Dichloroethylene",                                60.0, "H2C2Cl2",   "C(=C\Cl)\Cl", ["MrCleden", "Mr_Cleden", "(Z)1,2-Dichloroethylene", "(Z)1,2-dichloroethylene"])
]



#
#propene
#C3H4
#
#i-butane
#
#C4F8O
#
#HFO-1225yeZ
#n-butane
#HFO-1225yeE
#
#
#CFC-114a
#HFC-161
#NOVEC-4710
#C5F12
#neopentane
#propyne
#2-buteneE
#1-3-butadiene
