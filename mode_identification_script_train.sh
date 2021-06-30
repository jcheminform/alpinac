#!/bin/bash

# Initial Software, Myriam Guillevic and Aurore Guillevic,
# Copyright Empa and Inria, 2019 - 2021.

# This file is part of ALPINAC.

# ALPINAC is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# ALPINAC is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.

# You should have received a copy of the GNU Lesser General Public License
# along with ALPINAC.  If not, see <https://www.gnu.org/licenses/>.


PREFIX=data/nontarget_screening/fragments
mkdir --parents data/nontarget_screening/formulas/
mkdir --parents data/nontarget_screening/mass_spectra/

for FILE in `ls $PREFIX/*train*.txt`
do
    python3 mode_identification.py $FILE
done

exit 0

# the above loop corresponds to:

python3 mode_identification.py $PREFIX/200602.1550.std.7.frag.train_C2Cl4.txt
python3 mode_identification.py $PREFIX/200602.1550.std.7.frag.train_C2H6.txt
python3 mode_identification.py $PREFIX/200602.1550.std.7.frag.train_C2HCl3.txt # --target-elements CHCl --target-formula C2HCl3 # examples of options
python3 mode_identification.py $PREFIX/200602.1550.std.7.frag.train_C3H8.txt
python3 mode_identification.py $PREFIX/200602.1550.std.7.frag.train_C6F14.txt
python3 mode_identification.py $PREFIX/200602.1550.std.7.frag.train_C6H6.txt
python3 mode_identification.py $PREFIX/200602.1550.std.7.frag.train_C7H8.txt
python3 mode_identification.py $PREFIX/200602.1550.std.7.frag.train_CCl4.txt --verbose-knapsack
python3 mode_identification.py $PREFIX/200602.1550.std.7.frag.train_CFC113.txt
python3 mode_identification.py $PREFIX/200602.1550.std.7.frag.train_CFC114.txt
python3 mode_identification.py $PREFIX/200602.1550.std.7.frag.train_CFC115.txt
python3 mode_identification.py $PREFIX/200602.1550.std.7.frag.train_CFC11.txt
python3 mode_identification.py $PREFIX/200602.1550.std.7.frag.train_CFC12.txt
python3 mode_identification.py $PREFIX/200602.1550.std.7.frag.train_CFC13.txt
python3 mode_identification.py $PREFIX/200602.1550.std.7.frag.train_CH2Br2.txt
python3 mode_identification.py $PREFIX/200602.1550.std.7.frag.train_CH2Cl2.txt
python3 mode_identification.py $PREFIX/200602.1550.std.7.frag.train_CH3Br.txt
python3 mode_identification.py $PREFIX/200602.1550.std.7.frag.train_CH3Cl.txt
python3 mode_identification.py $PREFIX/200602.1550.std.7.frag.train_CH3I.txt
python3 mode_identification.py $PREFIX/200602.1550.std.7.frag.train_CHCl3.txt
python3 mode_identification.py $PREFIX/200602.1550.std.7.frag.train_COS.txt
python3 mode_identification.py $PREFIX/200602.1550.std.7.frag.train_Halon1211.txt
python3 mode_identification.py $PREFIX/200602.1550.std.7.frag.train_Halon1301.txt
python3 mode_identification.py $PREFIX/200602.1550.std.7.frag.train_Halon2402.txt
python3 mode_identification.py $PREFIX/200602.1550.std.7.frag.train_HCFC124.txt
python3 mode_identification.py $PREFIX/200602.1550.std.7.frag.train_HCFC141b.txt
python3 mode_identification.py $PREFIX/200602.1550.std.7.frag.train_HCFC142b.txt
python3 mode_identification.py $PREFIX/200602.1550.std.7.frag.train_HCFC22.txt
python3 mode_identification.py $PREFIX/200602.1550.std.7.frag.train_PFC116.txt
python3 mode_identification.py $PREFIX/200602.1550.std.7.frag.train_PFC218.txt
python3 mode_identification.py $PREFIX/200602.1550.std.7.frag.train_PFCc318.txt
python3 mode_identification.py $PREFIX/200602.1550.std.7.frag.train_SF5CF3.txt
python3 mode_identification.py $PREFIX/200602.1550.std.7.frag.train_SF6.txt
python3 mode_identification.py $PREFIX/200602.1550.std.7.frag.train_SO2F2.txt
python3 mode_identification.py $PREFIX/200714.2200.cal6L.7.frag.train_CF4.txt
python3 mode_identification.py $PREFIX/200714.2200.cal6L.7.frag.train_NF3.txt
