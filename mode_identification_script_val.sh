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

for FILE in `ls $PREFIX/*val*.txt`
do
    python3 mode_identification.py $FILE
done

exit 0

# the above loop corresponds to:

python3 mode_identification.py $PREFIX/190424.0810.tank.3.frag.val_HFC134a.txt
python3 mode_identification.py $PREFIX/190424.0810.tank.3.frag.val_HFC143a.txt
python3 mode_identification.py $PREFIX/190424.0810.tank.3.frag.val_HFC152.txt
python3 mode_identification.py $PREFIX/190424.0810.tank.3.frag.val_HFC227ea.txt
python3 mode_identification.py $PREFIX/190424.0810.tank.3.frag.val_HFC236cb.txt
python3 mode_identification.py $PREFIX/190424.0810.tank.3.frag.val_HFC236ea.txt
python3 mode_identification.py $PREFIX/190424.0810.tank.3.frag.val_HFC245fa.txt
python3 mode_identification.py $PREFIX/190424.0810.tank.3.frag.val_HFC41.txt
python3 mode_identification.py $PREFIX/190424.0810.tank.3.frag.val_HFC4310mee.txt
python3 mode_identification.py $PREFIX/190425.1731.tank.3.frag.val_HFC125.txt
python3 mode_identification.py $PREFIX/190425.1731.tank.3.frag.val_HFC134.txt
python3 mode_identification.py $PREFIX/190425.1731.tank.3.frag.val_HFC143.txt
python3 mode_identification.py $PREFIX/190425.1731.tank.3.frag.val_HFC152a.txt
python3 mode_identification.py $PREFIX/190425.1731.tank.3.frag.val_HFC236fa.txt
python3 mode_identification.py $PREFIX/190425.1731.tank.3.frag.val_HFC23.txt
python3 mode_identification.py $PREFIX/190425.1731.tank.3.frag.val_HFC245ca.txt
python3 mode_identification.py $PREFIX/190425.1731.tank.3.frag.val_HFC32.txt
python3 mode_identification.py $PREFIX/190425.1731.tank.3.frag.val_HFC365mfc.txt
python3 mode_identification.py $PREFIX/200602.1550.std.7.frag.val_HCFO1233zdE.txt
python3 mode_identification.py $PREFIX/200602.1550.std.7.frag.val_HFO1234yf.txt
python3 mode_identification.py $PREFIX/200602.1550.std.7.frag.val_HFO1234zeE.txt
