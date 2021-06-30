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
class RunType():
    _run_type_tofwerk_tof = 0
    _run_type_NIST = 1
    _run_type_no_mass_cal = 2
    _run_type_unit_mass = 3
    _str_ = ["tofwerk_tof", "NIST", "no_mass_cal"]
    def __init__(self):
        self._run_type = None
    def set_tofwerk_tof(self):
        self._run_type = self._run_type_tofwerk_tof
    def set_NIST(self):
        self._run_type = self._run_type_NIST
    def set_no_mass_cal(self):
        self._run_type = self._run_type_no_mass_cal
    def set_unit_mass(self):
        self._run_type = self._run_type_unit_mass

    def is_tofwerk_tof(self):
        return self._run_type == self._run_type_tofwerk_tof
    def is_NIST(self):
        return self._run_type == self._run_type_NIST
    def is_no_mass_cal(self):
        return self._run_type == self._run_type_no_mass_cal
    def is_unit_mass(self):
        return self._run_type == self._run_type_unit_mass

    def __repr__(self):
        if self._run_type is None:
            return "Undefined"
        return self._str_[self._run_type]
