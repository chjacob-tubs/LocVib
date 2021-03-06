# -*- coding: utf-8 -*-
#
# This file is part of the
# LocVib 1.1 suite of tools for the analysis for vibrational spectra.
# Copyright (C) 2009-2018 by Christoph R. Jacob and others.
#
#    LocVib is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    LocVib is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with LocVib.  If not, see <http://www.gnu.org/licenses/>.
#
# In scientific publications using the LocVib tools, please cite:
#   Ch. R. Jacob, J. Chem. Phys 130 (2009), 084106.
# 
# The most recent version of LocVib is available at
#   http://www.christophjacob.eu/software
"""
 Constants used in program.
"""

pi = 3.141592653589793

# speed of light in atomic units
cvel = 137.0359895
cvel_ms = 2.99792458e08

# conversion from bohr to Angstrom
Bohr_in_Angstrom = 0.5291772108
Bohr_in_Meter = Bohr_in_Angstrom * 1.0e-10

Avogadro = 6.02214199e23

amu_in_kg = 1.0e-3/Avogadro

Hartree_in_Joule = 4.35974381e-18
eV_in_Joule = 1.6021765654e-19


au_in_Debye =  2.54177
Debye_in_Cm = 3.33564e-30

epsilon0 = 8.854187817e-12  # in SI units
h_SI = 6.62606957e-34 # in SI units

me_in_amu = 5.4857990943e-4   # mass of electron in amu
atu_in_s = 2.41888432650516e-17  # atomic time unit in seconds

cm_in_au = atu_in_s * (2.0*pi*1e2*cvel_ms)   # cm-1 -> au

intfactor = 2.5066413842056297  # factor to calculate integral absorption coefficient having freq in [cm-1] and dipole moment in [Debye]
