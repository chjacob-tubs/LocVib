# -*- coding: utf-8 -*-
#
# This file is part of the
# LocVib 1.3 suite of tools for the analysis for vibrational spectra.
# Copyright (C) 2009-2023 by Christoph R. Jacob and others.
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
 Localizing normal modes
"""
import math
import numpy
import copy

from .Modes import VibModes
from . import Constants

class LocVib (object) :
    """
    Local Modes (LocVib) class.

    :Attributes: = Parameters.

    Parameters
    ----------  
    startmodes : VibTools.Modes class
        modes to be localized.  
    nmodes : int
        number of modes.
    natoms : int
        number of atoms in molecule.
    loctype : String
        Localization criterion: 
        PM (Pipek and Mezey); B (Boys)
    coords : ndarray 
        coordinates 
    transmat : ndarray (nmodes x nmodes)
        (localization) transformation matrix.
    """

    def __init__ (self, modes, loctype='PM') :
        """
        LocVib constructor.
        For more details see the class description/docstring
        """
        self.startmodes = modes 
        # TODO: name "modes" is used twice - risk of confusion -> calc_p,calc_ab
        # here: modes = VibTools.Modes class
        # elsewhere: modes = a set of normal modes
        self.nmodes = self.startmodes.nmodes
        self.natoms = self.startmodes.natoms
        self.loctype = loctype

        if self.loctype.startswith('B') :
            self.coords = self.startmodes.mol.coordinates

        self.locmodes = copy.copy(self.startmodes)
        self.transmat = numpy.identity(self.nmodes)
        # transmat is the transpose of the matrix U in the paper

        self.subsets = None

    def calc_p (self, modes) :
        """
        calculates the localization measure chi.

        Paper:
        *Ch. R. Jacob, J. Chem. Phys 130 (2009), 084106*

        Parameters
        ----------
        modes : ndarray (3*natoms-6,3*natoms)
            normal modes.

        Returns
        -------
        p : tba ?
        tba ?
        """
        squared_modes = modes**2
        # LocVib paper eq. 6
        c = squared_modes[:,0::3] + squared_modes[:,1::3] + squared_modes[:,2::3]
                
        if self.loctype == 'PM' : # Pipek + Mezey criterion
            p = numpy.linalg.norm(c)**2 # LocVibpaper eq. 5
        elif self.loctype.startswith('B') : # Boys criterion
            p_mode = numpy.zeros((3,self.nmodes))# LocVib paper eq. 9
            for idir in range(3):
                p_mode[idir] = (self.coords[:,idir] * c).sum(axis=1)
            p = (p_mode**2).sum()

        return p

    def calc_ab (self, modes, i, j):
        """
        TASK ??
        
        Parameters
        ----------
        modes : ndarray (3*natoms-6,3*natoms)
            normal modes.
        i,j : int, int
            index of modes array.

        Returns
        -------
        a,b : tba ?
        tba ?
        """
        # Pipek Mezey Ref-Paper? Eq.14 or 29???
        squared_modes = modes**2

        ci = squared_modes[i, 0::3] + squared_modes[i, 1::3] + squared_modes[i, 2::3]
        cj = squared_modes[j, 0::3] + squared_modes[j, 1::3] + squared_modes[j, 2::3]
        cij = modes[i, 0::3] * modes[j, 0::3] + modes[i, 1::3] * modes[j, 1::3] + modes[i, 2::3] * modes[j, 2::3]

        a = 0.0
        b = 0.0
            
        if self.loctype.startswith('PM') :
            a = (cij**2 - 0.25*(ci - cj)**2).sum()
            b = (cij*(ci - cj)).sum()

        elif self.loctype.startswith('B') :
            pi  = numpy.zeros((3,))
            pj  = numpy.zeros((3,))
            pij = numpy.zeros((3,))

            for idir in range(3):
                pi[idir] = (self.coords[:,idir] * ci).sum()
                pj[idir] = (self.coords[:,idir] * cj).sum()
                pij[idir] = (self.coords[:,idir] * cij).sum()

            a = numpy.dot(pij,pij) - 0.25*numpy.dot(pi-pj, pi-pj)
            b = numpy.dot(pij, pi-pj)

        return a, b
            
    def rotate (self, modes, i, j) :
        """
        rotates normal modes.

        Parameters
        ----------
        modes : ndarray (3*natoms-6,3*natoms)
            normal modes.
        i,j : int
            index of modes array.

        Returns
        -------
        rotated modes : tba ?
        tba ?
        alpha : tba ?
        tba ?
        """
        a, b = self.calc_ab(modes, i, j)

        alpha = math.acos(-a/math.sqrt(a*a+b*b))

        if math.sin(alpha) - (b/math.sqrt(a*a+b*b)) > 1e-6 :
            alpha = -alpha + 2.0*math.pi
        if (alpha < -math.pi) :
            alpha = alpha + 2.0*math.pi
        if (alpha >= math.pi) :
            alpha = alpha - 2.0*math.pi
        alpha = alpha / 4.0

        rotated_modes = modes.copy()
        rotated_modes[i,:] =   math.cos(alpha)*modes[i,:] + math.sin(alpha)*modes[j,:] 
        rotated_modes[j,:] = - math.sin(alpha)*modes[i,:] + math.cos(alpha)*modes[j,:] 

        return rotated_modes, alpha

    def try_localize (self, subset=None, thresh=1e-6, thresh2=1e-4, printing=False) :
        """
        maximizes the localization measure (procedere).
   
        Parameters
        ----------
        subset : ndarray(? x ?)
            subset of normal modes. or list of indicies?
        thresh : float
            first optimization threshold.
        thresh2 : float
            second optimization threshold. 
        printing : Bool
            If True: print out optimization informations.


        Returns
        -------
        transmat : tba ?
        tba ?
        del_p : ?
        tba ?
        """
        if subset is None :
            ss = list(range(self.nmodes))
        else:
            ss = subset

        loc_modes = self.locmodes.modes_mw[ss].copy()
        start_p = self.calc_p(loc_modes)

        transmat = numpy.identity(len(ss))

        def rotmat(i, j, alpha) :
            rmat = numpy.identity(len(ss))
            rmat[i,i] =  math.cos(alpha)
            rmat[i,j] =  math.sin(alpha)
            rmat[j,i] = -math.sin(alpha)
            rmat[j,j] =  math.cos(alpha)
            return rmat
        
        err  = 1e10
        err2 = 1e10
        isweep = 0
        while (err > thresh) or (err2 > thresh2) :
            isweep += 1

            err2 = 0.0
            if isweep == 1:
                old_p = start_p
            else:
                old_p = self.calc_p(loc_modes)
            for i in range(len(ss)) :
                for j in range(i+1, len(ss)) :
                    loc_modes, alpha = self.rotate(loc_modes, i, j)

                    transmat = numpy.dot(rotmat(i,j,alpha), transmat)
                    
                    err2 += abs(alpha)
                    
            p = self.calc_p(loc_modes)
            err = p - old_p

            if printing :
                print(" Normal mode localization: Cycle %3i    p: %8.3f   change: %10.7f  %10.5f " % \
                      (isweep, p, err, err2))

        del_p = p - start_p
        return transmat, del_p

    def localize (self, subset=None, thresh=1e-6, thresh2=1e-4, printing=True) :
        """
        maximizes the localization measure 
        and set the (localization) transformation matrix in LocVib class.
   
        Parameters
        ----------
        subset : narray(? x ?)
            subset of normal modes.
        thresh : float
            first optimization threshold.
        thresh2 : float
            second optimization threshold. 
        printing : Bool
            If True: print out optimization informations. 
        """
        if subset is None :
            ss = list(range(self.nmodes))
        else:
            ss = subset

        transmat, del_p = self.try_localize(subset, thresh, thresh2, printing)

        full_transmat = numpy.identity(self.nmodes)
        full_transmat[numpy.ix_(ss,ss)] = transmat
        self.set_transmat(numpy.dot(full_transmat, self.transmat))

    def localize_subsets(self, subsets, printing=True):
        """
        Task ?        

        Parameters
        ----------
        subsets : list of integer lists
            subsets = [[0,1],[2],[0,1,2]]
        """
        for subset in subsets:
            self.localize(subset, printing=printing)

        self.sort_by_freqs()
        self.subsets = subsets

    def localize_automatic_subsets(self, maxerr):
        """
        Automatic localization assignments of subsets.

        Parameters
        ----------
        maxerr : tba?
            tba ?
        """
        auto_assignment = AutomaticAssignment(self)
        self.subsets = auto_assignment.automatic_subsets(maxerr)
        self.sort_by_freqs()

    def try_localize_vcisdiff(self, subset=None, thresh=1e-6, thresh2=1e-4) :
        """
        tries to localize vibrational configuration intersection singles.
        *See: Panek, Hoeske, Jacob, J. Chem. Phys. 150, 054107 (2019)
        (doi: 10.1063/1.5083186)* 

        Parameters
        ----------
        subset : tba?
            tba?
        thresh : tba?
            tba?
        thresh2 : tba?
            tba?

        Returns
        -------
        vcis : tba
        tba
        del_p : tba
        tba
        """
        if subset is None :
            ss = list(range(self.nmodes))
        else:
            ss = subset

        transmat, del_p = self.try_localize(subset, thresh, thresh2, printing=False)
        tmat = numpy.dot(transmat, self.transmat[numpy.ix_(ss,ss)])

        diag2 = numpy.diag(self.startmodes.freqs[ss]**2)
        hmat = numpy.dot(numpy.dot(tmat, diag2), tmat.transpose())
        locfreqs = numpy.sqrt(numpy.diagonal(hmat))

        vcis_mat = hmat[:,:] / (2.0*numpy.sqrt(numpy.outer(locfreqs[:], locfreqs[:])))
        vcis_mat[numpy.diag_indices_from(vcis_mat)] = locfreqs

        ev, evecs = numpy.linalg.eigh(vcis_mat)
        vcis_maxdiff = numpy.max(numpy.abs(ev - numpy.sort(self.startmodes.freqs[ss])))

        return vcis_maxdiff, del_p

    def set_transmat (self, tmat) :
        """
        Sets the transformation matrix.

        Parameters
        ----------
        tmat : tba
            tba
        """
        self.transmat = tmat
        self.locmodes = self.startmodes.transform(tmat)

    def get_couplingmat (self, hessian=False) :
        """
        TEXT.

        Parameters
        ----------
        hessian : bool
            tba

        Returns
        -------
        cmat : ndarray
           Coupling matrix.
        """
        # coupling matrix is defined as in the paper:
        #   eigenvectors are rows of U / columns of U^t = transmat
        #   eigenvectors give normal mode in basis of localized modes
        if not hessian :
            diag = numpy.diag(self.startmodes.freqs)
        else :
            freqs = self.startmodes.freqs.copy() * 1e2 * Constants.cvel_ms
           # now freq is nu[1/s] 
            omega = 2.0 * math.pi * freqs
            ## now omega is omega[1/s]
            omega = omega**2 # matrix will be in [1/s2]
            omega = omega/(Constants.Hartree_in_Joule/(Constants.amu_in_kg*Constants.Bohr_in_Meter**2))  # 
            ### 

            # TEST PAWEL - correct conversion cm-1 -> au
            freqs = self.startmodes.freqs.copy()
            freqs = freqs * Constants.atu_in_s * (2.0*math.pi*1e2*Constants.cvel_ms)
            omega = freqs**2
            # END TEST PAWEL

            diag =  numpy.diag(omega)
        print('Obtaining coupling matrix in [a.u.]')
        # first normal modes    
        cmat = numpy.dot(numpy.dot(self.transmat, diag), self.transmat.transpose())

        # For consistency: shouldn't this function return values always in the same units?
        # now if hessian=False, returns [cm-1] otherwise returns [hartree]
        return cmat

    def get_vcismat (self) :
        """
        TEXT

        Returns
        -------
        vcis_mat : tba
        tba
        """
        diag = numpy.diag(self.startmodes.freqs**2)
        hmat = numpy.dot(numpy.dot(self.transmat, diag), self.transmat.transpose())

        locfreqs = numpy.sqrt(hmat.diagonal())

        vcis_mat = hmat[:,:] / (2.0*numpy.sqrt(numpy.outer(locfreqs[:], locfreqs[:])))
        vcis_mat[numpy.diag_indices_from(vcis_mat)] = locfreqs

        # VCIS matrix in cm-1
        return vcis_mat


    def sort_by_residue (self, pdb_mol=None) :
        """
        TEXT

        Parameters
        ----------
        pdb_mol : tba
            tba
        """
        if pdb_mol:
            print('Using external molecule definition')
            sortmat = self.locmodes.sortmat_by_residue(external_molecule=pdb_mol)
        else:
            sortmat = self.locmodes.sortmat_by_residue()
        tmat = numpy.dot(sortmat, self.transmat)
        self.set_transmat(tmat)


    def sort_by_groups (self, groups) :
        """
        TEXT

        Parameters
        ----------
        groups : tba
            tba
        """
        sortmat = self.locmodes.sortmat_by_groups(groups)
        tmat = numpy.dot(sortmat, self.transmat)
        self.set_transmat(tmat)


    def sort_by_freqs (self) :
        """
        TEXT
        """
        sortmat = self.locmodes.sortmat_by_freqs()
        tmat = numpy.dot(sortmat, self.transmat)
        self.set_transmat(tmat)

    def adjust_signs (self) :
        """
        TEXT
        """
        for imode in range(1,self.nmodes) :
            cmat = self.get_couplingmat()
            if cmat[imode-1,imode] < 0.0 :
                tmat = self.transmat
                tmat[imode,:] = -tmat[imode,:]
                self.set_transmat(tmat)

    def invert_signs (self, nums) :
        """
        TEXT

        Parameters
        ----------
        nums : tba
            tba
        """
        tmat = self.transmat
        tmat[nums,:] = -tmat[nums,:]
        self.set_transmat(tmat)

    def flip_modes (self, m1, m2) :
        """
        TEXT

        Parameters
        ----------
        m1 : tba
            tba
        m2 : tba
            tba
        """
        sortmat = numpy.identity(self.nmodes)
        sortmat[m1,m1] = 0.0
        sortmat[m2,m2] = 0.0
        sortmat[m1,m2] = 1.0
        sortmat[m2,m1] = 1.0

        tmat = numpy.dot(sortmat, self.transmat)
        self.set_transmat(tmat)

class AutomaticAssignment(object):
    """ TEXT ... """

    def __init__(self, lv) :
        """
        AutomaticAssignment constructor.
        For more details see the class description/docstring
        """
        self.lv = lv
        nmodes = self.lv.nmodes

        self.errmat = numpy.zeros((nmodes, nmodes))
        self.pmat = numpy.zeros((nmodes, nmodes))
        self.diffmat = numpy.zeros((nmodes, nmodes))

        self.subsets = [[m] for m in range(nmodes)]

        for i in range(len(self.subsets)) :
            for j in range(i+1, len(self.subsets)) :
                self.calc_ij(i, j)

    def calc_ij(self, i, j):
        """
        TEXT 
        
        Parameters
        ----------
        i,j : int
            tba
        """
        s1 = self.subsets[i]
        s2 = self.subsets[j]

        vciserr, p = self.lv.try_localize_vcisdiff(s1+s2)

        self.errmat[numpy.ix_(s1,s1)] = 0.0
        self.errmat[numpy.ix_(s2,s2)] = 0.0
        self.errmat[numpy.ix_(s1,s2)] = 0.0
        self.errmat[numpy.ix_(s2,s1)] = 0.0

        self.errmat[s1[0], s2[0]] = vciserr
        self.errmat[s2[0], s1[0]] = vciserr

        self.pmat[numpy.ix_(s1,s1)] = 0.0
        self.pmat[numpy.ix_(s2,s2)] = 0.0
        self.pmat[numpy.ix_(s1,s2)] = 0.0
        self.pmat[numpy.ix_(s2,s1)] = 0.0

        self.pmat[s1[0], s2[0]] = p
        self.pmat[s2[0], s1[0]] = p

        diff1 = abs(numpy.max(self.lv.locmodes.freqs[s1]) - numpy.min(self.lv.locmodes.freqs[s2]))
        diff2 = abs(numpy.min(self.lv.locmodes.freqs[s1]) - numpy.max(self.lv.locmodes.freqs[s2]))

        self.diffmat[numpy.ix_(s1,s1)] = 0.0
        self.diffmat[numpy.ix_(s2,s2)] = 0.0
        self.diffmat[numpy.ix_(s1,s2)] = 0.0
        self.diffmat[numpy.ix_(s2,s1)] = 0.0

        self.diffmat[s1[0], s2[0]] = max(diff1, diff2)
        self.diffmat[s2[0], s1[0]] = max(diff1, diff2)

    def find_maxp(self, maxerr, maxdiff):
        """
        TEXT

        Parameters
        ----------
        maxerr : tba
            tba
        maxdiff : tba
            tba

        Returns
        -------
        ind_maxp : tba
        tba
        maxp : tba
        tba
        """
        pmat_masked = numpy.where(self.errmat < maxerr, self.pmat, 0.0)
        pmat_masked = numpy.where(self.diffmat < maxdiff, pmat_masked, 0.0)

        ind_maxp_modes = numpy.unravel_index(numpy.argmax(pmat_masked), self.pmat.shape)
        maxp = numpy.max(pmat_masked)

        for i, ss in enumerate(self.subsets):
            if ind_maxp_modes[0] in ss:
                i_maxp = i
            if ind_maxp_modes[1] in ss:
                j_maxp = i
        ind_maxp = (i_maxp, j_maxp)

        return ind_maxp, maxp

    def update_subset(self, i):
        """
        TEXT

        Parameters
        ----------
        i : int
            tba
        """
        self.lv.localize(self.subsets[i], printing=False)
        for j in range(len(self.subsets)) :
            if not (i == j) :
                self.calc_ij(i, j)

    def automatic_subsets(self, maxerr) :
        """
        TEXT

        Parameters
        ----------
        maxerr : tba
            tba

        Returns
        -------
        subsets : tba
        tba
        """
        # find minimal frequency difference
        mask = numpy.ones(self.diffmat.shape, dtype=bool)
        numpy.fill_diagonal(mask, 0)
        m = self.diffmat[mask].min()

        # maxdiff starts from the minimal frequency difference and is doubled until
        # it is larger than 1e4
        for maxdiff in [m * 2**n for n in range(int(math.log(1e4/m, 2))+2)]:
            ind_maxp, maxp = self.find_maxp(maxerr, maxdiff)

            while (maxp > 0.01) :

                i, j = ind_maxp[0], ind_maxp[1]

                self.subsets[i] = sorted(self.subsets[i] + self.subsets[j])
                del self.subsets[j]

                self.update_subset(i)

                ind_maxp, maxp = self.find_maxp(maxerr, maxdiff)

        return self.subsets
