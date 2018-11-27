# -*- coding: utf-8 -*-
#
# This file is part of the
# LocVib 1.0 suite of tools for the analysis for vibrational spectra.
# Copyright (C) 2009-2012 by Christoph R. Jacob and others.
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
#   http://www.christophjacob.eu/locvib.php
"""
 Localizing normal modes
"""
import math
import numpy
import copy

from Modes import VibModes
import Constants

class LocVib (object) :

    def __init__ (self, modes, loctype='PM') :
        self.startmodes = modes
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
        squared_modes = modes**2

        c = squared_modes[:,0::3] + squared_modes[:,1::3] + squared_modes[:,2::3]
                
        if self.loctype == 'PM' :
            p = numpy.linalg.norm(c)**2
        elif self.loctype.startswith('B') :
            p_mode = numpy.zeros((3,self.nmodes))
            for idir in range(3):
                p_mode[idir] = (self.coords[:,idir] * c).sum(axis=1)
            p = (p_mode**2).sum()

        return p

    def calc_ab (self, modes, i, j) :
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
        if subset is None :
            ss = range(self.nmodes)
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
                print " Normal mode localization: Cycle %3i    p: %8.3f   change: %10.7f  %10.5f " % \
                      (isweep, p, err, err2)

        del_p = p - start_p
        return transmat, del_p

    def localize (self, subset=None, thresh=1e-6, thresh2=1e-4, printing=True) :
        if subset is None :
            ss = range(self.nmodes)
        else:
            ss = subset

        transmat, del_p = self.try_localize(subset, thresh, thresh2, printing)

        full_transmat = numpy.identity(self.nmodes)
        full_transmat[numpy.ix_(ss,ss)] = transmat
        self.set_transmat(numpy.dot(full_transmat, self.transmat))

    def localize_subsets(self, subsets, printing=True):
        for subset in subsets:
            self.localize(subset, printing=printing)

        self.sort_by_freqs()
        self.subsets = subsets

    def localize_automatic_subsets(self, maxerr):
        auto_assignment = AutomaticAssignment(self)
        self.subsets = auto_assignment.automatic_subsets(maxerr)
        self.sort_by_freqs()

    def try_localize_vcisdiff(self, subset=None, thresh=1e-6, thresh2=1e-4) :
        if subset is None :
            ss = range(self.nmodes)
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
        self.transmat = tmat
        self.locmodes = self.startmodes.transform(tmat)

    def get_couplingmat (self, hessian=False) :
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
        print 'Obtaining coupling matrix in [a.u.]'
        # first normal modes    
        cmat = numpy.dot(numpy.dot(self.transmat, diag), self.transmat.transpose())

        # For consistency: shouldn't this function return values always in the same units?
        # now if hessian=False, returns [cm-1] otherwise returns [hartree]
        return cmat

    def get_vcismat (self) :
        diag = numpy.diag(self.startmodes.freqs**2)
        hmat = numpy.dot(numpy.dot(self.transmat, diag), self.transmat.transpose())

        locfreqs = numpy.sqrt(hmat.diagonal())

        vcis_mat = hmat[:,:] / (2.0*numpy.sqrt(numpy.outer(locfreqs[:], locfreqs[:])))
        vcis_mat[numpy.diag_indices_from(vcis_mat)] = locfreqs

        # VCIS matrix in cm-1
        return vcis_mat

    def sort_by_residue (self) :
        sortmat = self.locmodes.sortmat_by_residue()
        tmat = numpy.dot(sortmat, self.transmat)
        self.set_transmat(tmat)

    def sort_by_groups (self, groups) :
        sortmat = self.locmodes.sortmat_by_groups(groups)
        tmat = numpy.dot(sortmat, self.transmat)
        self.set_transmat(tmat)

    def sort_by_freqs (self) :
        sortmat = self.locmodes.sortmat_by_freqs()
        tmat = numpy.dot(sortmat, self.transmat)
        self.set_transmat(tmat)

    def adjust_signs (self) :
        for imode in range(1,self.nmodes) :
            cmat = self.get_couplingmat()
            if cmat[imode-1,imode] < 0.0 :
                tmat = self.transmat
                tmat[imode,:] = -tmat[imode,:]
                self.set_transmat(tmat)

    def invert_signs (self, nums) :
        tmat = self.transmat
        tmat[nums,:] = -tmat[nums,:]
        self.set_transmat(tmat)

    def flip_modes (self, m1, m2) :
         sortmat = numpy.identity(self.nmodes)
         sortmat[m1,m1] = 0.0
         sortmat[m2,m2] = 0.0
         sortmat[m1,m2] = 1.0
         sortmat[m2,m1] = 1.0

         tmat = numpy.dot(sortmat, self.transmat)
         self.set_transmat(tmat)

class AutomaticAssignment(object):

    def __init__(self, lv) :
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
        self.lv.localize(self.subsets[i], printing=False)
        for j in range(len(self.subsets)) :
            if not (i == j) :
                self.calc_ij(i, j)

    def automatic_subsets(self, maxerr) :

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
