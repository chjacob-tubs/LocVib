
import math
import numpy

import Constants

class VibModes (object) :

    def __init__ (self, nmodes, mol) :
        self.nmodes = nmodes

        self.mol = mol
        self.natoms = mol.natoms

        self.freqs    = numpy.arange(1.0, nmodes+0.1, 1.0)
        
        self.modes_mw = numpy.zeros((nmodes, 3*self.natoms))
        self.modes_c  = numpy.zeros((nmodes, 3*self.natoms))

    def set_modes_mw (self, modes_mw)  :
        self.modes_mw[:,:] = modes_mw
        self.update_modes_c()

    def set_modes_c (self, modes_c)  :
        self.modes_c[:,:] = modes_c
        self.update_modes_mw()
        self.normalize_modes()

    def set_freqs (self, freqs) :
        self.freqs = freqs
        
    def get_modes_c_norm (self) :
        norm_modes = self.modes_c.copy()
        for imode in range(self.nmodes) :
            nc = math.sqrt((self.modes_c[imode,:]**2).sum())
            norm_modes[imode,:] = (1.0/nc) * norm_modes[imode,:]
        return norm_modes

    def normalize_modes (self) :
        for imode in range(self.nmodes) :
            nc = math.sqrt((self.modes_mw[imode,:]**2).sum())
            if abs(nc) > 0.0 :
                self.modes_mw[imode,:] = (1.0/nc) * self.modes_mw[imode,:]
        self.update_modes_c()

    modes_c_norm = property(get_modes_c_norm)
    
    def update_modes_c (self) :
        for i, atmass in enumerate(self.mol.atmasses) :
            self.modes_c[:,3*i:3*i+3] = self.modes_mw[:,3*i:3*i+3] / math.sqrt(atmass)

    def update_modes_mw (self) :
        for i, atmass in enumerate(self.mol.atmasses) :
            self.modes_mw[:,3*i:3*i+3] = self.modes_c[:,3*i:3*i+3] * math.sqrt(atmass)

    def get_subset (self, modelist) :
        subset = VibModes(len(modelist), self.mol)
        subset.set_modes_mw(self.modes_mw[modelist])
        subset.set_freqs(self.freqs[modelist])
        return subset

    def get_range (self, minfreq, maxfreq) :
        modelist = numpy.where((self.freqs > minfreq) & (self.freqs < maxfreq))[0]
        return self.get_subset(modelist)
            
    def write_g98out (self, massweighted=False, normalize=True, filename='g98.out') :

        if massweighted :
            modes = self.modes_mw
        else:
            if normalize :
                modes = self.modes_c_norm
            else:
                modes = self.modes_c

        if self.nmodes % 3 == 0 :
            temp_modes = modes
        else:
            temp_modes = numpy.zeros((3*(self.nmodes/3)+3,3*self.natoms))
            temp_modes[:self.nmodes,:] = modes[:,:]

        temp_freqs = numpy.zeros((temp_modes.shape[0],))
        temp_freqs[:self.nmodes] = self.freqs
        
        f = file(filename, 'w')
        f.write(' Entering Gaussian System \n')
        f.write(' *********************************************\n')
        f.write(' Gaussian 98:\n')
        f.write(' frequency fake output\n')
        f.write(' *********************************************\n')
        f.write('                         Standard orientation:\n')
        f.write(' --------------------------------------------------------------------\n')
        f.write('  Center     Atomic     Atomic              Coordinates (Angstroms)  \n')
        f.write('  Number     Number      Type              X           Y           Z \n')
        f.write(' --------------------------------------------------------------------\n')

        atnums = self.mol.atnums
        coords = self.mol.coordinates

        for iatom in range(self.natoms) :
            f.write(' %4i       %4i             0     %11.6f %11.6f %11.6f \n' %
                    (iatom+1, atnums[iatom], coords[iatom,0], coords[iatom,1], coords[iatom,2]))

        f.write(' --------------------------------------------------------------------\n')
        f.write('     1 basis functions        1 primitive gaussians \n')
        f.write('     1 alpha electrons        1 beta electrons\n')
        f.write('\n')
        f.write(' Harmonic frequencies (cm**-1), IR intensities (KM/Mole), \n')
        f.write(' Raman scattering activities (A**4/amu), Raman depolarization ratios, \n')
        f.write(' reduced masses (AMU), force constants (mDyne/A) and normal coordinates: \n')

        for i in range(0,temp_modes.shape[0],3) :
            f.write('                   %4i                   %4i                  %4i \n' %
                    (i+1, i+2, i+3))
            f.write('                     a                      a                      a  \n')
            f.write(' Frequencies -- %10.4f             %10.4f             %10.4f \n' %
                    (temp_freqs[i], temp_freqs[i+1], temp_freqs[i+2]))
            f.write(' Red. masses -- %10.4f             %10.4f             %10.4f \n' % (0.0, 0.0, 0.0))
            f.write(' Frc consts  -- %10.4f             %10.4f             %10.4f \n' % (0.0, 0.0, 0.0))
            f.write(' IR Inten    -- %10.4f             %10.4f             %10.4f \n' % (0.0, 0.0, 0.0))
            f.write(' Raman Activ --     0.0000                 0.0000                 0.0000 \n')
            f.write(' Depolar     --     0.0000                 0.0000                 0.0000 \n')
            f.write(' Atom AN      X      Y      Z        X      Y      Z        X      Y      Z \n')

            for iatom in range(self.natoms) :
                atnum = atnums[iatom]
                f.write('%4i %3i   %6.2f %6.2f %6.2f   %6.2f %6.2f %6.2f   %6.2f %6.2f %6.2f \n' %
                        (iatom+1, atnum,
                         temp_modes[i,3*iatom], temp_modes[i,3*iatom+1], temp_modes[i,3*iatom+2],
                         temp_modes[i+1,3*iatom], temp_modes[i+1,3*iatom+1], temp_modes[i+1,3*iatom+2],
                         temp_modes[i+2,3*iatom], temp_modes[i+2,3*iatom+1], temp_modes[i+2,3*iatom+2],))

        f.close()

        
    def get_composition (self, groups) :
        types = numpy.zeros((len(groups), self.nmodes))
        squared_modes = self.modes_mw**2
        for itype, attype in enumerate(groups) :
            for at in attype :
                types[itype,:] += squared_modes[:,3*at] 
                types[itype,:] += squared_modes[:,3*at+1] 
                types[itype,:] += squared_modes[:,3*at+2] 
        return types

    def print_composition (self, groupnames, comp, labels=None) :
        if labels==None :
            lab = ["%4i %8.1f   " % (i, self.freqs[i]) for i in range(comp.shape[1])]
        else:
            lab = labels
        
        print (" "*len(lab[0])),
        for t in groupnames :
            print "%4s  " % t,
        print

        for i in range(comp.shape[1]) :
            print lab[i],
            for t in range(comp.shape[0]) :
                print "%4.1f  " % (comp[t,i]*100.0),
            print

        print
        print "min: " + " "*(len(lab[0])-5),
        for t in range(comp.shape[0]) :
            print "%4.1f  " % (comp[t,:].min()*100.0),
        print

        print "max: " + " "*(len(lab[0])-5),
        for t in range(comp.shape[0]) :
            print "%4.1f  " % (comp[t,:].max()*100.0),
        print
        print

    def print_residue_composition(self, labels=None) :
        groups, groupnames = self.mol.residue_groups()
        comp = self.get_composition(groups)
        self.print_composition(groupnames, comp, labels)

    def print_attype_composition(self, labels=None) :
        groups, groupnames = self.mol.attype_groups()
        comp = self.get_composition(groups)
        self.print_composition(groupnames, comp, labels)

    def print_attype2_composition(self, labels=None) :
        groups, groupnames = self.mol.attype_groups_2()
        comp = self.get_composition(groups)
        self.print_composition(groupnames, comp, labels)

    def print_atom_composition(self, labels=None) :
        groups, groupnames = self.mol.atom_groups()
        comp = self.get_composition(groups)
        self.print_composition(groupnames, comp, labels)

    def transform (self, tmat) :
        tmodes = VibModes(self.nmodes, self.mol)
        tmodes.set_modes_mw(numpy.dot(tmat, self.modes_mw))

        cmat = numpy.dot(numpy.dot(tmat, numpy.diag(self.freqs)), tmat.transpose())

        tmodes.set_freqs(cmat.diagonal())

        return tmodes
    
    def sortmat_by_residue (self) :
        types = self.get_composition(self.mol.residue_groups()[0])

        max_res = []
        for i in range(self.nmodes) :
            maxind = -1
            for t in range(len(types)):
                if (types[t][i] > 0.4) :
                    maxind = t
            max_res.append((i,maxind))

        max_res.sort(lambda x,y: -cmp(self.freqs[x[0]], self.freqs[y[0]]))
        max_res.sort(lambda x,y: cmp(x[1],y[1]))

        sortmat = numpy.zeros((self.nmodes, self.nmodes))
        for i in range(self.nmodes) :
            sortmat[i, max_res[i][0]] = 1.0
            
        return sortmat

    def get_fragment_modes (self, atomlist) :
        frag = self.mol.get_fragment(atomlist)
        
        indices = [0]*3*len(atomlist)
        for i, a in enumerate(atomlist) :
            indices[3*i]   = 3*a
            indices[3*i+1] = 3*a+1
            indices[3*i+2] = 3*a+2

        modes = VibModes(self.nmodes, frag)
        modes.set_freqs(self.freqs)

        modes.set_modes_mw(self.modes_mw[:,indices])
        
        return modes
        
    def overlap (self, other) :
        ov = numpy.zeros((self.nmodes,))

        for i in range(self.nmodes) :
            ov[i] = numpy.dot(self.modes_mw[i], other.modes_mw[i])

        ov = ov**2

        return ov

    def center (self, imode) :
        sq_mode = self.modes_mw[imode,:]**2

        c = numpy.zeros((self.natoms,))
        for iatom in range(self.natoms) :
            c[iatom] = sq_mode[3*iatom] + sq_mode[3*iatom+1] + sq_mode[3*iatom+2]

        cen = numpy.zeros((3,))
        coords = self.mol.coordinates
        
        for iatom in range(self.natoms) :
            cen[0] += coords[iatom,0] * c[iatom]
            cen[1] += coords[iatom,1] * c[iatom]
            cen[2] += coords[iatom,2] * c[iatom]

        return cen

    def distance (self, imode, jmode) :
        dist = self.center(imode) - self.center(jmode)
        return math.sqrt(numpy.dot(dist,dist))
    
    def get_tdc_couplingmat (self, dipole) :
        tdcmat = numpy.diag(self.freqs)

        for i in range(self.nmodes) :
            dipi_mag = math.sqrt(numpy.dot(dipole[i], dipole[i]))
            dipi_vec = dipole[i] / dipi_mag

            #dipole is in  a.u./(amu^-0.5 * A)
            print i, dipi_mag * (Constants.au_in_Debye/Constants.Bohr_in_Angstrom)
        
            for j in range(i+1, self.nmodes) :
                dipj_mag = math.sqrt(numpy.dot(dipole[j], dipole[j]))
                dipj_vec = dipole[j] / dipj_mag

                dist_vec = self.center(i) - self.center(j)
                dist_mag = math.sqrt(numpy.dot(dist_vec,dist_vec))
                dist_vec = dist_vec / dist_mag
                dist_mag = dist_mag / Constants.Bohr_in_Angstrom
                #dist_mag is now in bohr

                x = numpy.dot(dipi_vec,dipj_vec) - \
                    3.0*numpy.dot(dipi_vec,dist_vec)*numpy.dot(dipj_vec,dist_vec)

                tdcmat[i,j] = (dipi_mag*dipj_mag / dist_mag**3) * x 

        # convert to Omega matrix in cm-1
        fact = (Constants.Hartree_in_Joule/(Constants.amu_in_kg*Constants.Bohr_in_Meter**2))
        fact = fact / (Constants.cvel_ms*1e2)**2
        for i in range(self.nmodes) :
            for j in range(i+1, self.nmodes) :
                tdcmat[i,j] = tdcmat[i,j] * fact
                tdcmat[i,j] = tdcmat[i,j] / (4*math.pi**2 * (tdcmat[i,i]+tdcmat[j,j]))
                tdcmat[j,i] = tdcmat[i,j]
 
        return tdcmat
    
        
