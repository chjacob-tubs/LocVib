#
# vibrational analysis for SNF
#

import numpy 
import math

from Constants import *

from Molecule import VibToolsMolecule
from Modes    import VibModes

class SNFRestartFile (object) :

    def __init__ (self) :
        self.n_per_line = None
        
        self.nvar   = None
        self.natoms = None
        self.iaus   = None
        self.displ  = None

        self.intonly = False
        self.nmodes  = 0

    def read (self, filename='restart') :
        f = file(filename, 'r')
        
        f.readline()  # first line: some path

        format = f.readline()  # second line: format
        format = format.replace('d', 'D')
        self.n_per_line = int(format[1:format.find('D')])
        
        line = f.readline() # third line
        line = [int(i) for i in line.split()]

        self.nvar   = line[4]
        self.natoms = self.nvar / 3
        self.iaus   = line[5]
        self.nmodes = line[6]
        if self.nmodes > 0 :
            self.intonly = True

        # number of sets of calculations 
        # normal run: one set are the 6 Cartesian displacements
        # intonly mode: one set are the two displacements per normal mode
        if self.intonly :
            self.ncalcsets = self.nmodes
        else :
            self.ncalcsets = self.natoms

        self.displ = float(f.readline().replace('D', 'E'))
        f.readline()  # fourth line: some other path
        f.readline()  # fifth line: and yet another path

        for i in xrange(self.ncalcsets) :
            line = f.readline()
            line = [int(i) for i in line.split()]
            for i in line[1:] :
                if not (i == 1) :
                    raise Exception('Not all points finished in SNF (step %i)' % line[0])
            f.readline()

        # FIXME: This only works with the ROA SNF version.
        #        It should be rewritten in a more general way that
        #        automatically detects which sections are present .
            
        self.dipole   = self._read_section(f, 3)
        if not self.intonly :
            self.gradient = self._read_section(f, self.nvar)
        self.pollen   = self._read_section(f, 8)
        self.polvel   = self._read_section(f, 8)
        self.gtenlen  = self._read_section(f, 9)
        self.gtenvel  = self._read_section(f, 9)
        self.aten     = self._read_section(f, 27)

        f.close()

    def _read_section(self, f, ncomp) :

        if self.intonly :
            nread = ncomp * self.iaus 
        else :
            nread = ncomp * self.iaus * 3

        result = numpy.zeros((self.ncalcsets, nread))
        
        nlines = nread / self.n_per_line
        if nlines*self.n_per_line < nread:
            nlines = nlines + 1

        for i in xrange(self.ncalcsets) :
            f.readline()
            n = 0
            for j in xrange(nlines) :
                line = f.readline().replace('D', 'E').split()
                for fl in line :
                    result[i,n] = float(fl)
                    n += 1

        result.shape = (self.ncalcsets, nread / ncomp, ncomp)
        return result

    
class SNFOutputFile (object) :

    def __init__ (self, filename='snf.out') :
        self.modes = None
        self.lwl = None
        self.filename = filename
        self.intonly = False

    def read (self, mol) :
        f = file(self.filename, 'r')
        lines = f.readlines()
        f.close()

        for l in lines :
            if 'intensities-only-mode' in l :
                self.intonly = True

        first = 0
        last  = 0
        for i, l in enumerate(lines) :
            if (first == 0) and ('root no.' in l) :
                first = i
            if 'Generate fake outputs for' in l :
                last = i
                break
            if 'W A R N I N G' in l :
                last = i-1
                break
            if 'Raman Optical Activity properties for freq.' in l :
                self.lwl = float(l[46:57])

        if not self.intonly :
            natoms = mol.natoms
            self.modes = VibModes(3*natoms-6, mol)
                
            lines = lines[first:last]
            lines = [l for l in lines if not len(l.strip()) == 0]
            lines = [l for l in lines if not l.startswith('1')]

            normalmodes = numpy.zeros((3*natoms-6, 3*natoms))
            freqs = numpy.zeros((3*natoms-6))
        
            ncol = len(lines[2][10:].split())

            icol = -ncol
            for i, l in enumerate(lines) :
                irow = i % (3*natoms+2) - 2
            
                if irow == -2:
                    if not ('root' in l) :
                        raise Exception('wrong number of lines for normal modes')
                    icol += ncol
                else:
                    line = [float(fl) for fl in l[10:].split()]

                    if irow == -1:
                        freqs[icol:icol+ncol] = line
                    else:
                        normalmodes[icol:icol+ncol,irow] = line
                
            self.modes.set_modes_mw(normalmodes)
            self.modes.set_freqs(freqs)

class SNFControlFile (object) :

    def __init__ (self) :
        self.modes = None 

    def read (self, mol) :
        f = file('snf_control', 'r')
        lines = f.readlines()
        f.close()

        for i, l in enumerate(lines) :
             if l.startswith('$nummodes') :
                 nmodes = int(l.split()[-1])
                 start = i+1

        natoms = mol.natoms

        self.modes = VibModes(nmodes, mol)

        normalmodes = numpy.zeros((nmodes, 3*natoms))
        freqs = numpy.zeros((nmodes))
       
        for imode in range(nmodes) :
            freqs[imode] = float(lines[start+imode*natoms].split()[0])
            for i in range(natoms) :
                x, y, z = lines[start+imode*natoms+i].split()[-3:]
                normalmodes[imode,3*i:3*i+3] = [x, y, z]

        self.modes.set_modes_c(normalmodes)
        self.modes.set_freqs(freqs)
        
class SNFResults (object) :

    def __init__ (self, outname='snf.out', restartname='restart', coordfile='coord') :
        self.mol          = VibToolsMolecule()
        self.snfoutput    = SNFOutputFile(filename=outname) 
        self.snfcontrol   = SNFControlFile() 
        self.restartfile  = SNFRestartFile()
        self.restartname  = restartname
        self.coordfile    = coordfile

        self.intonly      = False
    
    def _get_modes (self) :
        if self.snfoutput.intonly :
            return self.snfcontrol.modes
        else :
            return self.snfoutput.modes
    
    modes = property(_get_modes)
    freqs = property(lambda self: self.modes.freqs)
    lwl   = property(lambda self: self.snfoutput.lwl)
        
    def read (self) :
        self.mol.read_from_coord(filename=self.coordfile)

        self.snfoutput.read(self.mol)
        self.intonly = self.snfoutput.intonly
        if self.snfoutput.intonly :
            self.snfcontrol.read(self.mol)

        self.restartfile.read(filename=self.restartname)

    def get_mw_normalmodes (self) :
        return self.modes.modes_mw

    def get_c_normalmodes (self) :
        return self.modes.modes_c

    def get_tensor_mean (self, tens, ncomp=None) :
        tensor = eval('self.restartfile' + '.' + tens)

        ncalcsets = tensor.shape[0]
        nval = tensor.shape[1] / 2
        if (ncomp == None) :
            ncomp = tensor.shape[2]

        displ = self.restartfile.displ

        mean = numpy.zeros((ncalcsets, nval, ncomp))

        for i in range(nval) :
            mean[:,i,:] = (tensor[:,2*i,:ncomp] + tensor[:,2*i+1,:ncomp]) / 2.0 

        return mean

    def get_tensor_deriv_c (self, tens, ncomp=None) :

        if hasattr(self, tens+'_deriv_c') :
            return eval('self.'+tens+'_deriv_c')

        if not (self.restartfile.iaus == 2) :
            raise Exception('Differentiation only implemented for iaus = 2')

        if self.intonly :
            raise Exception('Cartesian derivatives not available in intensities only mode')

        tensor = eval('self.restartfile' + '.' + tens)

        natoms = tensor.shape[0]
        if (ncomp == None) :
            ncomp = tensor.shape[2]

        displ = self.restartfile.displ

        deriv = numpy.zeros((natoms, 3, ncomp))

        for i in range(3) :
            deriv[:,i,:] = (tensor[:,2*i,:ncomp] - tensor[:,2*i+1,:ncomp]) / (2.0 * displ)

        setattr(self, tens+'_deriv_c', deriv)

        return deriv

    def get_tensor_deriv_nm (self, tens, ncomp=None, modes=None) :

        if modes==None :
            if hasattr(self, tens+'_deriv_nm') :
                return eval('self.'+tens+'_deriv_nm')
        else :
            raise Exception('Argument modes of get_tensor_deriv_nm must be None in intensities-only mode')

        if self.intonly :
            tensor = eval('self.restartfile' + '.' + tens)

            nmodes = tensor.shape[0]
            if (ncomp == None) :
                ncomp = tensor.shape[2]

            displ = self.restartfile.displ

            deriv_nm = numpy.zeros((nmodes, ncomp))
            deriv_nm[:,:] = (tensor[:,0,:ncomp] - tensor[:,1,:ncomp]) / (2.0 * displ)

            modes = self.get_c_normalmodes()
            for i in range(nmodes):
                deriv_nm[i,:] = deriv_nm[i,:] * math.sqrt(sum(modes[i,:]**2))

        else :
            deriv_c  = self.get_tensor_deriv_c (tens, ncomp)

            ncomp = deriv_c.shape[2]

            if modes==None :
                modes_c = self.get_c_normalmodes()
            else:
                modes_c = modes.modes_c

            nmodes  = modes_c.shape[0]
            natoms  = modes_c.shape[1] / 3

            deriv_nm = numpy.zeros((nmodes, ncomp))

            for icomp in range(ncomp) :
                deriv_nm[:, icomp] = numpy.dot(modes_c,deriv_c[:,:,icomp].reshape((3*natoms)))

        if modes==None :
            setattr(self, tens+'_deriv_nm', deriv_nm)

        return deriv_nm

    def get_ir_intensity (self, modes=None) :
        mu = self.get_tensor_deriv_nm('dipole', modes=modes)

        irint = mu[:,0]*mu[:,0] + mu[:,1]*mu[:,1] + mu[:,2]*mu[:,2]

        # convert (d\mu/dR)^2 from  (au^2/amu) to (C^2 m^2 / kg)
        irint = irint * ( (au_in_Debye/Bohr_in_Meter)**2 * (Debye_in_Cm)**2 * (1.0/amu_in_kg) )

        # multiply by [ 1.0/(4\pi\epsilon0) * (N_A * \pi)/cvel ] (Eq 13 in SNF paper)
        irint = irint * (1.0/12.0) * (1.0/epsilon0) * Avogadro/cvel_ms 

        # divide by cvel (Eq 14 in SNF paper)
        irint = irint / (cvel_ms * 1000.0)

        # result is IR absoption in km/mol
        return irint

    def get_a2_invariant (self, gauge='len', modes=None) :
        pol    = self.get_tensor_deriv_nm('pol'+gauge, ncomp=6, modes=modes)
        nmodes = pol.shape[0]

        a2 = (1.0/3.0) * (pol[:,0] + pol[:,3] + pol[:,5])
        a2 = (a2**2)*(Bohr_in_Angstrom**4)
        return a2

    def get_g2_invariant (self, modes=None) :
        pol    = self.get_tensor_deriv_nm('pollen', ncomp=6, modes=modes)
        nmodes = pol.shape[0]

        g2 = (1.0/2.0) * (   (pol[:,0] - pol[:,3])**2
                           + (pol[:,3] - pol[:,5])**2
                           + (pol[:,5] - pol[:,0])**2
                           + 6.0 * pol[:,1]**2  
                           + 6.0 * pol[:,2]**2  
                           + 6.0 * pol[:,4]**2  
                         )
        g2 = g2 * (Bohr_in_Angstrom**4)
        return g2

    def get_raman_int (self, modes=None) :
        return 45.0*self.get_a2_invariant(modes=modes) + 7.0*self.get_g2_invariant(modes=modes)

    def get_aG_invariant (self, gauge='len', modes=None) :
        # for consistency with SNF alpha is always in length repr
        pol     = self.get_tensor_deriv_nm('pollen', ncomp=6, modes=modes)
        gten    = self.get_tensor_deriv_nm('gten'+gauge, modes=modes)

        aG = (1.0/9.0) * (pol[:,0] + pol[:,3] + pol[:,5])   \
                       * (gten[:,0] + gten[:,4] + gten[:,8])
        aG = aG*(Bohr_in_Angstrom**4) * (1 / cvel) * 1e6

        return aG

    def get_bG_invariant (self, gauge='len', modes=None) :
        pol     = self.get_tensor_deriv_nm('pol'+gauge, ncomp=6, modes=modes)
        gten    = self.get_tensor_deriv_nm('gten'+gauge, modes=modes)

        bG = 0.5 * (3*pol[:,0]*gten[:,0] - pol[:,0]*gten[:,0] +
                    3*pol[:,1]*gten[:,1] - pol[:,0]*gten[:,4] +
                    3*pol[:,2]*gten[:,2] - pol[:,0]*gten[:,8] +
                    3*pol[:,1]*gten[:,3] - pol[:,3]*gten[:,0] +
                    3*pol[:,3]*gten[:,4] - pol[:,3]*gten[:,4] +
                    3*pol[:,4]*gten[:,5] - pol[:,3]*gten[:,8] +
                    3*pol[:,2]*gten[:,6] - pol[:,5]*gten[:,0] +
                    3*pol[:,4]*gten[:,7] - pol[:,5]*gten[:,4] +
                    3*pol[:,5]*gten[:,8] - pol[:,5]*gten[:,8])
        bG = bG*(Bohr_in_Angstrom**4) * (1 / cvel) * 1e6

        return bG

    def get_bA_invariant (self, modes=None) :
        # always using length repr
        pol     = self.get_tensor_deriv_nm('pollen', ncomp=6, modes=modes)
        aten    = self.get_tensor_deriv_nm('aten', modes=modes)

        bA = 0.5 * self.lwl * (  (pol[:,3]-pol[:,0])*aten[:,11]
                               + (pol[:,0]-pol[:,5])*aten[:,6]
                               + (pol[:,5]-pol[:,3])*aten[:,15]
                               + pol[:,1]*(aten[:,19]-aten[:,20]+aten[:,8]-aten[:,14]) 
                               + pol[:,2]*(aten[:,25]-aten[:,21]+aten[:,3]-aten[:,4]) 
                               + pol[:,4]*(aten[:,10]-aten[:,24]+aten[:,12]-aten[:,5])
                              )
        bA = bA*(Bohr_in_Angstrom**4) * (1 / cvel) * 1e6

        return bA

    def get_backscattering_int (self, modes=None) :
         return 1e-6*96.0*(self.get_bG_invariant(gauge='vel', modes=modes) + \
                           (1.0/3.0)*self.get_bA_invariant(modes=modes))

    def check_consistency (self) :
        print
        print "Checking consistency of SNF restart file "
        print "  Large numbers mean that the change of the polarizabilities deviates from linear behavior."
        print

        def check (ten) :
            if ten.startswith('pol') :
                mean_pol = self.get_tensor_mean(ten, 6)
            else:
                mean_pol = self.get_tensor_mean(ten)
            ncalcsets = mean_pol.shape[0]
            nval      = mean_pol.shape[1]
            ncomp     = mean_pol.shape[2]

            mean_pol = mean_pol.reshape((ncalcsets*nval,ncomp))
            for icomp in range(ncomp) :
                print " %14.8f  min: %3i  max:  %3i" % \
                      (mean_pol[:,icomp].max()-mean_pol[:,icomp].min(),
                       mean_pol[:,icomp].argmin(), mean_pol[:,icomp].argmax() )
            print

            wrong = {} 

            for i in range(ncalcsets) :
                if (max(abs(mean_pol[3*i:3*i+3,0] - mean_pol[0,0])) >0.001) :
                    wrong[i+1] = [ abs(mean_pol[3*i,0] - mean_pol[0,0]) > 0.001, 
                                   abs(mean_pol[3*i+1,0] - mean_pol[0,0]) > 0.001, 
                                   abs(mean_pol[3*i+2,0] - mean_pol[0,0]) > 0.001]
                    print i+1, mean_pol[3*i:3*i+3,0] - mean_pol[0,0]

            return wrong

        print " Dipole moment "
        return check('dipole')

        return

        print " Polarizability (length) "
        check('pollen')

        print " Polarizability (vel) "
        check('polvel')

        print " G-tensor (len) "
        check('gtenlen')

        print " G-tensor (vel) "
        check('gtenvel')

        print " A-tensor "
        check('aten')
