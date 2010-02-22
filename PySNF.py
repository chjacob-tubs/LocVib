#
# vibrational analysis for SNF
#

import numpy 
import math

import Constants

from Molecule import VibToolsMolecule
from Modes    import VibModes

class Result (object):
    pass

class SNFRestartFile (object) :

    def __init__ (self) :
        self.n_per_line = None
        
        self.nvar   = None
        self.natoms = None
        self.iaus   = None
        self.displ  = None

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
        
        self.displ = float(f.readline().replace('D', 'E'))
        f.readline()  # fourth line: some other path
        f.readline()  # fifth line: and yet another path

        for i in xrange(self.natoms) :
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
        self.gradient = self._read_section(f, self.nvar)
        self.pollen   = self._read_section(f, 8)
        self.polvel   = self._read_section(f, 8)
        self.gtenlen  = self._read_section(f, 9)
        self.gtenvel  = self._read_section(f, 9)
        self.aten     = self._read_section(f, 27)

        f.close()

    def _read_section(self, f, ncomp) :

        nread = ncomp * self.iaus * 3

        result = numpy.zeros((self.natoms, nread))
        
        nlines = nread / self.n_per_line
        if nlines*self.n_per_line < nread:
            nlines = nlines + 1

        for i in xrange(self.natoms) :
            f.readline()
            n = 0
            for j in xrange(nlines) :
                line = f.readline().replace('D', 'E').split()
                for fl in line :
                    result[i,n] = float(fl)
                    n += 1

        result.shape = (self.natoms, self.iaus * 3, ncomp)
        return result

    
class SNFOutputFile (object) :

    def __init__ (self, filename='snf.out') :
        self.modes = None
        self.lwl = None
        self.filename = filename

    def read (self, mol) :
        natoms = mol.natoms
        self.modes = VibModes(3*natoms-6, mol)

        f = file(self.filename, 'r')
        lines = f.readlines()
        f.close()

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

        
class SNFResults (object) :

    def __init__ (self, outname='snf.out', restartname='restart', coordfile='coord') :
        self.mol          = VibToolsMolecule()
        self.snfoutput    = SNFOutputFile(filename=outname) 
        self.restartfile  = SNFRestartFile()

        self.restartname  = restartname
        self.coordfile    = coordfile
        
    freqs = property(lambda self: self.snfoutput.modes.freqs)
    modes = property(lambda self: self.snfoutput.modes)
    lwl   = property(lambda self: self.snfoutput.lwl)
        
    def read (self) :
        self.mol.read_from_coord(filename=self.coordfile)
        self.snfoutput.read(self.mol)
        self.restartfile.read(filename=self.restartname)

    def get_mw_normalmodes (self) :
        return self.modes.modes_mw

    def get_c_normalmodes (self) :
        return self.modes.modes_c

    def get_tensor_mean (self, tens, ncomp=None) :
        tensor = eval('self.restartfile' + '.' + tens)

        natoms = tensor.shape[0]
        if (ncomp == None) :
            ncomp = tensor.shape[2]

        displ = self.restartfile.displ

        mean = numpy.zeros((natoms, 3, ncomp))

        for i in range(3) :
            mean[:,i,:] = (tensor[:,2*i,:ncomp] + tensor[:,2*i+1,:ncomp]) / 2.0 

        return mean

    def get_tensor_deriv_c (self, tens, ncomp=None) :

        if hasattr(self, tens+'_deriv_c') :
            return eval('self.'+tens+'_deriv_c')

        if not (self.restartfile.iaus == 2) :
            raise Exception('Differentiation only implemented for iaus = 2')

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

        # FIXME: convert to absorption (in km/mol); scale factor stolen from SNF
        irint = irint * 863.865928384
        return irint

    def get_a2_invariant (self, gauge='len', modes=None) :
        pol    = self.get_tensor_deriv_nm('pol'+gauge, ncomp=6, modes=modes)
        nmodes = pol.shape[0]

        a2 = (1.0/3.0) * (pol[:,0] + pol[:,3] + pol[:,5])
        a2 = (a2**2)*(Constants.Bohr_in_Angstrom**4)
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
        g2 = g2 * (Constants.Bohr_in_Angstrom**4)
        return g2

    def get_raman_int (self, modes=None) :
        return 45.0*self.get_a2_invariant(modes=modes) + 7.0*self.get_g2_invariant(modes=modes)

    def get_aG_invariant (self, gauge='len', modes=None) :
        # for consistency with SNF alpha is always in length repr
        pol     = self.get_tensor_deriv_nm('pollen', ncomp=6, modes=modes)
        gten    = self.get_tensor_deriv_nm('gten'+gauge, modes=modes)

        aG = (1.0/9.0) * (pol[:,0] + pol[:,3] + pol[:,5])   \
                       * (gten[:,0] + gten[:,4] + gten[:,8])
        aG = aG*(Constants.Bohr_in_Angstrom**4) * (1 / Constants.cvel) * 1e6

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
        bG = bG*(Constants.Bohr_in_Angstrom**4) * (1 / Constants.cvel) * 1e6

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
        bA = bA*(Constants.Bohr_in_Angstrom**4) * (1 / Constants.cvel) * 1e6

        return bA

    def get_backscattering_int (self, modes=None) :
         return 1e-6*96.0*(self.get_bG_invariant(gauge='vel', modes=modes) + \
                           (1.0/3.0)*self.get_bA_invariant(modes=modes))

    def check_consistency (self) :
        res = SNFResults()
        res.read()

        print
        print "Checking consistency of SNF restart file "
        print "  Large numbers mean that the change of the polarizabilities deviates from linear behavior."
        print

        def check (ten) :
            if ten.startswith('pol') :
                mean_pol = res.get_tensor_mean(ten, 6)
            else:
                mean_pol = res.get_tensor_mean(ten)
            natoms = mean_pol.shape[0]
            ncomp = mean_pol.shape[2]

            mean_pol = mean_pol.reshape((natoms*3,ncomp))
            for icomp in range(ncomp) :
                print " %14.8f  min: %3i  max:  %3i" % \
                      (mean_pol[:,icomp].max()-mean_pol[:,icomp].min(),
                       mean_pol[:,icomp].argmin(), mean_pol[:,icomp].argmax() )
            print

        print " Dipole moment "
        check('dipole')

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
