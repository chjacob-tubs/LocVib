
from Molecule import VibToolsMolecule
from Modes    import VibModes

class AKIRAResults (object) :

    def __init__ (self, iterationsfile='akira_iterations.out', coordfile='coord') :
        self.mol          = VibToolsMolecule()

        self.iterationsfile = iterationsfile
        self.coordfile      = coordfile

        self.modes  = None
        self.irints = None

    freqs = property(lambda self: self.modes.freqs)

    def read (self) :
        import numpy, re

        self.mol.read_from_coord(filename=self.coordfile)

        #read akira iterations
        f = open(self.iterationsfile)
        lines = f.readlines()
        f.close()
        for i, l in enumerate(lines) :
            if re.search('Final results from module Davidson', l) :
                start = i
        nbas = int(lines[start+3].split()[-1])
        lines = lines[start+8:start+8+nbas]
        modes = []
        for l in lines : 
            spl = l.split()
            if spl[10] == 'YES' :
                modes.append([float(spl[2]),float(spl[8])])

        self.modes = VibModes(len(modes), self.mol)
        self.modes.set_freqs(numpy.array([m[0] for m in modes]))
        self.irints = numpy.array([m[1] for m in modes])

    def get_ir_intensity (self) :
        return self.irints        
