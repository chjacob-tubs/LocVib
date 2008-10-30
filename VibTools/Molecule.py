
import openbabel
import numpy

class AbstractMolecule (object) :

    def __init__ (self) :
        self.natoms      = None
        self.atmasses    = None
        self.atnums      = None
        self.coordinates = None
    
    def get_fragment (self, atomlist) :
        frag = AbstractMolecule()

        frag.natoms      = len(atomlist)
        frag.atmasses    = self.atmasses[atomlist]
        frag.atnums      = self.atnums[atomlist]
        frag.coordinates = self.coordinates[atomlist,:]

        return frag

class VibToolsMolecule (AbstractMolecule) :
    
    def __init__ (self) :
        self.obmol = openbabel.OBMol()

    def read_from_coord (self, filename='coord') :
        obconv = openbabel.OBConversion()
        obconv.SetInAndOutFormats("tmol", "tmol")
        
        if not obconv.ReadFile(self.obmol, filename) :
            raise Exception('Error reading file')

    def get_natoms (self) :
        return self.obmol.NumAtoms()
    natoms = property(get_natoms)
        
    def get_atmasses (self) :
        atmasses = numpy.zeros((self.natoms,))
        for at in openbabel.OBMolAtomIter(self.obmol) :
            i = at.GetIdx()-1
            atmasses[i] = at.GetExactMass()
        return atmasses
    atmasses = property(get_atmasses)

    def get_atnums (self) :
        atnums = numpy.zeros((self.natoms,))
        for at in openbabel.OBMolAtomIter(self.obmol) :
            i = at.GetIdx()-1
            atnums[i] = at.GetAtomicNum()
        return atnums
    atnums = property(get_atnums)

    def get_coordinates (self) :
        coords = numpy.zeros((self.natoms,3))
        for at in openbabel.OBMolAtomIter(self.obmol) :
            i = at.GetIdx()-1
            coords[i,0] = at.GetX()
            coords[i,1] = at.GetY()
            coords[i,2] = at.GetZ()
        return coords

    coordinates = property(get_coordinates)

    ## the following methods rely on Openbabel
    
    def residue_groups (self) :
        groupnames = []
        groups     = []

        # force residue perception
        self.obmol.GetAtom(1).GetResidue()

        for r in openbabel.OBResidueIter(self.obmol) :
            resgroup = []
            for at in openbabel.OBResidueAtomIter(r) :
                resgroup.append(at.GetIdx()-1)
            groupnames.append(str(r.GetNum()))
            groups.append(resgroup)

        return groups, groupnames

    def attype_groups (self) :
        groupnames = ['N', 'H', 'C', 'O',  'CA', 'HA', 'CB', 'HB', 'OXT', 'HXT']
        groups     = []

        for k in groupnames :
            groups.append([])
        for at in openbabel.OBMolAtomIter(self.obmol) :
            attype = at.GetResidue().GetAtomID(at)
            attype = attype.strip()

            if attype in ['1HB', '2HB', '3HB', 'HB1', 'HB2', 'HB3'] :
                attype = 'HB'
            if attype in ['1H', '2H', 'H 1', 'H 2'] :
                attype = 'HXT'

            ind = groupnames.index(attype)
            groups[ind].append(at.GetIdx()-1)

        return groups, groupnames

    def attype_groups_2 (self) :
        groupnames = ['NH', 'CO', 'CHA', 'CHB']
        groups     = []

        for k in groupnames :
            groups.append([])
        for at in openbabel.OBMolAtomIter(self.obmol) :
            attype = at.GetResidue().GetAtomID(at)
            attype = attype.strip()

            if attype in ['1HB', '2HB', '3HB', 'HB1', 'HB2', 'HB3'] :
                attype = 'HB'
            if attype in ['1H', '2H', 'H 1', 'H 2'] :
                attype = 'HXT'

            if attype in ['N', 'H'] :
                attype = 'NH'
            elif attype in ['C', 'O', 'OXT', 'HXT'] :
                attype = 'CO'
            elif attype in ['CA', 'HA'] :
                attype = 'CHA'
            elif attype in ['CB', 'HB'] :
                attype = 'CHB'

            ind = groupnames.index(attype)
            groups[ind].append(at.GetIdx()-1)

        return groups, groupnames

    def amide1_groups (self, res) :
        groupnames = ['CO', 'NH', 'NH2', 'NH3', 'CAH', 'CA0', 'CA2', 'COo', 'CBa', 'R']
        groups     = []

        for k in groupnames :
            groups.append([])
        for at in openbabel.OBMolAtomIter(self.obmol) :
            attype = at.GetResidue().GetAtomID(at)
            attype = attype.strip()

            rnum = at.GetResidue().GetNum()

            if (rnum==res) and (attype in ['C','O']) :
                typ = 'CO'
            elif (rnum==res+1) and (attype in ['N','H']) :
                typ = 'NH'
            elif (rnum==res+1) and (attype in ['CA','HA']) :
                typ ='CAH'
            elif (rnum==res+2) and attype in ['N','H'] :
                typ = 'NH2'
            elif (rnum==res+3) and attype in ['N','H'] :
                typ = 'NH3'
            elif (rnum==res) and attype in ['CA','HA'] :
                typ = 'CA0'
            elif (rnum==res+2) and attype in ['CA','HA'] :
                typ = 'CA2'
            elif attype in ['C','O'] :
                typ = 'COo'
            elif attype in ['CB','1HB','2HB','3HB'] :
                typ = 'CBa'
            elif attype in ['CB','HB'] :
                typ = 'CBo'
            else :
                typ = 'R'

            ind = groupnames.index(typ)
            groups[ind].append(at.GetIdx()-1)

        return groups, groupnames

    def amide2_groups (self, res) :

        groupnames = [ 'CO-1', 'NH', 'CHAB', 'CO', 'NH+1', 'R']
        groups     = []

        for k in groupnames :
            groups.append([])
        for at in openbabel.OBMolAtomIter(self.obmol) :
            attype = at.GetResidue().GetAtomID(at)
            attype = attype.strip()

            rnum = at.GetResidue().GetNum()

            if attype in ['1HB', '2HB', '3HB', 'HB1', 'HB2', 'HB3'] :
                attype = 'HB'

            if rnum == res :
                if attype in ['N', 'H'] :
                    attype = 'NH'
                elif attype in ['C', 'O'] :
                    attype = 'CO'
                elif attype in ['CA', 'HA', 'CB', 'HB'] :
                    attype = 'CHAB'
            elif rnum == res + 1 :
                if attype in ['N', 'H'] :
                    attype = 'NH+1'
                else :
                    attype = 'R'
            elif rnum == res -1 :
                if attype in ['C', 'O'] :
                    attype = 'CO-1'
                else :
                    attype = 'R'
            else :
                attype = 'R'

            ind = groupnames.index(attype)
            groups[ind].append(at.GetIdx()-1)

        return groups, groupnames



