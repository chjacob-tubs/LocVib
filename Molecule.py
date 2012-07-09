
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

        z = zip(groupnames, groups)
        z.sort(lambda x,y: cmp(int(x[0]), int(y[0])))
        groupnames = [x[0] for x in z]
        groups     = [x[1] for x in z]

        return groups, groupnames

    def attype_groups (self) :
        groupnames = ['N', 'H', 'C', 'O',  'CA', 'HA', 'CB', 'HB', 'OXT', 'HXT', 'H2O']
        groups     = []

        for k in groupnames :
            groups.append([])
        for at in openbabel.OBMolAtomIter(self.obmol) :
            attype = at.GetResidue().GetAtomID(at)
            attype = attype.strip()

            if attype in ['1HB', '2HB', '3HB', 'HB1', 'HB2', 'HB3'] :
                attype = 'HB'
            if attype in ['1H', '2H', 'H 1', 'H 2', 'H1', 'H2'] :
                attype = 'HXT'

            if at.GetResidue().GetName() == 'HOH' :
	        attype = 'H2O'

            ind = groupnames.index(attype)
            groups[ind].append(at.GetIdx()-1)

        return groups, groupnames

    def attype_groups_2 (self) :
        groupnames = ['NH', 'CO', 'CA', 'HA', 'CHB']
        groups     = []

        for k in groupnames :
            groups.append([])
        for at in openbabel.OBMolAtomIter(self.obmol) :
            attype = at.GetResidue().GetAtomID(at)
            attype = attype.strip()

            if attype in ['1HB', '2HB', '3HB', 'HB1', 'HB2', 'HB3'] :
                attype = 'HB'
            if attype in ['1H', '2H', 'H 1', 'H 2', 'H1', 'H2'] :
                attype = 'HXT'

            if attype in ['N', 'H', 'HXT'] :
                attype = 'NH'
            elif attype in ['C', 'O', 'OXT'] :
                attype = 'CO'
            elif attype in ['CB', 'HB'] :
                attype = 'CHB'

            ind = groupnames.index(attype)
            groups[ind].append(at.GetIdx()-1)

        return groups, groupnames

    def attype_groups_3 (self) :
        groupnames = ['NH', 'CO', 'CA', 'CHB1', 'CHB2', 'CHG1', 'CHD1']
        groups     = []

        for k in groupnames :
            groups.append([])
        for at in openbabel.OBMolAtomIter(self.obmol) :
            attype = at.GetResidue().GetAtomID(at)
            attype = attype.strip()

            if attype.startswith('HB1') : 
                attype = 'HB1'
            elif attype.startswith('HB2') : 
                attype = 'HB2'
            elif attype.startswith('HG1') : 
                attype = 'HG1'
            elif attype.startswith('HD1') : 
                attype = 'HD1'
            if attype in ['1H', '2H', 'H 1', 'H 2', 'H1', 'H2'] :
                attype = 'HXT'

            if attype in ['N', 'H', 'HXT'] :
                attype = 'NH'
            elif attype in ['C', 'O', 'OXT'] :
                attype = 'CO'
            elif attype in ['CB1', 'HB1'] :
                attype = 'CHB1'
            elif attype in ['CB2', 'HB2'] :
                attype = 'CHB2'
            elif attype in ['CG1', 'HG1'] :
                attype = 'CHG1'
            elif attype in ['CD1', 'HD1'] :
                attype = 'CHD1'

            ind = groupnames.index(attype)
            groups[ind].append(at.GetIdx()-1)

        return groups, groupnames

    def attype_groups_7B (self) :
        groupnames = ['NH', 'CO', 'CHA', 'CHB1', 'CHB2', 'CHG1', 'CHD1', 'Ring', 'Term']
        groups     = []

        for k in groupnames :
            groups.append([])
        for at in openbabel.OBMolAtomIter(self.obmol) :
            attype = at.GetResidue().GetAtomID(at)
            attype = attype.strip()
            resname = at.GetResidue().GetName()

            if attype.startswith('HA') : 
                attype = 'CHA'
            elif attype.startswith('HB1') : 
                attype = 'HB1'
            elif attype.startswith('HB2') : 
                attype = 'HB2'
            elif attype.startswith('HG1') : 
                attype = 'HG1'
            elif attype.startswith('HD1') : 
                attype = 'HD1'
            if attype in ['1H', '2H', 'H 1', 'H 2', 'H1', 'H2'] :
                attype = 'HXT'

            if attype in ['N', 'H', 'HXT'] :
                attype = 'NH'
            elif attype in ['C', 'O', 'OXT'] :
                attype = 'CO'
            elif attype in ['CA'] :
                attype = 'CHA'
            elif attype in ['CB', 'CB1', 'HB1'] :
                attype = 'CHB1'
            elif attype in ['CB2', 'HB2'] :
                attype = 'CHB2'
            elif attype in ['CG1', 'HG1'] :
                attype = 'CHG1'
            elif attype in ['CD1', 'HD1'] :
                attype = 'CHD1'

            if resname == 'BLA' :
                if not attype in ['CHA', 'CHB1', 'CHB2', 'NH', 'CO'] :
                   attype = 'Ring'
            if resname == 'UNK' :
                attype = 'Term'

            ind = groupnames.index(attype)
            groups[ind].append(at.GetIdx()-1)

        return groups, groupnames


    def atom_groups (self) :
        group_dict = {}

        pse = openbabel.OBElementTable()
        
        for at in openbabel.OBMolAtomIter(self.obmol) :
            atname = pse.GetSymbol(at.GetAtomicNum())
            if atname in group_dict :
                group_dict[atname].append(at.GetIdx()-1)
            else :
                group_dict[atname] = [at.GetIdx()-1]

        groupnames = []
        groups = []

        for g in group_dict :
            groupnames.append(g)
            groups.append(group_dict[g])

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

        groupnames = [ 'CO-1', 'NH', 'CO', 'NH+1', 'CO+1', 'NH+2', 'CH01', 'R']
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
                    attype = 'CH01'
            elif rnum == res + 1 :
                if attype in ['N', 'H'] :
                    attype = 'NH+1'
                elif attype in ['C', 'O'] :
                    attype = 'CO+1'
                elif attype in ['CA', 'HA', 'CB', 'HB'] :
                    attype = 'CH01'
                else :
                    attype = 'R'
            elif rnum == res + 2 :
                if attype in ['N', 'H'] :
                    attype = 'NH+2'
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


    def amide3_groups (self, res) :

        groupnames = [ 'Am-1/0', 'CA0', 'Am0/1', 'CA+1', 'Am1/2', 'CA+2', 'CB', 'R']
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
                    attype = 'Am-1/0'
                elif attype in ['C', 'O'] :
                    attype = 'Am0/1'
                elif attype in ['CA', 'HA'] :
                    attype = 'CA0'
                elif attype in ['CB', 'HB'] :
                    attype = 'CB'
            elif rnum == res + 1 :
                if attype in ['N', 'H'] :
                    attype = 'Am0/1'
                elif attype in ['C', 'O'] :
                    attype = 'Am1/2'
                elif attype in ['CA', 'HA'] :
                    attype = 'CA+1'
                elif attype in ['CB', 'HB'] :
                    attype = 'CB'
            elif rnum == res + 2 :
                if attype in ['N', 'H'] :
                    attype = 'Am1/2'
                elif attype in ['CA', 'HA'] :
                    attype = 'CA+2'
                elif attype in ['CB', 'HB'] :
                    attype = 'CB'
                else :
                    attype = 'R'
            elif rnum == res - 1 :
                if attype in ['C', 'O'] :
                    attype = 'Am-1/0'
                else :
                    attype = 'R'
            else :
                attype = 'R'

            try:
                ind = groupnames.index(attype)
            except ValueError :
                print attype
                raise
            groups[ind].append(at.GetIdx()-1)

        return groups, groupnames

    def skelCN_groups (self, res) :

        groupnames = ['NCA0', 'HA0', 'CO0',
                      'H+1', 'NCA+1', 'HA+1', 'CHB+1', 'CO+1', 
                      'H+2', 'NCA+2',
                      'R']
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
                if attype in ['C', 'O'] :
                    attype = 'CO0'
                elif attype in ['N','CA'] :
                    attype = 'NCA0'
                elif attype in ['HA'] :
                    attype = 'HA0'
                else:
                    attype = 'R'
            elif rnum == res + 1 :
                if attype in ['H'] :
                    attype = 'H+1'
                elif attype in ['C', 'O'] :
                    attype = 'CO+1'
                elif attype in ['N', 'CA'] :
                    attype = 'NCA+1'
                elif attype in ['HA'] :
                    attype = 'HA+1'
                elif attype in ['CB', 'HB'] :
                    attype = 'CHB+1'
            elif rnum == res + 2 :
                if attype in ['H'] :
                    attype = 'H+2'
                elif attype in ['N','CA'] :
                    attype = 'NCA+2'
                else:
                    attype = 'R'
            else :
                attype = 'R'

            try:
                ind = groupnames.index(attype)
            except ValueError :
                print attype
                raise
            groups[ind].append(at.GetIdx()-1)

        return groups, groupnames


    def ch3_groups (self, res) :

        groupnames = ['CO0', 'NH+1', 'CA+1', 'CB+1', 'CO+1', 'NH+2', 'R']
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
                if attype in ['C', 'O'] :
                    attype = 'CO0'
                else :
                    attype = 'R'
            elif rnum == res + 1 :
                if attype in ['N', 'H'] :
                    attype = 'NH+1'
                elif attype in ['C', 'O'] :
                    attype = 'CO+1'
                elif attype in ['CA', 'HA'] :
                    attype = 'CA+1'
                elif attype in ['CB', 'HB'] :
                    attype = 'CB+1'
            elif rnum == res + 2 :
                if attype in ['N', 'H'] :
                    attype = 'NH+2'
                else :
                    attype = 'R'
            else :
                attype = 'R'

            try:
                ind = groupnames.index(attype)
            except ValueError :
                print attype
                raise
            groups[ind].append(at.GetIdx()-1)

        return groups, groupnames




