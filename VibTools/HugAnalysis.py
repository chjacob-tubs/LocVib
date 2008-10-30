
import numpy

import Constants
import Modes

class HugAnalysis (object) :

    def __init__ (self, res, tensor, scale=1.0) :
        self.natoms = res.modes.natoms
        
        if tensor == 'ROA' :
            self.tensor_decomposed_c = self.get_backint_decomposed_c(res)
        elif tensor == 'Raman' :
            self.tensor_decomposed_c = self.get_ramanint_decomposed_c(res)
        elif tensor == 'IR' :
            self.tensor_decomposed_c = self.get_irint_decomposed_c(res)
        elif tensor in ['a2', 'g2', 'aG', 'bG', 'bA'] :
            self.tensor_decomposed_c = eval('self.get_'+tensor+'_decomposed_c(res)' )
        elif tensor == 'id' :
            self.tensor_decomposed_c = numpy.ones((self.natoms*3, self.natoms*3))
            
        self.tensor_decomposed_c = scale * self.tensor_decomposed_c

    def get_irint_decomposed_c (self, res) :
        dip  = res.get_tensor_deriv_c('dipole')
        dip  = dip.reshape((self.natoms*3, 3))

        mu = (numpy.outer(dip[:,0],dip[:,0])  + 
              numpy.outer(dip[:,1],dip[:,1])  +
              numpy.outer(dip[:,2],dip[:,2]))

        # FIXME: convert to absorption (in km/mol); scale factor stolen from SNF
        mu = mu * 863.865928384
        return mu
        
    def get_a2_decomposed_c (self, res, gauge='len') :
        pol  = res.get_tensor_deriv_c('pol'+gauge, 6)

        temp = (1.0/3.0)*(pol[:,:,0] + pol[:,:,3] + pol[:,:,5])
        temp = temp.reshape((self.natoms*3,))

        a2 = numpy.outer(temp[:], temp[:])
        a2 = a2*(Constants.Bohr_in_Angstrom**4)
        
        return a2

    def get_g2_decomposed_c (self, res) :
        pol     = res.get_tensor_deriv_c('pollen', 6)

        pol  = pol.reshape((self.natoms*3,6))

        g2 = 3.0*(  numpy.outer(pol[:,0], pol[:,0]) 
                  + numpy.outer(pol[:,1], pol[:,1])
                  + numpy.outer(pol[:,2], pol[:,2])
                  + numpy.outer(pol[:,1], pol[:,1])
                  + numpy.outer(pol[:,3], pol[:,3])
                  + numpy.outer(pol[:,4], pol[:,4])
                  + numpy.outer(pol[:,2], pol[:,2])
                  + numpy.outer(pol[:,4], pol[:,4])
                  + numpy.outer(pol[:,5], pol[:,5]))
        g2 = g2 - numpy.outer(pol[:,0]+pol[:,3]+pol[:,5], pol[:,0]+pol[:,3]+pol[:,5])

        g2 = 0.5*g2*(Constants.Bohr_in_Angstrom**4) 

        return g2

    def get_ramanint_decomposed_c (self, res) :
        return 45.0*self.get_a2_decomposed_c(res) + 7.0*self.get_g2_decomposed_c(res)

    def get_aG_decomposed_c (self, res, gauge='len') :
        # for consistency with SNF alpha is always in length repr
        pol     = res.get_tensor_deriv_c('pollen', 6)
        gten    = res.get_tensor_deriv_c('gten'+gauge)

        temp1 = (1.0/3.0)*(pol[:,:,0]  + pol[:,:,3]  + pol[:,:,5])
        temp1 = temp1.reshape((self.natoms*3,))
        temp2 = (1.0/3.0)*(gten[:,:,0] + gten[:,:,4] + gten[:,:,8]) 
        temp2 = temp2.reshape((self.natoms*3,))

        aG = numpy.outer(temp1[:], temp2[:])
        aG = aG*(Constants.Bohr_in_Angstrom**4) * (1/Constants.cvel) * 1e6

        return aG

    def get_bG_decomposed_c (self, res, gauge='len') :
        pol     = res.get_tensor_deriv_c('pol'+gauge, 6)
        gten    = res.get_tensor_deriv_c('gten'+gauge)

        pol  = pol.reshape((self.natoms*3,6))
        gten = gten.reshape((self.natoms*3,9))

        bG = 3.0*(  numpy.outer(pol[:,0], gten[:,0]) 
                  + numpy.outer(pol[:,1], gten[:,1])
                  + numpy.outer(pol[:,2], gten[:,2])
                  + numpy.outer(pol[:,1], gten[:,3])
                  + numpy.outer(pol[:,3], gten[:,4])
                  + numpy.outer(pol[:,4], gten[:,5])
                  + numpy.outer(pol[:,2], gten[:,6])
                  + numpy.outer(pol[:,4], gten[:,7])
                  + numpy.outer(pol[:,5], gten[:,8]))
        bG = bG - numpy.outer(pol[:,0]+pol[:,3]+pol[:,5], gten[:,0] + gten[:,4] + gten[:,8] )

        bG = 0.5*bG*(Constants.Bohr_in_Angstrom**4) * (1/Constants.cvel) * 1e6

        return bG

    def get_bA_decomposed_c (self, res) :
        pol     = res.get_tensor_deriv_c('pollen', 6)
        aten    = res.get_tensor_deriv_c('aten')

        pol  = pol.reshape((self.natoms*3,6))
        aten = aten.reshape((self.natoms*3,27))

        bA = (  numpy.outer(pol[:,3]-pol[:,0], aten[:,11])
              + numpy.outer(pol[:,0]-pol[:,5], aten[:,6] )
              + numpy.outer(pol[:,5]-pol[:,3], aten[:,15])
              + numpy.outer(pol[:,1], aten[:,19]-aten[:,20]+aten[:,8]-aten[:,14]) 
              + numpy.outer(pol[:,2], aten[:,25]-aten[:,21]+aten[:,3]-aten[:,4]) 
              + numpy.outer(pol[:,4], aten[:,10]-aten[:,24]+aten[:,12]-aten[:,5])
             )
        
        bA =  0.5 * res.lwl * bA * (Constants.Bohr_in_Angstrom**4) * (1/Constants.cvel) * 1e6

        return bA

    def get_backint_decomposed_c (self, res) :
         return 1e-6*96.0*(self.get_bG_decomposed_c(res, gauge='vel')+(1.0/3.0)*self.get_bA_decomposed_c(res))
      
    def project_on_modes (self, modes, nummode=None) :
        inv_decomposed_nm = numpy.zeros((self.natoms, self.natoms))

        if nummode is None :
            modelist = range(modes.nmodes)
        else:
            modelist = [nummode]

        for imode in modelist :
            # projection on normal modes
            temp = self.tensor_decomposed_c * numpy.outer(modes.modes_c[imode], modes.modes_c[imode])

            # now sum x, y, z-components (in both direction)
            temp = temp.reshape((self.natoms,3,self.natoms,3))
            temp = temp[:,:,:,0] + temp[:,:,:,1] + temp[:,:,:,2]
            temp = temp[:,0,:] + temp[:,1,:] + temp[:,2,:]

            inv_decomposed_nm += temp

        return inv_decomposed_nm
    
    def sum_groups (self, inv_decomposed_nm, groups) :
        ngroups = len(groups)
        inv_groups = numpy.zeros((ngroups,ngroups))

        for ig in range(ngroups) :
            for jg in range(ngroups) :
                for i in groups[ig] :
                    for j in groups[jg] :
                        inv_groups[ig,jg] += inv_decomposed_nm[i,j]

        inv_groups = inv_groups + inv_groups.transpose()

        for i in range(ngroups) :
            inv_groups[i,i] = inv_groups[i,i] / 2
            inv_groups[i,:i] = 0.0

        return inv_groups

    def get_group_coupling_matrix(self, groups, modes, num_mode=None) :
        inv_decomposed_nm = self.project_on_modes(modes, num_mode)
        inv_groups = self.sum_groups(inv_decomposed_nm, groups)  
        return inv_groups
        
    def print_gcm(self, inv_groups, groupnames) :
        for n in groupnames :
            print ("%6s " % n),
        print
        for i in range(inv_groups.shape[0]) :
            print " "*8*i,
            for j in range(i,inv_groups.shape[0]) :
                print "%6.1f " % inv_groups[i,j],
            print

    def print_group_coupling_matrix(self, groups, groupnames, modes, num_mode=None) :
        inv_groups = self.get_group_coupling_matrix(groups, modes, num_mode)
        print
        print "Total intensity: ", inv_groups.sum()
        print
        self.print_gcm(inv_groups, groupnames)

            
##############################################################################################
        
    def artmode_analysis (self, inv, cmodes) :
        nmodes = cmodes.shape[0]

        inv_decomposed_am = numpy.zeros((nmodes, nmodes))

        for imode in range(nmodes) :
            for jmode in range(nmodes) :
                temp = inv * numpy.outer(cmodes[imode],cmodes[jmode])

                inv_decomposed_am[imode, jmode] = temp.sum()

        inv_decomposed_am = inv_decomposed_am + inv_decomposed_am.transpose()

        for i in range(nmodes) :
            inv_decomposed_am[i,i] = inv_decomposed_am[i,i] / 2
            inv_decomposed_am[i,i+1:] = 0.0

        return inv_decomposed_am


