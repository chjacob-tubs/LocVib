import numpy as np
import VibTools as vt
import pytest





@pytest.fixture
def HugAna_Ala10_SNFsetup():
    path = 'test_data/Ala10/'
    res = vt.SNFResults(outname = path+'snf2.out', 
                        restartname = path+'restart2', 
                        coordfile = path+'coord2')
    res.read()
    a = range(0,len(res.modes.freqs))
    modes = res.modes.get_subset(a)
    HA = vt.HugAnalysis(res = res, tensor = 'IR')
    return res, modes, HA

def test_get_irint_decomposed_c(HugAna_Ala10_SNFsetup):
    res, modes, HA = HugAna_Ala10_SNFsetup
    mu = HA.get_irint_decomposed_c(res)
    refs = [10.273816239327385, 1.904448118426103]
    np.testing.assert_almost_equal(refs[0], mu[0,0],decimal=5)
    np.testing.assert_almost_equal(refs[1], mu.sum(),decimal=5)


def test_get_a2_decomposed_c(HugAna_Ala10_SNFsetup):
    res, modes, HA = HugAna_Ala10_SNFsetup
    a2 = HA.get_a2_decomposed_c(res)
    refs = [0.00109442551837648, 0.0004560748234168077]
    np.testing.assert_almost_equal(refs[0], a2[0,0],decimal=5)
    np.testing.assert_almost_equal(refs[1], a2.sum(),decimal=5)


def test_get_g2_decomposed_c(HugAna_Ala10_SNFsetup):
    res, modes, HA = HugAna_Ala10_SNFsetup
    g2 = HA.get_g2_decomposed_c(res)
    refs = [6.2416879144229, 0.001381634299067791]
    np.testing.assert_almost_equal(refs[0], g2[0,0],decimal=5)
    np.testing.assert_almost_equal(refs[1], g2.sum(),decimal=5)



def test_get_ramanint_decomposed_c(HugAna_Ala10_SNFsetup):
    res, modes, HA = HugAna_Ala10_SNFsetup
    mat = HA.get_ramanint_decomposed_c(res)
    refs = [43.741064549287245, 0.03019480714634426]
    np.testing.assert_almost_equal(refs[0], mat[0,0],decimal=5)
    np.testing.assert_almost_equal(refs[1], mat.sum(),decimal=5)



def test_get_aG_decomposed_c(HugAna_Ala10_SNFsetup):
    res, modes, HA = HugAna_Ala10_SNFsetup
    aG = HA.get_aG_decomposed_c(res)
    refs = [-0.8987664257501431, 0.9875102980278143]
    np.testing.assert_almost_equal(refs[0], aG[0,0],decimal=5)
    np.testing.assert_almost_equal(refs[1], aG.sum(),decimal=5)



def test_get_bG_decomposed_c(HugAna_Ala10_SNFsetup):
    res, modes, HA = HugAna_Ala10_SNFsetup
    bG = HA.get_bG_decomposed_c(res)
    refs = [2.845489405438741, -114.0731165362522]
    np.testing.assert_almost_equal(refs[0], bG[0,0],decimal=5)
    np.testing.assert_almost_equal(refs[1], bG.sum(),decimal=5)


def test_get_bA_decomposed_c(HugAna_Ala10_SNFsetup):
    res, modes, HA = HugAna_Ala10_SNFsetup
    bA = HA.get_bA_decomposed_c(res)
    refs = [-316.28640980489865, -111.27232171519427]
    np.testing.assert_almost_equal(refs[0], bA[0,0],decimal=5)
    np.testing.assert_almost_equal(refs[1], bA.sum(),decimal=5)



def test_get_backint_decomposed_c(HugAna_Ala10_SNFsetup):
    res, modes, HA = HugAna_Ala10_SNFsetup
    mat = HA.get_backint_decomposed_c(res)
    refs = [0.04793855158983615, -0.03944180412956566]
    np.testing.assert_almost_equal (refs[0], mat[0,0],decimal=5)
    np.testing.assert_almost_equal (refs[1], mat.sum(),decimal=5)



def test_project_on_modes(HugAna_Ala10_SNFsetup):
    res, modes, HA = HugAna_Ala10_SNFsetup
    refs = [-0.005341916285333833, -0.8928122953072393, 11.041252272527917,
            7498.640016759868, 1.781130199229294, 9080.512040619393, 
            0.00680756606311379, 90.97141390024858, 1.5335588285268298, 
            486.4180558926688, -1.005521192317628, -1328.517213065059, 
            -48.36665583186399, -6355.249047300887, -50.0646808902805, 
             2092.500764405779, 0.24440386699116035, 122.54646003899143]
    i = 0
    tens_list = ['ROA', 'Raman', 'IR', 'a2', 'g2', 'aG', 'bG', 'bA', 'id']
    for tensor in tens_list[:]:
        HA = vt.HugAnalysis(res = res, tensor = tensor)
        mat = HA.project_on_modes(modes)
        np.testing.assert_almost_equal(mat[0,0],refs[i],decimal=5)
        np.testing.assert_almost_equal(mat.sum(),refs[i+1],decimal=5)
        i += 2

def test_sum_groups(HugAna_Ala10_SNFsetup):
    res, modes, HA = HugAna_Ala10_SNFsetup
    # unter der Bedingung, dass project_on_modes funktioniert.
    g,gn = modes.mol.attype_groups()
    refs =[-0.07173694551076885, -0.8928122953072388, 85.44838132217726, 7498.640016759864, 827.3440296986114, 9080.512040619396, 0.8194948806523819, 90.97141390024862, 6.93873024183142, 486.41805589266903, -56.84293088816632, -1328.5172130650585, -726.6758454245762, -6355.24904730089, 339.2649533590012, 2092.5007644057764, 1.7167598464905927, 122.54646003899144]
    i = 0
    tens_list = ['ROA', 'Raman', 'IR', 'a2', 'g2', 'aG', 'bG', 'bA', 'id']
    for tensor in tens_list[:]:
        HA = vt.HugAnalysis(res = res, tensor = tensor)
        inv_decomposed_nm = HA.project_on_modes(modes)
        inv_groups = HA.sum_groups(inv_decomposed_nm,g)
        np.testing.assert_almost_equal(inv_groups[0,0],refs[i])
        np.testing.assert_almost_equal(inv_groups.sum(),refs[i+1])
        i += 2


def test_get_group_coupling_matrix(HugAna_Ala10_SNFsetup):
    res, modes, HA = HugAna_Ala10_SNFsetup
    # unter der Bedingung, dass sum_groups funktioniert.
    g,gn = modes.mol.attype_groups()
    tens_list = ['ROA', 'Raman', 'IR', 'a2', 'g2', 'aG', 'bG', 'bA', 'id']
    for tensor in tens_list[:]:
        HA = vt.HugAnalysis(res = res, tensor = tensor)
        inv_decomposed_nm = HA.project_on_modes(modes)
        inv_groups = HA.sum_groups(inv_decomposed_nm,g)
        inv_groups2 = HA.get_group_coupling_matrix(g,modes)
        np.testing.assert_almost_equal(inv_groups.sum(), inv_groups2.sum())



