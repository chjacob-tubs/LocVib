########
Examples
########

See `/LocVib/example/` directory for some examples of typical runs.


.. code-block:: console

   /LocVib/example/
   ├── 1_composition.py
   ├── 2_locmodes.py
   ├── 3_couplings.py
   ├── Ala10
   │   ├── coord
   │   ├── restart
   │   └── snf.out

As a typical example of a vibrational spectra calculation in terms of localized modes, we take the results of an **SNF** calculation (https://reiher.ethz.ch/software/snf.html) of an *Ala10* molecule (harmonic approximation).

Composition
===========

Here we see the script for the first part of the example `1_composition.py`:

.. code-block:: python

   import VibTools
   
   from numpy import *
   
   def main ():
       res = VibTools.SNFResults(outname='Ala10/snf.out',
                                 restartname='Ala10/restart',
                                 coordfile='Ala10/coord')
       res.read()
   
       # change here to range of modes which are possibly of interest
       ml = range(140,251)
   
       modes = res.modes.get_subset(ml)
   
       irint = res.get_ir_intensity(modes)
       ri    = res.get_raman_int(modes)
       bi    = res.get_backscattering_int(modes)
   
       labels = ["%i %6.1f %10.4f %10.4f %10.4f    " % (i, f, ir, r, b)
                 for i, f, ir, r, b in zip(ml, modes.freqs, irint, ri, bi)]
   
       modes.print_attype2_composition(labels)
   
   main()

When we run the script we get:

   >>> /LocVib/example/python3 1_composition.py

.. code-block:: console
   
                                                      NH      CO      CA      HA     CHB   
   140 1003.4     0.8322     0.4861     0.0018       8.4     4.6    23.3    23.2    40.6   
   141 1010.1     2.1909     0.2414     0.0018       9.8     5.5    22.6    21.8    40.3   
   142 1018.1     1.2079     0.3277     0.0017      11.1     6.7    21.8    20.0    40.3   
   143 1025.9     8.8684     0.8042    -0.0039      11.7     8.1    20.9    18.0    41.3   
   144 1032.4     4.1386     0.2587     0.0158      12.0     9.4    19.3    16.7    42.5   
   145 1038.2    28.0831     2.2072     0.0140      12.6     9.1    17.8    16.6    43.8   
   146 1044.3    41.6849     1.4998    -0.0598      13.9     7.6    17.3    16.5    44.7   
   147 1049.1    10.7513     0.7006    -0.0094      14.5     6.7    16.9    16.1    45.8   
   148 1052.4     0.8579     0.2134     0.0001      13.6     6.4    16.6    15.7    47.8   
   149 1054.3     2.6367     0.2466     0.0017       8.9     7.5    16.7    15.3    51.6   
   .
   .
   .
   238 1506.4   217.9075     7.7164    -0.0585      75.1    16.3     2.8     2.7     3.2   
   239 1508.5   407.3909     6.8705    -0.0143      76.4    18.1     2.3     1.1     2.2   
   240 1617.6    53.4093     8.3601    -0.0152      98.4     0.3     0.8     0.2     0.3   
   241 1644.0    22.0081     1.8999     0.0294       3.0    94.1     0.9     1.5     0.5   
   242 1650.0   322.4753    17.5984     0.0226       3.3    94.2     0.6     1.4     0.5   
   243 1651.1  1361.5743   112.2062     0.0108       5.1    93.3     0.2     0.8     0.5   
   244 1656.1   232.1278     2.5841    -0.0374       3.5    94.3     0.5     1.3     0.4   
   245 1657.2    41.6081     1.2601    -0.0010       4.4    93.9     0.3     1.0     0.4   
   246 1663.7    86.5379     1.8897    -0.0026       3.3    94.5     0.6     1.1     0.4   
   247 1668.3   247.8201     8.6080     0.0088       2.4    96.2     0.5     0.4     0.5   
   248 1670.0   290.1721    24.3792     0.0023       4.0    94.4     0.3     0.9     0.4   
   249 1677.3   217.4265     9.0573     0.0112       3.0    95.6     0.4     0.5     0.4   
   250 1736.0   326.1756    11.9512     0.0017       2.5    96.7     0.4     0.0     0.3   
   
   min:                                              0.3     0.2     0.1     0.0     0.3   
   max:                                             98.4    96.7    40.1    74.1    98.8 


locmodes
========

Here we see the script for the first part of the example `2_locmodes.py`:

.. code-block:: python

   import VibTools

   from numpy import *
   
   def main () :
       res = VibTools.SNFResults(outname='Ala10/snf.out',
                                 restartname='Ala10/restart',
                                 coordfile='Ala10/coord')
       res.read()
   
       # numbers of amide I normal modes (found with composition.py)
       ml = range(241,250)
   
       modes = res.modes.get_subset(ml)
   
       backint = res.get_backscattering_int(modes=modes)
   
       lv = VibTools.LocVib(modes, 'PM')
       lv.localize()
       lv.sort_by_residue()
       lv.adjust_signs()
   
       print()
       print("Composition of localized modes: ")
       lv.locmodes.print_attype2_composition()
       lv.locmodes.print_residue_composition()
   
       # write localized modes to g98-style file
       lv.locmodes.write_g98out(filename="locmodes-amide1.out")
   
       backint_loc = res.get_backscattering_int(modes=lv.locmodes)
   
       print()
       print("ROA backscattering intensities of normal modes and localized modes")
   
       print()
       print("  normal modes        |  localized modes ")
       print("     freq       int   |    freq       int")
   
       n = 0
       for i, f, bi, f_loc, bi_loc in zip(ml, modes.freqs, backint*1e3,
                                              lv.locmodes.freqs, backint_loc*1e3) :
           n = n+1
           print(r"%i  %6.1f  %8.2f | %2i %6.1f %8.2f" % (i, f, bi, n, f_loc, bi_loc))
   
       print()
       print("Total ROA Int    : ", backint.sum() * 1e3)
   
   main()

When we run the script we get:

   >>> /LocVib/example/python3 2_locmodes.py

.. code-block:: console
   
    Normal mode localization: Cycle   1    p:    4.501   change:  2.3145195     6.41573 
    Normal mode localization: Cycle   2    p:    4.565   change:  0.0635835     0.57591 
    Normal mode localization: Cycle   3    p:    4.565   change:  0.0000638     0.00852 
    Normal mode localization: Cycle   4    p:    4.565   change:  0.0000000     0.00000 
   Obtaining coupling matrix in [a.u.]
   Obtaining coupling matrix in [a.u.]
   Obtaining coupling matrix in [a.u.]
   Obtaining coupling matrix in [a.u.]
   Obtaining coupling matrix in [a.u.]
   Obtaining coupling matrix in [a.u.]
   Obtaining coupling matrix in [a.u.]
   Obtaining coupling matrix in [a.u.]
   
   Composition of localized modes: 
                       NH      CO      CA      HA     CHB   
      0   1663.7      3.6    94.5     0.6     0.8     0.5   
      1   1667.1      3.4    94.5     0.5     1.2     0.4   
      2   1654.1      3.9    93.9     0.5     1.2     0.5   
      3   1651.6      3.8    94.0     0.5     1.2     0.5   
      4   1651.7      3.9    94.0     0.5     1.2     0.5   
      5   1652.2      3.8    94.1     0.5     1.1     0.4   
      6   1652.2      4.2    93.6     0.5     1.2     0.5   
      7   1676.9      2.9    95.7     0.4     0.5     0.4   
      8   1668.2      2.6    96.0     0.5     0.4     0.5   
   
   min:               2.6    93.6     0.4     0.4     0.4   
   max:               4.2    96.0     0.6     1.2     0.5   
   
                        1      10       2       3       4       5       6       7       8       9   
      0   1663.7     95.5     0.0     3.6     0.5     0.3     0.0     0.0     0.0     0.0     0.0   
      1   1667.1      0.0     0.0    95.9     3.2     0.5     0.3     0.0     0.0     0.0     0.0   
      2   1654.1      0.0     0.0     0.0    95.2     3.8     0.6     0.3     0.0     0.0     0.0   
      3   1651.6      0.0     0.0     0.0     0.0    95.3     3.6     0.7     0.3     0.1     0.0   
      4   1651.7      0.0     0.0     0.0     0.0     0.0    95.3     3.6     0.6     0.4     0.1   
      5   1652.2      0.0     0.1     0.0     0.0     0.0     0.0    95.4     3.6     0.6     0.4   
      6   1652.2      0.0     0.7     0.0     0.0     0.0     0.0     0.0    94.9     3.7     0.8   
      7   1676.9      0.0     0.4     0.0     0.0     0.0     0.0     0.0     0.0    96.5     3.1   
      8   1668.2      0.0     3.2     0.0     0.0     0.0     0.0     0.0     0.0     0.0    96.8   
   
   min:               0.0     0.0     0.0     0.0     0.0     0.0     0.0     0.0     0.0     0.0   
   max:              95.5     3.2    95.9    95.2    95.3    95.3    95.4    94.9    96.5    96.8   
   
   
   ROA backscattering intensities of normal modes and localized modes
   
     normal modes        |  localized modes 
        freq       int   |    freq       int
   241  1644.0     29.38 |  1 1663.7     4.81
   242  1650.0     22.59 |  2 1667.1    12.55
   243  1651.1     10.84 |  3 1654.1     0.54
   244  1656.1    -37.40 |  4 1651.6    -0.25
   245  1657.2     -0.96 |  5 1651.7     5.95
   246  1663.7     -2.65 |  6 1652.2    -5.54
   247  1668.3      8.81 |  7 1652.2     1.44
   248  1670.0      2.27 |  8 1676.9    13.08
   249  1677.3     11.15 |  9 1668.2    11.46
   
   Total ROA Int    :  44.02970684750486   


couplings
=========

Here we see the script for the first part of the example `3_couplings.py`:

.. code-block:: python
   
   import VibTools

   from numpy import *
   
   def print_mat (a) :
       print
       for k in a :
           for l in k :
               print("%6.1f" % l,)
           print()
   
   def main () :
       res = VibTools.SNFResults(outname='Ala10/snf.out',
                                 restartname='Ala10/restart',
                                 coordfile='Ala10/coord')
       res.read()
   
       plots = []
   
       # numbers of the amide I modes
       ml = range(241,250)
   
       modes = res.modes.get_subset(ml)
   
       lv = VibTools.LocVib(modes, 'PM')
       lv.localize()
       lv.sort_by_residue()
       lv.adjust_signs()
   
       # vibrational coupling constants
       cmat = lv.get_couplingmat()
   
       print()
       print("Vibrational coupling matrix [in cm-1]")
       print_mat(cmat)
   
       # intensity coupling constants
       lma = VibTools.LocModeAnalysis(res, 'ROA', lv.locmodes)
       intcmat = lma.get_intensity_coupling_matrix()
   
       print()
       print("Intensity coupling matrix for ROA backscattering")
       print_mat(1000.0*intcmat)
   
   main()


When we run the script we get:

   >>> /LocVib/example/python3 couplings.py

.. code-block:: console
   
    Normal mode localization: Cycle   1    p:    4.501   change:  2.3145195     6.41573 
    Normal mode localization: Cycle   2    p:    4.565   change:  0.0635835     0.57591 
    Normal mode localization: Cycle   3    p:    4.565   change:  0.0000638     0.00852 
    Normal mode localization: Cycle   4    p:    4.565   change:  0.0000000     0.00000 
   Obtaining coupling matrix in [a.u.]
   Obtaining coupling matrix in [a.u.]
   Obtaining coupling matrix in [a.u.]
   Obtaining coupling matrix in [a.u.]
   Obtaining coupling matrix in [a.u.]
   Obtaining coupling matrix in [a.u.]
   Obtaining coupling matrix in [a.u.]
   Obtaining coupling matrix in [a.u.]
   Obtaining coupling matrix in [a.u.]
   
   Vibrational coupling matrix [in cm-1]
   1663.7
      3.6
     -3.7
     -0.8
     -0.7
     -0.5
     -0.3
     -0.2
     -0.2
   
     .
     .
     .
      
     -0.2
     -0.2
     -0.4
     -0.5
     -0.3
     -2.5
      1.0
   1676.9
      1.0
   
     -0.2
     -0.2
     -0.1
     -0.4
     -0.7
     -0.6
     -1.4
      1.0
   1668.2
   
   
   Intensity coupling matrix for ROA backscattering
      4.8
     -9.8
     42.1
    -24.6
     -1.4
     76.0
   -109.7
    -68.2
     57.6
   
      0.0
     12.6
      1.8
     19.1
     11.2
    -27.2
     83.0
    -92.6
    -21.5
   
    .
    .
    .
   
      0.0
      0.0
      0.0
      0.0
      0.0
      0.0
      1.4
    -75.7
     40.1
   
      0.0
      0.0
      0.0
      0.0
      0.0
      0.0
      0.0
     13.1
     -3.1
   
      0.0
      0.0
      0.0
      0.0
      0.0
      0.0
      0.0
      0.0
     11.5
