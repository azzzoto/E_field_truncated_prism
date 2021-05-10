# E_field_truncated_prism
Computation of the Electric field and potential of a truncated prism using Jacobi's method.
The output of the program is divided in four .dat files depending whether the reduced algorithm is used (that works thanks to the symmetries of the prism) or not.

E and V can be plotted using the dat files, for example (using Gnuplot):


$ splot 'vtot.dat'


$ plot 'etot.dat' u vec
