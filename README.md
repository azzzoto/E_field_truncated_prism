# E_field_truncated_prism
Computation of the Electric field and potential of a truncated prism using Jacobi's method.
The output of the program is divided in four .dat files depending whether the reduced algorithm is used (that works thanks to the symmetries of the prism) or not.

E and V can be plotted using the dat files, for example (using Gnuplot):


$ splot 'vtot.dat'


$ plot 'etot.dat' u vec

![efield](https://user-images.githubusercontent.com/58179857/117648216-3e157c80-b18e-11eb-92f7-a90a81bce9fe.png)
![potenziale](https://user-images.githubusercontent.com/58179857/117648235-44a3f400-b18e-11eb-9e63-7a9e6e062a14.png)
