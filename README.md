# MATLAB-SDOT

MATLAB functions for solving semi-discrete optimal transport problems in 3D: 
* ``SDOT_damped_Newton.m`` computes the optimal transport cost between the Lebesgue measure and a discrete measure on a 3D box with respect to the quadratic cost or the periodic quadratic cost, using the damped Newton method from
    * Jun Kitagawa, Quentin Mérigot, Boris Thibert, Convergence of a Newton algorithm for semi-discrete optimal transport. J. Eur. Math. Soc. 21 (2019), no. 9, pp. 2603–2651.
* ``SDOT_fminuc.m`` solves the same problem using the MATLAB function ``fminunc`` (slower). 
* ``kantorovich.m`` computes the Kantorovich function and its gradient and Hessian.

## Dependencies ##

To use these MATLAB functions you must first install
* [MATLAB-Voro](https://github.com/smr29git/MATLAB-Voro)

## Examples ##

See the MATLAB live script ``Examples_MATLAB_SDOT.mlx`` and the corresponding PDF file ``Examples_MATLAB_SDOT.pdf``. We have tested the code with target measusres with up to 100,000 Dirac masses.

## Licence ##

See LICENCE.md

## Software limitations ##

* quadratic cost
* Lebesgue measure

We plan to address some of these limitations in future updates.

## Related software ##

* Quentin's code
* Geogram: <https://github.com/BrunoLevy/geogram>
* Peyre's code

## Main contributors ##

* [Steve Roper](https://www.gla.ac.uk/schools/mathematicsstatistics/staff/stevenroper/#), University of Glasgow
* [David Bourne](http://www.macs.hw.ac.uk/~db92/), Heriot-Watt University and the Maxwell Institute for Mathematical Sciences
* Mason Pearce, Heriot-Watt University and the Maxwell Institute for Mathematical Sciences
