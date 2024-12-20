# MATLAB-SDOT

MATLAB functions for solving semi-discrete optimal transport problems in 2D & 3D.
* ``SDOT_damped_Newton.m`` computes the optimal transport cost between the Lebesgue measure and a discrete measure on a 3D box with respect to the quadratic cost or the periodic quadratic cost, using the damped Newton method from
    * Jun Kitagawa, Quentin Mérigot, Boris Thibert, Convergence of a Newton algorithm for semi-discrete optimal transport. *J. Eur. Math. Soc.* 21 (2019), no. 9, pp. 2603–2651.
* ``SDOT_fminunc.m`` solves the same problem using the MATLAB function ``fminunc`` (slower). 
* ``kantorovich.m`` computes the Kantorovich function and its gradient and Hessian.

## Getting started ##

To use these MATLAB functions you must first install
* [MATLAB-Voro](https://github.com/smr29git/MATLAB-Voro)

## Examples ##

See the MATLAB live script ``Examples_MATLAB_SDOT.mlx`` and the corresponding PDF file ``Examples_MATLAB_SDOT.pdf``. We have tested the code with target measures with up to 100,000 Dirac masses.

## Licence ##

See LICENCE.md

## Software limitations ##

* The source measure is the Lebesgue measure on a cuboid.
* The support of the discrete target measure must be contained in the support of the source measure.
* The transport cost is either the quadratic cost or the periodic quadratic cost.

We plan to address some of these limitations in future updates.

## Related software ##

* [pysdot](https://github.com/sd-ot/pysdot) by Quentin Mérigot & Hugo Leclerc
* [Geogram](https://github.com/BrunoLevy/geogram) by Bruno Lévy, [documentation](https://brunolevy.github.io/geogram/dir_dfcc9fc6d69b9d57f9f159e89cabbae9.html) 
* [sdot](https://github.com/nyorem/sdot) by Jocelyn Meyron
* [MongeAmpere](https://github.com/mrgt/MongeAmpere) and [PyMongeAmpere](https://github.com/mrgt/PyMongeAmpere) by Quentin Mérigot
* [Numerical Tours](https://nbviewer.org/github/gpeyre/numerical-tours/blob/master/matlab/optimaltransp_7_semidiscrete.ipynb)  by Gabriel Peyré

## Main contributors ##

* [Steve Roper](https://www.gla.ac.uk/schools/mathematicsstatistics/staff/stevenroper/#), University of Glasgow
* [David Bourne](http://www.macs.hw.ac.uk/~db92/), Heriot-Watt University and the Maxwell Institute for Mathematical Sciences
* Mason Pearce, Heriot-Watt University and the Maxwell Institute for Mathematical Sciences
