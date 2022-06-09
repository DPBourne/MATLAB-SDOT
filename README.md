# MATLAB-SDOT

MATLAB functions for solving semi-discrete optimal transport problems in 3D.
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

* The code is limited to 3D (2D code coming soon).
* The source measure is the Lebesgue measure on a cuboid (rather than a general absolutely continuous measure).
* The transport cost is either the quadratic cost or the periodic quadratic cost (optimal transport on a triply-periodic box).

We plan to address some of these limitations in future updates.

## Related software ##

* [MongeAmpere](https://github.com/mrgt/MongeAmpere) and [PyMongeAmpere](https://github.com/mrgt/PyMongeAmpere) by Quentin Mérigot
* [sdot by Jocelyn Meyron](https://github.com/nyorem/sdot)
* [Geogram by Bruno Lévy](https://github.com/BrunoLevy/geogram), [documentation](https://brunolevy.github.io/geogram/dir_dfcc9fc6d69b9d57f9f159e89cabbae9.html)
* [Numerical Tours by Gabriel Peyré](https://nbviewer.org/github/gpeyre/numerical-tours/blob/master/matlab/optimaltransp_7_semidiscrete.ipynb)

## Main contributors ##

* [Steve Roper](https://www.gla.ac.uk/schools/mathematicsstatistics/staff/stevenroper/#), University of Glasgow
* [David Bourne](http://www.macs.hw.ac.uk/~db92/), Heriot-Watt University and the Maxwell Institute for Mathematical Sciences
* Mason Pearce, Heriot-Watt University and the Maxwell Institute for Mathematical Sciences
