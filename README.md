## Reduction types of genus-3 curves in a special stratum of their moduli space
### Complementary material

In this repository we give extra material to complement the paper
> _Reduction types of genus-3 curves in a special stratum of their moduli space_, I. Bouw, N. Coppola, P. Kiliçer, S. Kunzweiler,
E. Lorenzo García and A. Somoza.

#### Proof of Proposition 3.4 in characteristic 0.
The file `InvariantsGenerateDO.txt` is an explicit proof that the invariants

<img src="https://render.githubusercontent.com/render/math?math=I_3 = ABC">

<img src="https://render.githubusercontent.com/render/math?math=I_3'= A(a^2-4BC) + B(b^2-4AC) + C(c^2-4AB)">

<img src="https://render.githubusercontent.com/render/math?math=I_3''= -4ABC + Aa^2 + Bb^2 + Cc^2 - abc">

<img src="https://render.githubusercontent.com/render/math?math=I_6 = (a^2-4BC)(b^2-4AC)(c^2-4AB)">

associated to a generic plane quartic in the strata <img src="https://render.githubusercontent.com/render/math?math=M_{3,V_4}"> given by the equation

<img src="https://render.githubusercontent.com/render/math?math=Ax^4 + By^4 + Cz^4 + ay^2z^2 + bx^2z^2 + cx^2y^2 = 0">

generate all the Dixmier-Ohno invariants, since it shows how one can write the latter in terms of the former.


#### Implementation of the reduction type
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/NirvanaC93/Invariants-Special-Strata-Genus-3-Curves/master?filepath=.%2Fexample.ipynb)

The given SageMath file implements a class with methods related to the results in the paper. See the [binder demo](https://mybinder.org/v2/gh/NirvanaC93/Invariants-Special-Strata-Genus-3-Curves/master?filepath=.%2Fexample.ipynb) for details.
