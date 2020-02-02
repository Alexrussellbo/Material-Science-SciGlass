# Modify the defination of distance bewtween the data points in the subset and the compound.

In this version, we still utilize the Euclidean distance to measure the similarity

e.g. there are 3 compounds located in the four compounds system (oxide_A, B, C, D)

		   |  oxide_A  |  oxide_B  |  oxide_C  |  oxide_D  |
Compound 1 |    a1%    |    b1%    |    c1%    |    d1%    |
Compound 2 |    a2%    |    b2%    |  c2% = 0% |    d2%    |
Compound 3 |    a3%    |  b3% = 0% |    c3%    |  d3% = 0% |

**NOTE:** $a_i$, $b_i$, $c_i$ and $d_i$ except c2, b3 and d3 are all none-zero.

However, there is one data point in one of the sub-dataset (at least one zero percent) which is as follows

		   |  oxide_A  |  oxide_B  |  oxide_C  |  oxide_D  |
    Data   |    aa%    |    bb%    |    cc%    |     0     |


The formula of these three distances:

dist_compound1-data = $\sqrt{(a1-aa)^2+(b1-bb)^2+(c1-cc)^2+(d1-0)^2}/100$
dist_compound2-data = $\sqrt{(a2-aa)^2+(b2-bb)^2+(0-cc)^2+(d2-0)^2}/100$
dist_compound3-data = $\sqrt{(a3-aa)^2+(0-bb)^2+(c3-cc)^2+(0-0)^2}/100$

All the four terms are included in these formula when calculating the Euclidean distance, which means each data point can obtain 3 distance no matter whether it has zero percent oxide component or not.










