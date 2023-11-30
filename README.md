# optcu
Cluster expansion of O-Pt/Cu(111) with **CELL**

This repository contains python files to run the example cluster expansion in Section 3 of the paper

[*CELL: a Python package for cluster expansion with a focus on complex alloys* Santiago Rigamonti, Maria Troppenz, Martin Kuban, Axel Huebner, and Claudia Draxl](https://arxiv.org/abs/2310.18223)

[**CELL** documentation](https://sol.physik.hu-berlin.de/cell/)

## Instructions

You can reproduce the example of Section 3 of the paper by running the script in file `optcu.py`. 

To do so, execute

```
python optcu.py task
```

here, task is an integer number which indicates what task to run. 

The available tasks are:

- 1: Create parent lattice (Listing 1)
- 2: Create super cell (Listing 2)
- 3: Create structures set (Listing 3)
- 4: Compute total energy, non relaxed (Listing 4)
- 5: Compute reference energies
- 6: Compute adsorption energy, relaxed (Listing 5)
- 7: Consider other example properties
- 8: Visualization: Property versus concentration (Listing 6 and Figure 8)
- 9: Create clusters pool (Listing 7)
- 10: Create input matrix X and vector of properties P (Listing 8)
- 11: Create CE model with ridge regression (Listing 9)
- 12: Create CE model with subset selection (Listing 10)
- 13: Visualization: Property versus concentration (Figure 13)
- 14: Visualization: Optimization by cluster selector (Figure 14)
- 15: Create CE model with LASSO (Listing 11)
- 16: Visualization: Optimization by cluster selector (Figure 15)

The execution of these tasks correspond to the listings and figures on the paper as indicated in parenthesis.

Figures are generated as PNG files in the same folder where you run the script.

You can open the file `optcu.py` with a Python code editor to see the code. 
The Listings corresponding to the paper are indicated by comment strings in the code. 

