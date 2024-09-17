Code for Sakamoto et al. (2024).

## Folder "simulation_continent_island"
This folder contains the codes for continent-island simulation.

### How to run
The c++ code can be compiled by the following command:
```
g++ *.cpp -Wall -Wextra -std=c++17 -O3 -o XXX.out
```

For execution, you need to pass eight parameters:
```
./XXX.out X1 X2 X3 X4 X5 X6 X7 X8
```
where  
``X1`` (int): size of island population  
``X2`` (int): size of continent population  
``X3`` (double): mutation rate in the background region  
``X4`` (double): migration rate from continent to island  
``X5`` (double): phenotypic optimum in the island population (optimum in the continent is fixed as 0)  
``X6`` (double): effect size at the focal locus  
``X7`` (double): standard deviation of effect size of mutations in the background region  
``X8`` (double): selection strength

Instead, you can compile and run the simulation using "run_simu.py":
```
python3 run_simu.py X
```
where ``X`` is int.
This python code refers to the Xth row of "parameter_list_full.csv" and pass its parameter set to the c++ simulation.

### Output
The simulation outputs two files.

"mean_var.txt": 11 columns  
time  
mean fitness in the island  
variance of fitness in the island  
mean phenotype in the island  
variance of phenotype in the island  
derived allele frequency at the focal locus in the island  
mean fitness in the continent  
variance of fitness in the continent  
mean phenotype in the continent  
variance of phenotype in the continent  
derived allele frequency at the focal locus in the continent

"bin_contribution.txt": 1002 columns  
time (1st column)  
difference of mean genotypic value at the $i$-th background locus ($i+1$-th column, $1\leq i \leq 1000$)  
difference of mean genotypic value at the focal locus (1002th column)

## Folder "simulation_symmetric_island"
The usage is almost identical as the folder "simulation_continent_island".
Only difference is that symmetric migration is assumed in this version. 
For convinience, we still use 'continent' and 'island' to call each population although there is no essential difference between the populations. For more details of the usage, see the section for Folder "simulation_continent_island.

## Folder "simulation_continent_island_linkage"
The usage is almost identical as the folder "simulation_continent_island".
Only difference is that the focal locus is located in the middle of the background region. 
For more details of the usage, see the section for Folder "simulation_continent_island.

## Folder "boundary"
This code calculates the threshold mutation rate using theory.

### How to run
First, specify the parameters in the header of "output_boundary.cpp".
Then, compile and run the codes by following commands:
```
g++ output_boundary.cpp -Wall -Wextra -std=c++17 -O3 -o XXX.out
./XXX.out
```

### Output
The program outputs one file "boundary.txt":
migration rate (first column)
threshold mutation rate (second column)

## Folder "determine_z2_p2"
This code calculates the equilibrium $\bar{z}_2$ and $p_2$ based on the theory.

### How to run
First, specify the parameters in the header of "calculate_equilibrium.cpp".
Then, compile and run the codes by following commands:
```
g++ calculate_equilibrium.cpp -Wall -Wextra -std=c++17 -O3 -o XXX.out
./XXX.out
```

### Output
The program outputs one final file "equilibrium_value.txt": 8 columns  
mutation rate  
$\bar{z}_2$ given no divergence at the focal locus  
Variance of $z_1$  
Variance of $z_2$ given no divergence at the focal locus  
$\bar{z}_2$ given stable divergence at the focal locus (when no stable equilibrium, it returns nan)  
$p_2$ given stable divergence at the focal locus (when no stable equilibrium, it returns nan)  
$\bar{y}_2$ given stable divergence at the focal locus (when no stable equilibrium, it returns nan)  
Variance of $z_2$ given stable divergence at the focal locus (when no stable equilibrium, it returns nan)  

The file "regi_beta.txt" is a intermediate file.
