# ASE tutorials - Molecular Dynamics 
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/UCL-DDMD/ASE-Tutorials-/main?urlpath=lab)


In these tutorials, we'll run some molecular dynamics (MD) simulation on bulk metals. 

## Exercise 1: 

In this example, you'll see a demonstration of running NVE MD simulations on Cu. You will then need to run the example 
yourself for FCC Al, Cu, Au, Ni, Ag, Pt and Pd. You can find the lattice constant (a) for each one in a table below:



| Metal     | Unit Cell | Lattice constant (Å) |
|-----------|-----------|----------------------|
| Aluminum  | FCC       | 4.046                |
| Copper    | FCC       | 3.597                |
| Gold      | FCC       | 4.065                |
| Nickel    | FCC       | 3.499                |
| Silver    | FCC       | 4.079                |
| Platinum  | FCC       | 3.912                |
| Palladium | FCC       | 3.859                |


Some questions to think about:

1) What would happen with the simulation if we changed the convergence criteria (fmax)? 
2) The lattice constants after the geometry optimisations are different to the experiment. Why? 
3) Why do we perform MD simulations on supercells? What would happen if we did MD on the unit cell?
4) If you wre to change the Maxwell-Boltzmann distribution how would this change the simulation? 
5) How can you improve the accuracy of your calculation, how would this effect the computational overhead? 




## Exercise 2: 

Post-processing is an important step after performing any calculation (especially MD). 
Here you'll extract the scientifically relevant information from your simulation and analyse your data. 

You'll see an example of plotting the simulation temperature against the timestep. Have a go at plotting
the temperature for all your systems and then the total, kinetic, and potential energy on the same graph. 

Here are some questions to consider:

1) For the plot of temperature, what do you see in the early stages of the simulation? What does this mean?
2) Looking at the plots for energy, what can you infer from this? 
3) How would this change in an NVT simulation? 





## Exercise 3: 

Here you'll see a plot of the 2D potential energy surface of Cu-Cu bond length from the MD trajectory. 
You'll also see a scatter plot of the Cu-Cu bonds and a heatmap. 
Try plotting these for the other metal systems you ran simulations on in Exercise 1. 

Here are some questions to think about: 

1) What do the different contours in the 2D PES correspond to? What does this tell you about the Cu-Cu bonds?
2) Similarly, what information can you deduce from the scatter plot? And what does this say about the nature of Cu-Cu bonds?

From your answers to 1) and 2) try simulating your systems as BCC instead of FCC
Hint: You'll have to work out a new lattice constant, a

3) First, run the simulation without a geometry optimisation. What packing efficiency do you get? Is this as expected? 
4) Secondly, now run the MD with a geometry optimisation, does the packing efficiency change? What does this tell us? 














