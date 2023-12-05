# ASE tutorials - Quantum Espresso

In this tutorial, you'll learn how to run DFT calculations with ASE. The exercise will take you through running DFT on Si, 
where you'll build a model using ASE, perform a geometry optimisation, and work out the electronic density of states and band structure for Si. 

## Getting started 

Firstly, login to Aristotle via:

`ssh  ucxxxxx@aristotle.rc.ucl.ac.uk`

where ucxxxxx is your UCL username. 

In your home directory on Aristotle, clone the tutorial GitHub repository through: 

`git clone https://github.com/UCL-DDMD/ASE-Tutorials.git`

To transfer files between your local machine and Aristotle run:

`scp -r ucxxxxx@aristotle.rc.ucl.ac.uk:path/to/your/file Desktop`


## Exercise 4:

In this exercise, you'll learn about the importance of converging parameters when running simulations. 
Here you are given a Python script (kpoints.py) to sample K-points through single point energy calculations on Si. 
After running the calculation, use the script plot_kpoints.py to plot your results. Which K-point sampling mesh is converged? 

run the exercise with:

`bash submission.script`

## Exercise 5:

Now that you've successfully chosen a k-grid setting, replace `KPOINT` with your converged value and  run a geometry optimisation 
for bulk Si. 
In this script, we create a unit cell of Si with `ase.build.bulk`, we then perform a geometry optimisation on the nuclear coordinates and the unit cell. 


## Exercise 6: 

Once you've completed your optimisation, use this script to obtain an electronic band structure of Si. 
This time instead of re-running the geometry optimisation, we read in the last (converged) image from the output trajectory like so:
`'Opt.traj@-1'`. 
We will then calculate the band structure for Si with ASE. 

## Exercise 7: 

Here, we'll look at the electronic density of states for our system with ASE. 


## Some things to think about: 

Here are some questions to think about whilst performing your calculation:

1) What would happen if we used an unconverged K-grid mesh for our electronic calculations? 
2) How would only optimising the coordinates of the atoms affect our results? 
3) What information does the band structure provide us? 
4) What information does the density of states provide us? 
5) What value is the fermi level from our plots? Does this make sense?
6) Does our DOS looks like our reference plots? What could we do to improve our method?

![alt text](https://github.com/UCL-DDMD/ASE-Tutorials/blob/main/Si/Reference_DOS.png)
*Reference DOS (PBE with QE): https://pranabdas.github.io/espresso/hands-on/dos* 
