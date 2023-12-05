# ASE tutorials - Quantum Espresso

In this tutorial, you'll learn how to run DFT calculations with ASE. The exercise will take you through running DFT on Si, 
where you'll build a model using ASE, perform a geometry optimisation, and work out the electronic density of states and band structure for Si. 

## Getting started 

In your home directory on Aristotle, clone the tutorial GitHub repository through: 

`git clone https://github.com/UCL-DDMD/ASE-Tutorials.git`

To transfer files between your local machine and Aristotle run:

`scp -r ucxxxxx@aristotle.rc.ucl.ac.uk:path/to/your/file Desktop`

where ucxxxxx is your UCL username. 

## Exercise 4:

In this exercise, you'll learn about the importance of converging parameters when running simulations. 
Here you are given a Python script (kpoints.py) to sample K-points through single point energy calculations on Si. 
After running the calculation, use the script plot_kpoints.py to plot your results. Which K-point sampling mesh is converged? 

## Exercise 5:

Now that you've successfully chosen a k-grid setting, run a geometry optimisation for bulk Si. 


## Exercise 6: 

Once you've completed your optimisation, use this script to obtain an electronic band structure of Si. 

## Exercise 7: 

Here, we'll look at the electronic density of states for our system. 

