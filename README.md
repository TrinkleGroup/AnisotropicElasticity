# Anisotropic elasticity solution for dislocations

Source codes and sample input files for generating initial dislocation geometries from the anisotropic elasticity solution. 

The anisotropic elasticity solution is evaluated by numerically integrating the integral expression derived using the sextic formalism. For the derivations, see for example, Chapter 12 of R. W. Balluffi, "Introduction to Elasticity Theory for Crystal Defects", Cambridge University Press, 2012. 

Most of the source codes here were originally written by D. R. Trinkle.

The source codes appended "-ref" have been modified slightly by A. M. Z. Tan. Use these if you need to evaluate the displacement field self-consistently, which is the correct way to set up any dislocation that has edge character.

The source code appended "_Lauren" has been modified slightly by L. Smith. It correctly treats atoms which cross the cut plane when displaced. (It might be hardcoded for the Ni edge dislocation that she looked at though...)

See the sample input files in the "examples" sub-directory for usage examples.
