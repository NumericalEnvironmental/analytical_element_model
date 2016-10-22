# analytical_element_model
This python script employs a variation on the analytic element method to simulate steady-state groundwater flow in a two-dimensional aquifer characterized by point sources/sinks, line source/sinks, faults of varying permeability, and regional flow. The steady-state model assumes a leaky aquifer is juxtaposed against an overlying aquitard of some specified thickness and vertical hydraulic conductivity. Line sources and faults are modeled using numerical integration (with a doublet element to address flow around fault segments). A particle tracking routine, constructed based upon a differential equation solver, is added as a post-processor.

The following tab-delimited input files are required (assisted with a tkinter-based GUI):

* base.txt - spatial scale and gridding (for output)
* elements.txt - definition of wells, line sources, and/or faults
* hydro.txt - hydraulic properties of aquifer and aquitard
* particlex.txt - initial positions of particles used in particle tracking routine
* track.txt - basic switches for particle tracking (on/off, forward/reverse, etc.)

More background information can be found here: https://numericalenvironmental.wordpress.com/2016/07/11/analytic-element-modeling/

Email me with questions at walt.mcnab@gmail.com. 
