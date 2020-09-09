Author: Simon Isphording
Date: 8-7-2019
------------------------------------------------------------
This folder contains two versions of the model: one models cells with the apical side only and the other models full cells.
In both cases, main.py contains everything needed for a simulation, with dependencies on:
* vertex_objects.py
* transitions.py
* parameters.py
* voronoi.py (in the apical only model)
For the initial configuration, the apical only model currently uses voronoi.py, while the apical-basal model uses the scripts in the initial configuration folder.
For setting the parameters, you either change the default parameters in parameters.py, or use a .pickle file including all parameters. For an example you can check hollandia_simulation.py
Hollandia_simulation.py and hollandia_analysis.py are scripts used to run simulations on the supercomputer.
------------------------------------------------------------

For any questions regarding the code, please contact: simonisphording@gmail.com