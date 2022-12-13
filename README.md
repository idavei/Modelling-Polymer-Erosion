
# Modelling Polymer Erosion

The files '2Dmodel.py' and '3Dmodel.py' are Python programs which were written to first, model the innate spatial structure of a specific polymer, and then to simulate erosion of that polymer in the presence of surrounding water. The 3D model is a natural extension to the 2D model.

The specific polymer modelled is a combination of 20 parts 1,3-bis(carboxy-phenoxy propane) (CPP) and 80 parts sebacic acid (SA), which together form p(CPP-SA) 20:80. This polymer's composition is notable for possessing both easy-to-erode (amorphous) regions and hard-to-erode (crystalline) regions. My aim was to represent the crystalline regions within the polymer in a novel manner, that is, closely mirroring physical observations of the polymer in empirical studies.

To run the programs locally, you may clone the repository and execute the scripts using the command line or a Python IDE of your choice. There are no arguments to be passed to the command line and model parameters must be altered within the script. The output from either program should be a single figure portraying an instantaneous snapshot of the eroded polymer grid after a certain amount of time has elapsed in the simulation.

The colour scheme in generated figures should be intepreted as follows:  
Green: Crystalline  
Red: Amorphous  
Blue: Water  
Maroon: Boundary

The model parameters which can be altered for meaningful outcomes are located at the following lines in each script:

2D model  
Line 715: Dimension for one edge of the square grid  
Line 783: Time after which the simulation terminates

Note: The 2D model also contains a concept model for _acid-driven self-inhibition_, which refers to the slowing down of erosion due to build-up of acidic by-products in eroded areas of the polymer.

3D model    
Line 376: Dimension for one edge of the cubic grid  
Line 411: Time after which the simulation terminates