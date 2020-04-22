# PQFDC

ExpectileCurves.R file includes the codes to obtain expectile curves using B-spline and P-spline fits. 
MquantileCurves.R file includes the codes to obtain M-quantile curves using B-spline and P-spline fits. 

Simulation_Example1_Exfclust.R file includes the codes to generate simulation dataset as in Example 1 and apply Exfclust methods. It requires that R functions defined in ExpectileCurves.R are attached. 

Simulation_Example2_Exfclust.R file includes the codes to generate simulation dataset as in Example 2 and apply Exfclust methods. It requires that R functions defined in ExpectileCurves.R are attached. 

Simulation_Example2_Mqfclust.R file includes the codes to apply Mqfclust methods. It requires that R functions defined in MquantileCurves.R and dataset generated in Simulation_Example_Exfclust.R are attached. 

The results of the simulation code are displayed in the paper "Pseudo-quantile based Functional Data Clustering" by Joonpyo Kim and Hee-Seok Oh. The paper was accepted for publication by the Journal of Multivariate Analysis: https://doi.org/10.1016/j.jmva.2020.104626
