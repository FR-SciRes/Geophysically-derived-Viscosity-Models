# Geophysically-derived-Viscosity-Models
The project aims to estimate upper mantle viscosity using seismic and magnetotelluric (MT) observations. Mainly, the user should have either seismic data only or both seismic and MT data to use the codes here. The idea is that temperatures are inferred from seismic data and water contents are inferred from MT data using MATE (see https://doi.org/10.1029/2021JB023824 for MATE reference). Both temperature and water content are then used in calculating viscosity. Thus, constraining temperatures and water contents from observations would also constrain viscosity estimates.
 
A. Below are the main codes used in this project.
(A.1) Vs_to_T_Conversion_HeFESTo.m --> converts seismic velocities to temperatures (and composition) using the theoretical calculations of HeFESTo (see https://doi.org/10.1029/2021JB023824 for reference on HeFESTo)
(A.2) TotalStrainRate_HK03_fH2O.m --> a function that calculates the total strain rate (= diffusion strain rate + dislocation strain rate + DisGBS strain rate) for certain stress, temperature, water content, grain size and pressure. This function has to be called when calculating viscosity, where viscosity equals stress over strain rate (a user has to write a code for this).
NOTE: These codes have important comments for the user, particularly the inputs. 

B. The following files have to be downloaded together with the codes.
(B.1) InterpolatedVelocities_allP&T_basalt.56Qcor; InterpolatedVelocities_allP&T_pyrolite.56Qcor; InterpolatedVelocities_allP&T_harzburgite.56Qcor; InterpolatedVelocities_allP&T_olivine.56Qcor

C. Step-by-step process in using the codes.
(C.1) Prepare seismic and MT profiles separately, which contain seismic velocity or electrical conductivity value (from MT) for each depth value. For MT data, refer to MATE as to its instruction and file structure. For seismic data, the user needs to have a file with depths and corresponding velocities (in columns). If the user has velocity ranges, then the file should have both the lower and upper velocity bounds, and potentially average velocities. 
(C.2) Run Vs_to_T_Conversion_HeFESTo.m in MATLAB, where the velocities in step 1 (named as Vs_lo and Vs_up for lower and upper velocity bounds, and Vs_ave for average velocities in my code) will be extracted, and then be converted to temperatures and phases (composition). Check the code for instructions and comments.
(C.3) The temperatures and phases together with electrical conductivities will be the inputs for MATE tool (see reference) to get water content values.
(C.4) Create a short code to calculate viscosity (=stress/total strain rate) using the function TotalStrainRate_HK03_fH2O(). See TotalStrainRate_HK03_fH2O.m for function inputs and outputs. The inputs are temperature, water content, pressure, stress and partial melt.

D. Important citations to include when using the codes.

Özaydin, S., & Selway, K. (2020). MATE: An analysis tool for the interpretation of magnetotelluric models of the mantle. Geochemistry, Geophysics, Geosystems, 21(9), 1–26. https://doi.org/10.1029/2020GC009126

Ramirez, F. D. C., Selway, K., Conrad, C. P., & Lithgow-Bertelloni, C. (2022). Constraining upper mantle viscosity using temperature and water content inferred from seismic and magnetotelluric data. Journal of Geophysical Research: Solid Earth, 127, e2021JB023824. https://doi.org/10.1029/2021JB023824

Stixrude, L., & Lithgow-Bertelloni, C. (2005b). Thermodynamics of mantle minerals - I. Physical properties. Geophysical Journal International, 162(2), 610–632. https://doi.org/10.1111/j.1365-246X.2005.02642.x

Stixrude, L., & Lithgow-Bertelloni, C. (2011). Thermodynamics of mantle minerals - II. Phase equilibria. Geophysical Journal International, 184(3), 1180–1213. https://doi.org/10.1111/j.1365-246X.2010.04890.x
