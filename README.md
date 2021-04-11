# Computer Model Emulation with High-Dimensional Functional Output in Large-Scale Observing System Uncertainy Experiments
 

# General description
 This file shows the entire pipelines to implement the NNGP-based emulator
 in combination with the active subspace approach for dimension reduction. 
 The details of the methods can be found in the paper titled "Computer Model
 Emulation with High-Dimensional Functional Output in Large-Scale Observing
 System Uncertainty Experiments" by  Ma et al. Technometrics (2021) (https://doi.org/10.1080/00401706.2021.1895890).

# Instructions for the overall picture 
The entire codes have three parts.
1. In the first part, R codes are used to implement functional principal
 component analysis (FPCA) so that the functional output from OCO-2 forward 
 model can be represented as a linear combination of basis functions and 
 corresponding FPC scores. 
2. In the second part, R codes are used to find the active subspace from 
 the space formed by the state vectors. The key quantity here is the 
 projection matrix that projects a high-dimensional space into a low-dimensional
 space.
3. In the third part, MATLAB codes are used to train the NNGP emulator based
 on a set of training data, in which the active variables and viewing geometry
 are treated as input variables, and FPC scores are treated as output. 


# Note on the data
The data from the full-physics forward model can be downloaded at the NASA Goddard Earth Science Data and Information Services Center (GES DISC https://disc.gsfc.nasa.gov/datasets?page=1\&keywords=OCO-2). 
 The dataset for the reduced order model

# Note on the code
1. We only provide the code to perform FPCA, active subspace, and NNGP emulation.  
 Detailed instrations for these three components can be found in the file main.R.
2. As performing active subspace requires downloading the full dataset that is too 
 large to upload, we only provide relevant code to implement active subspace, and 
 the resulting active variables. 
3. The NNGP emulator is provided as a general form that implements the independent emulator
 for all the bands, and separable emulator for the WCO2 and SCO2. 
4. For O2 band, we peform independent emulation with three different files: 
 * O21_AS_nugget_IndEmulation.m. This file contains the emulation code for the 
 first FPC score for the O2 band.
 * O22_AS_nugget_IndEmulation.m. This file contains the emulation code for the 
 second FPC score for the O2 band.   
 * O23_AS_nugget_IndEmulation.m. This file contains the emulation code for the 
 third FPC score for the O2 band.
5. For the WCO2 and SCO2 bands, we can use independent emulator and separable emulator.
 For independent emulator, the corresponding files are WCO2_AS_nugget_IndEmulation.m
 and SCO2_AS_nugget_IndEmulation.m; For Separable emulator, the corresponding files 
 are CO2_AS_nugget_SepEmulation.m

# What you can expect from the code?
 The results for NNGP emulators in Table 2 and Table 3 can be reproduced with the code along with Figures 3 to 6. 

