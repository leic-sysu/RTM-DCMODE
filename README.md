# RTM-DCMODE
This is the code of the Decomposition-based Constrained Multi-Objective Differential Evolution (DCMODE) algorithm for shallow water bathymetry retrieval.

Instructions:

1. Use mainDCMODE.m to execute the DCMODE algorithm. For detailed explanations of the input parameters, please refer to the comments in the code.
2. Two Landsat-8 OLI images covering the Zhaoshu Island and the North Island were provided for testing the algorithm in the file "testdata". For more detailed information of the Landsat-8 images, please refer to the manuscript entitled "Advancing Multispectral Image-Derived Physics-Based Bathymetry: Multi-Objective Evolutionary Computation for Shallow Water Depth Retrieval."
3. We also provided the in situ sonar depth for each area. Please refer to the file "validationdata".
4. If you want to use this algorithm in other study areas, remember to replace the benthic endmember spectra provided in testdata/Parameters.xlsx with the measured substrate spectra in the specific study region.
5. This code is used for retrieving bathymetry from Landsat-8 OLI images. If you want to apply this algorithm to other multispectral image, such as Sentinel-2A/B, you should modify the spectral response function.
