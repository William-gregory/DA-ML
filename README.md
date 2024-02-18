# Learning sea ice model errors from data assimilation increments

This repository is part of the larger [M2LInES](https://m2lines.github.io) project. M2LInES involves developing climate model parameterizations using machine learning, in order to improve phyiscs and reduce systematic model biases. Here we use Convolutional Neural Networks to derive a mapping from model state variables to sea ice concentration analysis increments from an ice-ocean data assimilation experiment. The model which this study has been applied to is an ice-ocean configuration of the Geophysical Fluid Dynamics Laboratory (GFDL) Seamless system for Prediction and EArth system Research (SPEAR) model.

### Offline learning

The `offline_learning` folder contains code and data relating to the article [Gregory et al., 2023](https://agupubs.onlinelibrary.wiley.com/doi/pdfdirect/10.1029/2023MS003757), with an example jupyter notebook of how the CNN was trained and model selection performed by carrying out 5-fold cross-validation tests.

### Online ice-ocean

The `online_iceocean` folder contains example scripts of how to implement the trained CNN into SPEAR ice-ocean simulations, by updating the sea ice restart files, as a way to correct short-forecasts. This methodology is outlined in the article [Gregory et al., 2024](https://doi.org/10.1029/2023GL106776).

### Online forPy

Example scripts of implementing the CNN into SIS2, via the Forpy Fortran-Python interface. This approach allows the CNN to be called at the model timestep, rather than the approach above, which relies on updating the model restart files. Using Forpy to implement a CNN into MOM6 is outlined in [Zhang et al., 2023](https://agupubs.onlinelibrary.wiley.com/doi/pdfdirect/10.1029/2023MS003697).
