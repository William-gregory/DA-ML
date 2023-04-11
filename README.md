# Deep learning of systematic sea ice model errors from data assimilation increments

This repository contains code and data relating to work as part of the [M2LInES](https://m2lines.github.io) project, which involves developing climate model parameterizations using machine learning, in order to reduce systematic model biases. Here we use Convolutional Neural Networks to derive a mapping from model state variables to sea ice concentration analysis increments from an ice-ocean data assimilation experiment (see [Gregory et al., 2023](https://doi.org/10.48550/arXiv.2304.03832)). The model which this study has been applied to is the Geophysical Fluid Dynamics Laboratory (GFDL) Seamless system for Prediction and EArth system Research (SPEAR) model.

### Offline learning

In the offline_learning folder is an example jupyter notebook of how the CNN was trained and model selection performed by carrying out 5-fold cross-validation tests

### Online ice-ocean

Example scripts of how to implement the trained CNN into SPEAR ice-ocean simulations as a way to correct short-forecasts.
