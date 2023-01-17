# Deep learning of systematic sea ice model errors from data assimilation increments

This repository contains code and data relating to work as part of the [M2LInES](https://m2lines.github.io) project, which involves developing climate model parameterizations using machine learning, in order to reduce systematic model biases. Here we use Convolutional Neural Networks to derive a mapping from model state variables to sea ice concentration analysis increments from an ice-ocean data assimilation. The model which this study has been applied to is the Geophysical Fluid Dynamics Laboratory (GFDL) Seamless system for Prediction and EArth system Research (SPEAR) model. A manuscript detailing this work is in preparation, and will be submitted to the Journal of Advancements in Model Earth Systems (JAMES) in due course.

### Network Training

For details on the data pre-processing, network architectures, and cross-validation testing, please see the code in the `DA-ML.ipynb` notebook. Dependencies for the code include:

`PyTorch`\
`Xarray`\
`Cartopy`

The weights for each network are saved in the CNN_weights folder. The weights labelled 'CV' are those generated from cross-validation experiments to assess the network's generalization performance. Those labelled 'allsamples' have been trained using all data points, and should be used for subsequent analysis.
