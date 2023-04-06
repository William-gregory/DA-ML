### Network Training

For details on the data pre-processing, network architectures, and cross-validation testing, please see the code in the `DA-ML.ipynb` notebook. Dependencies for the code include:

`PyTorch`\
`Xarray`\
`Cartopy`

The optimized weights for each network are saved in the `CNN_weights` folder. Weights labelled 'CV' are those generated from cross-validation experiments to assess the network's generalization performance. Those labelled 'allsamples' have been trained using all data points, and should be used for subsequent analysis.

### Data

The data used to train the networks, as outlined in the paper, are hosted on a permanent ftp server: `ftp://sftp.gfdl.noaa.gov/perm/William.Gregory/`. The inputs to the network are saved as `seaice_DA-ML_inputs_1982-2017.nc`, while the outputs (increments) are `seaice_DA-ML_outputs_1982-2017.nc`. These can be directly downloaded with e.g., `wget` (see jupyter notebook).
