### Network Training

For details on the data pre-processing, network architectures, and cross-validation testing, please see the code in the `DA-ML.ipynb` notebook. Dependencies for the code include:

`PyTorch`\
`Xarray`\
`Cartopy`

The notebook shows how the 5-fold cross validation analysis was performed in order to do model selection. The python script shows how the model was trained using all data to derive a "final" network.

### Data

The data used to train the networks, as outlined in the paper, are hosted on a permanent ftp server: `ftp://sftp.gfdl.noaa.gov/perm/William.Gregory/`. The inputs to the network are saved as `seaice_DA-ML_inputs_1982-2017.nc`, while the outputs (increments) are `seaice_DA-ML_outputs_1982-2017.nc`. These can be directly downloaded with e.g., `wget` (see jupyter notebook).
