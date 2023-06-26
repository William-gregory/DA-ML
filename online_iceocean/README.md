### Implementation of CNN into SPEAR ice-ocean simulations

The CNN correction is applied back into SPEAR ice-ocean simulations as a bias correction tool. In other words, we run a short simulation/forecast of the model, e.g., for 5 days, at which point the CNN updates the restart files as a post-processing step, and then the model continues from this updated state for the next 5 days. Note that the scripts within this repository are for implementing a few different variations of CNN and/or DA correction schemes:

`DAML_G23.csh` - train a CNN (based on Gregory et al., 2023) using all available data between 1982-2017.

`DAML_IFA.csh` - compute the mean SIC increment for each day of the year, over 1982-2017, and apply this as a "mean error correction" at each grid point location.

`DAML_OPTp1.csh` - train a CNN (as in G23) using all available data between 1982-2012 (this is then followed by a subsequent update with DA - see below).

`DAML_OPTp2.csh` - based on a final optimized network which has been tuned based on the increments produced from the DAML_OPTp1.csh simulation.

The post-processing step to implement a CNN or IFA into an xml is done in the following way:

    <freInclude name="OM4_postprocess">
        <postProcess>
            <csh><![CDATA[ 
         
             # implementing DA-based CNN correction scheme:
             cd $work/RESTART
             if (! -e coupler.res) then
                 echo model still running
             else
                 cp DAML_G23.csh .
                 csh DAML_G23.csh >& log_CNN.out
                 rm log_CNN.out
                 rm DAML_G23.csh
             endif
             
             # rest of post-processing script follows as normal
             
`DAML_G23.csh` is shown here as an example, so just swap this out for another (e.g., `DAML_IFA.csh` or `DAML_OPTp2.csh`). To run a 2-step correction approach where first the CNN updates the state, and then we do DA on top, just add a call to the DART script below, as follows:

    <freInclude name="OM4_postprocess">
        <postProcess>
            <csh><![CDATA[ 
         
             cd $work/RESTART
             if (! -e coupler.res) then
                 echo model still running
             else
                 # FIRST DO CNN
                 cp DAML_OPTp1.csh .
                 csh DAML_OPTp1.csh >& log_CNN.out
                 rm log_CNN.out
                 rm DAML_OPTp1.csh

                 # NOW DO DA
                 cp data_assimilation_DART.csh .
                 csh data_assimilation_DART.csh >&log.out || exit 12
                 wait
             endif

             # rest of post-processing script follows as normal

To then run the model as a series of short forecasts the exectuable needs to be adjusted (this is done *after* running frerun, where the exectuable then gets saved in /lustre/f2/dev/.../scripts/run/). An example of which parts of the executable need to be changed in order to run a 5-year simulation, with the CNN correction being applied every 5 days is shown below:

    set -r segmentsPerSimulation = 365 #365 segements each running 5 days = 1825 days
    set segmentsPerPPCall = 0 
    set -r segmentsPerJob = 40
    set -r monthslist = ( 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 )
    set -r dayslist = ( 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 )
    set -r hourslist = ( 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 )
