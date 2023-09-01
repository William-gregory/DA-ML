### Implementation of CNN into SPEAR ice-ocean simulations

The CNN correction is applied back into SPEAR ice-ocean simulations as a bias correction tool. In other words, we run a short simulation/forecast of the model, e.g., for 5 days, at which point the CNN updates the restart files as a post-processing step, and then the model continues from this updated state for the next 5 days. Note that the `CNN_SIS2.py' code will need to be changed depending on what CNN weights you want to use

The post-processing step to implement a CNN into an xml is done in the following way:

    <freInclude name="OM4_postprocess">
        <postProcess>
            <csh><![CDATA[ 
         
             # implementing DA-based CNN correction scheme:
             cd $work/RESTART
             if (! -e coupler.res) then
                 echo model still running
             else
                 cp DAML.csh .
                 csh DAML.csh >& log_CNN.out
             endif
             
             # rest of post-processing script follows as normal
             
To run a 2-step correction approach where first the CNN updates the state, and then we do DA on top, just add a call to the DART script below, as follows:

    <freInclude name="OM4_postprocess">
        <postProcess>
            <csh><![CDATA[ 
         
             cd $work/RESTART
             if (! -e coupler.res) then
                 echo model still running
             else
                 # FIRST DO CNN
                 cp DAML.csh .
                 csh DAML.csh >& log_CNN.out

                 # NOW DO DA
                 cp data_assimilation_DART.csh .
                 csh data_assimilation_DART.csh >& log.out || exit 12
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
