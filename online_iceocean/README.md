### Implementation ###

CNN correction is applied back into SPEAR ice-ocean simulations as a bias correction tool. This is implemented as post-processing to a short model forecast. E.g. in an xml script, this would be implemented in the following way:

<freInclude name="OM4_postprocess">
      <postProcess>
         <csh><![CDATA[ 
         
         #Implementing DA-based CNN correction scheme:                                                                                                                               
         cd $work/RESTART                                                                                                                                                              
         if (! -e coupler.res) then                                                                                                                                                    
             echo model still running                                                                                                                                                  
         else                                                                                                                                                                          
             cp DAML.csh .                                                                                                                                                             
             csh DAML.csh >& log_CNN.out                                                                                                                                               
             rm log_CNN.out                                                                                                                                                            
             rm DAML.csh                                                                                                                                                               
         endif
