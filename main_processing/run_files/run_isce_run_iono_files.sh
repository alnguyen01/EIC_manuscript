#!/bin/bash

# This runs the run_files as default and outputs the log files for ISCE v.2.6.4.

run.py -i run_17_subband_and_resamp -p 8 >& run_17_subband_and_resamp.log
run.py -i run_18_generateIgram_ion -p 8 >& run_18_generateIgram_ion.log
run.py -i run_19_mergeBurstsIon -p 8 >& run_19_mergeBurstsIon.log
run.py -i run_20_unwrap_ion -p 4 >& run_20_unwrap_ion.log
run.py -i run_21_look_ion -p 8 >& run_21_look_ion.log
run.py -i run_22_computeIon -p 8 >& run_22_computeIon.log
run.py -i run_23_filtIon -p 8 >& run_23_filtIon.log
run.py -i run_24_invertIon -p 8 >& run_24_invertIon.log
