#!/bin/bash

# This runs the run_files as default and outputs the log files for ISCE v.2.6.4.

run.py -i run_01_unpack_topo_reference -p 8 >& run_01_unpack_topo_reference.log
run.py -i run_02_unpack_secondary_slc -p 8 >& run_02_unpack_secondary_slc.log
run.py -i run_03_average_baseline -p 8 >& run_03_average_baseline.log
run.py -i run_04_extract_burst_overlaps -p 8 >& run_04_extract_burst_overlaps.log
run.py -i run_05_overlap_geo2rdr -p 8 >& run_05_overlap_geo2rdr.log
run.py -i run_06_overlap_resample -p 8 >& run_06_overlap_resample.log
run.py -i run_07_pairs_misreg -p 8 >& run_07_pairs_misreg.log
run.py -i run_08_timeseries_misreg -p 8 >& run_08_timeseries_misreg.log
run.py -i run_09_fullBurst_geo2rdr -p 8 >& run_09_fullBurst_geo2rdr.log
run.py -i run_10_fullBurst_resample -p 8 >& run_10_fullBurst_resample.log
run.py -i run_11_extract_stack_valid_region -p 8 >& run_11_extract_stack_valid_region.log
run.py -i run_12_merge_reference_secondary_slc -p 8 >& run_12_merge_reference_secondary_slc.log
run.py -i run_13_generate_burst_igram -p 8 >& run_13_generate_burst_igram.log
run.py -i run_14_merge_burst_igram -p 8 >& run_14_merge_burst_igram.log
run.py -i run_15_filter_coherence -p 8 >& run_15_filter_coherence.log
run.py -i run_16_unwrap -p 8 >& run_16_unwrap.log
