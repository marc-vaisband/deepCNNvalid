#!/bin/bash
# Will download in- and out-of-facility data. For detailed requirements,
# consult the comments in the shell scripts mentioned below.

# Note that this will download roughly 1TB of sequencing data!

IN_FACILITY_PATH="./In_facility_mutation_calling_pipeline_mouse.sh"
KOTANI_PATH="./Out_of_facility_mutation_calling_pipeline_mouse.sh"
GENERATE_DATA_PATH="./generateData.py"

. $IN_FACILITY_PATH
. $KOTANI_PATH

python $GENERATE_DATA_PATH