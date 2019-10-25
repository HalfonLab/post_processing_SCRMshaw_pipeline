

################################################################################												
   								  			
	`````````SCRMshaw post processing pipeline`````````````
				Halfon Lab			
  		     Date: Sept 2018		
  		      										
################################################################################

INTRODUCTION:
This script is written to do post processing steps on the SCRMshaw multiple offsets output including but not limited to use of peaks calling algorithm.


	1. INPUT
	2. USAGE
	3. PARAMETERS
	4. OUTPUT


1. INPUT:

The only file required for this script is the concatenated output(set of predictions) of SCRMshaw-HD (10 bp offset) and the number of predictions to take if more than 5000 (default).

2.  USAGE:

Pipeline requires MACs2 installed on the system.It will also check the status of MACs2 installation or will terminate with the error.
Following basic modules are required to run this script.Please make sure these modules have already been properly installed and are recognizable.
pybedtools, statistics, scipy, numpy, pandas, csv, subprocess, shutil and itertools.\
For running it on cluster, simply import the followings:

module load python/anaconda \
module load pybedtools/0.8.0 \
module load MACS2 

Following is an example of command line execution.
>python postProcessingScrmshawPipeline.py  -so SCRMshawConcatenatedOutputFile  -num 5000

3. PARAMETERS:

	-so	<str>	The output file for concatenation shell scripySCRMshaw\
	-num	<int>	Number of predictions to processs

4. OUTPUT:

Depending on user's input, if the input consists of multiple training sets the pipeline will create peaks file for each individual training set.
And a temporary folder which will have the files created during execution

