Due to the nature of HAVAC being based on an accelerator card, and comparing our results to the gold-standard SSV implementation that lives inside the HMMER codebase, generating the benchmarks used for this project is not a simple task. This readme outlines the steps that must be taken to reproduce this data.

# Hardware Installation
As outline in the main project's readme, the first step towards generating these benchmarks is installing an Alveo U50 data accelerator card, and the corresponding firmware and libraries. Please consult the AMD documentation currently located at https://docs.xilinx.com/v/u/en-US/ug1370-u50-installation. This also outlines the process for installing the XRT framework on compatible Linux systems.

Once the hardware is installed along with the XRT framework, and once the installation has been validated, the HAVAC and HMMER softwares must be built. please follow the main readme's instructions for building the HAVAC driver code.

Install the HAVAC driver. from the projects main directory:
	% git submodule init
	% git submodule update
	% cmake .
	% make
	% make install


# HMMER software, and modification to test the SSV algorithm
Unfortunately, HMMER does not have a native way to run only the SSV algorithm. In order to test against the reference implementation, we need to modify the HMMER code slightly to avoid running the downstream pipeline steps.

To begin, clone HMMER from its repository at https://github.com/EddyRivasLab/hmmer.

in the src directory, open p7_pipeline.c. After the call to p7_SSVFilter_longtarget() on line 1549, add a return statement. This bypasses all computation after the initial SSV step of the pipeline. 

Once this modification is complete, run    
	% ./configure --prefix /your/install/path   # replace /your/install/path with what you want, obv 
	% make
	
this will generate the nhmmer binaries inside the HMMER/src directory. The hmmbuild tool will also be generated in the HMMER/src directory. this will be used to generate the .hmm file for the full HMMER database


# Datasets
HAVAC is benchmarked against HMMER using the RFAM database of RNA families against human chromosome 22. 
A fasta of chromosome 22 should be used to reproduce the results from the publication. An example location where a user could download chromsome 22 is from https://hgdownload.soe.ucsc.edu/goldenPath/bosTau2/chromosomes/chr22.fa.gz.


To generate the HMM file containing all of RFAM, download the database in stockholm format at
	https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.seed.gz
extract the .gz file, and use the hmmbuild tool located inside HMMER/src as follows:
	% hmmbuild rfam.hmm Rfam.seed


Once the hmm file containing all models in the RFAM database is generated, the database must be split into sections of various lengths, starting with the lowest accession ID (RF00000), including families in ascending order of accession ID up to specified databases sizes. A Python script is included to perform this operation. In this (benchmarks) directory, run the script
	% hmmDbByLength.py path/to/rfam.hmm
	
to generate the benchmarking model databases. For your convenience, this script also generates a file containing the actual lengths of the model databases generated. The script concatenates the models up to the requested lengths, so the actual lengths will be usually be slightly over the requested length. This file generates the list in the form of a line that can be copied into a python script for visualization.

# Generating the Benchmark Results
## HAVAC

To generate the benchmark results,  enter this (benchmarks) directory and make the benchmark
	% make
	
Then, run the HAVAC benchmark executable refCompare.out, with the following arguments
	1. xclbin src (found by default in device/bin/havac.xclbin)
	2. chr22 fasta src
	3. RFAM section file src (generated from the hmmDbByLength script above)

This benchmark executable prints out the time taken by the various stages of the HAVAC driver. This should be run with /usr/bin/time to time the entire executable, e.g.,
	% /usr/bin/time refCompare.out ../device/bin/havac.xclbin chr22.fa rfam_10000.hmm
	
## HMMER
To test the HMMER timings, use /usr/bin/time on the nhmmer executable found in the HMMER/src directory, e.g.,
	% /usr/bin/time HMMER/src/nhmmer --watson --cpu 32 rfam_10000.hmm chr22.fa
	
the --watson flag forces nhmmer to only compare against one strand and not the reverse-compliment.
	
# Visualization script
A visualization script is included to generate the benchmark figure from the publication. runtime_table.py will generate figure 7 from the publication. It is up to the user to insert the timing data from their own run of the benchmarks into the script.
