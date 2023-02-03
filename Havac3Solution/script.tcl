############################################################
## This file is generated automatically by Vitis HLS.
## Please DO NOT edit it.
## Copyright 1986-2021 Xilinx, Inc. All Rights Reserved.
############################################################
open_project Havac3
set_top HavacKernelTopLevel
add_files Havac3/device/PublicDefines.h
add_files Havac3/device/HavacHls.hpp
add_files Havac3/device/HavacHls.cpp
add_files -tb Havac3/FastaVector/src/FastaVector.c -cflags "-Wno-unknown-pragmas" -csimflags "-Wno-unknown-pragmas"
add_files -tb Havac3/FastaVector/src/FastaVectorMetadataVector.c -cflags "-Wno-unknown-pragmas" -csimflags "-Wno-unknown-pragmas"
add_files -tb Havac3/FastaVector/src/FastaVectorString.c -cflags "-Wno-unknown-pragmas" -csimflags "-Wno-unknown-pragmas"
add_files -tb Havac3/PhmmReprojection/PhmmReprojection.cpp -cflags "-Wno-unknown-pragmas" -csimflags "-Wno-unknown-pragmas"
add_files -tb Havac3/test/softSsv/SoftSsv.cpp -cflags "-Wno-unknown-pragmas" -csimflags "-Wno-unknown-pragmas"
add_files -tb Havac3/test/byCellComparator/byCellComparator.cpp -cflags "-Wno-unknown-pragmas" -csimflags "-Wno-unknown-pragmas"
add_files -tb Havac3/test/generator/hmmSeqGenerator.cpp -cflags "-Wno-unknown-pragmas" -csimflags "-Wno-unknown-pragmas"
add_files -tb Havac3/test/generator/hmmSeqGenerator.h -cflags "-Wno-unknown-pragmas" -csimflags "-Wno-unknown-pragmas"
add_files -tb Havac3/P7HmmReader/src/p7HmmReader.c -cflags "-Wno-unknown-pragmas" -csimflags "-Wno-unknown-pragmas"
add_files -tb Havac3/P7HmmReader/src/p7HmmReaderLog.c -cflags "-Wno-unknown-pragmas" -csimflags "-Wno-unknown-pragmas"
add_files -tb Havac3/P7HmmReader/src/p7ProfileHmm.c -cflags "-Wno-unknown-pragmas" -csimflags "-Wno-unknown-pragmas"
add_files -tb Havac3/test/softwareTestbench.cpp -cflags "-Wno-unknown-pragmas" -csimflags "-Wno-unknown-pragmas"
add_files -tb Havac3/test/softwareTestbench.hpp -cflags "-Wno-unknown-pragmas" -csimflags "-Wno-unknown-pragmas"
open_solution "Havac3Solution" -flow_target vivado
set_part {xcu50-fsvh2104-2-e}
create_clock -period 2.65 -name default
source "./Havac3/Havac3Solution/directives.tcl"
csim_design -clean
csynth_design
cosim_design -trace_level all -argv {-v}
export_design -format ip_catalog
