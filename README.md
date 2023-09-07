# HAVAC
Hardware Accellerated Viterbi Additional Coprocessor (HAVAC) is an FPGA-accellerated implementation of the Single-segment Ungapped Viterbi algorithm for use in nucleotide sequence homology search. The HAVAC device design is implemented with Vitis HLS; the design code is provided, but should not need to be used unless you want to generate an implementation for a different accelerator card. 

HAVAC is designed to be implemented on the Xilinx Alveo U50 Datacenter accelerator card; The FPGA configuration binary (.xclbin file) is provided in the device/bin directory. This design implements 12288 cell processing elements (PEs) at 145MHz and is able to achieve matrix calculation speeds up to 1,739 billion cells updated per second (GCUPS).


## Getting Started
To get started on the Alveo U50 Datacenter accelerator card, consult the Xilinx Alveo U50 Data Center Accelerator Card Data Sheet (at time of writing, found at https://www.xilinx.com/content/dam/xilinx/support/documents/data_sheets/ds965-u50.pdf). Once installed, generate the driver executable, and HAVAC can be invoked with the provided .xclbin file.

### Building the Software Driver

Building the software driver requires the following
* GCC C++ compiler
* Cmake Version 3.10 or higher
* Xilinx Runtime (XRT) framework

To install XRT, consult Xilinx UG1939 (https://docs.xilinx.com/r/en-US/ug1393-vitis-application-acceleration/Installing-Xilinx-Runtime-and-Platforms)

After cloning the repo, the user should initialize and update the project submodules.

```
git clone https://github.com/TravisWheelerLab/HAVAC.git
git submodule init
git submodule update
```

Now the project is ready to build and install.
```
cmake .
make
make install
```

### Searching with HAVAC

To use the HAVAC driver, create an object of the Havac type with your device index (usually 0), the p-value to search with, and the src for the .xclbin file.
```
uint32_t deviceIndex = 0; //this may be a different index if you have multiple FPGA cards.
float pValue = 0.02; //this is the p-value used by nhmmer's SSV filter
std::string xclbinSrc = "../device/bin/havac.xclbin";
Havac *havacDriver = new Havac(deviceIndex, pValue, xclbinSrc); 
```

This will configure the accelerator device, and set the p-value used to project pHMM data for correct use on the device. Then, load the database sequences and the model collection to the device. 
```
havacDriver->loadSequence("path/to/sequences.fasta");
havacDriver->loadPhmm("path/to/models.hmm);
```

this will load all sequences in the fasta and all models in the .hmm file to the device for search. Note that HAVAC supports HMMER3 model files.

HAVAC can then be invoked to search the sequence collection with all the models. This can be done synchronously
```
havacDriver->runHardwareClient()
```

to force the driver to wait until the hardware is done. Alternatively, the hardware can be invoked asynchronously:
```
havacDriver->runHardwareClientAsync();
```

While running asynchronously, the state of the device can be accesed with 
```
enum havac_cmd_state currentState = havacDriver->currentHardwareState();
//   HAVAC_CMD_STATE_NEW = 1, 
//   HAVAC_CMD_STATE_QUEUED = 2,
//   HAVAC_CMD_STATE_RUNNING = 3,
//   HAVAC_CMD_STATE_COMPLETED = 4,
//   HAVAC_CMD_STATE_ERROR = 5,
//   HAVAC_CMD_STATE_ABORT = 6,
//   HAVAC_CMD_STATE_SUBMITTED = 7,
//   HAVAC_CMD_STATEIMEOUT = 8,
//   HAVAC_CMD_STATE_NORESPONSE = 9
```

While running asynchronously, the abort function may be invoked to stop processing.
```
havacDriver->abortHardwareClient();
```

To wait until the hardware is finished, invoke the async wait function. Note that must be invoked before reading the hit reports when running asynchronously.
```
havacDriver->waitHardwareClientAsync();
```

After the hardware is finished with its computation, the hits can be recovered via 
```
vector<HavacHit> reportedHits = havacDriver->getHitsFromFinishedRun()
```

This function resolves the hits recovered from the accelerator device to the index of the sequences and models in their corresponding files, and identifies the local position within those sequences and models. Here's what the hit reports will look like:
```
class HavacHit {
public:
  HavacHit(const uint64_t sequencePosition, const uint32_t sequenceIndex, const uint32_t phmmPosition, const uint32_t phmmIndex);
  uint64_t sequencePosition;
  uint32_t sequenceIndex;
  uint32_t phmmPosition;
  uint32_t phmmIndex;
  std::string toString();
};
```