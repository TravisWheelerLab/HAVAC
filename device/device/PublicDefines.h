#ifndef HAVAC_HLS_PUBLIC_DEFINES_H
#define HAVAC_HLS_PUBLIC_DEFINES_H

// this file contains definitions that may be useful for host code and test code, as well as device code.
//alternatively could be 64 for a 16B hit report

//#define HAVAC_TESTING
//#define HAVAC_PER_CELL_DATA_TESTING

#ifdef HAVAC_TESTING
#define NUM_CELL_GROUPS 64 //next 128
#define CELLS_PER_GROUP 4//4	//next 80
#define HAVAC_MAX_SUPPORTED_PHMM_LENGTH
#else
#define NUM_CELL_GROUPS 128		//128
#define CELLS_PER_GROUP 48		//48
#define HAVAC_MAX_SUPPORTED_PHMM_LENGTH
#endif
#define NUM_CELL_PROCESSORS (NUM_CELL_GROUPS *CELLS_PER_GROUP)


#define TEST_NUM_SEQUENCE_SEGMENTS 2
//deprecated, since I'm not going to support Amino Acids and you cant make me
#define IS_NUCLEOTIDE_IMPLEMENTATION

#ifdef IS_NUCLEOTIDE_IMPLEMENTATION
#define SCORES_PER_PHMM_VECTOR 4

#else
#define SCORES_PER_PHMM_VECTOR 20
#endif

#endif
