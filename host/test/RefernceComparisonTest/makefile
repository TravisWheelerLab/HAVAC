CXX = g++
CXXFLAGS = -std=c++17 -I/opt/xilinx/xrt/include -I../../../generator
LIBS = -L../build/ -L../../../build/ -lhavac -L/opt/xilinx/xrt/lib -lxrt_coreutil -luuid -lfastavector -lP7HmmReader

SRCS =  ../../../test/generator/hmmSeqGenerator.cpp ../../../PhmmReprojection/PhmmReprojection.cpp ReferenceComparisonTest.cpp Ssv.cpp
OBJS = $(SRCS:.cpp=.o)
EXEC = refCompare.out

all: $(EXEC)

$(EXEC): $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) $(LIBS) -o $(EXEC)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(EXEC)