CXX = g++
CXXFLAGS = -std=c++17 -I/opt/xilinx/xrt/include -I../../../generator -Wall -Wextra -O3
CXXFLAGS_DEBUG = -ggdb -fsanitize=address -fno-omit-frame-pointer -static-libstdc++ -static-libasan -lrt -g -O0
LIBS =  -L../build/ -L../../../build/ -lhavac -L/opt/xilinx/xrt/lib -lxrt_coreutil -luuid -lfastavector \
	-lP7HmmReader -Wl,-rpath,/opt/xilinx/xrt/lib -L/usr/local/lib

SRCS = ../../../test/generator/hmmSeqGenerator.cpp  \
	ReferenceComparisonTest.cpp ../Ssv.cpp

OBJS = $(SRCS:.cpp=.o)
EXEC = refCompare.out


all: $(EXEC)

$(EXEC): $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) $(LIBS) -o $(EXEC)

debug: CXXFLAGS += $(CXXFLAGS_DEBUG)
debug: clean $(EXEC) 

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@


clean:
	rm -f $(OBJS) $(EXEC)