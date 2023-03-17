CXX = g++
CXXFLAGS = -std=c++20
#havac static library should be in the host build directory
LIBS = -L../build -lhavac

SRCS = ReferenceComparisonTest.cpp Ssv.cpp
OBJS = $(SRCS:.cpp=.o)
EXEC = refCompare.out

all: $(EXEC)

$(EXEC): $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) $(LIBS) -o $(EXEC)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(EXEC)