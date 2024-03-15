# Define the compiler
CXX=g++

# Define any compile-time flags
CXXFLAGS=-Wall -g -std=c++20

# Define the directories to include
INCLUDES=-I. -Iboundary -Istate

# Define the C++ source files
SOURCES=$(wildcard *.cpp) \
        $(wildcard boundary/*.cpp) \
        $(wildcard state/*.cpp)

# Define the object files
OBJS=$(SOURCES:.cpp=.o)

# Define the executable file 
MAIN=lp

.PHONY: clean

all:    $(MAIN)
	@echo Program has been compiled

$(MAIN): $(OBJS) 
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $(MAIN) $(OBJS)

# Pattern rule for object files, taking into account source files
# in subdirectories
%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

clean:
	$(RM) *.o boundary/*.o state/*.o *~ $(MAIN)

