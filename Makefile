DEBUG ?= 0
ifeq ($(DEBUG), 0)
        CXXFLAGS +=-DNDEBUG -O3 -g
else ifeq ($(DEBUG), 2)
        CXXFLAGS +=-DNDEBUG -O3
else
        CXXFLAGS +=-DDEBUG -O0 -g
endif

all: orbit fullperm nine-point-graphs

orbit: orbit.cc
	g++ $(CXXFLAGS) -o orbit orbit.cc

fullperm: fullperm.cc
	g++ $(CXXFLAGS) -o fullperm fullperm.cc

nine-point-graphs: nine-point-graphs.cc
	g++ $(CXXFLAGS) -o nine-point-graphs nine-point-graphs.cc

clean:
	rm -f orbit fullperm nine-point-graphs
