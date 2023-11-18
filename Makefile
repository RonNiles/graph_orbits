all: orbit fullperm nine-point-graphs

orbit: orbit.cc
	g++ -O3 -o orbit orbit.cc

fullperm: fullperm.cc
	g++ -O3 -o fullperm fullperm.cc

nine-point-graphs: nine-point-graphs.cc
	g++ -O3 -o nine-point-graphs nine-point-graphs.cc

clean:
	rm orbit fullperm nine-point-graphs
