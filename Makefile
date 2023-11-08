all: orbit fullperm

orbit: orbit.cc
	g++ -O3 -o orbit orbit.cc

fullperm: fullperm.cc
	g++ -O3 -o fullperm fullperm.cc

clean:
	rm orbit fullperm
