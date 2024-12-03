
CXXFLAGS ?= -O3 -pipe 

.PHONY: analysis

RootLibs:=$(shell root-config --libs)
RootInclude:=$(shell root-config --incdir)

analysis: analysis.o
	g++ -std=c++17 -O3 -o $@ $^ $(RootLibDir) $(RootLibs) -I./* 

analysis.o: analysis.cpp
	g++ -std=c++17 -O3 -c $< -I $(RootInclude) -I./*
clean:
	rm -rf analysis *.o
	echo Clean done
