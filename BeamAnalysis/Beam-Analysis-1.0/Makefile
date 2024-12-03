
CXXFLAGS ?= -O3 -pipe 

.PHONY: all

RootLibs:=$(shell root-config --libs)
RootInclude:=$(shell root-config --incdir)

all: analyse_tuple 

analyse_tuple: utils.o buffer.o analyse_tuple.o 
	g++ -std=c++17 -O3 -o $@ $^ $(RootLibDir) $(RootLibs) -I./* 

buffer.o: buffer.cpp
	g++ -std=c++17 -O3 -c $< -I $(RootInclude) -I./*
utils.o: utils.cpp
	g++ -std=c++17 -O3 -c $< -I $(RootInclude) -I./*
analyse_tuple.o: analyse_tuple.cpp
	g++ -std=c++17 -O3 -c $< -I $(RootInclude) -I./*

clean:
	rm -rf analyse_tuple *.o
	echo Clean done
