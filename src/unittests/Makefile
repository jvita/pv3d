INCLUDE = /home/jvita/scripts/pv3d/src/
OBJS = Pv3d.o Grain.o Tools.o Lammps.o
TESTS = ToolsTest.cpp
VPATH = ../src

CC = g++
DEBUG = -g
STD = -std=c++11
CFLAGS = -Wall -c $(DEBUG) $(STD)
LFLAGS = -Wall $(DEBUG) $(STD)

test_suite: $(TESTS) $(OBJS)
	export CPATH=$(INCLUDE)
	echo ${CPATH}
	g++ -Wall -std=c++11 $(TESTS) -lUnitTest++ -o $(OBJS)test_suite
