# Something weird happened when edited Tools.cpp and did 'make tests'
OBJS = $(patsubst %.cpp, obj/%.o, $(wildcard *.cpp))
TEST_OBJS = $(patsubst %.cpp, %.o, $(wildcard unittests/*.cpp)) \
			$(patsubst %.o, obj/%.o, Tools.o Lammps.o Grain.o)

CC = g++
DEBUG = -g
STD = -std=c++11

CFLAGS = -Wall -c $(DEBUG) $(STD)
LFLAGS = -Wall $(DEBUG) $(STD)

INCLUDE = -I ./

all: pv3d tests

pv3d: $(OBJS)
	$(CC) $(LFLAGS) $(OBJS) -o pv3d

tests: $(OBJS) $(TEST_OBJS)
	$(CC) $(LFLAGS) $(TEST_OBJS) -lUnitTest++ -o tests

unittests/%.o: unittests/%.cpp
	$(CC) $(CFLAGS) $(INCLUDE) $< -lUnitTest++ -o $@

# $@ is target, $< is first dependency, $^ is all dependencies
obj/%.o: %.cpp | obj
	$(CC) $(CFLAGS) $< -o $@

obj:
	mkdir -p obj

clean:
	rm -rf obj/
	rm -f unittests/*.o
	rm -f pv3d
	rm -f tests
