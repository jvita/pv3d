OBJS = pv3d.o Grain.o Tools.o
CC = g++
DEBUG = -g
STD = -std=c++11
CFLAGS = -Wall -c $(DEBUG) $(STD)
LFLAGS = -Wall $(DEBUG) $(STD)

all: pv3d $(OBJS)

pv3d: $(OBJS)
	$(CC) $(LFLAGS) $(OBJS) -o pv3d

pv3d.o: pv3d.cpp Grain.cpp Tools.cpp Grain.h Tools.h
	$(CC) $(CFLAGS) pv3d.cpp

Grain.o: Grain.cpp Tools.cpp Tools.h
	$(CC) $(CFLAGS) Grain.cpp

Tools.o: Tools.cpp Tools.h define.h
	$(CC) $(CFLAGS) Tools.cpp
