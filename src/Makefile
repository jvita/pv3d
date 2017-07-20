OBJS = Pv3d.o Grain.o Tools.o Lammps.o
VPATH = /src
CC = g++
DEBUG = -g
STD = -std=c++11
CFLAGS = -Wall -c $(DEBUG) $(STD)
LFLAGS = -Wall $(DEBUG) $(STD)

pv3d: $(OBJS)
	$(CC) $(LFLAGS) $(OBJS) -o pv3d

Pv3d.o: Pv3d.cpp Pv3d.h Grain.cpp Tools.cpp Grain.h Tools.h define.h
	$(CC) $(CFLAGS) Pv3d.cpp

Grain.o: Grain.cpp Tools.cpp Tools.h define.h
	$(CC) $(CFLAGS) Grain.cpp

Tools.o: Tools.cpp Tools.h define.h
	$(CC) $(CFLAGS) Tools.cpp

Lammps.o: Lammps.cpp Lammps.h
	$(CC) $(CFLAGS) Lammps.cpp
