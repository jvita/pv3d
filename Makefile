pv3d: pv3d.o Tools.o
	g++ -Wall -std=c++11 pv3d.o Tools.o -o pv3d

pv3d.o: pv3d.cpp Tools.cpp Tools.h
	g++ -Wall -c std=c++11 pv3d.cpp

Tools.o: Tools.cpp Tools.h define.h
	g++ -Wall -c -std=c++11 Tools.cpp
