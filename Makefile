pv3d: pv3d.o Tools.o Grain.o
	g++ -Wall -std=c++11 pv3d.o Tools.o Grain.o -o pv3d

pv3d.o: pv3d.cpp Tools.cpp Tools.h Grain.cpp Grain.h
	g++ -Wall -c -std=c++11 pv3d.cpp

Grain.o: Grain.cpp Tools.cpp Tools.h
	g++ -Wall -c -std=c++11 Grain.cpp

Tools.o: Tools.cpp Tools.h define.h
	g++ -Wall -c -std=c++11 Tools.cpp

clean:
	\rm *.o pv3d
