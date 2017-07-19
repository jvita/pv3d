all:
	cd src/ && $(MAKE)
	mv src/pv3d .
	cd tests/ && $(MAKE)

clean:
	rm pv3d
	cd src/ && \rm *.o
	cd tests/ && \rm tests
