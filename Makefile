# TODO: fix dependencies, update at proper time
all:
	cd src/ && $(MAKE)
	cd tests/ && $(MAKE)

clean:
	rm src/pv3d
	rm tests/test_suite
	cd src/ && \rm *.o
