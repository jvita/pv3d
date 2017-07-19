# TODO: fix dependencies, update at proper time
all:
	cd src/ && $(MAKE)
	cd tests/ && $(MAKE)

clean:
	rm pv3d
	rm test_suite
	cd src/ && \rm *.o
