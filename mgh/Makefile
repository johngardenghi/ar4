all: mgh.o

mgh.o: mgh.f08 set_precision.o
	$(FC) $(FFLAGS) -c -o $@ $<

set_precision.o: set_precision.f08
	$(FC) $(FFLAGS) -c -o $@ $<

clean:
	rm -f mgh.o mgh.mod
	rm -f set_precision.o set_precision.mod
