VPATH := $(SRC)

COMMONOBJ := arp.o arp_model.o arp_types.o

all: $(COMMONOBJ)
	mkdir -p $(INC)
	mkdir -p $(LIB)
	cp arp.mod $(INC)
	$(AR) rcs $(LIB)/libarp.a $(COMMONOBJ)

clean:
	rm -f *.mod *.o
	rm -f $(LIB)/libarp.a
	rm -f $(INC)/arp.mod

arp.o: arp.f08 arp_model.o arp_types.o
	$(FC) $(FFLAGS) -c -o $@ $<

arp_model.o: arp_model.f08 arp_types.o
	$(FC) $(FFLAGS) -c -o $@ $<

arp_types.o: arp_types.f08
	$(FC) $(FFLAGS) -c -o $@ $<
