LIBRARY_DIR=/home/ITensor/
include $(LIBRARY_DIR)/this_dir.mk
include $(LIBRARY_DIR)/options.mk
#LINUXWSLIBS   = -L/usr/include/X11/lib -lxview -lolgx -lX11 -lm -lf2c -lgfortran
#Define Flags ----------

TENSOR_HEADERS=$(PREFIX)/itensor/all.h
CCFLAGS= -I. $(ITENSOR_INCLUDEFLAGS) $(CPPFLAGS) $(OPTIMIZATIONS)
CCGFLAGS= -I. $(ITENSOR_INCLUDEFLAGS) $(DEBUGFLAGS)
LIBFLAGS=-L'$(ITENSOR_LIBDIR)' $(ITENSOR_LIBFLAGS)
LIBGFLAGS=-L'$(ITENSOR_LIBDIR)' $(ITENSOR_LIBGFLAGS)

#Rules ------------------

%.o: %.cc $(ITENSOR_LIBS) $(TENSOR_HEADERS)
	$(CCCOM) -c $(CCFLAGS) -o $@ $<

.debug_objs/%.o: %.cc $(ITENSOR_GLIBS) $(TENSOR_HEADERS)
	$(CCCOM) -c $(CCGFLAGS) -o $@ $<

#Targets -----------------

build: dmrg dmrg_table dmrgj1j2 exthubbard trg ctmrg mixedspin hubbard_2d 2channeldmrgthird hubbard_2d_conserve_momentum helloitensor 2channeldmrg 2channeldmrgsecond 1channeldmrg 2channeldmrgfourth 2channeldmrgfifth 2channeldmrgsixth 2channeldmrgseventh 2channeldmrgeighth try2 2channeldmrgninth 2channeldmrgtenth 2channeldmrgeleventh 2channeldmrgtwelveth 2channeldmrgthirteenth 2channeldmrgfourteenth 2channeldmrgfifteenth 2channeldmrgsixteenth dmrggrover flatband1 kitaev kitaev2 kitaev3 lieblattice try1 liebhubbard liebhubbard2 liebhubbard3 liebhubbard4 liebhubbard5 liebhubbard6 liebhubbard7 liebhuabbard9 liebkitaevread liebkitaevread2 lieblattice2 liebhubbrad10 liebhubbard11 kitaevhubbard

debug: dmrg-g dmrg_table-g dmrgj1j2-g exthubbard-g trg-g ctmrg-g mixedspin-g hubbard_2d-g hubbard_2d_conserve_momentum-g flatband1

all: dmrg dmrg_table dmrgj1j2 exthubbard trg ctmrg mixedspin hubbard_2d hubbard_2d_conserve_momentum helloitensor 2channeldmrgthird 2channeldmrg 2channeldmrgsecond 1channeldmrg 2channeldmrgfourth 2channeldmrgfifth 2channeldmrgsixth 2channeldmrgseventh 2channeldmrgeighth try2 2channeldmrgninth 2channeldmrgtenth 2channeldmrgeleventh 2channeldmrgtwelveth 2channeldmrgthirteenth 2channeldmrgfourteenth 2channeldmrgfifteenth 2channeldmrgsixteenth dmrggrover flatband1 kitaev kitaev2 kitaev3 lieblattice try1 liebhubbard liebhubbard2 liebhubbard3 liebhubbard4 liebhubbard5 liebhubbard6 liebhubbard7 liebhubbard9 liebkitaevread liebkitaevread2 lieblattice2 liebhubbard10 liebhubbard11 kitaevhubbard

dmrg: dmrg.o $(ITENSOR_LIBS) $(TENSOR_HEADERS)
	$(CCCOM) $(CCFLAGS) dmrg.o -o dmrg $(LIBFLAGS)

liebkitaevread: liebkitaevread.o $(ITENSOR_LIBS) $(TENSOR_HEADERS)
	$(CCCOM) $(CCFLAGS) liebkitaevread.o -o liebkitaevread $(LIBFLAGS)

liebkitaevread2: liebkitaevread2.o $(ITENSOR_LIBS) $(TENSOR_HEADERS)
	$(CCCOM) $(CCFLAGS) liebkitaevread2.o -o liebkitaevread2 $(LIBFLAGS)

kitaevhubbard: kitaevhubbard.o $(ITENSOR_LIBS) $(TENSOR_HEADERS)
	$(CCCOM) $(CCFLAGS) kitaevhubbard.o -o kitaevhubbard $(LIBFLAGS)

liebhubbard: liebhubbard.o $(ITENSOR_LIBS) $(TENSOR_HEADERS)
	$(CCCOM) $(CCFLAGS) liebhubbard.o -o liebhubbard $(LIBFLAGS)

liebhubbard2: liebhubbard2.o $(ITENSOR_LIBS) $(TENSOR_HEADERS)
	$(CCCOM) $(CCFLAGS) liebhubbard2.o -o liebhubbard2 $(LIBFLAGS)

liebhubbard3: liebhubbard3.o $(ITENSOR_LIBS) $(TENSOR_HEADERS)
	$(CCCOM) $(CCFLAGS) liebhubbard3.o -o liebhubbard3 $(LIBFLAGS)

liebhubbard4: liebhubbard4.o $(ITENSOR_LIBS) $(TENSOR_HEADERS)
	$(CCCOM) $(CCFLAGS) liebhubbard4.o -o liebhubbard4 $(LIBFLAGS)

liebhubbard5: liebhubbard5.o $(ITENSOR_LIBS) $(TENSOR_HEADERS)
	$(CCCOM) $(CCFLAGS) liebhubbard5.o -o liebhubbard5 $(LIBFLAGS)

liebhubbard6: liebhubbard6.o $(ITENSOR_LIBS) $(TENSOR_HEADERS)
	$(CCCOM) $(CCFLAGS) liebhubbard6.o -o liebhubbard6 $(LIBFLAGS)

liebhubbard7: liebhubbard7.o $(ITENSOR_LIBS) $(TENSOR_HEADERS)
	$(CCCOM) $(CCFLAGS) liebhubbard7.o -o liebhubbard7 $(LIBFLAGS)
	
liebhubbard9: liebhubbard9.o $(ITENSOR_LIBS) $(TENSOR_HEADERS)
	$(CCCOM) $(CCFLAGS) liebhubbard9.o -o liebhubbard9 $(LIBFLAGS)	

liebhubbard10: liebhubbard10.o $(ITENSOR_LIBS) $(TENSOR_HEADERS)
	$(CCCOM) $(CCFLAGS) liebhubbard10.o -o liebhubbard10 $(LIBFLAGS)	

liebhubbard11: liebhubbard11.o $(ITENSOR_LIBS) $(TENSOR_HEADERS)
	$(CCCOM) $(CCFLAGS) liebhubbard11.o -o liebhubbard11 $(LIBFLAGS)	

try1: try1.o $(ITENSOR_LIBS) $(TENSOR_HEADERS)
	$(CCCOM) $(CCFLAGS) try1.o -o try1 $(LIBFLAGS)

kitaev: kitaev.o $(ITENSOR_LIBS) $(TENSOR_HEADERS)
	$(CCCOM) $(CCFLAGS) kitaev.o -o kitaev $(LIBFLAGS)

kitaev2: kitaev2.o $(ITENSOR_LIBS) $(TENSOR_HEADERS)
	$(CCCOM) $(CCFLAGS) kitaev2.o -o kitaev2 $(LIBFLAGS)

kitaev3: kitaev3.o $(ITENSOR_LIBS) $(TENSOR_HEADERS)
	$(CCCOM) $(CCFLAGS) kitaev3.o -o kitaev3 $(LIBFLAGS)

lieblattice: lieblattice.o $(ITENSOR_LIBS) $(TENSOR_HEADERS)
	$(CCCOM) $(CCFLAGS) lieblattice.o -o lieblattice $(LIBFLAGS)

lieblattice2: lieblattice2.o $(ITENSOR_LIBS) $(TENSOR_HEADERS)
	$(CCCOM) $(CCFLAGS) lieblattice2.o -o lieblattice2 $(LIBFLAGS)

1channeldmrg: 1channeldmrg.o $(ITENSOR_LIBS) $(TENSOR_HEADERS)
	$(CCCOM) $(CCFLAGS) 1channeldmrg.o -o 1channeldmrg $(LIBFLAGS)

helloitensor: helloitensor.o $(ITENSOR_LIBS) $(TENSOR_HEADERS)
	$(CCCOM) $(CCFLAGS) helloitensor.o -o helloitensor $(LIBFLAGS)

2channeldmrg: 2channeldmrg.o $(ITENSOR_LIBS) $(TENSOR_HEADERS)
	$(CCCOM) $(CCFLAGS) 2channeldmrg.o -o 2channeldmrg $(LIBFLAGS\
		)

2channeldmrgthird: 2channeldmrgthird.o $(ITENSOR_LIBS) $(TENSOR_HEADERS)
	$(CCCOM) $(CCFLAGS) 2channeldmrgthird.o -o 2channeldmrgthird $(LIBFLAGS)

2channeldmrgthird-g: mkdebugdir .debug_objs/2channeldmrgthird.o $(ITENSOR_GLIBS) $(TENSOR_HEADERS)
	$(CCCOM) $(CCGLAGS) .debug_objs/2channeldmrgthird.o -o 2channeldmrgthird-g $(LIBGFLAGS)

2channeldmrgfourth: 2channeldmrgfourth.o $(ITENSOR_LIBS) $(TENSOR_HEADERS)
	$(CCCOM) $(CCFLAGS) 2channeldmrgfourth.o -o 2channeldmrgfourth $(LIBFLAGS)

2channeldmrgfourth-g: mkdebugdir .debug_objs/2channeldmrgfifth.o $(ITENSOR_GLIBS) $(TENSOR_HEADERS)
	$(CCCOM) $(CCGLAGS) .debug_objs/2channeldmrgfifth.o -o 2channeldmrgfifth-g $(LIBGFLAGS)

2channeldmrgfifth: 2channeldmrgfifth.o $(ITENSOR_LIBS) $(TENSOR_HEADERS)
	$(CCCOM) $(CCFLAGS) 2channeldmrgfifth.o -o 2channeldmrgfifth $(LIBFLAGS)

2channeldmrgfifth-g: mkdebugdir .debug_objs/2channeldmrgfifth.o $(ITENSOR_GLIBS) $(TENSOR_HEADERS)
	$(CCCOM) $(CCGLAGS) .debug_objs/2channeldmrgfifth.o -o 2channeldmrgfifth-g $(LIBGFLAGS)

2channeldmrgsixth: 2channeldmrgsixth.o $(ITENSOR_LIBS) $(TENSOR_HEADERS)
	$(CCCOM) $(CCFLAGS) 2channeldmrgsixth.o -o 2channeldmrgsixth $(LIBFLAGS)

2channeldmrgsixth-g: mkdebugdir .debug_objs/2channeldmrgsixth.o $(ITENSOR_GLIBS) $(TENSOR_HEADERS)
	$(CCCOM) $(CCGLAGS) .debug_objs/2channeldmrgsixth.o -o 2channeldmrgsixth-g $(LIBGFLAGS)

2channeldmrgseventh: 2channeldmrgseventh.o $(ITENSOR_LIBS) $(TENSOR_HEADERS)
	$(CCCOM) $(CCFLAGS) 2channeldmrgseventh.o -o 2channeldmrgseventh $(LIBFLAGS)

2channeldmrgseventh-g: mkdebugdir .debug_objs/2channeldmrgseventh.o $(ITENSOR_GLIBS) $(TENSOR_HEADERS)
	$(CCCOM) $(CCGLAGS) .debug_objs/2channeldmrgseventh.o -o 2channeldmrgseventh-g $(LIBGFLAGS)

2channeldmrgeighth: 2channeldmrgeighth.o $(ITENSOR_LIBS) $(TENSOR_HEADERS)
	$(CCCOM) $(CCFLAGS) 2channeldmrgeighth.o -o 2channeldmrgeighth $(LIBFLAGS)

2channeldmrgeighth-g: mkdebugdir .debug_objs/2channeldmrgeighth.o $(ITENSOR_GLIBS) $(TENSOR_HEADERS)
	$(CCCOM) $(CCGLAGS) .debug_objs/2channeldmrgeighth.o -o 2channeldmrgeighth-g $(LIBGFLAGS)

2channeldmrgninth: 2channeldmrgninth.o $(ITENSOR_LIBS) $(TENSOR_HEADERS)
	$(CCCOM) $(CCFLAGS) 2channeldmrgninth.o -o 2channeldmrgninth $(LIBFLAGS)

2channeldmrgninth-g: mkdebugdir .debug_objs/2channeldmrgninth.o $(ITENSOR_GLIBS) $(TENSOR_HEADERS)
	$(CCCOM) $(CCGLAGS) .debug_objs/2channeldmrgninth.o -o 2channeldmrgninth-g $(LIBGFLAGS)

2channeldmrgtenth: 2channeldmrgtenth.o $(ITENSOR_LIBS) $(TENSOR_HEADERS)
	$(CCCOM) $(CCFLAGS) 2channeldmrgtenth.o -o 2channeldmrgtenth $(LIBFLAGS)

2channeldmrgtenth-g: mkdebugdir .debug_objs/2channeldmrgtenth.o $(ITENSOR_GLIBS) $(TENSOR_HEADERS)
	$(CCCOM) $(CCGLAGS) .debug_objs/2channeldmrgtenth.o -o 2channeldmrgtenth-g $(LIBGFLAGS)

2channeldmrgeleventh: 2channeldmrgeleventh.o $(ITENSOR_LIBS) $(TENSOR_HEADERS)
	$(CCCOM) $(CCFLAGS) 2channeldmrgeleventh.o -o 2channeldmrgeleventh $(LIBFLAGS)

2channeldmrgeleventh-g: mkdebugdir .debug_objs/2channeldmrgeleventh.o $(ITENSOR_GLIBS) $(TENSOR_HEADERS)
	$(CCCOM) $(CCGLAGS) .debug_objs/2channeldmrgeleventh.o -o 2channeldmrgeleventh-g $(LIBGFLAGS)

2channeldmrgtwelveth: 2channeldmrgtwelveth.o $(ITENSOR_LIBS) $(TENSOR_HEADERS)
	$(CCCOM) $(CCFLAGS) 2channeldmrgtwelveth.o -o 2channeldmrgtwelveth $(LIBFLAGS)

2channeldmrgtwelveth-g: mkdebugdir .debug_objs/2channeldmrgtwelveth.o $(ITENSOR_GLIBS) $(TENSOR_HEADERS)
	$(CCCOM) $(CCGLAGS) .debug_objs/2channeldmrgtwelveth.o -o 2channeldmrgtwelveth-g $(LIBGFLAGS)

2channeldmrgthirteenth: 2channeldmrgthirteenth.o $(ITENSOR_LIBS) $(TENSOR_HEADERS)
	$(CCCOM) $(CCFLAGS) 2channeldmrgthirteenth.o -o 2channeldmrgthirteenth $(LIBFLAGS)

2channeldmrgthirteenth-g: mkdebugdir .debug_objs/2channeldmrgthirteenth.o $(ITENSOR_GLIBS) $(TENSOR_HEADERS)
	$(CCCOM) $(CCGLAGS) .debug_objs/2channeldmrgthirteenth.o -o 2channeldmrgthirteenth-g $(LIBGFLAGS)

2channeldmrgfourteenth: 2channeldmrgfourteenth.o $(ITENSOR_LIBS) $(TENSOR_HEADERS)
	$(CCCOM) $(CCFLAGS) 2channeldmrgfourteenth.o -o 2channeldmrgfourteenth $(LIBFLAGS)

2channeldmrgfourteenth-g: mkdebugdir .debug_objs/2channeldmrgfourteenth.o $(ITENSOR_GLIBS) $(TENSOR_HEADERS)
	$(CCCOM) $(CCGLAGS) .debug_objs/2channeldmrgfourteenth.o -o 2channeldmrgfourteenth-g $(LIBGFLAGS)

2channeldmrgfifteenth: 2channeldmrgfifteenth.o $(ITENSOR_LIBS) $(TENSOR_HEADERS)
	$(CCCOM) $(CCFLAGS) 2channeldmrgfifteenth.o -o 2channeldmrgfifteenth $(LIBFLAGS)

2channeldmrgfifteenth-g: mkdebugdir .debug_objs/2channeldmrgfifteenth.o $(ITENSOR_GLIBS) $(TENSOR_HEADERS)
	$(CCCOM) $(CCGLAGS) .debug_objs/2channeldmrgfifteenth.o -o 2channeldmrgfifteenth-g $(LIBGFLAGS)

2channeldmrgsixteenth: 2channeldmrgsixteenth.o $(ITENSOR_LIBS) $(TENSOR_HEADERS)
	$(CCCOM) $(CCFLAGS) 2channeldmrgsixteenth.o -o 2channeldmrgsixteenth $(LIBFLAGS)

2channeldmrgsixteenth-g: mkdebugdir .debug_objs/2channeldmrgsixteenth.o $(ITENSOR_GLIBS) $(TENSOR_HEADERS)
	$(CCCOM) $(CCGLAGS) .debug_objs/2channeldmrgsixteenth.o -o 2channeldmrgsixteenth-g $(LIBGFLAGS)

dmrggrover: dmrggrover.o $(ITENSOR_LIBS) $(TENSOR_HEADERS)
	$(CCCOM) $(CCFLAGS) dmrggrover.o -o dmrggrover $(LIBFLAGS)

dmrggrover-g: mkdebugdir .debug_objs/dmrggrover.o $(ITENSOR_GLIBS) $(TENSOR_HEADERS)
	$(CCCOM) $(CCGLAGS) .debug_objs/dmrggrover.o -o dmrggrover-g $(LIBGFLAGS)

try2: try2.o $(ITENSOR_LIBS) $(TENSOR_HEADERS)
	$(CCCOM) $(CCFLAGS) try2.o -o try2 $(LIBFLAGS)

try2-g: mkdebugdir .debug_objs/try2.o $(ITENSOR_GLIBS) $(TENSOR_HEADERS)
	$(CCCOM) $(CCGLAGS) .debug_objs/try2.o -o try2-g $(LIBGFLAGS)

2channeldmrgsecond: 2channeldmrgsecond.o $(ITENSOR_LIBS) $(TENSOR_HEADERS)
	$(CCCOM) $(CCFLAGS) 2channeldmrgsecond.o -o 2channeldmrgsecond $(LIBFLAGS)

dmrg-g: mkdebugdir .debug_objs/dmrg.o $(ITENSOR_GLIBS) $(TENSOR_HEADERS)
	$(CCCOM) $(CCGFLAGS) .debug_objs/dmrg.o -o dmrg-g $(LIBGFLAGS)

mixedspin: mixedspin.o $(ITENSOR_LIBS) $(TENSOR_HEADERS)
	$(CCCOM) $(CCFLAGS) mixedspin.o -o mixedspin $(LIBFLAGS)

mixedspin-g: mkdebugdir .debug_objs/mixedspin.o $(ITENSOR_GLIBS) $(TENSOR_HEADERS)
	$(CCCOM) $(CCGFLAGS) .debug_objs/mixedspin.o -o mixedspin-g $(LIBGFLAGS)

dmrg_table: dmrg_table.o $(ITENSOR_LIBS) $(TENSOR_HEADERS)
	$(CCCOM) $(CCFLAGS) dmrg_table.o -o dmrg_table $(LIBFLAGS)

dmrg_table-g: mkdebugdir .debug_objs/dmrg_table.o $(ITENSOR_GLIBS) $(TENSOR_HEADERS)
	$(CCCOM) $(CCGFLAGS) .debug_objs/dmrg_table.o -o dmrg_table-g $(LIBGFLAGS)

dmrgj1j2: dmrgj1j2.o $(ITENSOR_LIBS) $(TENSOR_HEADERS)
	$(CCCOM) $(CCFLAGS) dmrgj1j2.o -o dmrgj1j2 $(LIBFLAGS)

dmrgj1j2-g: mkdebugdir .debug_objs/dmrgj1j2.o $(ITENSOR_LIBS) $(TENSOR_HEADERS)
	$(CCCOM) $(CCFLAGS) .debug_objs/dmrgj1j2.o -o dmrgj1j2-g $(LIBGFLAGS)

exthubbard: exthubbard.o $(ITENSOR_LIBS) $(TENSOR_HEADERS)
	$(CCCOM) $(CCFLAGS) exthubbard.o -o exthubbard $(LIBFLAGS)

exthubbard-g: mkdebugdir .debug_objs/exthubbard.o $(ITENSOR_LIBS) $(TENSOR_HEADERS)
	$(CCCOM) $(CCFLAGS) .debug_objs/exthubbard.o -o exthubbard-g $(LIBGFLAGS)

trg: trg.o $(ITENSOR_LIBS) $(TENSOR_HEADERS)
	$(CCCOM) $(CCFLAGS) trg.o -o trg $(LIBFLAGS)

trg-g: mkdebugdir .debug_objs/trg.o $(ITENSOR_GLIBS) $(TENSOR_HEADERS)
	$(CCCOM) $(CCGFLAGS) .debug_objs/trg.o -o trg-g $(LIBGFLAGS)

ctmrg: ctmrg.o $(ITENSOR_LIBS) $(TENSOR_HEADERS)
	$(CCCOM) $(CCFLAGS) ctmrg.o -o ctmrg $(LIBFLAGS)

ctmrg-g: mkdebugdir .debug_objs/ctmrg.o $(ITENSOR_GLIBS) $(TENSOR_HEADERS)
	$(CCCOM) $(CCGFLAGS) .debug_objs/ctmrg.o -o ctmrg-g $(LIBGFLAGS)

hubbard_2d: hubbard_2d.o $(ITENSOR_LIBS) $(TENSOR_HEADERS)
	$(CCCOM) $(CCFLAGS) hubbard_2d.o -o hubbard_2d $(LIBFLAGS)

hubbard_2d-g: mkdebugdir .debug_objs/hubbard_2d.o $(ITENSOR_GLIBS) $(TENSOR_HEADERS)
	$(CCCOM) $(CCGFLAGS) .debug_objs/hubbard_2d.o -o hubbard_2d-g $(LIBGFLAGS)

hubbard_2d_conserve_momentum: hubbard_2d_conserve_momentum.o $(ITENSOR_LIBS) $(TENSOR_HEADERS)
	$(CCCOM) $(CCFLAGS) hubbard_2d_conserve_momentum.o -o hubbard_2d_conserve_momentum $(LIBFLAGS)

hubbard_2d_conserve_momentum-g: mkdebugdir .debug_objs/hubbard_2d_conserve_momentum.o $(ITENSOR_GLIBS) $(TENSOR_HEADERS)
	$(CCCOM) $(CCGFLAGS) .debug_objs/hubbard_2d_conserve_momentum.o -o hubbard_2d_conserve_momentum-g $(LIBGFLAGS)

flatband1: flatband1.o $(ITENSOR_LIBS) $(TENSOR_HEADERS)
	$(CCCOM) $(CCFLAGS) flatband1.o -o flatband1 $(LIBFLAGS)

mkdebugdir:
	mkdir -p .debug_objs

clean:
	@rm -fr *.o .debug_objs dmrg dmrg-g helloitensor 2channeldmrg 2channeldmrgthird 1channeldmrg 2channeldmrgfourth 2channeldmrgfifth 2channeldmrgsixth 2channeldmrgseventh try2 2channeldmrgninth 2channeldmrgtenth 2channeldmrgeleventh 2channeldmrgtwelveth 2channeldmrgthirteenth 2channeldmrgfourteenth 2channeldmrgfifteenth 2channeldmrgsixteenth flatband1 kitaev kitaev2 kitaev3 lieblattice lieblattice2 try1 liebhubbard liebhubbard2 liebhubbard3 liebhubbard4 liebhubbard5 liebhubbard6 liebhubbard7 liebhubbard9 liebkitaevread liebkitaevread2 dmrggrover liebhubbard10 liebhubbard11 kitaevhubbard\
	dmrg_table dmrg_table-g dmrgj1j2 dmrgj1j2-g exthubbard exthubbard-g 2channeldmrgsecond\
    mixedspin mixedspin-g trg trg-g ctmrg ctmrg-g hubbard_2d hubbard_2d-g hubbard_2d_conserve_momentum hubbard_2d_conserve_momentum-g
