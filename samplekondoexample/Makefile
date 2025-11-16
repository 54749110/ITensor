LIBRARY_DIR=/home/ITensor/
include $(LIBRARY_DIR)/this_dir.mk
include $(LIBRARY_DIR)/options.mk

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

build: dmrg_kondo dmrg_kondo2channel dmrg_kondo2channelv2 dmrg_kondo3channel dmrg_kondo4channel

debug: dmrg_kondo-g kondo2channel-g kondo2channelv2-g kondo3channel-g kondo4channel-g

all: dmrg_kondo kondo2channel kondo2channelv2 kondo3channel kondo4channel

dmrg_kondo: dmrg_kondo.o $(ITENSOR_LIBS) $(TENSOR_HEADERS)
	$(CCCOM) $(CCFLAGS) dmrg_kondo.o -o dmrg_kondo $(LIBFLAGS)

dmrg_kondo-g: mkdebugdir .debug_objs/dmrg_kondo.o $(ITENSOR_GLIBS) $(TENSOR_HEADERS)
	$(CCCOM) $(CCGFLAGS) .debug_objs/dmrg_kondo.o -o dmrg_kondo-g $(LIBGFLAGS)

dmrg_kondo2channel: dmrg_kondo2channel.o $(ITENSOR_LIBS) $(TENSOR_HEADERS)
	$(CCCOM) $(CCFLAGS) dmrg_kondo2channel.o -o dmrg_kondo2channel $(LIBFLAGS)

dmrg_kondo2channel-g: mkdebugdir .debug_objs/dmrg_kondo2channel.o $(ITENSOR_GLIBS) $(TENSOR_HEADERS)

dmrg_kondo2channelv2: dmrg_kondo2channelv2.o $(ITENSOR_LIBS) $(TENSOR_HEADERS)
	$(CCCOM) $(CCFLAGS) dmrg_kondo2channelv2.o -o dmrg_kondo2channelv2 $(LIBFLAGS)

dmrg_kondo2channelv2-g: mkdebugdir .debug_objs/dmrg_kondo2channelv2.o $(ITENSOR_GLIBS) $(TENSOR_HEADERS)

dmrg_kondo3channel: dmrg_kondo3channel.o $(ITENSOR_LIBS) $(TENSOR_HEADERS)
	$(CCCOM) $(CCFLAGS) dmrg_kondo3channel.o -o dmrg_kondo3channel $(LIBFLAGS)

dmrg_kondo3channel-g: mkdebugdir .debug_objs/dmrg_kondo3channel.o $(ITENSOR_GLIBS) $(TENSOR_HEADERS)

dmrg_kondo4channel: dmrg_kondo4channel.o $(ITENSOR_LIBS) $(TENSOR_HEADERS)
	$(CCCOM) $(CCFLAGS) dmrg_kondo4channel.o -o dmrg_kondo4channel $(LIBFLAGS)

dmrg_kondo4channel-g: mkdebugdir .debug_objs/dmrg_kondo4channel.o $(ITENSOR_GLIBS) $(TENSOR_HEADERS)

mkdebugdir:
	mkdir -p .debug_objs

clean:
	@rm -fr *.o .debug_objs dmrg_kondo dmrg_kondo-g \
	dmrg_kondo2channel dmrg_kondo2channel-g \
	dmrg_kondo2channelv2 dmrg_kondo2channelv2-g \
	dmrg_kondo3channel dmrg_kondo3channel-g \
	dmrg_kondo4channel dmrg_kondo4channel-g
