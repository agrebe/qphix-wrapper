QDP = /work/01112/tg803910/stampede2/install/qdpxx-scalar
QPHIX = /work/01112/tg803910/stampede2/install/qphix-scalar
QDPCONFIG = $(QDP)/bin/qdp++-config


CC = icpc
CFLAGS = -g -qopenmp $(shell $(QDPCONFIG) --cxxflags) -I$(QPHIX)/include/
#CFLAGS = -g -O2 -I$(QDP)/include -I/usr/include/libxml2 -I/work/01112/tg803910/stampede2/install/qmp/include -qopenmp -I$(QPHIX)/include/ -Wall
LIBS = -lqphix_solver -lqphix_codegen $(shell $(QDPCONFIG) --libs)
LDFLAGS = -L$(QPHIX)/lib $(shell $(QDPCONFIG) --ldflags) -L.

clover-test: main.o libqphix-wrapper.a
	$(CC) $(CFLAGS) $^ -o $@ -lqphix-wrapper $(LIBS) $(LDFLAGS)

%.o: %.cc
	$(CC) $(CFLAGS) $< -c -o $@

libqphix-wrapper.a: mesfield.o dslashm_w.o setup.o sources.o invert.o
	ar rcs $@ $^

clean:
	rm -f *.o *.a clover-test
