QDPCONFIG = /work/01112/tg803910/stampede2/install/qdpxx/bin/qdp++-config

CC = mpiicc
CFLAGS = -g -qopenmp $(shell $(QDPCONFIG) --cxxflags) -I/work/01112/tg803910/stampede2/install/qphix/include/
#CFLAGS = -g -O2 -I/work/01112/tg803910/stampede2/install/qdpxx/include -I/usr/include/libxml2 -I/work/01112/tg803910/stampede2/install/qmp/include -qopenmp -I/work/01112/tg803910/stampede2/install/qphix/include/ -Wall
LIBS = -lqphix_solver -lqphix_codegen $(shell $(QDPCONFIG) --libs)
LDFLAGS = -L/work/01112/tg803910/stampede2/install/qphix/lib $(shell $(QDPCONFIG) --ldflags)

clover-test: mesfield.o main.o dslashm_w.o setup.o sources.o invert.o
	$(CC) $(CFLAGS) $^ -o $@ $(LIBS) $(LDFLAGS)

%.o: %.cc
	$(CC) $(CFLAGS) $< -c -o $@

clean:
	rm -f *.o clover-test
