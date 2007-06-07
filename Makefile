DEBUG_FLAG = -g
TEST_WARNING_FLAG = -Wall
WARNING_FLAG = -Wall
PRNGDIR = ${HOME}/Software/prng-3.0.2

test : main_test.c libgraph.a
	gcc -o test main_test.c -I${PRNGDIR}/include/ \
	-L${PRNGDIR}/lib/ -L./ -lprng -lgraph -lm \
	${DEBUG_FLAG} ${TEST_WARNING_FLAG}

libgraph.a :  graph.o tools.o modules.o
	ar rc libgraph.a graph.o tools.o modules.o
	ranlib libgraph.a

tools.o : tools.c tools.h
	gcc -c tools.c ${DEBUG_FLAG} ${WARNING_FLAG}

graph.o : graph.c graph.h tools.h
	gcc -c graph.c -I${PRNGDIR}/include ${DEBUG_FLAG} ${WARNING_FLAG}

modules.o : modules.c modules.h tools.h
	gcc -c modules.c -I${PRNGDIR}/include ${DEBUG_FLAG} ${WARNING_FLAG}

clean:
	rm -f *.o
