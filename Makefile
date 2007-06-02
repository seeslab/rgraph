PRNGDIR=${HOME}/Software/prng-3.0.2

test : main_test.c libgraph.a
	gcc -Wall -o test main_test.c -I${PRNGDIR}/include/ \
	-L${PRNGDIR}/lib/ -L./ -lprng -lgraph -lm -g

libgraph.a :  graph.o tools.o modules.o
	ar rc libgraph.a graph.o tools.o modules.o
	ranlib libgraph.a

tools.o : tools.c tools.h
	gcc -c tools.c

graph.o : graph.c graph.h tools.h
	gcc -c graph.c -I${PRNGDIR}/include

modules.o : modules.c modules.h tools.h
	gcc -c modules.c -I${PRNGDIR}/include -Wall

clean:
	rm -f *.o
