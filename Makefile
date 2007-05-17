

PRNGDIR=/home/rguimera/Software/prng-3.0.2

test : main_test.c libgraph.a
	gcc -o test main_test.c -I${PRNGDIR}/include/ \
	-L${PRNGDIR}/lib/ -L./ -lprng -lgraph -lm

libgraph.a :  graph.o tools.o
	ar rc libgraph.a graph.o tools.o
	ranlib libgraph.a

tools.o : tools.c tools.h
	gcc -c tools.c

graph.o : graph.c graph.h tools.h
	gcc -c graph.c -I /home/rguimera/Software/prng-3.0.2/include/
