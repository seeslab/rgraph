DEBUG_FLAG = -g
TEST_WARNING_FLAG = -Wall
WARNING_FLAG = -Wall
PRNGDIR = ${HOME}/Software/prng-3.0.2

test : main_test.c librgraph.a
	gcc -o test main_test.c -I"${PRNGDIR}/include/" \
	-L"${PRNGDIR}/lib/" -L./ -lprng -lrgraph -lm \
	${DEBUG_FLAG} ${TEST_WARNING_FLAG}

librgraph.a :  tools.o graph.o models.o modules.o 
	ar rc librgraph.a tools.o graph.o models.o modules.o
	ranlib librgraph.a

tools.o : tools.c tools.h
	gcc -c tools.c ${DEBUG_FLAG} ${WARNING_FLAG}

graph.o : graph.c graph.h tools.h
	gcc -c graph.c -I"${PRNGDIR}/include" ${DEBUG_FLAG} ${WARNING_FLAG}

models.o : models.c models.h tools.h graph.h
	gcc -c models.c -I"${PRNGDIR}/include" ${DEBUG_FLAG} ${WARNING_FLAG}

modules.o : modules.c modules.h tools.h graph.h
	gcc -c modules.c -I"${PRNGDIR}/include" ${DEBUG_FLAG} ${WARNING_FLAG}

clean:
	rm -f *.o
