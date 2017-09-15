makeRDF.o: makeRDF.cpp makeRDF.h snapshot.h
main.o: main.cpp snapshot.h makeRDF.h
snapshot.o: snapshot.cpp snapshot.h

clean:
	rm *.o
