MAJOR=1
MINOR=0
PATCH=4


build: main.cpp makeRDF_omp.cpp makeRDF.h kh_math_fourier.h kh_math_fourier.cpp snapshot.cpp makeRDF.h snapshot.h 
	g++ main.cpp makeRDF_omp.cpp kh_math_fourier.cpp  snapshot.cpp -O2 -fopenmp -o lmprdf_$(MAJOR).$(MINOR).$(PATCH).out

x3d: main.cpp makeRDF_omp.cpp makeRDF.h kh_math_fourier.h kh_math_fourier.cpp snapshot.cpp makeRDF.h snapshot.h 
	g++ main.cpp makeRDF_omp.cpp kh_math_fourier.cpp  snapshot.cpp -O2 -fopenmp -o lmprdf_$(MAJOR).$(MINOR).$(PATCH).out
