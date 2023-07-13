MAJOR=1
MINOR=1
PATCH=0
HOTFIX=.1



build: main.cpp makeRDF_omp.cpp makeRDF.h kh_math_fourier.h kh_math_fourier.cpp snapshot.cpp makeRDF.h snapshot.h  common.c common.h  hydro_math.cpp hydro_math.h
	g++ -std=c++17 main.cpp makeRDF_omp.cpp kh_math_fourier.cpp  snapshot.cpp common.c hydro_math.cpp -O2 -fopenmp -o lmprdf_$(MAJOR).$(MINOR).$(PATCH)$(HOTFIX).out

x3d: main.cpp makeRDF_omp.cpp makeRDF.h kh_math_fourier.h kh_math_fourier.cpp snapshot.cpp makeRDF.h snapshot.h 
	g++ main.cpp makeRDF_omp.cpp kh_math_fourier.cpp  snapshot.cpp -O2 -fopenmp -o lmprdf_$(MAJOR).$(MINOR).$(PATCH).out
