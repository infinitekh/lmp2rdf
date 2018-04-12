lmprdf_omp.out:
	g++ main.cpp makeRDF_omp.cpp kh_math_fourier.cpp  snapshot.cpp -O2 -fopenmp -o lmprdf_omp.out

