gcc -o exe openmp.c -lm -llapacke -llapack -lblas -fopenmp
./exe #>&out&
