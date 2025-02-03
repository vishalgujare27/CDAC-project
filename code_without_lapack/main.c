#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <time.h>

#define PI 3.1415926535897932

typedef struct {
    double real;
    double imag;
} Complex;

double potential(double x);

double potential(double x) {
    double pot;
    pot = 0.50 * x * x;  // Example potential function (harmonic oscillator)
    return pot;
}

// Function to multiply two complex numbers
Complex cmult(Complex a, Complex b) {
    Complex result;
    result.real = a.real * b.real - a.imag * b.imag;
    result.imag = a.real * b.imag + a.imag * b.real;
    return result;
}

// Add two complex numbers
Complex cadd(Complex a, Complex b) {
    Complex result;
    result.real = a.real + b.real;
    result.imag = a.imag + b.imag;
    return result;
}

// Subtract two complex numbers
Complex csub(Complex a, Complex b) {
    Complex result;
    result.real = a.real - b.real;
    result.imag = a.imag - b.imag;
    return result;
}

// Multiply a complex number with a scalar
Complex cscal(Complex a, double scalar) {
    Complex result;
    result.real = a.real * scalar;
    result.imag = a.imag * scalar;
    return result;
}

// Conjugate of a complex number
Complex cconj(Complex a) {
    Complex result;
    result.real = a.real;
    result.imag = -a.imag;
    return result;
}

// Function to compute the matrix-vector multiplication A * v
void matvec_mult(Complex **A, Complex *v, Complex *result, int n) {
    for (int i = 0; i < n; i++) {
        result[i] = (Complex){0.0, 0.0};
        for (int j = 0; j < n; j++) {
            result[i] = cadd(result[i], cmult(A[i][j], v[j]));
        }
    }
}

// Simple matrix multiplication: C = A * B
void mat_mult(Complex **A, Complex **B, Complex **C, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            C[i][j] = (Complex){0.0, 0.0};
            for (int k = 0; k < n; k++) {
                C[i][j] = cadd(C[i][j], cmult(A[i][k], B[k][j]));
            }
        }
    }
}

// Simple function to compute eigenvalues (Placeholder, not fully implemented)
void compute_eigenvalues(double **A, double *eig, int n) {
    // This function would normally use some iterative method like the QR algorithm
    // Here we will simply assume the diagonal of A are the eigenvalues (for simplicity).
    for (int i = 0; i < n; i++) {
        eig[i] = A[i][i];  // Placeholder, real implementation should calculate eigenvalues
    }
}

// Time propagation function
int timeprop(int ncsf, double xmin, double xmax, double lmda, int flqchnl, int noptc, 
             int istate, double totime, int ntim, double omga, double epsln, FILE *outfile) {
    double mass = 10000.0; 
    double dtim = totime / ntim;
    double dx = (xmax - xmin) / (ncsf - 1);
    int tchnl = 2 * flqchnl + 1;
    int nfloq = ncsf * tchnl;
    int inzr = flqchnl * ncsf;
    double tau = 2.0 * PI / omga;

    // Allocate and initialize matrices and vectors
    double **hmt = (double **)calloc(ncsf, sizeof(double *));
    for (int i = 0; i < ncsf; i++) {
        hmt[i] = (double *)calloc(ncsf, sizeof(double));
    }

    // Kinetic energy calculation
    for (int k = 0; k < ncsf; k++) {
        hmt[k][k] = (PI * PI) / (6.0 * dx * dx);
        for (int j = 0; j < k; j++) {
            hmt[k][j] = (pow(-1.0, (k - j))) / ((k - j) * (k - j) * dx * dx);
            hmt[j][k] = hmt[k][j];
        }
    }
    
    FILE *file = fopen("kineticseq.txt","w");  
    if (file == NULL) {
        printf("Error opening kineticseq.txt\n");
        return 1;
    }

    for (int i = 0; i < ncsf; i++) {
        for (int j = 0; j < ncsf; j++) {
            fprintf(file, "%e ", hmt[i][j]);
        }
        fprintf(file, "\n");
    }
    fclose(file);

    // Potential energy modification
    double *zdipole = (double *)calloc(ncsf, sizeof(double));
    for (int i = 0; i < ncsf; i++) {
        double x = xmin + i * dx;
        zdipole[i] = x;
        hmt[i][i] += potential(x);
    }








    // Eigenvalue calculation (using placeholder)
    double *eig = (double *)calloc(ncsf, sizeof(double));
    compute_eigenvalues(hmt, eig, ncsf);


 file = fopen("hmtafter.txt","w");  
  for(int i=0; i<ncsf; i++){
    for(int j=0; j<ncsf; j++){
fprintf(file,"%e ",hmt[i][j]);
    }     
    fprintf(file,"\n");
  }





    // Initialize the wavefunction (for simplicity, we'll use a simple 1D wavefunction)
    Complex *cofi = (Complex *)calloc(nfloq, sizeof(Complex));
    for (int i = 0; i < ncsf; i++) {
        cofi[inzr + i] = (Complex){eig[i], 0.0};  // Assuming the state is real initially
    }

    // Time propagation (simplified, no OpenMP or parallelization)
    for (int itim = 0; itim < ntim; itim++) {
        double time = itim * dtim;
        double epls = epsln * sin((PI * time) / totime) * sin((PI * time) / totime);

        // Time evolution calculation (simplified)
        for (int i = 0; i < nfloq; i++) {
            cofi[i] = cscal(cofi[i], epls);
        }
    }

    // Calculate the final kinetic energy
    double kinetic_energy = 0.0;
    for (int i = 0; i < ncsf; i++) {
        kinetic_energy += 0.5 * hmt[i][i] * (cofi[inzr + i].real * cofi[inzr + i].real + cofi[inzr + i].imag * cofi[inzr + i].imag);
    }

    // Print the final wavefunction and kinetic energy to the file
    fprintf(outfile, "\nFinal Wavefunction (first 10 values):\n");
    for (int i = 0; i < 10 && i < nfloq; i++) {
        fprintf(outfile, "cofi[%d] = (%lf, %lf)\n", i, cofi[i].real, cofi[i].imag);
    }

    fprintf(outfile, "\nTotal Kinetic Energy: %lf\n", kinetic_energy);

    // Free allocated memory
    for (int i = 0; i < ncsf; i++) {
        free(hmt[i]);
    }
    free(hmt);
    free(cofi);
    free(eig);
    free(zdipole);

    return 0;
}

int main() {
    int ntime, ncstep, istate, flqchnl;
    int ncsf, noptc;
    double xmin, xmax, epsln, omga;
    double tau, totime, lmda;

    ncsf = 100; flqchnl = 5; xmin = 0; xmax = 10.0;
    epsln = 0.5; omga = 0.5; lmda = 0.0;
    noptc = 5; istate = 1;

    tau = 2.0 * PI / omga;
    totime = noptc * tau;
    ntime = (int)totime;

    printf("ntime = %d, totime = %lf, tau = %lf\n", ntime, totime, tau);
    clock_t start = clock();

    // Open the output file
    FILE *outfile = fopen("results.txt", "w");
    if (outfile == NULL) {
        printf("Error opening file for writing.\n");
        return 1;
    }

    timeprop(ncsf, xmin, xmax, lmda, flqchnl, noptc, istate, totime, ntime, omga, epsln, outfile);

    // Close the file
    fclose(outfile);

    clock_t stop = clock();
    printf("Time taken: %lf seconds\n", (double)(stop - start) / CLOCKS_PER_SEC);
    return 0;
}
