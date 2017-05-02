#include <stdlib.h>

#ifndef max
#define max(x,y) ((x)>(y)?(x):(y))
#endif

extern int dsvd(float **a, int m, int n, float *w, float **v);

// eps(x) for single input
static double eps(double x)
{
	double r;
	int exponent;

	if (x < 0.0) 
		x = -x;
	if (x != x)	// test if x is Not-a-Number
		return x;
	
	if (x <= 1.17549435E-38F) 
		r = 1.4013E-45F;
	else {
		r = frexp(x, &exponent);
		r = ldexp(1.0, exponent - 24);
    }
	return r;
}

// *   m = row dimension of a
// *   n = column dimension of a
// Output: X
void pinv(double *a, int m, int n, double *X)
{
	int i, j, k, r1, pos = 0;
	double t, norm;
	double *w, *v_mem;
	double **v, **ppa, **pX;
	double *mem, **ppmem;

	mem = (double *)malloc(sizeof(double) * (m + m * n));
	w     = mem;
	v_mem = mem + m;

	ppmem = (double **)malloc(sizeof(double *) * m * 3);
	v   = ppmem + pos, pos += m;
	ppa = ppmem + pos, pos += m;
	pX  = ppmem + pos, pos += m;

	for (i = 0; i < m; ++i) {
		v[i]   = &v_mem[i * n];
		ppa[i] = &a[i * n];
		pX[i]  = &X[i * n];
	}

	dsvd(ppa, m, n, w, v);

	// norm(w, inf)
	norm = (w[0] < 0) ? (-w[0]) : (w[0]);
	for (i = 1; i < m; ++i) {
		t = (w[i] < 0) ? (-w[i]) : (w[i]); 
		if (t > norm) norm = t;
	}
	
	// t = tolerance, any singular values less than a tolerance are treated as zero
	t = max(m, n) * eps(norm);
	for (i = 0, r1 = 0; i < m; ++i) {
		if (t > w[i])
			++r1;
	}	
	 
	for (i = 0; i < (m - r1); ++i) {
		for (j = 0; j < n; ++j)
			v[j][i] *= w[i];	// bsxfun(@times,V,s.')
	}

	for (i = 0; i < m; ++i) {
		for (j = 0; j < n; ++j) {
			double sum = 0;
			for (k = 0; k < (m - r1); ++k)
				sum += v[i][k] * ppa[j][k];
			pX[i][j] = sum;
		}
	}

	free(mem);
	free(ppmem);
}
