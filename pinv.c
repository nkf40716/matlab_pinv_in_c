#include <stdlib.h>

extern int dsvd(float **a, int m, int n, float *w, float **v);

// *   m = row dimension of a
// *   n = column dimension of a
// Output: X
void pinv(float *a, int m, int n, float *X)
{
	int i, j, k, pos = 0;
	float t;
	float *w, *v_mem;
	float **v, **ppa, **pX;
	float *mem, **ppmem;

	mem = (float *)malloc(sizeof(float) * (m + m * n));
	w     = mem;
	v_mem = mem + m;

	ppmem = (float **)malloc(sizeof(float *) * m * 3);
	v   = ppmem + pos, pos += m;
	ppa = ppmem + pos, pos += m;
	pX  = ppmem + pos, pos += m;

	for (i = 0; i < m; ++i) {
		v[i]   = &v_mem[i * n];
		ppa[i] = &a[i * n];
		pX[i]  = &X[i * n];
	}

	// Singular Value Decomposition
	dsvd(ppa, m, n, w, v);

	// swap for match the format in matlab
	for (i = 0, j = m-1; i < (m >> 1); ++i, --j) {
		t = w[i];
		w[i] = 1 / w[j];
		w[j] = 1 / t;		// s = 1./s
	}
	for (k = 0; k < n; ++k) {
		for (i = 0, j = m-1; i < (m >> 1); ++i, --j) {
			t = ppa[k][i];
			ppa[k][i] = ppa[k][j];
			ppa[k][j] = t;
			t = v[k][i];
			v[k][i] = v[k][j];
			v[k][j] = t;
		}
	}

// 	for (i = 0; i < m; ++i) {
// 		float wt = (w[i] < 0) ? (-w[i]) : (w[i]);
// 		if (wt > tol) tol = wt;		// find eps(norm(s,inf))
	 
	for (i = 0; i < m; ++i) {
		for (j = 0; j < n; ++j)
			v[j][i] *= w[i];	// bsxfun(@times,V,s.')
	}

	for (i = 0; i < m; ++i) {
		for (j = 0; j < n; ++j) {
			float sum = 0;
			for (k = 0; k < n; ++k)
				sum += v[i][k] * ppa[j][k];
			pX[i][j] = sum;
		}
	}

	free(mem);
	free(ppmem);
}
