#include <stdlib.h>

// Index matrix (by columns)
int mi(int i, int j, int d1, int d2)
{
	return j*d1 + i;
}

// Index 3-tensor (by columns, matrices ordered by last dim)
int ti(int i, int j, int k, int d1, int d2, int d3)
{
	return k*d1*d2 + j*d1 + i;
}

// Empirical cross-moment of order 2 between X size nxd and Y size n
void Moments_M2(double* X, double* Y, int* pn, int* pd, double* M2)
{
	int n=*pn, d=*pd;
	//double* M2 = (double*)calloc(d*d,sizeof(double));

	// M2 = E[Y*X^*2] - E[Y*e^*2] = E[Y (X^*2 - I)]
	for (int j=0; j<d; j++)
	{
		for (int i=0; i<n; i++)
		{
			M2[mi(j,j,d,d)] -= Y[i] / n;
			for (int k=0; k<d; k++)
				M2[mi(j,k,d,d)] += Y[i] * X[mi(i,j,n,d)]*X[mi(i,k,n,d)] / n;
		}
	}
}

// Empirical cross-moment of order 3 between X size nxd and Y size n
void Moments_M3(double* X, double* Y, int* pn, int* pd, double* M3)
{
	int n=*pn, d=*pd;
	//double* M3 = (double*)calloc(d*d*d,sizeof(double));

	// M3 = E[Y*X^*3] - E[Y*e*X*e] - E[Y*e*e*X] - E[Y*X*e*e]
	for (int j=0; j<d; j++)
	{
		for (int k=0; k<d; k++)
		{
			for (int i=0; i<n; i++)
			{
				double tensor_elt = Y[i]*X[mi(i,k,n,d)] / n;
				M3[ti(j,k,j,d,d,d)] -= tensor_elt;
				M3[ti(j,j,k,d,d,d)] -= tensor_elt;
				M3[ti(k,j,j,d,d,d)] -= tensor_elt;
				for (int o=0; o<d; o++)
					M3[ti(j,k,o,d,d,d)] += Y[i] * X[mi(i,j,n,d)]*X[mi(i,k,n,d)]*X[mi(i,o,n,d)] / n;
			}
		}
	}
}

// W = 1/N sum( t(g(Zi,theta)) g(Zi,theta) )
// with g(Zi, theta) = i-th contribution to all moments (size dim) - real moments
void Compute_Omega(double* X, double* Y, double* M, int* pn, int* pd, double* W)
{
	int n=*pn, d=*pd;
  //int dim = d+d*d+d*d*d
  //double* W = (double*)calloc(dim*dim,sizeof(double));
  double* g = (double*)malloc(dim * sizeof(double));
  for (int i=0; i<n; i++)
  {
    // Fill gi:
    for (int j=0; j<d; j++)
      g[j] = Y[i] * X[mi(i,j,n,d)] - M[i]
    for (int j=d; j<d+(d*d); j++)
    {
      int idx1 = (j-d) % d; //num row
      int idx2 = ((j-d) - idx1) / d; //num col
      g[j] = 0.0;
      if (idx1 == idx2)
        g[j] -= Y[i];
			g[j] += Y[i] * X[mi(i,idx1,n,d)]*X[mi(i,idx2,n,d)];
    }
    for (int j=d+d*d; j<dim; j++)
    {
      int idx1 = (j-d-d*d) % d; //num row
      int idx2 = ((j-d-d*d - idx1) / d) %d; //num col
      int idx3 = (((j-d-d*d - idx1) / d) - idx2) / d; //num "depth"
      g[j] = 0.0;
      double tensor_elt = Y[i]*X[mi(i,k,n,d)] / n;
      if (idx1 == idx2)
        g[j] -= Y[i] * X[mi(i,idx3,n,d)];
      if (idx1 == idx3)
        g[j] -= Y[i] * X[mi(i,idx2,n,d)];
      if (idx2 == idx3)
        g[j] -= Y[i] * X[mi(i,idx1,n,d)];
      g[j] += Y[i] * X[mi(i,idx1,n,d)]*X[mi(i,idx2,n,d)]*X[mi(i,idx3,n,d)];
    }
    // Add 1/n t(gi) %*% gi to W
    for (int j=0; j<dim; j++)
    {
      for (int k=0; k<dim; k++)
        W[j*dim+k] += g[j] * g[k] / n;
    }
  }
  free(g);
}
