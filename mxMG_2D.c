#include <math.h>
#include <matrix.h>
#include <mex.h>
#include "mg_2d.c"

/*
mg_2d_MATLAB_gateway(f, iter, r)
f     - f-vector
iter  - FMGV iterations
(k     - 2^k - 1 finest grid)
r     - 2^r - 1 roughest grid
*/

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

  double *f, *iter,  *r, *v, *k, gamma;
  size_t mrows, ncols;

  if(nrhs!=5) {
    mexErrMsgIdAndTxt( "MATLAB:mg_1d_MATLAB_gateway:invalidNumInputs",
      "Five inputs required.");
  }

  mrows = mxGetM(prhs[0]);
  ncols = mxGetN(prhs[0]);

  //double k = log2((double) sqrt(mrows)+1);

  plhs[0] = mxCreateDoubleMatrix((mwSize)mrows, (mwSize)ncols, mxREAL);
  
  f = mxGetPr(prhs[0]);
  iter = mxGetPr(prhs[1]);
  k = mxGetPr(prhs[2]);
  r = mxGetPr(prhs[3]);
  gamma = *mxGetPr(prhs[4]);

  v = mxGetPr(plhs[0]);

  if(r[0] >= k[0]) {
    mexErrMsgIdAndTxt( "MATLAB:mg_1d_MATLAB_gateway:InconsistentDataType",
      "Roughest grid must be smaller than finest.");
  }

  grid_t* g = initiate_grid((pow_t) k[0], (pow_t) r[0]);
  int i, j;
  
  // Set v_0
  for(i = 0; i<g->n; ++i){
    for(j = 0; j<g->n; ++j){
      *eget(g->g, g->c, j, i) = 0.0;
    }
  }
  // Set Boundary conditions
  for(i = 0; i<g->n; ++i){
    g->bc0[i] = 0;
    g->bc1[i] = 0;
    g->bc2[i] = 0;
    g->bc3[i] = 0;
  }

  g->bc0[g->n+1] = 0;
  g->bc0[g->n+2] = 0;
  g->bc3[g->n+1] = 0;
  g->bc3[g->n+2] = 0;

  g->f = (data_t*) f;

  for(i = 0; i<iter[0]; ++i){
    FMGV(g, gamma);
  }

  for(i=0; i<mrows; i++){
    v[i] = g->g[i];
  }

  free(g->bc0);
  free(g);
}