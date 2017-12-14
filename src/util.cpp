#include "EEGpal.h"

/***********************************************************************************/
/* meanNA()                                                                        */
/***********************************************************************************/
double meanNA(vec &y)
{
  int i, n, ny;
  double m;
  ny = y.n_elem;
  m  = 0.0;
  n  = 0;
  for(i=0; i<ny; i++){
    if(y(i) != 12345.12345){
      m += y(i);
      n ++;
    }
  }
  return m/n;
}


/***********************************************************************************/
/* longslice()                                                                     */
/***********************************************************************************/
mat longslice(cube &A, int r)
{
  int i, j;
  mat b(A.n_rows, A.n_slices); b.fill(0.0);
  for(i=0; i<A.n_rows; i++){
    for(j=0; j<A.n_slices; j++){
      b(i,j) = A(i,r,j);
    }
  }
  return b;
}
/***********************************************************************************/
/* flatslice()                                                                     */
/***********************************************************************************/
mat flatslice(cube &A, int i)
{
  int r, j;
  mat b(A.n_cols, A.n_slices); b.fill(0.0);
  for(j=0; j<A.n_cols; j++){
    for(r=0; r<A.n_slices; r++){
      b(j,r) = A(i,j,r);
    }
  }
  return b;
}


/***********************************************************************************/
/* trisolve()                                                                      */
/* Sole triangular system Lx = b                                                   */
/***********************************************************************************/
void trisolve(vec &solution, mat const &L, vec const &x, bool const &lower)
{
  int i,j;
  double sum;
  int n = x.n_elem;
  //
  if(lower){ // Lower (forward) ---------------------
    solution(0) = x(0)/L(0,0);
    for(i=1; i<n; i++){
      sum = 0.0;
      for(j=0; j<i; j++)  sum +=  L(i,j)*solution(j);
      solution(i) = (x(i) - sum)/L(i,i);
    }
  }
  if(!lower){ // Upper (backward) --------------------
    solution((n-1)) = x((n-1))/L((n-1),(n-1));
    for(i=(n-2); i>=0; i--) {
      sum = 0.0;
      for(j=i+1; j<n; j++) sum +=  L(i,j)*solution(j);
      solution(i) = (x(i) - sum)/L(i,i);
    }
  }
}
/***********************************************************************************/
/* BayesReg()                                                                      */
/* Sample x ~ N_p(Q^-1b, Q^-1)                                                     */
/***********************************************************************************/
vec BayesReg(vec const &b, mat const &Q)
{
  vec z, x, m;
  int p = b.n_elem;
  x = zeros<vec>(p); m = zeros<vec>(p);
  //
  mat U = chol(Q, "upper");       // Q = U'U : U is upper triang
  trisolve(x, trans(U), b, true); // U'(Um) = b; Um = x
  trisolve(m, U, x, false);       // m = Q^-1 b
  z = randn<vec>(p);              // z ~ N_p(0, I)
  trisolve(x, U, z, false);       // x ~ N_p(0, Q^-1)
  z = x + m;                      // z ~ N_p(b, Q^-1)
  return z;
}
/***********************************************************************************/
/* BayesRegL()                                                                     */
/* Sample x ~ N_p((LL')^-1b, (LL')^-1)                                             */
/***********************************************************************************/
vec BayesRegL(vec const &b, mat const &L)
{ // L is lower triang
  vec z, x, m;
  int p = b.n_elem;
  //
  x = zeros<vec>(p); m = zeros<vec>(p);
  //
  trisolve(x, L, b, true);            //
  trisolve(m, trans(L), x, false);    // m = Q^-1 b
  z = randn<vec>(p);                  // z ~ N_p(0, I)
  trisolve(x, trans(L), z, false);    // x ~ N_p(0, Q^-1)
  z = x + m;                          // z ~ N_p(b, Q^-1)
  return z;
}
/***********************************************************************************/
/* rWish()                                                                         */
/* Sample W ~ Wishart(S, v)                                                        */
/***********************************************************************************/
mat rWish(mat const &V, double v)
{
  int m = V.n_rows;
  mat T = zeros<mat>(m,m);
  // Barlett decomposition;
  for(int i = 0; i < m; i++) {
    T(i,i) = std::pow(rchisq(1,v-i)[0], 0.5);
  }
  for(int j = 0; j < m; j++) {
    for(int i = j+1; i < m; i++) {
      T(j,i) = rnorm(1)[0];
    }
  }
  mat C = T * chol(V);
  return trans(C) * C ;
}
/***********************************************************************************/
/* rMN()                                                                           */
/* Sample X ~ MN(m, S, V)                                                          */
/***********************************************************************************/
mat rMN(mat const &m, mat const &S, mat const& V)
{
  int p, k;
  mat z, Sc, Vc;
  //
  p  = S.n_cols;
  k  = V.n_cols;
  //
  Sc = chol(S, "lower");                    // Cholesky S = ScSc'
  Vc = chol(V, "upper");                    // Cholesky V = Vc'Vc
  //
  z  = randn<mat>(p, k);           // MN(0, I_p, I_k)
  return m + Sc*z*Vc;       // MN(m, Sc'Sc, Vc'Vc)
}
/***********************************************************************************/
/* cov2cor()                                                                       */
/* Transform a covariance matrix W into a correlation matrix R                     */
/***********************************************************************************/
mat cov2cor(mat S)
{
  int m = S.n_cols;
  mat D = eye<mat>(m, m);
  for(int i=0; i<m; i++) D(i,i) = pow(S(i,i), -0.5);
  return D * S * D;
}
/***********************************************************************************/
/* vec2mat()                                                                       */
/* column-wise matricization                                                       */
/***********************************************************************************/
mat vec2mat(vec const &v, int const& nr, int const& nc)
{
  int i, j, index;
  mat out(nr,nc);
  for(j=0; j<nc; j++){
    for(i=0; i<nr; i++){
      index = nr*j+i;
      out(i,j) = v[index];
    }
  }
  return out;
}


