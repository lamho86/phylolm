#include "phylolm.h"

void transbranchlengths_IvesGarland2010(int *Np, int *des, int *look, int *look2, double *distFromRoot, int *externalEdge, double *mu, double *pp, double *alphap, double *D, double *edgelength, double *diag){
  // Np:   number of edges --or pointer to
  // des:  vector of indices of the edge descendants
  // look: index of the ancestor   of edge i in 'distFromRoot'
  // look2 index of the descendant of edge i in 'distFromRoot'
  // externalEdge: 0/1 values to indicate which edges are external
  // mu:   mean at tips on linear scale of logistic model: mu = X*beta
  // pp:   mean of mu (arithmetic)
  // alpha phylogenetic correlation parameter
  // D:    extra length that was added to external edges to make tree ultrametric
  // edgelength: BM + edgelength provide the same covariance matrix
  // diag:       diagonal correction to a generalized 3-point structure.
	
  int N = *Np;
  double p = *pp;
  double alpha= *alphap;

  int i=0;
  int di=0;
  double d1=0;
  double d2=0;
  double psq=sqrt(1.0 - p)/sqrt( p );
  double ipsq=1/psq;
  double m=0;
  double m2=0;
  for (i = 0; i<N; i++) {
    di=des[i]-1;
    d1 = distFromRoot[look[i]-1];
    if (externalEdge[i]) {
      if (mu[di] < p) {m = mu[di]*psq;} 
      else {m = (1-mu[di])*ipsq;}
      m2=m*m;
      d2 = exp(-2.0 * alpha * D[di]) * mu[di] * (1-mu[di]) / m2;
      diag[di] = m*exp(alpha*D[di]);				
    } 
    else {
      d2 = distFromRoot[look2[i]-1];
    }
    edgelength[i] = d2 - d1;				
  }
}
