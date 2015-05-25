#include "phylolm.h"

double smartadd(double a, double b){
  // returns log(exp(a) + exp(b)), but avoids numerical issues: a + log(1+exp(b-a)) if a>b
  if (a>b) return(a + log(1+exp(b-a)));
  else return(b + log(1+exp(a-b)));
}

void logistreglikelihood (int *Ne, int *ntip, int *Ni, int *rootp, double *len, int *des, int *anc, int *y, double *mu, int *dkp, double *alphap, double *loglikelihood){

  int Nno  = *ntip + *Ni; // total # nodes: tips + internal (including root). should = Ne+1.
  int root = *rootp;      // root node ID (in R)
  root--;
  // len:  vector of edge lengths
  // des:  vector of descendents (one for each edge)
  // anc:  vector of ancestor nodes
  // y:    vector of binary response values
  // mu:   vector of predicted response: 1/(1+exp(-X*beta))
  int dk   = *dkp;     // # columns in X. Assumes intercept present, i.e. dk=1 equivalent to intercept only.
  double alpha = *alphap;  // phylogenetic correlation.
  // output: fixit vector of size 4+2dX+dX^2: logdet(V), 1'V^{-1}1,
  //         y'V^{-1}1 (dY values), y'V^{-1}y (dY^2 values),  
  //         X'V^{-1}1 (dX values), X'V^{-1}X (dX^2 values), X'V^{-1}y (dX*dY values)

  double* llk0 = (double*)calloc(Nno, sizeof(double));
  double* llk1 = (double*)calloc(Nno, sizeof(double));
  // llk0 = log-likelihood of data at descendants of a node, conditional on that node being in state 0
  // llk1 = same, but conditional on that node being in state 1.

  double meanp = 0.0;
  for (int i=0; i< *ntip; i++)
    meanp += mu[i];
  meanp /= *ntip;
  double meanq = 1 - meanp;
  // two-state Markov process with rates q_01=meanp*alpha, q_10=meanq*alpha.

  for (int i=0; i<Nno; i++){
    llk0[i] = 0.0;
    llk1[i] = 0.0;
  }

  // loop over all Ne edges. Assumes postorder traversal.
  for (int iedge=0; iedge < *Ne; iedge++){
    int di= des[iedge]-1; // C-index of descendant node. R-index: di+1. 
    int ai= anc[iedge]-1;

    if (di< *ntip){ // external edge: initialize llk0 and llk1 based on data y
      if (y[di]) 
	llk0[di] = -1.0/0.0; // y=1. This is to get -infinity = log(0)
      else
	llk1[di] = -1.0/0.0; // y=0
      if (dk > 1) { // X has covariates, in addition to intercept
	// et = min( mu[di]/meanp , (1-mu[di])/meanq);
	if ( mu[di] < meanp ){
	  double et = mu[di]/meanp; // if y=1, it's flipped to 0 with prob 1-et
	  llk1[di] = (y[di] ? log(et) : log(1-et));
	} else {
	  double et = (1-mu[di])/meanq; // if y=0, it's flipped to 1 with prob 1-et
	  llk0[di] = (y[di] ? log(1-et) : log(et));
	}
      }
    }
    // update likelihood at parent node, contribution from current edge
    double eta = exp(-len[iedge]*alpha);
    double peta = meanp * eta;
    double qeta = meanq * eta;
    llk0[ai] += smartadd(log(meanq + peta)+llk0[di] , log(meanp - peta)+llk1[di]);
    llk1[ai] += smartadd(log(meanq - qeta)+llk0[di] , log(meanp + qeta)+llk1[di]);
  }
  // assuming stationary distribution at the root:
  loglikelihood[0] = smartadd(log(meanq) + llk0[root] , log(meanp) + llk1[root]);
  free(llk0);
  free(llk1);
}
