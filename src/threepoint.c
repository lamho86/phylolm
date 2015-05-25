#include "phylolm.h"

void threepoint (int *Npo, int *npo, int *pNpo, int *dYpo, int *dXpo, int *rootpo, double *transa, double *transb, int *des, int *anc, double *y, double *X, double *output){

  int N = *Npo;     // number of edges, excluding root edge
  int n = *npo;     // number of tips
  int pN=*pNpo;     // number of internal nodes
  int npN = n + pN; //      total # nodes, including root node. N+1 if binary tree.
  int dY=*dYpo;     // number of columns in Y (dependent variables, usually dY=1)
  int dX=*dXpo;     // number of columns in X (independent variables)
  int r=*rootpo;    // root node ID (in R)
  r--;
  double rootEdge= *transa; // length of root edge
  // transb: vector of edge lengths
  // des:    vector of descendents (one for each edge)
  // anc:    vector of ancestor nodes
  // y:  vector of response values
  // X:  vectorized matrix X, where matrix model is: y = X*beta + error
  // output: vector of size 4+2dX+dX^2: logdet(V), 1'V^{-1}1,
  //         y'V^{-1}1 (dY values), y'V^{-1}y (dY^2 values),  
  //         X'V^{-1}1 (dX values), X'V^{-1}X (dX^2 values), X'V^{-1}y (dX*dY values)

  double* logd=(double*)calloc(npN, sizeof(double));
  // logd[i-1] = log|V_subtree| for subtree at *node* i, or parent edge of node i.
  double* vec11=(double*)calloc(npN, sizeof(double));
  double* yy=(double*)calloc(npN*dY*dY, sizeof(double));
  // first (Y'V^{-1}Y)[1,1] at all nodes, then (Y'V^{-1}Y)[2,1] at all nodes, etc.
  double* y1=(double*)calloc(npN*dY, sizeof(double));
  double* Xy=(double*)calloc(npN*dX*dY, sizeof(double));
  double* X1=(double*)calloc(npN*dX, sizeof(double));
  double* XX=(double*)calloc(npN*dX*dX, sizeof(double));
  int* zero =(int*)calloc(npN, sizeof(int));
  // zero[i-1] will be -1 if node i has no children edge of length 0
  //                    d if node i has exactly 1 child edge of length 0,
  //                      to child node d+1.
  // error returned if some node has 2 or more child edges of length 0
  for (int iedge=0; iedge<N+1; iedge++)
    zero[iedge] = -1;

  // loop over all N+1 edges, with root edge last. Assumes postorder traversal.
  for (int iedge=0; iedge < N+1; iedge++){
    double el, invel;  // edge length and its inverse
    int di, anci=0;
    if (iedge<N){      // non-root edge
      el=transb[iedge]; 
      di= des[iedge]-1;// subtree at node di+1. 
      anci= anc[iedge]-1;
    } else {           // iedge=N -> root edge
      el = rootEdge;
      di = r; // but anci meaningless
    }

    int iY1 = di; // index in y1 for current edge and column j in Y: di + npN*j
    int iX1 = di; //          X1                                  X: di + npN*j
    int iYY = di; // index in yy for columns k and j in Y: di + npN*k + npN*dY*j
    int iXY = di; // index in Xy for column k in X and j in Y: di + npN*k + npN*dX*j
    int iXX = di; // di + npN*k + npN*dX*j for XX: columns k and j in X.

    if (di<n){ // external edge
      if(el>0) invel = 1/el;// y1, X1 will actually hold ybar*p and Xbar*p.
      else     invel = 1.0; // in this case, we need the staight y'1, y'y, X'1, X'X etc.
      if (el>0){
	logd[di]=log(el);
	vec11[di]=invel;
      } else{
	if (zero[anci] >= 0) // anci already found with 1 child edge of length 0.
	  error("two or more sister external edges have length 0, V is singular\n (node %d in pruning-wise-ordered tree)", anci+1);
	else
	  zero[anci] = di; // which is >= 0
      }
      int jY  = di; // index in y for observations in species di, column j: di + n*j
      for (int j=0; j<dY; j++){
	y1[iY1] = y[jY]*invel;
	int kY = di; // index in y for observation in species di, column k: di + n*k
	for (int k=0; k<dY; k++){
	  yy[iYY] = y1[iY1] * y[kY];
	  iYY += npN;
	  kY  += n;
	}
	int kX = di; // index in X for observation in species di, column k: di + n*k
	for (int k=0; k<dX; k++){
	  Xy[iXY] = y1[iY1] * X[kX];
	  iXY += npN;
	  kX  += n;
	}
	iY1 += npN;
	jY  += n;
      }
      int jX  = di; // di + n*j   for X, species di, column j.
      for (int j=0; j<dX; j++){
	X1[iX1] = X[jX] * invel;
	int kX = di;
	for (int k=0; k<dX; k++){
	  XX[iXX] = X1[iX1] * X[kX];
	  iXX += npN;
	  kX  += n;
	}
	iX1 += npN;
	jX  += n;
      }
    }
    else{ // internal edge. contributions from descendants of di have already been collected.
      int goodchildren = 1;
      double ev, ev2; 
      if (zero[di] >= 0) { // exactly 1 child edge of length 0, to descendant d0=zero[di]
	if (el<=0)
	  error("One external edge and its parent both have length 0\n (node %d in pruning-wise-ordered tree). To avoid this situation,\n please make a polytomy or resolve it differently ",di+1);
	goodchildren = 0;
	// still assuming di has more than 1 child, with current vec11[di]>0.
      }
      if (goodchildren) {
	logd[di] += log(1+ el * vec11[di]);
	ev  = el/(1+ el * vec11[di]);
	ev2 =  1/(1+ el * vec11[di]);
	for (int j=0; j<dY; j++){
	  double tmp1 = ev * y1[iY1];
	  int kY1 = di; // index in ybar, column k: di + npN*k
	  for (int k=0; k<dY; k++){
	    yy[iYY] -= tmp1 * y1[kY1];
	    iYY += npN;
	    kY1 += npN;
	  }
	  int kX1 = di; // index in Xbar, column k: di + npN*k
	  for (int k=0; k<dX; k++){
	    Xy[iXY] -= tmp1 * X1[kX1];
	    iXY += npN;
	    kX1 += npN;
	  }
	  iY1 += npN;
	}
	for (int j=0; j<dX; j++){
	  double tmp1 = ev * X1[iX1];
	  int kX1 = di; // index in Xbar, column k: di + npN*k
	  for (int k=0; k<dX; k++){
	    XX[iXX] -= tmp1 * X1[kX1];
	    iXX += npN;
	    kX1 += npN;
	  }
	  iX1 += npN;
	}
      }
      else { // 1 bad child
	logd[di] += log(el);
	int d0 = zero[di];
	double fac = 1/el + vec11[di];
	// X'V^{-1}y = Xytilde + (1/t +ptilde_A) X'd0 yd0
	//            - ptilde_A X'd0 y1tilde - ptilde_A X1tilde' yd0
	// tilde: values collected from "normal" children
	// d0:    values from the special child with edge length 0
	iY1 = 0; iX1 = 0;
	for (int j=0; j<dY; j++){
	  double tmp1 = fac * y1[d0+iY1];
	  int kY1 = 0; // index in ybar, column k: (di or d0) + npN*k
	  for (int k=0; k<dY; k++){
	    yy[iYY] += tmp1 * y1[d0+kY1] - y1[d0+iY1]*y1[di+kY1] - y1[di+iY1]*y1[d0+kY1];
	    iYY += npN;
	    kY1 += npN;
	  }
	  int kX1 = 0; // index in Xbar, column k: (di or d0) + npN*k
	  for (int k=0; k<dX; k++){
	    Xy[iXY] += tmp1 * X1[d0+kX1] - y1[d0+iY1]*X1[di+kX1] - y1[di+iY1]*X1[d0+kX1];
	    iXY += npN;
	    kX1 += npN;
	  }
	  iY1 += npN;
	}
	for (int j=0; j<dX; j++){
	  double tmp1 = fac * X1[d0+iX1];
	  int kX1 = 0;
	  for (int k=0; k<dX; k++){
	    XX[iXX] += tmp1 * X1[d0+kX1] - X1[d0+iX1]*X1[di+kX1] - X1[di+iX1]*X1[d0+kX1];
	    iXX += npN;
	    kX1 += npN;
	  }
	  iX1 += npN;
	}
      }
      if (goodchildren) {
	// rescaling to get vec11 correct ("p" instead of "p_A"),
	// and Y1 = ybar*p, X1 = bar*p instead of ybar*p_A etc.
	iY1 = di;
	for (int j=0; j<dY; j++){
	  y1[iY1] *= ev2;
	  iY1 += npN;
	}
	iX1 = di;
	for (int j=0; j<dX; j++){
	  X1[iX1] *= ev2;
	  iX1 += npN;
	}
	vec11[di] *= ev2;
      }
      else {
	invel = 1/el;
	iY1 = 0;
	for (int j=0; j<dY; j++){
	  y1[di+iY1] = y1[zero[di]+iY1]*invel;
	  iY1 += npN;
	}
	iX1 = 0;
	for (int j=0; j<dX; j++){
	  X1[di+iX1] = X1[zero[di]+iX1]*invel;
	  iX1 += npN;
	}
	vec11[di] = invel;
      }
    }
    // next: collect contribution of iedge to its parent, i.e. to ancestor node anci.
    // *except* if root edge (anci meaningless) or if external edge of zero length
    if ((iedge < N) && ((di>=n) || (el>0))){
      logd[anci] += logd[di];
      iY1 = 0; iX1 = 0;
      iYY = 0; iXX = 0; iXY = 0;
      for (int j=0; j<dY; j++){
        y1[anci+iY1] += y1[di+iY1];
        for (int k=0; k<dY; k++){
          yy[anci + iYY] += yy[di + iYY];
          iYY += npN;
        }
        for (int k=0; k<dX; k++){
          Xy[anci + iXY] += Xy[di + iXY];
          iXY += npN;
        }
        iY1 += npN;
      }
      for (int j=0; j<dX; j++){
        X1[anci+iX1] += X1[di+iX1];
        for (int k=0; k<dX; k++){
          XX[anci + iXX] += XX[di + iXX];
          iXX += npN;
        }
        iX1 += npN;
      }
      vec11[anci] += vec11[di];	
    }
  }
  // combine all results at root edge into output
  output[0]=logd[r];
  output[1]=vec11[r];
  int p=2; // # parameters already logged
  int ikXY = r;
  for (int j=0; j<dY; j++){
    output[p+j] = y1[ikXY];
    ikXY += npN;
  }
  p += dY;
  ikXY = r;
  for (int j=0; j<dY*dY; j++){
    output[p+j] = yy[ikXY];
    ikXY += npN;
  }
  p += dY*dY;
  ikXY = r;
  for (int j=0; j<dX; j++){
    output[p+j] = X1[ikXY];
    ikXY += npN;
  }
  p += dX;
  ikXY = r;
  for (int j=0; j<dX*dX; j++){
    output[p+j] = XX[ikXY];
    ikXY += npN;
  }
  p += dX*dX;
  ikXY = r;
  for (int j=0; j<dX*dY; j++){
    output[p+j] = Xy[ikXY];
    ikXY += npN;
  }
  
  free(logd);
  free(vec11);
  free(y1);
  free(yy);
  free(X1);
  free(XX);
  free(Xy);
  free(zero);
}	

