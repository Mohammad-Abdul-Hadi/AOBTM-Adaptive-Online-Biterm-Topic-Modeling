#include <cstdlib>
#include <string.h>
#include <string>
#include <iostream>
#include <ctime>

#include "str_util.h"
#include "obtm.h"
#include "ibtm.h"

using namespace std;

int main(int argc, char* argv[]) {
  int i = 1;
  string method(argv[i++]); //method to be used, obtm or ibtm
  // common parameters
  int K = atoi(argv[i++]);          // topic number
  int W = atoi(argv[i++]);			    // vocabulary size
  double alpha = atof(argv[i++]);   // hyperparameters of p(z)
  double beta = atof(argv[i++]);    // hyperparameters of p(w|z)
  string input_dir(argv[i++]);
  int n_day = atoi(argv[i++]);	// 
  string res_dir(argv[i++]);
  int win = atoi(argv[i++]);
  printf("=== %s: end_day %d , K %d, W %d, window-size %d, alpha %f, beta %f ====\n", method.c_str(), n_day, K, W, win, alpha, beta);
 
  // model-sepcific parameters
  if (method == "obtm") {
      int n_iter = atoi(argv[i++]);
      double lam = atof(argv[i++]); //lambda
      printf("Model=obtm, n_iter=%d, lam=%f\n", n_iter, lam);
      OBTM mod(K, W, win, alpha, beta, n_iter, lam);
      mod.run(input_dir, n_day, res_dir);
  }

  else if (method == "ibtm") {
      int win = atoi(argv[i++]);
      int n_rej = atoi(argv[i++]);
      printf("Model=ibtm, win=%d, n_rej=%d\n", win, n_rej);
      IBTM mod(K, W, alpha, beta, win, n_rej);
      mod.run(input_dir, n_day, res_dir);
  }

  else {
      cout << "Unrecognized method:" << argv[1] << endl;
      exit(-1);
  }
}
