#include <iostream>
#include <fstream>
#include <string>
#include <ctime>
#include <cassert>

#include "doc.h"
#include "sampler.h"
#include "str_util.h"

#include "obtm.h"

OBTM::OBTM(int K, int W, int win, double a, double b, int n_iter, double l):
  K(K), W(W), win(win), n_iter(n_iter), lam(l) {
  srand(time(NULL));
  // In order to generate random-like numbers, 
  // srand is usually initialized to some distinctive runtime value, 
  // like the value returned by function time (declared in header <ctime>). 
  // This is distinctive enough for most trivial randomization needs.
	
  alpha.resize(K, a);
  beta.resize(K, W, b);
  beta_sum = beta.rowSum();

  nb_z.resize(K);
  nwz.resize(K, W);
}

/**
 *   input_dir    contains docs orangnized by days start from 0,
 *                as day0.doc_wids, day1.day_wids, day2.day_wids, ...
 *   n_day        the number of days (starting from 0)
 *   output_dir   stores the p(z) and p(w|z) results for each day
 */
void OBTM::run(string input_dir, int n_day, string res_dir) {  
  
// printf("\n\n\nWindow Size: ");
// cout << win << endl;
// printf("\n\n\n");

  for (int d = 0; d < n_day; ++d) {
      cout << "----- proc day " << d << " ------" << endl;
      if (d != 0 && d != (n_day - 1))
        prepare_day();
      else if (d == (n_day - 1))
        prepare_final_day();
      string pt = input_dir + str_util::itos(d) + ".txt";
      proc_day(pt);
      string dir = res_dir + "k" + str_util::itos(K) + ".day" + str_util::itos(d) + ".";
      save_res(dir);
  }

//  int i;
//  for(i = 0; i < phis.size(); ++i){  
//    printf("\n\n\n");
//    cout << phis[i].str() << endl;
//    printf("\n\n\n");
//  }
}

void OBTM::proc_day(string pt) {
  // Load biterms in current day
  load_docs(pt);
  // initalize topic assignments
  model_init();
  // perform Gibbs sampling on current day biterms
  Gibbs_sampling();  
}

// load biterms from docs
// input, each line is a doc with format: wid  wid  wid ...
void OBTM::load_docs(string dfile) {
  cout << "load docs: " << dfile << endl;
  ifstream rf( dfile.c_str() );
  if (!rf) {
    cout << "file not found:" << dfile << endl;
    exit(-1);
  }

  string line;
  while(getline(rf, line)) {
    Doc doc(line);
    doc.gen_biterms(bs);
  }
  cout << "n(biterms)=" << bs.size() << endl;
}

// random initialize
void OBTM::model_init() {
  // nb_z[i] == 0 at this point <-- tested
  for (vector<Biterm>::iterator b = bs.begin(); b != bs.end(); ++b) {
	int k = Sampler::uni_sample(K);

	assign_biterm_topic(*b, k);
  }
}

void OBTM::prepare_day() {
  // prepare the Dirichlet priors
  Pvec<double> t = nb_z.toDouble() * lam;
  alpha += t;
  beta += nwz.toDouble() * lam;
  beta_sum = beta.rowSum();

  // reset the states members
  bs.clear();
  nb_z.fill(0);
  nwz.fill(0);
}

void OBTM::prepare_final_day() {

//  printf("\n\n\nWindow Size: ");
//  cout << win << endl;
//  printf("\n\n\n");

  // prepare the Dirichlet priors
  Pvec<double> t = nb_z.toDouble() * lam;
  alpha += t;

  Pmat<double> weights;
  weights.resize(win, K, 0.00);

  Pmat<double> tempT_1;
  Pmat<double> tempT_2;
  Pmat<double> beta_new;

  beta_new.resize(K, W, 0.000);

  tempT_1 = phis[phis.size() - 1];

  //einstein summation
  int i, j, k, row = win - 1;
  double sum;
  for(i = phis.size() - win ; i < phis.size(); ++i, --row){
    tempT_2 = phis[i];
    //printf("\n\n\ntest = %.3f\n\n\n", tempT_2[0][0]);
    for(j = 0; j < K; ++j){
      sum = 0;
      for(k = 0; k < W; ++k){
          sum += tempT_1[j][k] * tempT_2[j][k];
      }
      weights[row][j] = sum;
    }
  }

  // taking exp for every element
  for(i = 0; i < win; ++i){
    for(j = 0; j < K; ++j){
      weights[i][j] = exp(weights[i][j]);
    }
  }

  // applying normalization for each topic
  for(i = 0; i < K; ++i){
    sum = 0;
    for(j = 0; j < win; ++j){
      sum  += weights[j][i];
    }
    for(j = 0; j < win; ++j){
      weights[j][i]  = weights[j][i]/sum;
    }
  }
  
//  printf("\n\n\n");
//  cout << weights.str() << endl;
//  printf("\n\n\n");

  row = win - 1;
  for(i = phis.size() - win; i < phis.size() ; ++i, --row){ //(phis.size()+1-win)
      tempT_2 = phis[i];
      for(j = 0; j < K; ++j){ //topic
        //cout << weights[(phis.size() - 1)-i][j] << endl;
        for(k = 0; k < W; ++k){ //word
          beta_new[j][k] += (weights[row][j]) * tempT_2[j][k];
        }
      }
  }

  beta_new += nwz.toDouble() * lam;
  beta = beta_new;
  //beta += nwz.toDouble() * lam;
  beta_sum = beta.rowSum();

//  printf("\n\n\n");
//  cout << beta.str() << endl;
//  printf("\n\n\n");
//  cout << beta_new.str() << endl;
//  printf("\n\n\n");

  // reset the states members
  bs.clear();
  nb_z.fill(0);
  nwz.fill(0);
}


void OBTM::Gibbs_sampling(){  
  
  // printf("\n\n\n");
  // cout << beta.str() << endl;
  // printf("\n\n\n");

  cout << "Begin iteration" << endl;
  for (int it = 1; it < n_iter + 1; ++it) {
    cout << "\riter " << it << '/' << n_iter;
    fflush(stdout);
    for (int b = 0; b < bs.size(); ++b) {
      update_biterm(bs[b]);
    }
  }
}

// sample procedure for ith biterm 
void OBTM::update_biterm(Biterm& bi) {
  if (bi.get_z() != -1) 
	reset_biterm_topic(bi);
  
  // compute p(z|b)
  Pvec<double> pz;
  compute_pz_b(bi, pz);

  // sample topic for biterm b
  int k = Sampler::mult_sample(pz.to_vector());
  assign_biterm_topic(bi, k);
}

// reset topic assignment of biterm i
void OBTM::reset_biterm_topic(Biterm& bi) {
  int w1 = bi.get_wi();
  int w2 = bi.get_wj();
  int k = bi.get_z();
  
  --nb_z[k];				// update proportion of biterms in topic K
  --nwz[k][w1];			// update proportion w1's occurrence times in topic K
  --nwz[k][w2];
  assert(nb_z[k]>=0 && nwz[k][w1]>=0 && nwz[k][w2]>=0);
  
  bi.reset_z();
}

// compute p(z|w_i, w_j)
void OBTM::compute_pz_b(Biterm& bi, Pvec<double>& pz) {
  pz.resize(K);
  
  int w1 = bi.get_wi();
  int w2 = bi.get_wj();
  
  double pw1k, pw2k;
  for (int k = 0; k < K; ++k) {	
	pw1k = (nwz[k][w1] + beta[k][w1]) / (2 * nb_z[k] + beta_sum[k]);
	pw2k = (nwz[k][w2] + beta[k][w2]) / (2 * nb_z[k] + 1 + beta_sum[k]);
	pz[k] = (nb_z[k] + alpha[k]) * pw1k * pw2k;
  }
  pz.normalize();
}

// assign topic k to biterm i
void OBTM::assign_biterm_topic(Biterm& bi, int k) {
  int w1 = bi.get_wi();
  int w2 = bi.get_wj();
  bi.set_z(k);
	
  ++nb_z[k];
  ++nwz[k][w1];
  ++nwz[k][w2];
}

void OBTM::save_res(string dir) {  
  string pt = dir + "pz";
  cout << "\nwrite p(z): " << pt << endl;
  save_pz(pt);
  
  string pt2 = dir + "pw_z";
  cout << "write p(w|z): " << pt2 << endl;
  save_pw_z(pt2);
}

// p(z) is determinated by the overall proportions
// of biterms in it
void OBTM::save_pz(string pt) {  
  Pvec<double> pz(K);	          // p(z) = theta
  for (int k = 0; k < K; k++) 
	pz[k] = (nb_z[k] + alpha[k]);
  pz.normalize();

  pz.write(pt);
}

void OBTM::save_pw_z(string pt) {  
  Pmat<double> pw_z = nwz.toDouble() + beta;   // p(w|z) = phi, size K * M
  pw_z.normr();
  phis.push_back(pw_z);
  pw_z.write(pt);
}
