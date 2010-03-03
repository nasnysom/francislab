/* NTE.cpp
 *
 * Written October 2009 by Marcello DiStasio
 * 
 * Computes the Normalized Transfer Entropy from one signal to
 * another
 */


//Includes for MSVC Express 2005 used in MATLAB
#include <iostream>
#include <iomanip>
#include <cstdio> // for writing results to a file
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>
#include <fstream>
#include <cmath>
#include <cstdlib>

#include "mex.h" // MATLAB
//#include "matrix.h" // MATLAB


using namespace std;

class NTE_Calculator {
  
public: 

  //Vector of binned spike times (or other binned data)
  vector<int> X_binned[2];

  //Past sums
  vector<int> X_P[2];
  //Future sums
  vector<int> X_F[2];

  // Size of the signal vectors (which must be all equal) after windowing by tau
  unsigned int SizeOfAllWindowed;

  //Unique elements in each (see findUniqueElements method)
  vector<int> X_F_1_uniques;
  vector<int> X_P_1_uniques;
  vector<int> X_P_0_uniques;

  int tau_bins_F;
  int tau_bins_P;

  // The joint and single match counts (one for every element in X_F_1_uniques, X_P_1_uniques, X_P_0_uniques)
  vector<int> klm_Matches, kl_Matches, lm_Matches, l_Matches;

  // the return variable TE_total for method TE() is the total transfer entropy calculated 
  double TE_total;

  //Constructor
  NTE_Calculator();

  //Methods
  void set_Xbinned(vector<int>& spikes_binned0, vector<int>& spikes_binned1);
  
  void set_tau_bins_P(int bins) { tau_bins_P = bins; }
  void set_tau_bins_F(int bins) { tau_bins_F = bins; }


  void findUniqueElements();
  vector<double> readFile(std::istream& is);
  double Cond_Entropy();
  double TE();
  double NTE();
  
};

//Constructor
inline NTE_Calculator::NTE_Calculator() { 

  SizeOfAllWindowed = 0;

  set_tau_bins_F(10); 
  set_tau_bins_P(10); 

}

//Setter methods
inline void NTE_Calculator::set_Xbinned(vector<int>& spikes_binned0, vector<int>& spikes_binned1) {

  X_binned[0].clear();
  X_binned[1].clear();

  X_binned[0] = spikes_binned0;
  X_binned[1] = spikes_binned1;

}

inline double NTE_Calculator::NTE() {

  const unsigned int shuff_reps = 100;

  // Calculate the Transfer Entropy of the signals
  double TE_real = TE();
  //  cout << "TE_real: " << TE_real << "\n";
  
  double CE = Cond_Entropy();
  //  cout << "CE: " << CE << "\n";

  // Shuffle the signals and calculate the transfer entropy shuff_reps times
  vector<double> shuff_TE;
  unsigned int i,j;
  for (i = 0; i < shuff_reps; i++) {
    random_shuffle( X_binned[0].begin(), X_binned[0].end() );
    //recalculate the TE and add it to the shuffled results vector
    shuff_TE.push_back(TE());
  }

  double TE_shuff_mean = 0;
  for (j=0; j < shuff_TE.size(); j++) {
    TE_shuff_mean = TE_shuff_mean + shuff_TE[j];
  }
  TE_shuff_mean = TE_shuff_mean / shuff_TE.size();
  
  return (TE_real - TE_shuff_mean) / CE;

}

inline double NTE_Calculator::Cond_Entropy( ) {
  
  // Calculates the conditional entropy
  // i.e. H(X_2_Future|X_2_Past)
  
  double CE_Summand;
  double CE_Total = 0;

  double P_l;
  double P_kl;
  
  //Note log() is the natural log, so to convert to base 2, we divide the result by log(2)
  
  for (unsigned int l = 0; l < X_P_1_uniques.size(); l++) {
    
    P_l = double(l_Matches[l]) / double(SizeOfAllWindowed);
    for (unsigned int k = 0; k < X_F_1_uniques.size(); k++) {
      
      P_kl = double(kl_Matches[(X_P_1_uniques.size()*k) + l]) / double(SizeOfAllWindowed);
      CE_Summand = 0.0;
      
      if (P_l <= 0) {
	CE_Summand = 0.0; 
      } else {
	if (P_kl <= 0) {
	  CE_Summand = 0.0; 
	}
	else {
	  CE_Summand = -(P_kl * log((P_kl)/(P_l))/log(double(2)));
	  
	}
	CE_Total = CE_Total + CE_Summand;
      }
    }
  }
  
  return CE_Total;
  
}


inline void NTE_Calculator::findUniqueElements() {
  
  unsigned int i;
  unsigned int j;
  bool match;

  X_F_1_uniques.clear();
  X_P_1_uniques.clear();
  X_P_0_uniques.clear();

  //-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
  for (i = 0; i < X_F[1].size(); i++) {
    match = 0;
    for (j = 0; j < X_F_1_uniques.size(); j++) {
      if (X_F[1][i] == X_F_1_uniques[j])
	match = 1;
    }
    if (!match)
      X_F_1_uniques.push_back(X_F[1][i]);
  }
  //----------------------------------------------
  for (i = 0; i < X_P[1].size(); i++) {
    match = 0;
    for (j = 0; j < X_P_1_uniques.size(); j++) {
      if (X_P[1][i] == X_P_1_uniques[j])
	match = 1;
    }
    if (!match)
      X_P_1_uniques.push_back(X_P[1][i]);
  }
  //----------------------------------------------
  for (i = 0; i < X_P[0].size(); i++) {
    match = 0;
    for (j = 0; j < X_P_0_uniques.size(); j++) {
      if (X_P[0][i] == X_P_0_uniques[j])
	match = 1;
    }
    if (!match)
      X_P_0_uniques.push_back(X_P[0][i]);
  }
  //-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

  /*
  cout << "Unique Elements of X_F[1]: \n\n";
  for (unsigned int y = 0; y < X_F_1_uniques.size(); y++) {
    cout << X_F_1_uniques[y] << "\n";
  }
  cout << "Unique Elements of X_P[1]: \n\n";
  for (unsigned int y = 0; y < X_P_1_uniques.size(); y++) {
    cout << X_P_1_uniques[y] << "\n";
  }
  cout << "Unique Elements of X_P[0]: \n\n";
  for (unsigned int y = 0; y < X_P_0_uniques.size(); y++) {
    cout << X_P_0_uniques[y] << "\n";
  }
  cout << ".,.-'-.,.,.-'-.,.,.-'-.,.\n\n";
  */

  return;

}

inline vector<double> NTE_Calculator::readFile(std::istream& is) {
  
  char chdummy;
  is >> std::ws >> chdummy >> std::ws; 
  std::vector<double> result;
  for(;;) {
    double number;
    is>>number;
    if (is.eof()) break;
    result.push_back(number);
    //	cout << number << '\n';
  }
  
  return result;
}


inline double NTE_Calculator::TE() {

  //return variable giving the Transfer Entropy between spike time signals X_binned[0] and [1]
  TE_total = 0.00;

  //  cout << "Size of vectors: " << X_binned[0].size() << " and " << X_binned[1].size() << "\n\n";

  if (X_binned[0].size() != X_binned[1].size()) {
    mexPrintf("Binned spike time vectors must be equal length");
    return -8;
  }

  //find the maximum value in the vec X_binned of spike times
  /*
  vector<int>::const_iterator X_max[2];
  int X_max_both;
  X_max[0]=max_element( X_binned[0].begin(), X_binned[0].end() );
  X_max[1]=max_element( X_binned[1].begin(), X_binned[1].end() );
  if (*X_max[0] > *X_max[1])
    X_max_both = *X_max[0];
  else
    X_max_both = *X_max[1];
  */
  //cout << "Vector 0 Max (before summing over tau-window): " << *X_max[0] << "\n";
  //cout << "Vector 1 Max (before summing over tau-window): " << *X_max[1] << "\n";

  int au = 0;
  int bin_sum[2] = {0,0};

  //Clear the Past and Future Windowed vectors
  X_P[0].clear();
  X_P[1].clear();
  X_F[0].clear();
  X_F[1].clear();

  
  // Sum Vectors 1 AND 2 in the Past over a window of size tau_bins_P
  // Iterate by tau_bins_P (to avoid smoothing)
  for (unsigned int k = tau_bins_P; k < (X_binned[0].size() - tau_bins_F); k = k + tau_bins_P) {
    bin_sum[0] = 0;     
    bin_sum[1] = 0;

    for (au = 0; au < tau_bins_P; au ++){
      bin_sum[0] = bin_sum[0] + X_binned[0][k-au];
      bin_sum[1] = bin_sum[1] + X_binned[1][k-au];
    }

    X_P[0].push_back(bin_sum[0]);
    X_P[1].push_back(bin_sum[1]);
  }


  // Sum Vector 2 in the future over a window of size tau_bins_F.  ONLY VECTOR 2 is used beyond here for the future.  
  // Iterate by tau_bins_P (to avoid smoothing)
  for (unsigned int k = tau_bins_P; k < (X_binned[0].size() - tau_bins_F); k = k + tau_bins_F) {
    //    bin_sum[0] = 0;
    bin_sum[1] = 0;

    for (au = 0; au < tau_bins_F; au ++){
      //bin_sum[0] = bin_sum[0] + X_binned[0][((k+1)+au)];
      bin_sum[1] = bin_sum[1] + X_binned[1][((k+1)+au)];
    }
    //X_F[0].push_back(bin_sum[0]);
    X_F[1].push_back(bin_sum[1]);
  }





  // Get a set of unique elements in each of X_F[1], X_P[1] ,X_P[0]
  // So that only probabilities for these must be summed, and all others
  // may be set to zero.
  findUniqueElements();


  /*
  // Reporting of summed signal vectors
  for (unsigned int i = 0; i < X_P[0].size(); i++) {
    cout << "X_P[0][" << i << "]: " << X_P[0][i] << "\n";
    cout << "X_P[1][" << i << "]: " << X_P[1][i] << "\n";
  }
  for (unsigned int i = 0; i < X_F[1].size(); i++) {
    cout << "X_F[0][" << i << "]: " << X_F[0][i] << "\n";
    cout << "X_F[1][" << i << "]: " << X_F[1][i] << "\n";
    } */
  
  
  /*---------------------------------------------
    DO CALC HERE. (Eqn (2) from Gourevitch)
    ---------------------------------------------- */
  

  SizeOfAllWindowed = 0;
  if ((X_F[1].size() != X_P[1].size()) | (X_F[1].size()  != X_P[0].size())) {
    cerr << "In TE(): X_F[1], X_P[1], X_P[0] are not the same size!";
    return -1;
  } else { SizeOfAllWindowed = X_F[1].size(); }
  
  //k iterates all possible values of X_F[1] (receiver signal)
  //l iterates all possible values of X_P[1] (receiver signal)
  //m iterates all possible values of X_P[0] (source signal)

  int matches = 0; // temporary counter for matches (ALWAYS SET TO ZERO BEFORE USING IT IN A LOOP)
  unsigned int cnt;
  unsigned int k,l,m;


  klm_Matches.clear();
  kl_Matches.clear();
  lm_Matches.clear();
  l_Matches.clear();

  // TE_Summand is the TE calculated on each loop -- the total TE is the sum of thesee
  double TE_Summand = 0.00;

  //------------------------------------------------
  for (k = 0; k < X_F_1_uniques.size(); k++) {

    //------------------------------------------------
    for (l = 0; l < X_P_1_uniques.size(); l++) {

      //__Only calculate the number of matches for l alone on the first pass through k.
      //__This will save computation time, and we can look up the result in the array at
      //__the time of the probability calculation within each loop of m, below
      if (k==0) {

	matches = 0;
	for (cnt = 0; cnt < SizeOfAllWindowed; cnt++) {
	  if (X_P_1_uniques[l] == X_P[1][cnt]) 
	    matches++;
	}  
	//cout << "Matches for (l)=(" << X_P_1_uniques[l] << ") : "  << matches<< "\n";
	l_Matches.push_back(matches);

      }

      //__Calculate the joint probability for k and l on every pass
      matches = 0;
      for (cnt = 0; cnt < SizeOfAllWindowed; cnt++) {
	if ( (X_F_1_uniques[k] == X_F[1][cnt]) & (X_P_1_uniques[l] == X_P[1][cnt]) )
	  matches++;
      }  
      //cout << "Matches for (k,l)=(" << X_F_1_uniques[k] << "," << X_P_1_uniques[l] << ") : "  << matches<< "\n";
      kl_Matches.push_back(matches);
      
      //------------------------------------------------  
      for (m = 0; m < X_P_0_uniques.size(); m++) {
	
	//__Only calculate the number of matches for l and m (but not k) on the first pass through k.
	//__This will save computation time, and we can look up the result in the array at
	//__the time of the probability calculation at each loop through m
	
	if (k==0) {

	  matches = 0;
	  for (cnt = 0; cnt < SizeOfAllWindowed; cnt++) {
	    if ( (X_P_1_uniques[l] == X_P[1][cnt]) & (X_P_0_uniques[m] == X_P[0][cnt]) )
	      matches++;
	  }  
	  //cout << "Matches for (l,m)=(" << X_P_1_uniques[l] << "," << X_P_0_uniques[m] << ") : "  << matches << "\n";
	  lm_Matches.push_back(matches);

	}

	// Looking for a match for all three k,l,m (on every pass)
	matches = 0;
	for (cnt = 0; cnt < SizeOfAllWindowed; cnt++) {
	  if ( (X_F_1_uniques[k] == X_F[1][cnt]) & (X_P_1_uniques[l] == X_P[1][cnt]) & (X_P_0_uniques[m] == X_P[0][cnt]) )
	    matches++;
	}  
	//	cout << "Matches for (k,l,m)=(" << X_F_1_uniques[k] << "," << X_P_1_uniques[l] << "," << X_P_0_uniques[m] << ") : "  << matches<< "\n";
	
	klm_Matches.push_back(matches); // This may be removed
	// matches will be used directly in the main calculation below to save on speed

	//****************************************************************************************************************************
	//*********                                                                                                        ***********
	//*********                               Main Probability Calculation HERE                                        ***********
	//*********                                                                            -Chelly                     ***********
	//****************************************************************************************************************************

	double P_l;
	double P_lm;
	double P_kl;
	double P_klm;	

	P_l = double(l_Matches[l]) / double(SizeOfAllWindowed);
	P_lm = double(lm_Matches[(X_P_0_uniques.size()*l) + m]) / double(SizeOfAllWindowed);
	//	cout << "lm Matches: " << lm_Matches[(X_P_0_uniques.size()*l) + m] << "\n";
	P_kl = double(kl_Matches[(X_P_1_uniques.size()*k) + l]) / double(SizeOfAllWindowed);
	//	cout << "kl Matches: " <<  kl_Matches[(X_P_1_uniques.size()*k) + l] << "\n";
	P_klm = double(matches) / double(SizeOfAllWindowed);

	//Note log() is the natural log, so to convert to base 2, we divide the result by log(2)
	double prob_product;
	TE_Summand = 0.0;

	double p = (P_lm * P_kl);
	if (p <= 0)
	  TE_Summand = 0.0;
	else {
	  prob_product = (P_klm * P_l)/(P_lm * P_kl);
	  if (prob_product <= 0)
	    TE_Summand = 0.0;
	  else {
	    TE_Summand = P_klm * log((P_klm * P_l)/(P_lm * P_kl))/log(double(2));
	  }
	}

	TE_total = TE_total + TE_Summand;

	// if (TE_Summand > 0) {
	//  cout << "(k,l,m) = (" << k << "," << l << "," << m << "):  TE_total: " << TE_total;
	//  cout << "    TE_Summand: " << TE_Summand << "\n";
	// }

	//*********                                                                                                        ***********
	//*********                                            End of Main Probability Calculation                         ***********
	//****************************************************************************************************************************

      } // end for m loop
    }
  }


  return TE_total;

}



void mexFunction(int nlhs, mxArray *plhs[ ], int nrhs, const mxArray *prhs[ ]) 
{

  double *InputArray[2]; //pointer to data in MATLAB input arrays for the signals
  double *InputTauBins; // pointer to array of taubins sizes to loop through
  double *OutputArray; //pointer to data in MATLAB output array

  // Variables to hold dimensions of input arrays to MATLAB function
  int r0, c0;
  int r1, c1;
  int r2, c2;
  
  // Vector with tau_bins sizes, to be read in from input
  vector<int> tau_bins;
  
  /* Input checking ***************************************/
  //
  
  /* check: only one input and one output argument */ 
  if (nrhs < 2) {
    mexPrintf("NTE : Normalized Transfer Entropy \n\tAn Approximation to the amount of information transfer between two spike trains\n");
    mexPrintf("Usage: \n\ntotal_ NTE = NTE(X1,X2)\n\n");
    mexPrintf("total_NTE is the normalized transfer entropy (unbiased based on tau by subtracting shuffled TE, and normalized by H(X2F|X2P) )\n");
    mexPrintf("X1,X2 are row vectors of non-negative spike/event arrival times.\n\n");
    mexErrMsgTxt("");
  }
  
  /* get the dimensions of the input and check them (r = rows, c = columns)*/
  r0 = mxGetM(prhs[0]);
  c0 = mxGetN(prhs[0]);
  r1 = mxGetM(prhs[1]);
  c1 = mxGetN(prhs[1]);
  if(r0!=1 || r1!=1 || c0<1 || c1<1) {
    mexPrintf("First 2 arguments must be row vectors of length > 0\n");
    mexErrMsgTxt("");
  }

  // If there is a third argument given, use it to specify an array of tau_bins to try
  // If not, set tau_bins to 10
  if (nrhs ==3) {
    InputTauBins = mxGetPr(prhs[2]);
    r2 = mxGetM(prhs[2]);
    c2 = mxGetN(prhs[2]);
    if (r2 != 1 || c2<1) {
      mexPrintf("Third argument must be a row vector of length > 0\n");
      mexErrMsgTxt("");
    }
    //Fill the tau_bins vector with the input array
    for (int ta = 0; ta < c2; ta++) {
      tau_bins.push_back(InputTauBins[ta]);
    }

  } else {
    mexPrintf("Using default tau: 10 bins\n\n");
    tau_bins.clear();
    tau_bins.push_back(10);
  }

  /*** End Input checking *********************************/
  
  
  /* Create the output arrays */
  plhs[0] = mxCreateDoubleMatrix(1,tau_bins.size(),mxREAL); // this is the array that actually gets returned
  OutputArray = mxGetPr(plhs[0]);
  
  /* Get the input */
  InputArray[0] = mxGetPr(prhs[0]);
  InputArray[1] = mxGetPr(prhs[1]);

  
  // Set input variables, as would be done by a main function
  // X[0] and X[1] are spike times
  vector<double> X[2];
  for (int i = 0; i < c0; i++) {
    X[0].push_back(InputArray[0][i]);
  }
  for (int i = 0; i < c1; i++) {
    X[1].push_back(InputArray[1][i]);
  }

  /*
  // Check to make sure we have the same number of arrival times (Not neccessary)
  int size[2];
  size[0] = X[0].size();
  size[1] = X[1].size();
  
  if (size[0] != size[1]) {
    mexPrintf("\n\n Inputs must have the same number of arrival times \n\n");
    //    cout << "\n\n Inputs must have the same number of arrival times \n\n";
    //    return 8;
  }
  */
  
  //find the maximum value in the input vector of spike times
  vector<double>::const_iterator largest_X0 =
    max_element( X[0].begin(), X[0].end() );
  vector<double>::const_iterator largest_X1 =
    max_element( X[1].begin(), X[1].end() );

  double max_both;
  if (*largest_X0 > *largest_X1)
    max_both = *largest_X0;
  else
    max_both = *largest_X1;


  //Create an NTE_Calculator object
  NTE_Calculator NTE_Calc;

  //bin the data (production version)
  const double binsize = 0.001;
  vector<int> spikes_binned[2];
  int spikes_in_bin[2];
  double edge;
  NTE_Calc.X_binned[0].clear();
  NTE_Calc.X_binned[1].clear();
  for (edge = 0.00; edge < max_both; edge = edge + binsize) {
    spikes_in_bin[0] = 0;
    spikes_in_bin[1] = 0;
    for (int i = 0; i < (int)X[0].size(); i++) {
      if ( (X[0][i] > edge) & (X[0][i] < edge + binsize) )
	spikes_in_bin[0]++;
      if ( (X[1][i] > edge) & (X[1][i] < edge + binsize) )
	spikes_in_bin[1]++;
    }
    spikes_binned[0].push_back(spikes_in_bin[0]);
    spikes_binned[1].push_back(spikes_in_bin[1]);
  }
  
  // Set the binned data vectors in the NTE_Calculator object
  NTE_Calc.set_Xbinned(spikes_binned[0], spikes_binned[1]);


  // Loop through the taubin sizes specified as the third input
  vector<double> myNTE;

  int tau;
  for (int i = 0; i < tau_bins.size(); i++) {
    //mexPrintf("Calculating NTE for tau: %d; ", tau_bins[i]);
    NTE_Calc.set_Xbinned(spikes_binned[0], spikes_binned[1]);
    tau = (int)tau_bins[i];
    NTE_Calc.set_tau_bins_P(tau);
    NTE_Calc.set_tau_bins_F(tau);
    myNTE.push_back(NTE_Calc.NTE());
    //mexPrintf("binsize=%f;  tauP = tauF = %d; NTE = %f \n\n", binsize, tau, myNTE[i]);
  }
  
  for (int q = 0; q < myNTE.size(); q++) {
    OutputArray[q] = myNTE[q]; // normalized transfer entropy
  }

}      


