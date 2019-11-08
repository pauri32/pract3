/// @file

#include <iostream>
#include <math.h>
#include <cmath>

#include "pitch_analyzer.h"

using namespace std;

/// Name space of UPC
namespace upc {
  void PitchAnalyzer::autocorrelation(const vector<float> &x, vector<float> &r) const {
    for (unsigned int l = 0; l < r.size(); ++l) {
  		/// \TODO Compute the autocorrelation r[l]
      for(unsigned int n = 0; n < x.size()-l; n++){
        r[l] += x[n]*x[n+l];
      }
      r[l] /= x.size();
    }
    if (r[0] == 0.0F){ //to avoid log() and divide zero
      r[0] = 1e-10;
    }
  }

  void PitchAnalyzer::set_window(Window win_type) {
    if (frameLen == 0)
      return;

    window.resize(frameLen);

    switch (win_type) {
    case HAMMING:
      /// \TODO Implement the Hamming window
      for(int n = 0; n < frameLen; n++){
        window[n] = 0.53836-0.46164*cos((2*M_PI*n)/(frameLen-1));
      }
      break;
    case RECT:
      for(int n = 0; n < frameLen; n++){
        window[n] = 1;
      }
    default:
      window.assign(frameLen, 1);
    }
  }

  void PitchAnalyzer::set_f0_range(float min_F0, float max_F0) {
    npitch_min = (unsigned int) samplingFreq/max_F0;
    if (npitch_min < 2)
      npitch_min = 2;  // samplingFreq/2

    npitch_max = 1 + (unsigned int) samplingFreq/min_F0;

    //frameLen should include at least 2*T0
    if (npitch_max > frameLen/2)
      npitch_max = frameLen/2;
  }

  bool PitchAnalyzer::unvoiced(float pot, float r1norm, float rmaxnorm) const {
    /// \TODO Implement a rule to decide whether the sound is voiced or not.
    /// * You can use the standard features (pot, r1norm, rmaxnorm),
    ///   or compute and use other ones.
    // cout << rmaxnorm << endl;
    cout << rmaxnorm << endl;
    if (pot > -55 && rmaxnorm > 100){
      //cout << rmaxnorm << endl;
      return false;
    }
    else{
      //cout << rmaxnorm << endl;
      return true;
    }
  }

  float PitchAnalyzer::compute_pitch(vector<float> & x) const {
    if (x.size() != frameLen)
      return -1.0F;

    //Window input frame
    for (unsigned int i=0; i<x.size(); ++i)
      x[i] *= window[i];

    vector<float> r(npitch_max);

    //Compute correlation
    float maxframe = *max_element(x.begin(),x.end());
    for (int i = 0; i < x.size(); i++) {
      if(abs(x[i])/maxframe < 0.1){
        x[i] = 0;
      }
    }
    autocorrelation(x, r);

    vector<float>::const_iterator iR = r.begin(), iRMax = iR;

    /// \TODO
	/// Find the lag of the maximum value of the autocorrelation away from the origin.<br>
	/// Choices to set the minimum value of the lag are:
	///    - The first negative value of the autocorrelation.
	///    - The lag corresponding to the maximum value of the pitch.
    ///	   .
	/// In either case, the lag should not exceed that of the minimum value of the pitch.
    unsigned int lag = iRMax-r.begin();
    int maxpos = 0;
    int max = samplingFreq/50;
    int min = samplingFreq/500;
    float maxval = 0;
    advance(iR,min);
    lag = max_element(r.begin()+min, r.begin()+max)-r.begin();
    for(int i = min; i < max-min; i++){
      if(maxval < x[i]) maxpos = i;
    }
    // cout << samplingFreq/maxpos << endl;
/*    for (unsigned int i = 0; i < x.size(); i++) {
      if(cross == false){
        if(r[i]<0)  cross = true;
      }
      else{
        if(maxpeak < r[i]){
          lag = lag + i;
          maxpeak = r[i];
        }
      }
    } */
    float pot = 10 * log10(r[0]);

    //You can print these (and other) features, look at them using wavesurfer
    //Based on that, implement a rule for unvoiced
    //change to #if 1 and compile
#if 0
    if (r[0] > 0.0F)
      cout << pot << '\t' << r[1]/r[0] << '\t' << r[lag]/r[0] << endl;
#endif

    if (unvoiced(pot, r[1]/r[0], r[lag]/r[0]))
      return 0;
    else
      return (float)samplingFreq/(float)lag;
  }
}
