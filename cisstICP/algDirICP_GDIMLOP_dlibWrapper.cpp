// ****************************************************************************
//
//    Copyright (c) 2014, Seth Billings, Russell Taylor, Johns Hopkins University
//    All rights reserved.
//
//    Redistribution and use in source and binary forms, with or without
//    modification, are permitted provided that the following conditions are
//    met:
//
//    1. Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//
//    2. Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//
//    3. Neither the name of the copyright holder nor the names of its
//    contributors may be used to endorse or promote products derived from
//    this software without specific prior written permission.
//
//    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
//    "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
//    LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
//    A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
//    HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
//    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
//    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
//    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
//    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
//    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
//    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//  
// ****************************************************************************
#include "algDirICP_GDIMLOP_dlibWrapper.h"
#include "algDirICP_GDIMLOP.h"

#include <assert.h>
#undef NDEBUG       // enable debug in release mode

// define to print iteration details
// #define DLIB_VERBOSE


//--- Globals ---//

namespace
{
  // Global variables
  //  (referenced from global functions)
  algDirICP_GDIMLOP *Kent_dlib = NULL;

  // Global functions
  //  (needed for function pointers)
  double fValue( const wrapper_dlib::dlib_vector &x_dlib )
  {
    //static vct6 x;
	static vctDynamicVector<double> x;
	x.SetSize(x_dlib.size());
    //x.Assign( x_dlib(0), x_dlib(1), x_dlib(2), 
    //          x_dlib(3), x_dlib(4), x_dlib(5) );
	for (int i = 0; i < x_dlib.size(); i++)
		/*if (i < 3)
		x[i] = 0;
		else*/
		x[i] = x_dlib(i);
 
    return Kent_dlib->CostFunctionValue( x );
  }

  algDirICP_GDIMLOP_dlibWrapper::dlib_vector fDerivative(const algDirICP_GDIMLOP_dlibWrapper::dlib_vector &x_dlib)
  {
 //   static vct6   x;
	//static vct6   g;
	static vctDynamicVector<double>   x;
	static vctDynamicVector<double>   g;
	static algDirICP_GDIMLOP_dlibWrapper::dlib_vector  g_dlib/*(6)*/;    // 6-element vector

	x.SetSize(x_dlib.size());
	g.SetSize(x_dlib.size());
	g_dlib.set_size(x_dlib.size());

    //x.Assign( x_dlib(0), x_dlib(1), x_dlib(2), 
    //          x_dlib(3), x_dlib(4), x_dlib(5) );
	for (int i = 0; i < x_dlib.size(); i++)
		x[i] = x_dlib(i);

    Kent_dlib->CostFunctionGradient( x, g );
    //g_dlib(0) = g[0];
    //g_dlib(1) = g[1];
    //g_dlib(2) = g[2];
    //g_dlib(3) = g[3];
    //g_dlib(4) = g[4];
    //g_dlib(5) = g[5];
	for (int i = 0; i < x_dlib.size(); i++)
		/*if (i < 3)
		g_dlib(i) = 0;
		else*/
		g_dlib(i) = g[i];

    return g_dlib;
  }
} // local namespace



//--- Non-Globals ---//

// Constructor
algDirICP_GDIMLOP_dlibWrapper::algDirICP_GDIMLOP_dlibWrapper(algDirICP_GDIMLOP *argAlg)  // algDirICP_GIMLOP *kent )
  : maxIter( 20 ), //maxIter(40),
  //tol_df( 1.0e-6 ),
  gradientNormThresh( 1.0e-3 )
{
	Kent_dlib = argAlg; // NULL;     // initialize global variable
  //Kent_dlib = kent;  // initialize global variable
}


// passing a pointer to the algorithm is necessary if this library is being
//  used with multiple Kent algorithms simultaneously (even if single threaded)
/*vct6*/vctDynamicVector<double> algDirICP_GDIMLOP_dlibWrapper::ComputeRegistration(const /*vct6*/vctDynamicVector<double> &x0, algDirICP_GDIMLOP *kent)
{
  // initialize global pointer to algorithm
  Kent_dlib = kent;

  int nComponents = x0.size();
  dlib_vector x_dlib(nComponents);	// n_trans+n_modes-element vector
  dlib_vector x_lower(nComponents);
  dlib_vector x_upper(nComponents);

  vctDynamicVector<double> x;
  x.SetSize(nComponents);

  //double f0 = Kent_gsl->CostFunctionValue( x0 );
  //double thresh_df = tol_df * f0;
  //double thresh_df = tol_df * Kent_gsl->k_sum;
  //std::cout << "thresh_df: " << thresh_df << std::endl;

  try
  {
    // Now we use the find_min() function to find the minimum point.  The first argument
    // to this routine is the search strategy we want to use.  The second argument is the 
    // stopping strategy.
    //   objective_delta_stop_strategy:  stop if df < threshold

    // The other arguments to find_min() are the function to be minimized, its derivative, 
    // then the starting point, and the last is an acceptable minimum value.  
    // That is, if the algorithm finds any inputs that give an output value less than
    // this then it will stop immediately.  Usually you supply a number smaller than 
    // the actual global minimum.

    //x_dlib(0) = x0[0];
    //x_dlib(1) = x0[1];
    //x_dlib(2) = x0[2];
    //x_dlib(3) = x0[3];
    //x_dlib(4) = x0[4];
    //x_dlib(5) = x0[5];

	int nTrans;
	if (Kent_dlib->bScale)
		nTrans = 7;
	else
		nTrans = 6;

	for (int i = 0; i < nComponents; i++)
	{
		x_dlib(i) = x0[i];

		if (i < nTrans)
		{
			x_lower(i) = -DBL_MAX;
			x_upper(i) = DBL_MAX;

			// should we constrain the scale as well?

			//if (i == 6) {
			// x_lower(i) = 0.9;
			// x_upper(i) = 1.1;
			//}
		}
		else
		{
			x_lower(i) = -Kent_dlib->spb;
			x_upper(i) = Kent_dlib->spb;
		}
	}

    //std::cout << "Difference between analytic derivative and numerical approximation of derivative: " 
    //  << dlib::length(dlib::derivative(fValue)(x_dlib) - fDerivative(x_dlib)) << std::endl;


    dlib::find_min_box_constrained( 
	  dlib::bfgs_search_strategy(),
#ifdef DLIB_VERBOSE
      //dlib::objective_delta_stop_strategy( Tol_df,maxIter ).be_verbose(),
      dlib::gradient_norm_stop_strategy( gradientNormThresh,maxIter ).be_verbose(),
#else
      //dlib::objective_delta_stop_strategy( Tol_df,maxIter ),
      dlib::gradient_norm_stop_strategy( gradientNormThresh, maxIter ),
#endif
	  fValue, fDerivative, //dlib::derivative(fValue), 
      x_dlib, x_lower, x_upper);

    //find_min_using_approximate_derivatives( bfgs_search_strategy(),
    //                                        objective_delta_stop_strategy(Tol_df).be_verbose(),
    //                                        dlib_fValue, 
    //                                        x0_dlib, -1.0);
  }
  catch (std::exception& e)
  {
    std::cout << "DLIB EXCEPTION: " << e.what() << std::endl;
    assert(0);
  }

  Kent_dlib = NULL;

  //return vct6( x_dlib(0), x_dlib(1), x_dlib(2), 
  //             x_dlib(3), x_dlib(4), x_dlib(5) );

  for (int i = 0; i < nComponents; i++)
  {
	  /*if (i < 3)
	  x[i] = 0;
	  else
	  */
	  x[i] = x_dlib(i);
	  //if (i > 5) {
	  // x[i] = std::min(x[i], 3.0);
	  // x[i] = std::max(x[i], -3.0);
	  //}
  }

  return x;
}