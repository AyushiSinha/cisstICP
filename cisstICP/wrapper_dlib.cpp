// ****************************************************************************
//
//    Copyright (c) 2014, Seth Billings, Ayushi Sinha Russell Taylor, Johns Hopkins University.
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
#include "wrapper_dlib.h"
#include "algDirICP_GIMLOP.h"

#include <assert.h>
#undef NDEBUG       // enable debug in release mode

// define to print iteration details
//#define DLIB_VERBOSE


//--- Globals ---//

namespace
{
  // Global variables
  //  (referenced from global functions)
  algDirICP_GIMLOP *Kent_dlib = NULL;

  // Global functions
  //  (needed for function pointers)
  double fValue( const wrapper_dlib::dlib_vector &x_dlib )
  {
    static vctDynamicVector<double> x;
	x.SetSize(x_dlib.size());
	for (int i = 0; i < x_dlib.size(); i++)
		x[i] = x_dlib(i);
 
    return Kent_dlib->CostFunctionValue( x );
  }

  wrapper_dlib::dlib_vector fDerivative( const wrapper_dlib::dlib_vector &x_dlib )
  {
    static vctDynamicVector<double>  x;
	static vctDynamicVector<double>  g;
    static wrapper_dlib::dlib_vector  g_dlib;    // 7-element vector

	x.SetSize(x_dlib.size());
	g.SetSize(x_dlib.size());
	g_dlib.set_size(x_dlib.size());

	for (int i = 0; i < x_dlib.size(); i++)
		x[i] = x_dlib(i);

    Kent_dlib->CostFunctionGradient( x, g );
	for (int i = 0; i < x_dlib.size(); i++)
		g_dlib(i) = g[i];

    return g_dlib;
  }
} // local namespace



//--- Non-Globals ---//

// Constructor
wrapper_dlib::wrapper_dlib()  // algDirICP_GIMLOP *kent )
  : maxIter(40), //maxIter( 20 ),
  //tol_df( 1.0e-6 ),
  gradientNormThresh( 1.0e-3 )
{
  Kent_dlib = NULL;     // initialize global variable
  //Kent_dlib = kent;  // initialize global variable
}


// passing a pointer to the algorithm is necessary if this library is being
//  used with multiple Kent algorithms simultaneously (even if single threaded)
/*vct6*/ vctDynamicVector<double> wrapper_dlib::ComputeRegistration(const vctDynamicVector<double> &x0, algDirICP_GIMLOP *kent)
{
  // initialize global pointer to algorithm
  Kent_dlib = kent;

  int nComponents = x0.size();
  dlib_vector x_dlib(nComponents);  // 7-element vector
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

	  int nRot = 3;
	  int nTrans = 6;
	  int nScale;
	  if (Kent_dlib->bScale) nScale = 7; else nScale = 6;

	  for (int i = 0; i < nComponents; i++)
	  {
		  x_dlib(i) = x0[i];

		  if (i < nRot)
		  {
			  x_lower(i) = x_dlib(i) - Kent_dlib->rb;
			  x_upper(i) = x_dlib(i) + Kent_dlib->rb;
		  }
		  else if (i >= nRot && i < nTrans)
		  {
			  x_lower(i) = x_dlib(i) - Kent_dlib->tb;
			  x_upper(i) = x_dlib(i) + Kent_dlib->tb;
		  }
		  else if (i >= nTrans && i < nScale)
		  {
			  x_lower(i) = 1 - Kent_dlib->sb;
			  x_upper(i) = 1 + Kent_dlib->sb;
		  }
	  }

#ifdef DLIB_VERIFY_DERIVATIVE
    std::cout << "Difference between analytic derivative and numerical approximation of derivative: " 
		<< dlib::length(dlib::derivative(fValue)(x_dlib)-fDerivative(x_dlib)) << std::endl;
#endif

    dlib::find_min_box_constrained( 
		dlib::bfgs_search_strategy(),
#ifdef DLIB_VERBOSE
      dlib::gradient_norm_stop_strategy( gradientNormThresh,maxIter ).be_verbose(),
#else
      dlib::gradient_norm_stop_strategy( gradientNormThresh,maxIter ),
#endif
      fValue, fDerivative,
	  x_dlib, x_lower, x_upper);
  }
  catch (std::exception& e)
  {
    std::cout << "DLIB EXCEPTION: " << e.what() << std::endl;
    assert(0);
  }

  Kent_dlib = NULL;  
  
  for (int i = 0; i < nComponents; i++)
  {
	  x[i] = x_dlib(i);
  }

  return x;
}