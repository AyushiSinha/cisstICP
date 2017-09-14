// ****************************************************************************
//
//    Copyright (c) 2014, Ayushi Sinha, Seth Billings, Russell Taylor, 
//	  Johns Hopkins University. 
//	  All rights reserved.
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
#include "algDirICP_DIMLOP_dlibWrapper.h"
#include "algDirICP_DIMLOP.h"

#include <assert.h>
#undef NDEBUG       // enable debug in release mode

// uncomment to activate debug modes
//#define DLIB_VERBOSE    
//#define DLIB_VERIFY_DERIVATIVE

//--- Globals ---//
namespace
{
  // Global variables
  //  (referenced from global functions)
  algDirICP_DIMLOP *alg = NULL;

  // Global functions
  //  (needed for function pointers)
  double fValue(const algDirICP_DIMLOP_dlibWrapper::dlib_vector &x_dlib)
  {
    //static vct7 x;
	  static vctDynamicVector<double> x;
	  x.SetSize(x_dlib.size());
   // x.Assign(x_dlib(0), x_dlib(1), x_dlib(2),
   //   x_dlib(3), x_dlib(4), x_dlib(5),
	  //x_dlib(6));

	  for (int i = 0; i < x_dlib.size(); i++)
		  /*if (i < 3)
			  x[i] = 0;
		  else*/
			  x[i] = x_dlib(i);

    return alg->CostFunctionValue(x);
  }

  algDirICP_DIMLOP_dlibWrapper::dlib_vector fDerivative(
	  const algDirICP_DIMLOP_dlibWrapper::dlib_vector &x_dlib)
  {
    //static vct7   x;
    //static vct7   g;
	  static vctDynamicVector<double> x;
	  static vctDynamicVector<double> g;
	//static algICP_DIMLP_dlibWrapper::dlib_vector  g_dlib(7);  // 7-element vector
	  static algDirICP_DIMLOP_dlibWrapper::dlib_vector  g_dlib;
	  x.SetSize(x_dlib.size());
	  g.SetSize(x_dlib.size());
	  g_dlib.set_size(x_dlib.size());

  //  x.Assign(x_dlib(0), x_dlib(1), x_dlib(2),
		//x_dlib(3), x_dlib(4), x_dlib(5),
		//x_dlib(6));
	  for (int i = 0; i < x_dlib.size(); i++)
		  x[i] = x_dlib(i);

    alg->CostFunctionGradient(x, g);
 //   g_dlib(0) = g[0];
 //   g_dlib(1) = g[1];
 //   g_dlib(2) = g[2];
 //   g_dlib(3) = g[3];
 //   g_dlib(4) = g[4];
	//g_dlib(5) = g[5];
	//g_dlib(6) = g[6];

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
algDirICP_DIMLOP_dlibWrapper::algDirICP_DIMLOP_dlibWrapper(algDirICP_DIMLOP *argAlg)
  : maxIter( 20 ),
  //tol_df( 1.0e-6 ),
  gradientNormThresh( 1.0e-3 )
{
  alg = argAlg;  // initialize global reference to algorithm object
}


//vct7 algICP_DIMLP_dlibWrapper::ComputeRegistration(const vct7 &x0)
vctDynamicVector<double> algDirICP_DIMLOP_dlibWrapper::ComputeRegistration(const vctDynamicVector<double> &x0)
{
  //dlib_vector x_dlib(7);  // 7-element vector
	int nComponents = x0.size();

	dlib_vector x_dlib(nComponents); // n_trans+n_modes-element vector
	dlib_vector x_lower(nComponents);
	dlib_vector x_upper(nComponents);

	vctDynamicVector<double> x;
	x.SetSize(nComponents);

  try
  {
    // Now we use the find_min_box_constrained() function to find the minimum point.  The first argument
    // to this routine is the search strategy we want to use.  The second argument is the 
    // stopping strategy:
    //   objective_delta_stop_strategy:  stop if df < threshold

    // The other arguments to find_min_box_constrained() are the function to be minimized, its derivative, 
    // then the starting point, and the last two are the lower and upper bounds of the constraints on the
	// parameters being optimized:
	//	 transformation paramters:	unconstrained
	//	 shape paramters:			constrained between -3.0 and 3.0 (default) or user specified limits

 //   x_dlib(0) = x0[0];
 //   x_dlib(1) = x0[1];
 //   x_dlib(2) = x0[2];
 //   x_dlib(3) = x0[3];
 //   x_dlib(4) = x0[4];
	//x_dlib(5) = x0[5];
	//x_dlib(6) = x0[6];

	  for (int i = 0; i < nComponents; i++)
	  {
		  x_dlib(i) = x0[i];

		  int nTrans;
		  if (alg->bScale)
			  nTrans = 7;
		  else
			  nTrans = 6;

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
			  x_lower(i) = -alg->spb;
			  x_upper(i) = alg->spb;
		  }
	  }

#ifdef DLIB_VERIFY_DERIVATIVE
    std::cout << "Difference between analytic derivative and numerical approximation of derivative: \n" 
		<< dlib::derivative(fValue)(x_dlib) << " - \n" << fDerivative(x_dlib) << " = \n"
      << dlib::derivative(fValue)(x_dlib) - fDerivative(x_dlib) << std::endl;
#endif


#ifdef DLIB_VERBOSE
    //dlib::find_min(
    //  dlib::bfgs_search_strategy(),
    //  //dlib::objective_delta_stop_strategy( tol_df,maxIter ).be_verbose(),
    //  dlib::gradient_norm_stop_strategy(gradientNormThresh, maxIter).be_verbose(),
    //  fValue, dlib::derivative(fValue), //fDerivative,
    //  x_dlib, -1.0);

	dlib::find_min_box_constrained(
		dlib::bfgs_search_strategy(),
		//dlib::objective_delta_stop_strategy( tol_df,maxIter ).be_verbose(),
		dlib::gradient_norm_stop_strategy(gradientNormThresh, maxIter).be_verbose(),
		fValue, dlib::derivative(fValue), //fDerivative,
		x_dlib, x_lower, x_upper);
#else 
	//dlib::find_min(
	//  dlib::bfgs_search_strategy(),
	//  //dlib::objective_delta_stop_strategy( tol_df,maxIter ),
	//  dlib::gradient_norm_stop_strategy(gradientNormThresh, maxIter),
	//  fValue, fDerivative, //dlib::derivative(fValue),
	//  x_dlib, -1.0);

	dlib::find_min_box_constrained(
		dlib::bfgs_search_strategy(),
		//dlib::objective_delta_stop_strategy( tol_df,maxIter ),
		dlib::gradient_norm_stop_strategy(gradientNormThresh, maxIter),
		fValue, fDerivative, //dlib::derivative(fValue),
		x_dlib, x_lower, x_upper);
#endif

    //dlib::find_min_using_approximate_derivatives(
    //  dlib::bfgs_search_strategy(),
    //  dlib::objective_delta_stop_strategy(Tol_df).be_verbose(),
    //  dlib_fValue, 
    //  x0_dlib, -1.0);
  }
  catch (std::exception& e)
  {
    std::cout << "DLIB EXCEPTION: " << e.what() << std::endl;
    assert(0);
  }

  //return vct7( x_dlib(0), x_dlib(1), x_dlib(2), 
		//	   x_dlib(3), x_dlib(4), x_dlib(5),
		//	   x_dlib(6));

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