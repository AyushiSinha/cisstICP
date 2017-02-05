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
#include "algICP_DIMLP_dlibWrapper.h"
#include "algICP_DIMLP.h"

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
	algICP_DIMLP *alg = NULL;

  // Global functions
  //  (needed for function pointers)
  double fValue(const algICP_DIMLP_dlibWrapper::dlib_vector &x_dlib)
  {
    //static vct7 x;
	  static vctDynamicVector<double> x;
	  x.SetSize(x_dlib.size());
   // x.Assign(x_dlib(0), x_dlib(1), x_dlib(2),
   //   x_dlib(3), x_dlib(4), x_dlib(5),
	  //x_dlib(6));

	  for (int i = 0; i < x_dlib.size(); i++)
		  x[i] = x_dlib(i);

    return alg->CostFunctionValue(x);
  }

  algICP_DIMLP_dlibWrapper::dlib_vector fDerivative(
	  const algICP_DIMLP_dlibWrapper::dlib_vector &x_dlib)
  {
    //static vct7   x;
    //static vct7   g;
	  static vctDynamicVector<double> x;
	  static vctDynamicVector<double> g;
	//static algICP_DIMLP_dlibWrapper::dlib_vector  g_dlib(7);  // 7-element vector
	  static algICP_DIMLP_dlibWrapper::dlib_vector  g_dlib;
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
		g_dlib(i) = g[i];

    return g_dlib;
  }
} // local namespace


//--- Non-Globals ---//

// Constructor
algICP_DIMLP_dlibWrapper::algICP_DIMLP_dlibWrapper(algICP_DIMLP *argAlg)
  : maxIter( 20 ),
  //tol_df( 1.0e-6 ),
  gradientNormThresh( 1.0e-3 )
{
  alg = argAlg;  // initialize global reference to algorithm object
}


//vct7 algICP_DIMLP_dlibWrapper::ComputeRegistration(const vct7 &x0)
vctDynamicVector<double> algICP_DIMLP_dlibWrapper::ComputeRegistration(const vctDynamicVector<double> &x0)
{
  //dlib_vector x_dlib(7);  // 7-element vector
	int nComponents = x0.size();
	dlib_vector x_dlib(nComponents); // 6+n_modes-element vector
	vctDynamicVector<double> x;
	x.SetSize(nComponents);

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
	  }

#if 0 //#ifdef DLIB_VERIFY_DERIVATIVE
    std::cout << "Difference between analytic derivative and numerical approximation of derivative: \n" 
		<< dlib::derivative(fValue)(x_dlib) << " - " << fDerivative(x_dlib) << " = "
      << dlib::derivative(fValue)(x_dlib) - fDerivative(x_dlib) << std::endl;
#endif


#ifdef DLIB_VERBOSE
    dlib::find_min(
      dlib::bfgs_search_strategy(),
      //dlib::objective_delta_stop_strategy( Tol_df,maxIter ).be_verbose(),
      dlib::gradient_norm_stop_strategy(gradientNormThresh, maxIter).be_verbose(),
      fValue, fDerivative,
      x_dlib, -1.0);
#else
    dlib::find_min( 
      dlib::bfgs_search_strategy(),
      //dlib::objective_delta_stop_strategy( Tol_df,maxIter ),
      dlib::gradient_norm_stop_strategy(gradientNormThresh, maxIter),
	  fValue, dlib::derivative(fValue),//fDerivative, 
	  x_dlib, -1.0);
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
	  x[i] = x_dlib(i);
  }

  return x;
}