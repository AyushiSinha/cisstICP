// ****************************************************************************
//
//    Copyright (c) 2014, Seth Billings, Ayushi Sinha, Russell Taylor, Johns Hopkins University.
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

#include <algorithm>
#include <limits>

#include "algDirICP_GIMLOP.h"
#include "DirPDTreeNode.h"
#include "utilities.h"

#include "mins.h"   // Numerical Recipes

#include <assert.h>
#undef NDEBUG       // enable assert in release mode

#define EPS  1e-14
#define ENABLE_PARALLELIZATION

// Constructor
algDirICP_GIMLOP::algDirICP_GIMLOP(
  DirPDTreeBase *pDirTree,
  vctDynamicVector<vct3> &samplePts,
  vctDynamicVector<vct3> &sampleNorms,
  const vctDynamicVector<double> &argK,
  const vctDynamicVector<double> &argE,
  const vctDynamicVector<vct3x2> &argL,
  const vctDynamicVector<vct3x3> &argM,
  double scale, bool bScale,
  PARAM_EST_TYPE paramEst)
  : algDirICP(pDirTree, samplePts, sampleNorms),
  algDirPDTree(pDirTree),
  dlib(),
  paramEstMethod(paramEst)
{
  SetSamples(samplePts, sampleNorms, argK, argE, argL, argM, scale, bScale, paramEst);
  //SetNoiseModel(argK, argE, argL, argM, paramEst);

#ifdef SAVE_MATCHES
  MatchIsotropic = false;
#endif

#ifdef ENABLE_DEBUG_FILE
  std::string debugFile = "KentDebug.txt";
  debugStream.open( debugFile.c_str() );
#endif
}

void algDirICP_GIMLOP::SetConstraints(double argRotbounds,
	double argTransbounds,
	double argScalebounds)
{
	rb = argRotbounds;
	tb = argTransbounds;
	sb = argScalebounds;
}

void algDirICP_GIMLOP::ComputeMatchStatistics(double &Avg, double &StdDev) //gotta do this right with mahalanobis distance
{
	double sumSqrMatchDist = 0.0;
	double sumMatchDist = 0.0;
	double sqrMahalDist;
	double sqrMatchDist;
	double matchAngle;
	double axisAngle1;
	double axisAngle2;

	totalSumSqrMahalDist = 0.0;
	sumSqrMahalDist = 0.0;
	double sumMahalDist = 0.0;
	totalSumSqrMatchAngle = 0.0;
	sumSqrMatchAngle = 0.0;

	nGoodSamples = 0;

	vct3 residual;

	// NOTE: if using a method with outlier rejection, it may be desirable to
	//       compute statistics on only the inliers
	for (unsigned int i = 0; i < nSamples; i++)
	{
		residual = matchPts[i] - (Freg * samplePts[i]) * sc;

		sqrMahalDist = residual*invM[i] * residual;
		totalSumSqrMahalDist += sqrMahalDist;

		matchAngle = acos(std::fmod(matchNorms[i] * (Freg.Rotation() * sampleNorms[i]), 2 * cmnPI));
		axisAngle1 = asin(std::fmod(R_L[i].Column(0) * matchNorms[i], 2 * cmnPI));
		axisAngle2 = asin(std::fmod(R_L[i].Column(1) * matchNorms[i], 2 * cmnPI));

		totalSumSqrMatchAngle += (k[i] * matchAngle * matchAngle) + ((k[i] - 2 * B[i]) * axisAngle1 * axisAngle1) + ((k[i] + 2 * B[i]) * axisAngle2 * axisAngle2);

		if (outlierFlags[i])	continue;	// skip outliers

		sumSqrMahalDist += sqrMahalDist;
		sumMahalDist += sqrt(sqrMahalDist);

		sqrMatchDist = residual.NormSquare();
		sumSqrMatchDist += sqrMatchDist;
		sumMatchDist += sqrt(sqrMatchDist);

		sumSqrMatchAngle += (k[i] * matchAngle * matchAngle) + ( (k[i] - 2*B[i]) * axisAngle1 * axisAngle1 ) + ( (k[i] + 2*B[i]) * axisAngle2 * axisAngle2 ) ;
		nGoodSamples++;
	}
	Avg = sumMahalDist / nGoodSamples;
	StdDev = sqrt((sumSqrMahalDist / nGoodSamples) + Avg*Avg);
}

void algDirICP_GIMLOP::PrintMatchStatistics(std::stringstream &tMsg)
{
	// For registration rejection purpose:
	tMsg << "\nSum square mahalanobis distance = " << totalSumSqrMahalDist << " over " << nSamples << " samples";
	tMsg << "\nSum square match angle = " << totalSumSqrMatchAngle << " over " << nSamples << " samples";
	tMsg << "\nSum square mahalanobis distance = " << sumSqrMahalDist << " over " << nGoodSamples << " inliers";
	tMsg << "\nSum square match angle = " << sumSqrMatchAngle  << " over " << nGoodSamples << " inliers\n";
}

double algDirICP_GIMLOP::ICP_EvaluateErrorFunction()
{
  // Return the negative log likelihood of the Kent
  //  and Gaussian distributions under the assumption
  //   of independence between the two distributions
  //
  //   Negative Log-Likelihood:
  //    -log[ C * exp( ... ) ]

  // Just return match error for now, shifted to be lower bounded at 0
  double Error = 0.0;
  for (unsigned int i = 0; i < nSamples; i++)
  {
	if (outlierFlags[i])	continue;	// skip outliers

    Error += MatchError(samplePtsXfmd[i], sampleNormsXfmd[i],
      matchPts[i], matchNorms[i],
      k[i], B[i], R_L[i], R_invM_Rt[i]);
  }

  return Error;
}

vctFrm3 algDirICP_GIMLOP::ICP_RegisterMatches()
{
  vctFrm3 dF;
  vctDynamicVector<double> x0;
  vctDynamicVector<double> x;

  x0.SetSize(nTrans);
  x.SetSize(nTrans);

  x0.SetAll(0.0);
  // x_prev must begin at a different value than x0
  x_prev.SetAll(std::numeric_limits<double>::max());

  //x = gsl.ComputeRegistration( x0 );
  x = dlib.ComputeRegistration(x0, this);

  // update transform
  vctFixedSizeVectorRef<double, 3, 1> alpha(x, 0);
  vctFixedSizeVectorRef<double, 3, 1> t(x, 3);
  double scale;
  if (bScale)
	  scale = x[6];
  dF.Rotation() = vctRot3(vctRodRot3(alpha));
  dF.Translation() = t;
  Freg = dF * Freg;
  if (bScale)
	  sc = scale;

  return Freg;
}

void algDirICP_GIMLOP::ICP_InitializeParameters(vctFrm3 &FGuess)
{
  // initialize base class
  algDirICP::ICP_InitializeParameters(FGuess);

  bFirstIter_Matches = true;
  nOutliers = 0;
  sigma2 = 0.0;

  // monitoring variables
  errFuncNormWeight = 0.0;
  errFuncPosWeight = 0.0;

  if (k.size() != nSamples || B.size() != nSamples || L.size() != nSamples ||
    R_L.size() != nSamples || invM.size() != nSamples || N.size() != nSamples || invN.size() != nSamples ||
    R_invM_Rt.size() != nSamples || N_Rt.size() != nSamples || inv_N_Rt.size() != nSamples ||
    Yp_t.size() != nSamples || Rat_Yp_RaXp_t.size() != nSamples || 
    invM_Rat_Yp_RaXp_t.size() != nSamples ||
    RaXn.size() != nSamples || RaRL.size() != nSamples || M.size() != nSamples ||
    k_msmt.size() != nSamples || E_msmt.size() != nSamples ||
    M_msmt.size() != nSamples || invM_msmt.size() != nSamples || N_msmt.size() != nSamples ||
    invN_msmt.size() != nSamples || Dmin_msmt.size() != nSamples || Emin_msmt.size() != nSamples)
  {
    std::cout << "ERROR: noise parameters not initialized for Kent algorithm" << std::endl;
    assert(0);
  }
  
  R_invM_Rt.SetAll(vct3x3::Eye());
  R_MsmtM_Rt.SetAll(vct3x3(0.0));
  outlierFlags.SetAll(0);

  // Initialize noise-model values
  switch (paramEstMethod)
  {
  case PARAMS_FIXED:
  {
    k = k_msmt;
    M = M_msmt;
    invM = invM_msmt;
    N = N_msmt;
    invN = invN_msmt;
    Dmin = Dmin_msmt;
    Emin = Emin_msmt;
    for (unsigned int i = 0; i < nSamples; i++)
    {
      // B = E * k/2
      B[i] = E_msmt[i] * k_msmt[i] / 2.0;
    }
    meanSigma2 = meanTraceM_msmt / 3.0;
    meanK = meanK_msmt;
    meanE = meanE_msmt;
    break;
  }
  case PARAMS_DYN_THRESH:
  {
    // use isotropic positions-only noise model for first iteration
    //  a. allows correcting for macro translational misalignment
    //  b. prevents giving too much weight to orientations in the first
    //     match when we don't yet know the match uncertainty parameters
    //  NOTE: could do an initial position-only match to estimate initial
    //        match uncertainty and then do a follow-up position-orientation match;
    //        how would this behave with large initial offset in translation?
    k.SetAll(0.0);
    B.SetAll(0.0);
    M.SetAll(vct3x3::Eye());
    invM.SetAll(vct3x3::Eye());
    N.SetAll(vct3x3::Eye());
    invN.SetAll(vct3x3::Eye());
    Dmin.SetAll(1.0);
    Emin.SetAll(1.0);
    meanSigma2 = 1.0;
    meanK = 0.0;
    meanE = 0.0;
    break;
  }
  case PARAMS_FIXED_1stITER_POS:
  {
    // use isotropic positions-only noise model for first iteration
    //  a. allows correcting for macro translational misalignment
    //  b. prevents giving too much weight to orientations in the first
    //     match when we don't yet know the match uncertainty parameters
    k.SetAll(0.0);
    B.SetAll(0.0);
    M.SetAll(vct3x3::Eye());
    invM.SetAll(vct3x3::Eye());
    N.SetAll(vct3x3::Eye());
    invN.SetAll(vct3x3::Eye());
    Dmin.SetAll(1.0);
    Emin.SetAll(1.0);
    meanSigma2 = 1.0;
    meanK = 0.0;
    meanE = 0.0;
    break;
  }
  default:
  {
    std::cout << std::endl
      << "ERROR: unrecognized parameter estimation setting!" << std::endl;
    assert(0);
  }
  }

  // initialize rotation dependent parameters to prepare for initial match
  UpdateNoiseModel_SamplesXfmd(FGuess);
}

void algDirICP_GIMLOP::ICP_UpdateParameters_PostMatch() // CHECK IF YOU NEED NORMAL INFO HERE - should not for sigma because that's from positional data
{
	// base class
	algDirICP::ICP_UpdateParameters_PostMatch();

	// compute sum of square distance of inliers
	sumSqrDist_Inliers = 0.0;
	double sumNormProducts_Inliers = 0.0;
	for (unsigned int s = 0; s < nSamples; s++)
	{
		residuals_PostMatch[s] = samplePtsXfmd[s] - matchPts[s];
		sqrDist_PostMatch[s] = residuals_PostMatch[s].NormSquare();

		if (outlierFlags[s])	continue;	// skip outliers

		sumSqrDist_Inliers += sqrDist_PostMatch[s];

		sumNormProducts_Inliers +=
			vctDotProduct(sampleNormsXfmd.Element(s), matchNorms.Element(s));
	}

	// update the match uncertainty factor
	sigma2 = sumSqrDist_Inliers / (nSamples - nOutliers);

	// apply max threshold
	if (sigma2 > sigma2Max)
		sigma2 = sigma2Max;

	//// update noise models of hte matches
	//for (unsigned int s = 0; s < nSamples; s++)
	//{
	//	// update target covariances
	//	Myi[s] = pDirTree->DatumCovPtr(matchDatums[s]);	// use pointer here for efficiency

	//	// target covariance with match uncertainty		// TODO
	//	//Myi_sigma2[s] = *Myi[s];	
	//	//Myi_sigma2[s].Element(0, 0) += sigma2;
	//	//Myi_sigma2[s].Element(1, 1) += sigma2;
	//	//Myi_sigma2[s].Element(2, 2) += sigma2;
	//}

	// update noise model
	UpdateNoiseModel(sumSqrDist_Inliers, sumNormProducts_Inliers);

	if (bFirstIter_Matches)
	{
		UpdateNoiseModel_SamplesXfmd(FGuess);
	}

	vctRot3 R(FGuess.Rotation());
	for (unsigned int s = 0; s < nSamples; s++)
	{
		R_MsmtM_Rt[s] = R*M_msmt[s] * R.Transpose();
	}

	bFirstIter_Matches = false;
}

void algDirICP_GIMLOP::ICP_UpdateParameters_PostRegister(vctFrm3 &Freg)
{
  // base class
	algDirICP::ICP_UpdateParameters_PostRegister(Freg);

	if (bScale)
		for (unsigned int s = 0; s < nSamples; s++)
			samplePtsXfmd.Element(s) = sc * samplePtsXfmd.Element(s); // move this to IMLOP later

  UpdateNoiseModel_DynamicEstimates();
  UpdateNoiseModel_SamplesXfmd(Freg);  
}

unsigned int algDirICP_GIMLOP::ICP_FilterMatches()
{
	//return 0;
	//
	// Filer Matches for Outliers
	//  
	// The Square Mahalanobis Distance of the matches follow a chi-square distribution
	//  with 3 degrees of freedom (1 DOF for each spatial dimension).
	//
	//  Detect outliers as:  Square Mahalanobis Distance > ChiSquare(c)
	//
	//  Note:  ChiSquare(0.95) = 7.81     (1.96 Std Dev)
	//         ChiSquare(0.975) = 9.35    (2.24 Std Dev)
	//         ChiSquare(0.99) = 11.34    (2.56 Std Dev)
	//         ChiSquare(0.9973) = 14.16  (3.0 Std Dev)     MATLAB: chi2inv(0.9973,3)
	//
	//  When an outlier is identified, increase the variance of its noise model
	//  such that residual for that match is considered to be w/in 1 standard 
	//  deviation of its mean. This will reduce the impact of this match error
	//  on the registration result.
	//

	double StdDevExpansionFactor = 3.0;    // std dev expansion factor
	double varExpansionFactor = StdDevExpansionFactor * StdDevExpansionFactor;

	double ThetaThresh = StdDevExpansionFactor * circSD;
	ThetaThresh = ThetaThresh > cmnPI_2 ? cmnPI_2 : ThetaThresh;
	double NormProductThresh = cos(ThetaThresh);

	nOutliers = 0;
	vct3x3 Mo, inv_Mo;
	double sqrMahalDist = 0.0;
	double normProduct = 0.0;

	for (unsigned int s = 0; s < nSamples; s++)
	{
		// compute outlier noise model based on mearurment noise and sigma2 only
		//  and not including the surface model covariance
		//
		// Note: the code below assumes that the covariance model of the target
		//       shape is comprised of only a surface model covariance with zero
		//       measurement noise; if this is not true, then the target measurement
		//       noise should be added to the outlier covariance test below as well
		//   
		//Mo = R_M_Rt.Element(s);
		//Mo.Element(0, 0) += sigma2;
		//Mo.Element(1, 1) += sigma2;
		//Mo.Element(2, 2) += sigma2;
		//Mo = R_Mxi_Rt.Element(s) + Myi_sigma2.Element(s);
		//Mo.Element(0, 0) += outlier_alpha;
		//Mo.Element(1, 1) += outlier_alpha;
		//Mo.Element(2, 2) += outlier_alpha;

		// compute Mahalanobis distance
		//ComputeCovInverse_NonIter(Mo, inv_Mo);
		sqrMahalDist = residuals_PostMatch.Element(s)*R_invM_Rt[s]*residuals_PostMatch.Element(s);
		normProduct = vctDotProduct(sampleNormsXfmd.Element(s), matchNorms.Element(s));

		// check if outlier
		if (sqrMahalDist > ChiSquareThresh)
		{ // an outlier
			nOutliers++;
			outlierFlags[s] = 1;
			// add isotropic outlier term to noise model for this match
			//  with magnitude of half the square match distance
			//double outlierScale = 0.5 * sqrDist_PostMatch.Element(s) * varExpansionFactor;
			//Myi_sigma2[s].Element(0, 0) += outlierScale;
			//Myi_sigma2[s].Element(1, 1) += outlierScale;
			//Myi_sigma2[s].Element(2, 2) += outlierScale;
			//R_M_Rt[s].Element(0, 0) += outlierScale;
			//R_M_Rt[s].Element(1, 1) += outlierScale;
			//R_M_Rt[s].Element(2, 2) += outlierScale;
		}
		else if (normProduct < NormProductThresh)
		{
			nOutliers++;
			outlierFlags[s] = 1;
		}
		else
		{
			outlierFlags[s] = 0;
		}
	}

	return nOutliers;
}

void algDirICP_GIMLOP::UpdateNoiseModel(double sumSqrDist, double sumNormProducts)
{
	// angular error of normal orientations
	ComputeCircErrorStatistics(sumNormProducts, Rnorm, circSD);

#ifdef TEST_STD_ICP
	k = 0.0;
#endif
}

void algDirICP_GIMLOP::UpdateNoiseModel_DynamicEstimates()
{
  switch (paramEstMethod)
  {
  case PARAMS_FIXED:
    break;
  case PARAMS_DYN_THRESH:
  { // modify noise model if measured uncertainty is greater than the defined
    //  measurement uncertainties

    // add isotropic variance to noise model to account for match uncertainty
    ComputeMatchUncertaintyEstimates();

    // update positional effective noise model
    if (traceM_est > meanTraceM_msmt)
    {
      // divide difference in trace estimate by dimensionality of the space
      //  to distribute increase in covariance equally along each dimension
      double sigma2_diff = (traceM_est - meanTraceM_msmt) / 3.0;
      vct3x3 isoVar;
      isoVar.SetAll(0.0);
      isoVar.Element(0, 0) = sigma2_diff;
      isoVar.Element(1, 1) = sigma2_diff;
      isoVar.Element(2, 2) = sigma2_diff;
      for (unsigned int i = 0; i < nSamples; i++)
      {
        // TODO: can speed this up by only recomputing diagonal matrix D
        //       and re-multiplying to get invM, N, etc.
        M[i] = M_msmt[i] + isoVar;
      }
      ComputeCovDecompositions( M, invM, N, invN, Dmin, Emin);      
    }
    meanSigma2 = traceM_est > meanTraceM_msmt ? traceM_est : meanTraceM_msmt;
    meanSigma2 /= 3.0;

    // update orientation effective noise model
    if (k_est < meanK_msmt)
    {
      double diff = meanK_msmt - k_est;
      // reduce eccentricity when expanding k to maintain similar
      //  physical size of anisotropy
      double E;
      for (unsigned int i = 0; i < nSamples; i++)
      {
        // NOTE: E = Emsmt * k_est / k_msmt
        k[i] = (k_msmt[i] - diff) > 0.0 ? k_msmt[i] - diff : 0.0;        
        E = E_msmt[i] * k[i] / k_msmt[i];
        B[i] = E * k[i] / 2.0;
      }
    }
    meanK = k_est < meanK_msmt ? k_est : meanK_msmt;
    meanE = k_est < meanK_msmt ? meanE_msmt * k_est / meanK_msmt : meanE_msmt;
    break;
  }
  case PARAMS_FIXED_1stITER_POS:
  {
    meanK = meanK_msmt;
    meanE = meanE_msmt;
    M = M_msmt;
    k = k_msmt;
    for (unsigned int i = 0; i < nSamples; i++)
    {
      B[i] = E_msmt[i] * k_msmt[i] / 2.0;
    }
    ComputeCovDecompositions(M, invM, N, invN, Dmin, Emin);
    break;
  }
  default:
  {
    std::cout << std::endl
      << "ERROR: unrecognized parameter estimation setting!" << std::endl;
    assert(0);
  }
  }
}

void algDirICP_GIMLOP::UpdateNoiseModel_SamplesXfmd(vctFrm3 &Freg)
{
  // update noise models of the transformed sample points
  vctRot3 R(Freg.Rotation());
  for (unsigned int s = 0; s < nSamples; s++)
  {
    R_invM_Rt[s] = R*invM[s] * R.Transpose();
    N_Rt[s] = N[s] * R.Transpose();
    inv_N_Rt[s] = R*invN[s];
    R_L[s] = R*L[s];
  }
}


// dynamic portion of noise model
void algDirICP_GIMLOP::ComputeMatchUncertaintyEstimates()
{
  double sumSqrDist = 0.0;
  double sumNormProducts = 0.0;
  for (unsigned int s = 0; s < nSamples; s++)
  {
	  if (outlierFlags[s])	continue;	// skip outliers

    sumSqrDist += (samplePtsXfmd.Element(s) - matchPts.Element(s)).NormSquare();
    sumNormProducts += vctDotProduct(sampleNormsXfmd.Element(s), matchNorms.Element(s));
  }

  // Gaussian Parameters
  //  NOTE: the trace estimate here is 3 times the sigma2 estimate
  traceM_est = sumSqrDist / nSamples;

  // Kent Parameters
  //
  //  MLE estimate for k is an implicit ratio of
  //   two Bessel functions => no analytical soln.
  //
  //  MLE for k:  set Ap(k) = R
  //              
  //    Ap(k) = Ip/2(k) / Ip/2-1(k)       <--- Bessel Functions
  //    p =  spatial dimension (i.e. 3)
  //    R =  Sum_i(dot(Ny,Nx))/N
  //
  //  Closed form approximation for k:  R(p-R^2)/(1-R^2)
  //

  //  angular error of orientation vectors
  // R must be >= 0 for a dist'n that is biased towards its mean direction
  //  (i.e. for a Fisher-like distn')
  // It is possible that all matches point in opposite directions
  //   => R could be negative from the calculation above.
  double Rnorm = sumNormProducts / nSamples;
  Rnorm = Rnorm < 0.0 ? 0.0 : Rnorm;
  //  angular error of positions (wrt ctr of rotation)
  double Rpos = ComputeRpos();
  //  effective angular error
  // only include positional error if it reduces confidence in orientation
  //  (i.e. if it reduces the value of k)
  double R;
  double wRpos = 0.5;
  if (Rpos < Rnorm)
  {
    R = wRpos*Rpos + (1.0 - wRpos)*Rnorm;
  }
  else
  {
    R = Rnorm;
  }
  double R2 = R*R;
  k_est = R*(3.0 - R2) / (1.0 - R2);    // approx for k MLE
}


double algDirICP_GIMLOP::ComputeRpos()
{

#define NMLZ_METHOD 1   // Normalization method

  // Compute angular match error of sample positions relative 
  //  to the center of rotation; this is to include angular error
  //  of position matches measured from center of rotation of
  //  Procrustes problem as added indicator of rotational uncertainty

  vctDynamicVectorRef<vct3> S(samplePtsXfmd);
  vctDynamicVectorRef<vct3> C(matchPts);
  unsigned int N = S.size();
  vct3 Smean = vctCentroid(S);
  vct3 Cmean = vctCentroid(C);
  vct3 Sp, Cp;
  double dotProducts = 0.0;
  double nmlzTerm = 0.0;
  for (unsigned int i = 0; i < N; i++)
  {
    Sp.DifferenceOf(S[i], Smean);
    Cp.DifferenceOf(C[i], Cmean);

    //  TODO: should normalization be a linear or square term?
    //
    //  Do not scale down to unit vectors, as points farther
    //    away should have greater effect on rotation.
    //  R calculation assumes unit vector normalization => need some
    //    form of normalization in the end.    
#if (NMLZ_METHOD == 1)
    // use square normalization
    dotProducts += vctDotProduct(Sp, Cp);
    nmlzTerm += Sp.Norm()*Cp.Norm();
#elif (NMLZ_METHOD == 2)
    // use linear normalization
    double avgLen = (Sp.Norm() + Cp.Norm()) / 2.0;
    dotProducts += avgLen*vctDotProduct(Sp.Normalized(), Cp.Normalized());
    nmlzTerm += avgLen;
#elif (NMLZ_METHOD == 3)
    // use unit vectors
    dotProducts += vctDotProduct(Sp.Normalized(), Cp.Normalized());
#else
    std::cout << "Error: normalization method unrecognized" << std::endl;
#endif
  }
#if (NMLZ_METHOD == 3)
  nmlzTerm = N;
#endif

  double Rpos = dotProducts / nmlzTerm;
  Rpos = Rpos < 0.0 ? 0.0 : Rpos;
  //double posCircSD = sqrt(-2*log(Rpos));

  return Rpos;
}


// PD Tree Methods

void algDirICP_GIMLOP::SamplePreMatch(unsigned int sampleIndex)
{
  sampleK = k[sampleIndex];
  sampleB = B[sampleIndex];
  sampleR_L = R_L[sampleIndex];
  sampleDmin = Dmin[sampleIndex];
  sampleEmin = Emin[sampleIndex];
  sampleR_InvM_Rt = R_invM_Rt[sampleIndex];
  sampleN_Rt = N_Rt[sampleIndex];
  sample_inv_N_Rt = inv_N_Rt[sampleIndex];

#ifdef SAVE_MATCHES
  if (MatchIsotropic)
  {
    sampleB = 0.0;
  }
#endif
}


// fast check if a node might contain a datum having smaller match error
//  than the error bound
int algDirICP_GIMLOP::NodeMightBeCloser(
  const vct3 &Xp, const vct3 &Xn,
  DirPDTreeNode const *node,
  double ErrorBound)
{

  // Check if point lies w/in search range of the bounding box for this node
  //
  // Enlarge node boundary relative to the error bound:
  //
  //  cost:  [k - k*Nclosest'*Xn - B((L2'*Xn)^2 - (L3'*Xn)^2)] + (1/2)*(v-closest)'*inv(M)*(v-closest)
  //   
  //  Orientation Bound:  Compute lowest possible orientation error by choosing the point on the
  //                      node orientation boundary that has minimal Kent error. The node orienation
  //                      boundary is formed by rotating the average node orientation in all directions
  //                      by the max angular deviation from the average orientation of the node.
  //
  //    If we know the average orientation (Navg) and the maximum angular deviation (dThetaMax)
  //     from the average normal for all datums in a node, then:
  //    
  //      maxSearchDist = sqrt( (CurrentBestMatchError - MinKentError)*2*sigma2 )
  //
  //  Note: Xn and L1 are the same value  
  // 

  // --- Compute Lower Bound on Orientation Error --- //

  double MinKentError;
  vct3 yNewton;

  // compute unconstrained "Newton step" point that intersects the node boundary
  //  Note: this is simply the boundary point that intersects the arc on the unit sphere
  //        that runs most directly from Navg to Xn
  // rotation axis: Navg -> Xn & Navg -> Nnewton
  vct3 axis = vctCrossProduct(node->Navg, Xn);
  // Navg angle: Navg -> Xn
  double dThetaNavg = acos(vctDotProduct(node->Navg, Xn));
  // if node boundary overlaps Xn or if angle bewteen Navg and Xn < 0.1 degrees
  //  then assume perfect orientation match 
  // Note: sin(0.1*PI/180.0)^2 = 3.0462e-06
  if (dThetaNavg <= node->dThetaMax || axis.NormSquare() < 3.0462e-06)
  {
    // NOTE: since we add an extra "k" to the orientation match error, 0 is the minimal
    //       possible match error for orientation
    MinKentError = 0.0;
  }
  else
  {
    axis.NormalizedSelf();
    // shortened "newton" step boundary point
    //  (rotate Navg towards Xn by the max angular deviation within the node)
    vctAxAnRot3 Rnewton(axis, node->dThetaMax);
    yNewton = vctRot3(Rnewton) * node->Navg;

    // evaluate error +/- 0.05 deg on either side of newton point
    //  as initial test points for the routine to bracket the global min
    // NOTE: we must use very small step size because it is possible to have
    //       two local minimums and the 2nd local min may be closer physically
    //       than the global local min; note as well that there is a local max
    //       between the shortened newton step and the 2nd local min => if we
    //       choose a small step size then the bracket routine will always find
    //       the global min.
    static const double x1 = 0.05*cmnPI / 180.0;
    static const double x2 = -0.05*cmnPI / 180.0;
    //double f1 = KentBoundaryMatchError( x1,Xn,yNewton,node->Navg,sampleK,sampleB,sampleL );
    //double f2 = KentBoundaryMatchError( x2,Xn,yNewton,node->Navg,sampleK,sampleB,sampleL );

    // Brent Method
    // bracket the global min
    Brent brent;
    brent.bracket(x1, x2, Xn, yNewton, node->Navg, sampleK, sampleB, sampleR_L);
    // find the global min
    MinKentError = brent.minimize(Xn, yNewton, node->Navg, sampleK, sampleB, sampleR_L);

    //// Golden Method
    //Golden golden;
    //golden.bracket(x1,x2,Xn,yNewton,node->Navg,sampleK,sampleB,sampleL);
    //// compute optimal boundary
    //MinKentError = golden.minimize(Xn,yNewton,node->Navg,sampleK,sampleB,sampleL);
  }


  // --- Positional Node Distance Test --- //  

#ifdef KENT_POS_ISOTROPIC
  vct3 Xp_node = node->F * Xp;        // transform point into local coordinate system of node
  //vct3 Xn_node = F.Rotation()*Xn;   // don't need this since average node orientation is calculated
  //  wrt the world reference frame, not local reference frame

  // Assume ErrorBound = MinKentError + (1/2)*(v-closest)'*inv(M)*(v-closest) 
  //                   = MinKentError + (1/2)*(v-closest)'*diag(Emin)*(v-closest)
  //  => MaxDist = sqrt[ (ErrorBound - MinKentError)*2.0/Emin ]
  double searchDist2 = (ErrorBound - MinKentError)*2.0/sampleEmin;
  if (searchDist2 < 0.0)
  { // orientation error alone dismisses this node
    return 0;
  }

  // Rather than comparing only the x-axis value, check all coordinate directions
  //  of the node bounding box to further refine whether this node may be within 
  //  the search range of this point. Using the node coordinate frame is still 
  //  useful in this context, because it ensures a node bounding box of minimum size.
  // Conceptually, this check places another bounding box of size search distance
  //  around the point and then checks if this bounding box intersects the bounding
  //  box of this node.
  return node->Bounds.Includes(Xp_node, sqrt(searchDist2));
#else
  // Assume ErrorBound = MinKentError + (1/2)*(v-closest)'*inv(M)*(v-closest) 
  //  => MaxSqrMahalanobisError = (ErrorBound - MinKentError)*2.0
  double MaxSqrMahalError = (ErrorBound - MinKentError)*2.0;
  // Ellipsoid / OBB Intersection Test
  //  test if the ellipsoid centered at the sample point and defined by the 
  //  level set of MaxSqrMahalError intersects the oriented bounding box of this node
  return IntersectionSolver.Test_Ellipsoid_OBB_Intersection(
    Xp, node->Bounds, node->F,
    MaxSqrMahalError,
    sampleN_Rt, sampleDmin);
#endif
}

// Helper Methods

void algDirICP_GIMLOP::SetSamples(
  const vctDynamicVector<vct3> &argSamplePts,
  const vctDynamicVector<vct3> &argSampleNorms,
  const vctDynamicVector<double> &argK,
  const vctDynamicVector<double> &argE,
  const vctDynamicVector<vct3x2> &argL,
  const vctDynamicVector<vct3x3> &argM,
  double argScale, bool argbScale,
  PARAM_EST_TYPE paramEst)
{
  // base class
  nSamples = argSamplePts.size();

  R_MsmtM_Rt.SetSize(nSamples);
  outlierFlags.SetSize(nSamples);
  residuals_PostMatch.SetSize(nSamples);
  sqrDist_PostMatch.SetSize(nSamples);

  // scale sample points (remove this from here when you move it to IMLOP)
  sc = argScale;
  bScale = argbScale;
  if (bScale) {
	  for (unsigned int i = 0; i < nSamples; i++)
		  samplePts[i] = argSamplePts[i] * sc;
	  nTrans = 7; // 3 for rotation, 3 for translation, 1 for scale
  }
  else
	  nTrans = 6; // 3 for rotation, 3 for translation
  std::cout << "# tranformation parameters: " << nTrans << std::endl;

  algDirICP::SetSamples(argSamplePts, argSampleNorms);
  
  SetNoiseModel(argK, argE, argL, argM, paramEst);

  x_prev.SetSize(nTrans);
}


void algDirICP_GIMLOP::SetNoiseModel(
  const vctDynamicVector<double> &argK,
  const vctDynamicVector<double> &argE,
  const vctDynamicVector<vct3x2> &argL,
  const vctDynamicVector<vct3x3> &argM,
  PARAM_EST_TYPE paramEst)
{
  if (nSamples != argK.size() || nSamples != argE.size() ||
    nSamples != argL.size() || nSamples != argM.size())
  {
    std::cout << "ERROR: noise model does not match sample size: " 
      << nSamples << std::endl;
    assert(0);
  }

  paramEstMethod = paramEst;

  // measurement noise model parameters
  L = argL;
  k_msmt = argK;
  E_msmt = argE; 
  M_msmt = argM;
  invM_msmt.SetSize(nSamples);
  N_msmt.SetSize(nSamples);
  invN_msmt.SetSize(nSamples);
  Dmin_msmt.SetSize(nSamples);
  Emin_msmt.SetSize(nSamples);

  ComputeCovDecompositions(
    M_msmt, invM_msmt, N_msmt, invN_msmt,
    Dmin_msmt, Emin_msmt);

  // TODO: only include non-zero values in the means
  //       (since some data may leave out orientation or position for some samples)
  meanK_msmt = k_msmt.SumOfElements() / nSamples;
  meanE_msmt = E_msmt.SumOfElements() / nSamples;

  //traceM_msmt.SetSize(nSamples);
  meanTraceM_msmt = 0.0;
  for (unsigned int i = 0; i < nSamples; i++)
  {
    meanTraceM_msmt += M_msmt[i].Trace();
  }
  meanTraceM_msmt /= nSamples;

  // effective noise model parameters
  k.SetSize(nSamples);
  B.SetSize(nSamples);
  M.SetSize(nSamples);
  invM.SetSize(nSamples);
  N.SetSize(nSamples);
  invN.SetSize(nSamples);
  Dmin.SetSize(nSamples);
  Emin.SetSize(nSamples);

  R_L.SetSize(nSamples);
  R_invM_Rt.SetSize(nSamples);
  N_Rt.SetSize(nSamples);
  inv_N_Rt.SetSize(nSamples);

  // optimizer variables
  Yp_t.SetSize(nSamples);
  Rat_Yp_RaXp_t.SetSize(nSamples);  
  invM_Rat_Yp_RaXp_t.SetSize(nSamples);
  RaXn.SetSize(nSamples);
  RaRL.SetSize(nSamples);
}


// Note: this function depends on the SamplePreMatch() function
//       to set the noise model of the current transformed sample
//       point before this function is called
void algDirICP_GIMLOP::ComputeCovDecompositions(
  const vctDynamicVector<vct3x3> &M, vctDynamicVector<vct3x3> &invM,
  vctDynamicVector<vct3x3> &N, vctDynamicVector<vct3x3> &invN,
  vctDoubleVec &Dmin, vctDoubleVec &Emin)
{
  if (M.size() != nSamples)
  {
    std::cout << "ERROR: cannot compute Gaussian noise model decomposition; "
      << "noise model has not been specified or no samples exist" << std::endl;
    assert(0);
  }

  // compute decomposition of noise covariances for each sample
  //static vctFixedSizeMatrix<double, 3, 3, VCT_COL_MAJOR> A;
  //static vctFixedSizeMatrix<double, 3, 3, VCT_COL_MAJOR> U;
  //static vctFixedSizeMatrix<double, 3, 3, VCT_COL_MAJOR> Vt;
  //static vct3 S;
  //static vct3 D;
  //static nmrSVDFixedSizeData<3, 3, VCT_COL_MAJOR>::VectorTypeWorkspace workspace;
  for (unsigned int i = 0; i < nSamples; i++)
  {
    // Compute Decomposition of Minv = N'*N
    //   Minv = N'*N = R*D^2*R' = U*S*V'   =>   R'=V', D=sqrt(S), R=U
    //   N = D*R'   invN = R*inv(D)
    ComputeCovDecomposition_NonIter(
      M[i], invM[i], N[i], invN[i], Dmin[i], Emin[i]);
  }
}

void algDirICP_GIMLOP::ComputeCovDecomposition_NonIter(
  const vct3x3 &M, vct3x3 &Minv, vct3x3 &N, vct3x3 &Ninv,
  double &Dmin, double &Emin)
  //double &det_M)
{
  // Compute eigen decomposition of M
  //   M = V*diag(S)*V'
  // eigen values in descending order
  // eigen vectors listed by column   (has determinant = 1)
  vct3    eigenValues;
  vct3x3  eigenVectors;
  ComputeCovEigenDecomposition_NonIter(M, eigenValues, eigenVectors);

  // Compute Minv
  //   Minv = V*diag(1/S)*V'
  static vctFixedSizeMatrix<double, 3, 3, VCT_COL_MAJOR> V_Sinv;
  static vct3 Sinv;
  Sinv[0] = 1.0 / eigenValues[0];
  Sinv[1] = 1.0 / eigenValues[1];
  Sinv[2] = 1.0 / eigenValues[2];
  V_Sinv.Column(0) = eigenVectors.Column(0)*Sinv[0];
  V_Sinv.Column(1) = eigenVectors.Column(1)*Sinv[1];
  V_Sinv.Column(2) = eigenVectors.Column(2)*Sinv[2];
  Minv.Assign(V_Sinv * eigenVectors.TransposeRef());

  // Compute Decomposition of Minv = N'*N
  //   Minv = R*D^2*R' = N'*N     M = R*Dinv^2*R' => R' = V', Dinv = sqrt(S)
  //   N = D*R'      Ninv = R*inv(D)
  //   Dmin  ~ sqrt of smallest eigenvalue of Minv
  //   Emin  ~ smallest eigenvalue of Minv => Emin = Dmin.^2
  vct3 Dinv(
    sqrt(eigenValues[0]),
    sqrt(eigenValues[1]),
    sqrt(eigenValues[2])
    );
  N.Row(0) = eigenVectors.Column(0) / Dinv[0];
  N.Row(1) = eigenVectors.Column(1) / Dinv[1];
  N.Row(2) = eigenVectors.Column(2) / Dinv[2];
  Ninv.Column(0) = eigenVectors.Column(0)*Dinv[0];
  Ninv.Column(1) = eigenVectors.Column(1)*Dinv[1];
  Ninv.Column(2) = eigenVectors.Column(2)*Dinv[2];

  Emin = 1.0 / eigenValues[0];
  Dmin = 1.0 / Dinv[0];

  // Compute determinant of M
  //det_M = eigenValues.ProductOfElements();
}


void algDirICP_GIMLOP::UpdateOptimizerCalculations(const /*vct6*/ vctDynamicVector<double> &x)
{
  a.Assign(x[0], x[1], x[2]);
  t.Assign(x[3], x[4], x[5]);
  if (bScale)
	  sc = x[6];

  // matrix for rotation increment
  Ra = vctRot3(vctRodRot3(a));

  Xp_xfm = samplePtsXfmd;
  vctDynamicVectorRef<vct3>   Xp_xfm(samplePtsXfmd);
  vctDynamicVectorRef<vct3>   Xn_xfm(sampleNormsXfmd);
  vctDynamicVectorRef<vct3>   Yp(matchPts);
  vctDynamicVectorRef<vct3>   Yn(matchNorms);

  unsigned int i;
#ifdef ENABLE_PARALLELIZATION
#pragma omp parallel for
#endif
  for (i = 0; i < nSamples; i++)
  {
    //--- orientation ---//
    RaXn[i] = Ra*Xn_xfm[i];
    RaRL[i] = Ra*R_L[i];

    //--- position---//
    Yp_t[i] = Yp[i] - t;
    Rat_Yp_RaXp_t[i] = Ra.TransposeRef() * Yp_t[i] - sc * Xp_xfm[i];
#ifdef KENT_POS_ISOTROPIC
    invM_Rat_Yp_RaXp_t[i] = Emin[i]*(Ra.TransposeRef() * Yp_t[i] - Xp_xfm[i]);
#else
    invM_Rat_Yp_RaXp_t[i] = R_invM_Rt[i] * Rat_Yp_RaXp_t[i];
#endif
  }
  x_prev = x;
}

double algDirICP_GIMLOP::CostFunctionValue(const vctDynamicVector<double> &x)
{
  // don't recompute these if already computed for gradient
  if (x.NotEqual(x_prev))
  {
    UpdateOptimizerCalculations(x);
  }

  double f = 0.0;
  unsigned int i;
#ifdef ENABLE_PARALLELIZATION
#pragma omp parallel for
#endif
  vctDynamicVectorRef<vct3>   Yn(matchNorms);
  for (i = 0; i < nSamples; i++)
  {
	if (outlierFlags[i])	continue;

    double major = RaRL[i].Column(0)*Yn[i];
    double minor = RaRL[i].Column(1)*Yn[i];
    // add extra k[i] * 1.0 to make cost function > 0
    f += k[i] * (1.0 - RaXn[i]*Yn[i]) - B[i] * (major*major - minor*minor) 
      + (Rat_Yp_RaXp_t[i] * invM_Rat_Yp_RaXp_t[i]) / 2.0;
  }

  return f;
}


void algDirICP_GIMLOP::CostFunctionGradient(const /*vct6*/ vctDynamicVector<double> &x, /*vct6*/ vctDynamicVector<double> &g)
{
  //vct3 alpha;         // alpha = a/norm(a)
  //double theta;       // theta = norm(a)
  //vct3x3 sk_alpha, AlphaAlphaT_I, aaT, theta_daat_adat;
  //double sTheta, cTheta, theta2, theta3, dTheta;
  vct3 gTheta;
  vctFixedSizeVector<vctRot3, 3> dRa;  // Jacobian of R(a) wrt ax,ay,az

  // don't recompute these if already computed for cost function value
  if (x.NotEqual(x_prev))
  {
    UpdateOptimizerCalculations(x);
  }

  ComputeRodriguesJacobians(a, dRa);

  // form the cost function gradient
  g.SetAll(0.0);
  vctFixedSizeVectorRef<double, 3, 1> ga(g, 0);
  vctFixedSizeVectorRef<double, 3, 1> gt(g, 3);
  vctFixedSizeVectorRef<double, 1, 1> gsc;
  if (bScale)
	  gsc.SetRef(g, 6);

  vctDynamicVectorRef<vct3>   Xp_xfm(samplePtsXfmd);
  vctDynamicVectorRef<vct3>   Xn_xfm(sampleNormsXfmd);
  vctDynamicVectorRef<vct3>   Yp(matchPts);
  vctDynamicVectorRef<vct3>   Yn(matchNorms);
  vct3x3 Jz_a;

  unsigned int j;
#ifdef ENABLE_PARALLELIZATION
#pragma omp parallel for
#endif
  for (j = 0; j < nSamples; j++)
  {
    // TODO: vectorize Kent term better
	  if (outlierFlags[j])	continue;

    //--- Kent term ---//   (orientations)        
    for (unsigned int i = 0; i < 3; i++)
    {
      //  rotational effect
      ga[i] += -k[j] * ((dRa[i] * Xn_xfm[j])*Yn[j])
        - 2.0*B[j] * ((dRa[i] * R_L[j].Column(0))*Yn[j] * (Yn[j] * RaRL[j].Column(0))
        - (dRa[i] * R_L[j].Column(1))*Yn[j] * (Yn[j] * RaRL[j].Column(1)));
    }

    //--- Gaussian term ---//   (positions)
    //  rotational effect
    Jz_a.Column(0) = dRa[0].TransposeRef() * Yp_t[j];
    Jz_a.Column(1) = dRa[1].TransposeRef() * Yp_t[j];
    Jz_a.Column(2) = dRa[2].TransposeRef() * Yp_t[j];
    ga += invM_Rat_Yp_RaXp_t[j] * Jz_a;
    // translational effect
    gt -= Ra * invM_Rat_Yp_RaXp_t[j];

	if (bScale)
		gsc += invM_Rat_Yp_RaXp_t[j] * (-Xp_xfm.Element(j));
  }
}

void algDirICP_GIMLOP::ReturnScale(double &scale)
{
	scale = sc;
}

double algDirICP_GIMLOP::MatchError(
  const vct3 &Xp, const vct3 &Xn,
  const vct3 &Yp, const vct3 &Yn,
  double k, double B, const vctFixedSizeMatrix<double, 3, 2> &L,
  const vct3x3 &invM)
{
  double major = L.Column(0)*Yn;
  double minor = L.Column(1)*Yn;

#ifdef KENT_POS_ISOTROPIC
  // add extra k term so that match error >= 0
  return k*(1 - Xn*Yn) - B*(major*major - minor*minor) + (Xp-Yp).NormSquare()*invM(0,0)/2.0;
#else
  // add extra k term so that match error >= 0  
  return k*(1 - Xn*Yn) - B*(major*major - minor*minor) + ((Xp - Yp)*invM*(Xp - Yp)) / 2.0;
#endif
}
