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

#include "algDirICP_IMLOP.h"
#include "DirPDTreeNode.h"
#include "RegisterP2P.h"
#include "utilities.h"

#define EPS  1e-12

void  algDirICP_IMLOP::SetNoiseModel(
  double initK, double initSigma2, double w_Rpos, bool dynamicallyEstParams)
{
  k_init = initK;
  sigma2_init = initSigma2;
  wRpos = w_Rpos;
  dynamicParamEst = dynamicallyEstParams;
}

void algDirICP_IMLOP::ComputeMatchStatistics(double &Avg, double &StdDev) //TODO: do this right with mahalanobis distance
{
	sumSqrMatchDist = 0.0;
	double sumMatchDist = 0.0;
	double sqrMatchDist;
	double matchAngle;
	sumSqrMatchAngle = 0.0;

	// NOTE: if using a method with outlier rejection, it may be desirable to
	//       compute statistics on only the inliers
	for (unsigned int i = 0; i < nGoodSamples; i++)
	{
		sqrMatchDist = (goodMatchPts[i] - Freg * goodSamplePts[i]).NormSquare();

		sumSqrMatchDist += sqrMatchDist;
		sumMatchDist += sqrt(sqrMatchDist);

		matchAngle = acos(std::fmod(goodMatchNorms[i] * (Freg.Rotation() * goodSampleNorms[i]), 2 * cmnPI));

		sumSqrMatchAngle += k_init * matchAngle * matchAngle;
	}
	Avg = sumMatchDist / (sigma2*nGoodSamples);
	StdDev = sqrt( (sumSqrMatchDist / (sigma2*nGoodSamples))+Avg*Avg );
}

void algDirICP_IMLOP::PrintMatchStatistics(std::stringstream &tMsg)
{
	// For registration rejection purpose:
	tMsg << "\nSum square mahalanobis distance = " << sumSqrMatchDist/sigma2 << " over " << nGoodSamples << " inliers";
	tMsg << "\nSum square match angle = " << sumSqrMatchAngle << " over " << nGoodSamples << " inliers\n";
}

void algDirICP_IMLOP::SetSamples(
  const vctDynamicVector<vct3> &argSamplePts,
  const vctDynamicVector<vct3> &argSampleNorms)
{
  // base class
  algDirICP::SetSamples(argSamplePts, argSampleNorms);

  //std::cout << "algDirICP_IMLOP::SetSamples()" << std::endl;

  unsigned int nSamples = samplePts.size();

  // allocate buffers
  goodSamplePtsBuf.SetSize(nSamples);
  goodSampleNormsBuf.SetSize(nSamples);
  goodMatchPtsBuf.SetSize(nSamples);
  goodMatchNormsBuf.SetSize(nSamples);
}

double algDirICP_IMLOP::ICP_EvaluateErrorFunction()
{
  //// Return the negative log likelihood of the vonMises-Fisher
  ////  and Gaussian distributions under the assumption
  ////   of independence between the two distributions
  ////
  ////   Negative Log-Likelihood:
  ////    -log[ C * exp( k*dot(Ny,Nx) - B*||Y-X||^2 ) ]
  ////        where C  =  product of normalizations terms
  ////                    for Fisher and 3D Gaussian distributions
  ////              C  =  [k/(2*PI*(e^k-e^-k))]*[1/(2*PI*sigma^2)^(3/2)]
  ////              B  =  1/(2*sigma^2)
  //double logC;
  //static const double log2PI = log(2 * cmnPI); // compute this once for efficiency

#ifdef TEST_STD_ICP
  //// Test Standard ICP Condition:
  ////   compute the log of the normalization constant C for
  ////   only a Gaussian distribution
  ////   (k=0 and orientations are not being considered)
  //logC = -(3.0/2.0)*log(2*cmnPI*sigma2);

  // include match errors of good samples and of thresholded outliers
  //  to improve smoothness of cost function
  return B*(pICP->gSumSqrDist_PostMatch + pICP->oSumSqrDist_PostMatch); //- nSamples*logC;
#endif

  vct3 residual;
  double sumSqrDist = 0.0;
  double sumNormProducts = 0.0;
  for (unsigned int s = 0; s < nGoodSamples; s++)
  {
    residual = goodSamplePts.Element(s) - goodMatchPts.Element(s);
    sumSqrDist += residual.NormSquare();

    sumNormProducts += vctDotProduct(goodSampleNorms.Element(s), goodMatchNorms.Element(s));
  }

  // include match errors of good samples and of thresholded outliers
  //  to improve smoothness of cost function
  // NOTE: adding an extra k*nSamples to the cost produces a cost function >= 0
  return B*(sumSqrDist)+k*(nSamples - sumNormProducts); // -logC*nSamples;
}

vctFrm3 algDirICP_IMLOP::ICP_RegisterMatches()
{
  RegisterP2P_Normals_vMFG(
    goodSamplePts,
    goodMatchPts,
    goodSampleNorms,
    goodMatchNorms,
    B, k, Freg);

  return Freg;
}


void algDirICP_IMLOP::ICP_InitializeParameters(vctFrm3 &FGuess)
{
  // initialize base class
  algDirICP::ICP_InitializeParameters(FGuess);

  k = k_init;
  sigma2 = sigma2_init;
  B = 1.0 / (2.0*sigma2_init);

#ifdef TEST_STD_ICP
  k = 0.0;
#endif

  // monitoring variables
  errFuncNormWeight = k;
  errFuncPosWeight = B;
}

void algDirICP_IMLOP::ICP_UpdateParameters_PostRegister(vctFrm3 &Freg)
{
  // base class
  algDirICP::ICP_UpdateParameters_PostRegister(Freg);

  double sumSqrDist_PostRegister = 0.0;  
  double sumNormProducts_PostRegister = 0.0;
  for (unsigned int s = 0; s < nSamples; s++)
  { 
    sumSqrDist_PostRegister += 
      (samplePtsXfmd.Element(s) - matchPts.Element(s)).NormSquare();
    sumNormProducts_PostRegister += 
      vctDotProduct(sampleNormsXfmd.Element(s), matchNorms.Element(s));
  }

  // update noise model
  UpdateNoiseModel(sumSqrDist_PostRegister, sumNormProducts_PostRegister);

  // update monitoring variables
  errFuncNormWeight = k;
  errFuncPosWeight = B;
}

unsigned int algDirICP_IMLOP::ICP_FilterMatches()
{

  // Method 1: No Outlier Detection
#if 1
  nGoodSamples = 0;
  nOutliers = 0;
  nPosOutliers = 0;
  nNormOutliers = 0;

  // use all samples
  for (unsigned int s = 0; s < nSamples; s++)
  { // copy all samples to buffer
    goodSamplePtsBuf.Element(s) = samplePts.Element(s);
    goodMatchPtsBuf.Element(s) = matchPts.Element(s);

    goodSampleNormsBuf.Element(s) = sampleNorms.Element(s);
    goodMatchNormsBuf.Element(s) = matchNorms.Element(s);
  }
  nGoodSamples = nSamples;

  // Non-destructively resize good sample reference vectors
  goodSamplePts.SetRef(goodSamplePtsBuf, 0, nGoodSamples);
  goodMatchPts.SetRef(goodMatchPtsBuf, 0, nGoodSamples);

  goodSampleNorms.SetRef(goodSampleNormsBuf, 0, nGoodSamples);
  goodMatchNorms.SetRef(goodMatchNormsBuf, 0, nGoodSamples);

  return nOutliers;
#endif 

  // Method 2: Chi-Square Outlier Test
  // TODO
}


// PD Tree Methods

int algDirICP_IMLOP::NodeMightBeCloser(
  const vct3 &v, const vct3 &n,
  DirPDTreeNode const *node,
  double ErrorBound)
{
  vct3 Fv = node->F*v;          // transform point into local coordinate system of node
  //vct3 Fn = F.Rotation()*n;   // don't need this since normal statistics for a node are calculated
  //  wrt the world reference frame, not local reference frame

  // Check if point lies w/in search range of the bounding box for this node
  //
  // Enlarge node boundary relative to the error bound:
  //
  //  cost:  k*(1-N'*Nclosest) + B*||v - closest||^2 
  //                                  (dist^2)
  //
  //  Improved Bound:  (Assume Best Match Normal is aligned by Max Angular Deviation
  //                    from the Mean Orientation of the Node)
  //
  //    If we know the avg normal (Navg) and the maximum angular deviation (dThetaMax)
  //     from the average normal for all triangles in a node, then:
  //     
  //     given: 0 < dThetaMax < 180 (deg)
  //     N'*Navg = cos(dThetaAvg)
  //      set ThetaC = dThetaAvg - dThetaMax    (assume avg normal tilted towards n by max deviation)
  //      if dThetaC < 0 then dThetaC = 0
  //     set N'*Nc = cos(dThetaC)    (this is the max possible norm product (lowest possible orienation error) 
  //                                  since N'*Nc <= cos(dThetaC) by reason of max deviation)
  //     
  //     =>  cost = k*(1-cos(dThetaC)) + B*dist^2
  //         maxSearchDist = sqrt([ErrorBound - k*(1-cos(dThetaC))]/B)
  //
  //     =>  search a node if the node boundary enlarged by
  //         maxSearchDist contains the point v
  //

  // Simple Bound
  //double searchDist2 = ErrorBound/posWeight;

  // Improved Bound
  double dThetaAvg = acos(n.DotProduct(node->Navg));
  double dThetaC = dThetaAvg - node->dThetaMax;
  dThetaC = dThetaC > 0.0 ? dThetaC : 0.0;   // enforce ThetaC >= 0
  double searchDist2 = (ErrorBound - k*(1 - cos(dThetaC))) / B;
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
  return node->Bounds.Includes(Fv, sqrt(searchDist2));


  // Search Method for Alternate Cost Function:
  // -----------------------------------------
  // Enlarge node boundary by the weighted normal distance
  //  to ensure all possible candidates are compared for cost function:
  //
  //  cost:  -k*N'*Nclosest + ||v - closest||^2 
  //
  //  Since the PD tree sorts only by position and not by the normal
  //  vector, the farthest positional search distance is obtained by assuming
  //  a triangle normal parallel to the point normal.  => cost = -k + dist^2
  //  
  //  Assuming DistBound is the smallest "cost" found so far, we have maximum
  //  search distance:  dist^2 = DistBound + k
  //                    dist = sqrt(DistBound + k)
  //
  //  NOTE: by this formulation, "cost" and "DistBound" may be negative, however,
  //        the search distance (after adding k) is always >= 0.
  //
}

// Helper Methods:

void algDirICP_IMLOP::UpdateNoiseModel(double sumSqrDist, double sumNormProducts)
{
  // Position Parameters
  // compute Gaussian variables
  //  B = 1/(2*sigma2)
  // divide sum of square distances by number of samples times 3 
  //  to get the estimate of variance because square distance is 
  //  a summation of three square Gaussian RVs, one for each coordinate axis
  double S2 = sumSqrDist / (3.0*nSamples);
  if (dynamicParamEst)
  {
    B = 1.0 / (2.0*S2);
    B = B > threshB ? threshB : B;  // threshold in case of perfect matches
    sigma2 = 1.0 / (2.0*B);
  }

  // Orientation Parameters
  // compute Fisher variables
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
  //   NOTE: circular standard deviation of theta
  //         circSD = sqrt(-2*ln(R))
  //            where theta is angle between matched vectors
  //            i.e. theta = acos(dot(Ny,Nx))
  //

  // angular error of normal orientations
  ComputeCircErrorStatistics(sumNormProducts, Rnorm, circSD);

  // angular error of positions (wrt ctr of rotation)
  Rpos = ComputeRpos();

  // effective angular error
  // only include positional error if it reduces the value of K
  if (Rpos < Rnorm)
  {
    R = wRpos*Rpos + (1.0 - wRpos)*Rnorm;
  }
  else
  {
    R = Rnorm;
  }

  double R2 = R*R;
  if (dynamicParamEst)
  {
    // protect from division by zero
    if (R2 >= 1.0)
    {
      k = threshK;  // set k to its max value
    }
    else
    {
      k = R*(3.0 - R2) / (1.0 - R2);  // approx for k
      k = k > threshK ? threshK : k;  // perfect match threshold
    }

    // reduce k by a factor
    k = k * k_factor;
  }

#ifdef TEST_STD_ICP
  k = 0.0;
#endif
}


double algDirICP_IMLOP::ComputeRpos()
{

#define NMLZ_METHOD 1   // Normalization method

  // NOTE: could improve efficiency by computing this value as an 
  //       optional output argument in the vMFG P2P Registration
  //       function

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

    //  Do not scale down to unit vectors, as points farther
    //    away should have greater effect on rotation.
    //  R calculation assumes unit vector normalization => need some
    //    form of normalization in the end.
    //  TODO: should normalization be a linear or square term?
#if (NMLZ_METHOD == 1)
    // use square normalization
    //  this is the best solution, as it works with the data directly
    //  i.e. the orientations do not have to be normalized prior to
    //       taking their dot product
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
  // 0 <= Rpos <= 1
  Rpos = Rpos < 0.0 ? 0.0 : Rpos;
  Rpos = Rpos > 1.0 ? 1.0 : Rpos;
  //double posCircSD = sqrt(-2*log(Rpos));
  return Rpos;
}
