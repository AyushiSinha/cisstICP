// ****************************************************************************
//
//    Copyright (c) 2014, Ayushi Sinha, Seth Billings, Russell Taylor, Johns Hopkins University
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

#include "algDirICP_VIMLOP.h"
#include "cisstICP.h"
#include "PDTreeNode.h"
#include "RegisterP2P.h"
#include "utilities.h"

#include <limits.h>

#define COMPUTE_ERROR_FUNCTION
//#define DEBUG_IMLP
//#define REMOVE_OUTLIERS


algDirICP_VIMLOP::algDirICP_VIMLOP(
  PDTree_Mesh *pTree,
  //DirPDTree2DBase *pDirTree,
  vctDynamicVector<vct3> &samplePts, 
  vctDynamicVector<vct3x3> &sampleCov,
  vctDynamicVector<vct3x3> &sampleMsmtCov, 
  vctDynamicVector<vct2> &sample2dPts,
  vctDynamicVector<vct2> &sample2dNorms,
  vctDynamicVector<vct2x2> &sample2dMsmtCov,
  vctDynamicVector<vct2x2> &sample2dNormCov,
  camera cam,
  double outlierChiSquareThreshold,
  double sigma2Max)
  : algICP_IMLP_Mesh(pTree, samplePts, sampleCov, sampleMsmtCov, outlierChiSquareThreshold, sigma2Max)//, 
  //alg2D_DirICP(pDirTree, sample2dPts, sample2dNorms),
  //alg2D_DirPDTree_vonMises(pDirTree)
{
	SetSamples(samplePts, sampleCov, sampleMsmtCov, sample2dPts, sample2dNorms, sample2dMsmtCov, sample2dNormCov);
}


void algDirICP_VIMLOP::ComputeMatchStatistics(double &Avg, double &StdDev)
{
  // base class
	algICP_IMLP_Mesh::ComputeMatchStatistics(Avg, StdDev);
	//alg2D_DirICP::ComputeCircErrorStatistics()
}

// TODO: change this so that covariances for measurement noise
//       and surface model are specified independently
//       rather than specifying measurement noise model
//       and the combination of measurement noise and surface models
// NOTE: Mxi includes both MsmtMxi (measurement noise model)
//       and planar surface approximation model
void algDirICP_VIMLOP::SetSamples(
  vctDynamicVector<vct3> &argSamplePts,
  vctDynamicVector<vct3x3> &argMxi,
  vctDynamicVector<vct3x3> &argMsmtMxi,
  vctDynamicVector<vct2> &argSample2dPts,
  vctDynamicVector<vct2> &argSample2dNorms,
  vctDynamicVector<vct2x2> &arg2dMsmtMxi,
  vctDynamicVector<vct2x2> &arg2dNormMxi)
{
	algICP_IMLP_Mesh::SetSamples(argSamplePts, argMxi, argMsmtMxi);
	//alg2D_DirICP::SetSamples(argSample2dPts, argSample2dNorms);

	// 2D
	//MsmtMxi = arg2dMsmtMxi; 
	//NormMxi = arg2dNormMxi;

	//outlierFlags.SetSize(alg2D_DirICP::nSamples);
}

//// TODO: change this so that covariances for measurement noise
////       and surface model are specified independently
////       rather than specifying measurement noise model
////       and the combination of measurement noise and surface models
//// NOTE: Mi includes both MsmtMi (measurement noise model)
////       and planar surface approximation model
//void algICP_IMLP::SetSampleCovariances(
//  vctDynamicVector<vct3x3> &Mi,
//  vctDynamicVector<vct3x3> &MsmtMi)
//{
//  if (Mi.size() != nSamples || MsmtMi.size() != nSamples)
//  {
//    std::cout << "ERROR: number of covariances matrices does not match number of samples" << std::endl;
//  }
//
//  Mxi.SetSize(nSamples);
//  eigMxi.SetSize(nSamples);
//  for (unsigned int s = 0; s<nSamples; s++)
//  {
//    Mxi[s] = Mi[s];
//    ComputeCovEigenValues_SVD(Mi[s], eigMxi[s]);
//  }
//
//  this->MsmtMxi = MsmtMi;
//}

void algDirICP_VIMLOP::ICP_InitializeParameters(vctFrm3 &FGuess)
{
	algICP_IMLP_Mesh::ICP_InitializeParameters(FGuess);
	//alg2D_DirICP::ICP_InitializeParameters(FGuess2d);
}

void algDirICP_VIMLOP::ICP_ComputeMatches()
{
	// Compute 3D matches
	algICP_IMLP_Mesh::ICP_ComputeMatches();

	// Compute 2D matches
	//alg2D_DirICP::ICP_ComputeMatches(); // Currently, single frame; will have to rewrite this to handle frames

	//outlierFlags.SetAll(0);
}

void algDirICP_VIMLOP::ICP_UpdateParameters_PostMatch()
{
  // base class
  algICP_IMLP_Mesh::ICP_UpdateParameters_PostMatch();
  //alg2D_DirICP::ICP_UpdateParameters_PostMatch(); 

  // 2D - change vars later if need be
  //sumSqrDist_Inliers = 0.0;

  //for (unsigned int s = 0; s < alg2D_DirICP::nSamples; s++)
  //{
	 // residuals_PostMatch.Element(s) = alg2D_DirICP::samplePtsXfmd.Element(s) - alg2D_DirICP::matchPts.Element(s);
	 // sqrDist_PostMatch.Element(s) = residuals_PostMatch.Element(s).NormSquare();

	 // if (!outlierFlags[s])
	 // {
		//  sumSqrDist_Inliers += sqrDist_PostMatch.Element(s);
	 // }
  //}

  //// update the match uncertainty factor
  //sigma2 = sumSqrDist_Inliers / (alg2D_DirICP::nSamples - alg2D_DirICP::nOutliers);

  //// apply max threshold
  //if (sigma2 > sigma2Max)
  //{
	 // sigma2 = sigma2Max;
  //}

  // update noise model of the matches
  // TODO
  // here...

  //if (bFirstIter_Matches)
  //{
	 // UpdateNoiseModel_SamplesXfmd(FGuess);
  //}

  // Outlier noise model
  //vctRot2 R(FGuess.Rotation()); // project FGuess
  //for (unsigned int s = 0; s < alg2D_DirICP::nSamples; s++)
  //{
	 // R_MsmtMxi_Rt[s] = (R.operator*(MsmtMxi[s])).operator*(R.Transpose());
  //}
  //bFirstIter_Matches = false;
}

void algDirICP_VIMLOP::ICP_UpdateParameters_PostRegister(vctFrm3 &Freg)
{
  // base class
  algICP_IMLP_Mesh::ICP_UpdateParameters_PostRegister(Freg);
  //alg2D_DirICP::ICP_UpdateParameters_PostRegister(Freg);
}

void algDirICP_VIMLOP::UpdateNoiseModel_SamplesXfmd(vctFrm3 &Freg)
{
  // update noise models of the transformed sample points
	algICP_IMLP::UpdateNoiseModel_SamplesXfmd(Freg);
}


double algDirICP_VIMLOP::ICP_EvaluateErrorFunction()
{
	return algICP_IMLP_Mesh::ICP_EvaluateErrorFunction();

	//vct2 residual;
	//double sumSqrDist = 0.0;

	////for (unsigned int s = 0; s < nGoodSamples; s++)
	////{
	////	residual = goodSamplePts.Element(s) - goodMatchPts.Element(s);
	////	sumSqrDist += residual.NormSquare();
	////}

	//for (unsigned int s = 0; s < alg2D_DirICP::nSamples; s++)
	//{
	//	residual = alg2D_DirICP::samplePts.Element(s) - alg2D_DirICP::matchPts.Element(s);
	//	sumSqrDist += residual.NormSquare();
	//}

	//return sqrt(sumSqrDist / alg2D_DirICP::nSamples);
}

//bool algDirICP_VIMLOP::ICP_Terminate(vctFrm3 &F)
//{
//  if (bTerminateAlgorithm)
//  {
//    // return registration for last iteration of non-increasing cost
//    F = Fdec;
//    return true;
//  }
//  else
//  {
//    return false;
//  }
//}

// TODO: change &F to Freg
//  even better: return registration rather than pass-by-reference
vctFrm3 algDirICP_VIMLOP::ICP_RegisterMatches()
{
	return algICP_IMLP_Mesh::ICP_RegisterMatches();
	//Freg = alg2D_DirICP::ICP_RegisterMatches();

	//return Freg; // TODO
	// TODO: implement the 2D version here
}

unsigned int algDirICP_VIMLOP::ICP_FilterMatches()
{
	return algICP_IMLP_Mesh::ICP_FilterMatches();
	//alg2D_DirICP::ICP_FilterMatches();

	//return 1; // TODO
	// TODO: implement the 2D version here
}



// PD Tree Methods

//void algDirICP_VIMLOP::SamplePreMatch(unsigned int sampleIndex)
//{
//	algICP_IMLP::SamplePreMatch(sampleIndex);
//	//alg2D_DirICP::SamplePreMatch(sampleIndex);
//
//	// 2D
//	// measurement noise
//	sample_RMxRt_sigma2 = R_Mxi_Rt[sampleIndex];
//
//	// add term for match uncertainty here to the RMxiRt
//	//  term rather than to each Myi term for improved efficiency
//	sample_RMxRt_sigma2.Element(0, 0) += sigma2;
//	sample_RMxRt_sigma2.Element(1, 1) += sigma2;
//
//	// compute eigenvalues
//	//  Note: R does not change eigenvalues of Mxi => these may be precomputed
//	//        also, adding isotropic terms adds same amount to each eigenvalue
//	vct2 &sampEig = eigMxi[sampleIndex];
//	sample_RMxRt_sigma2_Eig[0] = sampEig[0] + sigma2;
//	sample_RMxRt_sigma2_Eig[1] = sampEig[1] + sigma2;
//
//#ifdef DEBUG_IMLP
//	if (sampleIndex == 0)
//	{
//		std::cout << "SamplePreMatch():" << std::endl;
//		std::cout << "sigma2: " << sigma2 << std::endl;
//		std::cout << "R_Mx0_Rt: " << std::endl << R_Mxi_Rt[0] << std::endl;
//		std::cout << "sample_RMxRt_sigma2: " << std::endl << sample_RMxRt_sigma2 << std::endl;
//	}
//#endif
//}

// fast check if a datum might have smaller match error than error bound
int algDirICP_VIMLOP::DatumMightBeCloser(const vct3 &v,
                                                int datum,
                                                double ErrorBound)
{
	return algICP_IMLP_Mesh::DatumMightBeCloser(v, datum, ErrorBound);
}

int algDirICP_VIMLOP::DatumMightBeCloser(const vct2 &v, const vct2 &n,
	int datum,
	double ErrorBound)
{
	return true;
}

// fast check if a node might contain a datum having smaller match error
//  than the error bound
int algDirICP_VIMLOP::NodeMightBeCloser(const vct3 &v,
	PDTreeNode *node,
	double ErrorBound)
{
	return algICP_IMLP_Mesh::NodeMightBeCloser(v, node, ErrorBound);
}


// Helper Methods

// Note: this function depends on the SamplePreMatch() function
//       to set the noise model of the current transformed sample
//       point before this function is called
void algDirICP_VIMLOP::ComputeNodeMatchCov(PDTreeNode *node, DirPDTree2DNode *dirNode)
{
	algICP_IMLP_Mesh::ComputeNodeMatchCov(node);
  //// Note:  This function is called when searching a node that is using 
  ////        its own noise model rather than that of its parent node

  //// Compute the effective noise model for this node, assuming the noise 
  ////  model of the transformed sample point has already been computed

  //// noise model of transformed sample
  //M = sample_RMxRt_sigma2;
  //// add the effective My for this node
  //M.Element(0, 0) += node->EigMax; // change to dirNode
  //M.Element(1, 1) += node->EigMax;

  //// TODO: can this be done using only the eigen decomposition
  ////       of RMxRt
  //// Compute Decomposition of M
  ////   M = V*S*V'
  //vct2    eigenValues;
  //vct2x2  eigenVectors;  
  //ComputeCovEigenDecomposition_NonIter(M, eigenValues, eigenVectors);

  //// Compute Decomposition of Minv = N'*N
  ////   Minv = R*D^2*R' = N'*N     M = R*Dinv^2*R' => R' = V', Dinv = sqrt(S)
  ////   N = D*R'      Ninv = R*inv(D)
  //vct3 Dinv(
  //  sqrt(eigenValues[0]),
  //  sqrt(eigenValues[1])
  //  );
  //N.Row(0) = eigenVectors.Column(0) / Dinv[0];
  //N.Row(1) = eigenVectors.Column(1) / Dinv[1];
  ////N.Row(0) = eigenVectors.TransposeRef().Row(0) / Dinv[0];
  ////N.Row(1) = eigenVectors.TransposeRef().Row(1) / Dinv[1];
  ////N.Row(2) = eigenVectors.TransposeRef().Row(2) / Dinv[2];
  //Dmin = 1.0/Dinv[0]; // eigenvalues are arranged in order of decreasing magnitude
}


void algDirICP_VIMLOP::ComputeCovDecomposition_NonIter(const vct3x3 &M, vct3x3 &Minv, double &det_M)
{
	algICP_IMLP_Mesh::ComputeCovDecomposition_NonIter(M, Minv, det_M);
}

void algDirICP_VIMLOP::ComputeCovDecomposition_NonIter(const vct3x3 &M, vct3x3 &Minv, vct3x3 &N, vct3x3 &Ninv, double &det_M)
{
	algICP_IMLP_Mesh::ComputeCovDecomposition_NonIter(M, Minv, N, Ninv, det_M);
}


void algDirICP_VIMLOP::ComputeCovDecomposition_SVD(const vct3x3 &M, vct3x3 &Minv, double &det_M)
{
	algICP_IMLP_Mesh::ComputeCovDecomposition_SVD(M, Minv, det_M);
}

void algDirICP_VIMLOP::ComputeCovDecomposition_SVD(const vct3x3 &M, vct3x3 &Minv,
                                                      vct3x3 &N, vct3x3 &Ninv, double &det_M )
{
	algICP_IMLP_Mesh::ComputeCovDecomposition_SVD(M, Minv, N, Ninv, det_M);
}

// finds the point on this datum with lowest match error
//  and returns the match error and closest point
double algDirICP_VIMLOP::FindClosestPointOnDatum(
	const vct3 &point,
	vct3 &closest,
	int datum)
{
	return algICP_IMLP_Mesh::FindClosestPointOnDatum(point, closest, datum);
}

double algDirICP_VIMLOP::FindClosestPointOnDatum(
	const vct2 &point,
	const vct2 &norm,
	vct2 &closestP,
	vct2 &closestN,
	int datum)
{
	return 1.0;
}


