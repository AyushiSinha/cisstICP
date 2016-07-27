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

#include "algICP_DIMLP.h"
#include "cisstICP.h"
#include "PDTreeNode.h"
#include "RegisterP2P.h"
#include "utilities.h"

#include <limits.h>

#define COMPUTE_ERROR_FUNCTION

algICP_DIMLP::algICP_DIMLP(
	PDTree_Mesh *pTree,
	vctDynamicVector<vct3> &samplePts,
	vctDynamicVector<vct3x3> &sampleCov,
	vctDynamicVector<vct3x3> &sampleMsmtCov,
	vctDynamicVector<vct3> &meanShape,
	vctDynamicVector<vct3> &sampleModes,
	vctDynamicVector<double> &sampleModeWts,
	double outlierChiSquareThreshold,
	double sigma2Max)
	: algICP_IMLP(pTree, samplePts, sampleCov, sampleMsmtCov, outlierChiSquareThreshold, sigma2Max),
	dlib(this),
	pTree(pTree),
	pMesh(pTree->MeshP),
	TCPS(*(pTree->MeshP))
{
	SetSamples(samplePts, sampleCov, sampleMsmtCov, meanShape, sampleModes, sampleModeWts);
}

void algICP_DIMLP::ComputeMatchStatistics(double &Avg, double &StdDev)
{
	// return the average mahalanobis distance of the matches
	//  based on point noise models only (measurement and surface model covariances)
	//  i.e. do not include sigma2
	double sumSqrMahalDist = 0.0;
	double sumMahalDist = 0.0;
	double sqrMahalDist;
	int nGoodSamples = 0;
	vct3x3 M, Minv;
	vct3 residual;

	for (unsigned int i = 0; i < nSamples; i++)
	{
		if (outlierFlags[i])	continue;	// skip outliers

		residual = Tssm_matchPts[i] - Freg * samplePts[i];
		M = Freg.Rotation() * Mxi[i] * Freg.Rotation().Transpose() + *Myi[i];
		ComputeCovInverse_NonIter(M, Minv);
		sqrMahalDist = residual*Minv*residual;

		sumSqrMahalDist += sqrMahalDist;
		sumMahalDist += sqrt(sqrMahalDist);
		nGoodSamples++;
	}
	Avg = sumMahalDist / nGoodSamples;
	StdDev = (sumSqrMahalDist / nGoodSamples) + Avg*Avg;
}

void algICP_DIMLP::SetSamples(
	vctDynamicVector<vct3> &argSamplePts,
	vctDynamicVector<vct3x3> &argMxi,
	vctDynamicVector<vct3x3> &argMsmtMxi,
	vctDynamicVector<vct3> &argMeanShape,
	vctDynamicVector<vct3> &argSampleModes,
	vctDynamicVector<double> &argSampleModeWts)
{
	/*if (argSampleModes.size() != nSamples)
	{
		std::cout << "ERROR: number of sample mode matrices (" << argSampleModes.size() << ") does not match number of samples (" << nSamples << ")" << std::endl;
	}*/

	// base class
	algICP_IMLP::SetSamples(argSamplePts, argMxi, argMsmtMxi);

	meanShape = argMeanShape;
	sampleModes = argSampleModes;
	sampleModeWts = argSampleModeWts;

	Si = pTree->MeshP->Si;
	wi = pTree->MeshP->wi;

	Tssm_Y_t.resize(nSamples);
	Rat_Tssm_Y_t_x.resize(nSamples);
	invMx_Rat_Tssm_Y_t_x.resize(nSamples);

	// Moved to cisstMesh
	// 
	//wi.SetSize(sampleModes.size());
	//Si.SetSize(sampleModes.size());
	//for (unsigned int s = 0; s < sampleModes.size(); s++)
	//{
	//	for (unsigned int c = 0; c < 3; c++)
	//	{
	//		// wi = sqrt(lambda_i*mi)
	//		wi[s].Element(c) = sqrt(sampleModeWts.Element(0)*sampleModes[s].Element(c));
	//	}
	//}
	//// si = wi*(V-meanV)
	//Si = wi * (pTree->MeshP->vertices - meanShape);
	
	eta = vct3(0.33, 0.33, 0.33);
}

void algICP_DIMLP::ICP_InitializeParameters(vctFrm3 &FGuess)
{
	// initialize base class 
	algICP_IMLP::ICP_InitializeParameters(FGuess);
	this->FGuess = FGuess;

	//Tssm_Y_t.SetAll(vct3(0.0));
	//Rat_Tssm_Y_t_x.SetAll(vct3(0.0));
	//invMx_Rat_Tssm_Y_t_x.SetAll(vct3(0.0));
}

void algICP_DIMLP::ICP_UpdateParameters_PostMatch()
{
	// base class
	//algICP_IMLP::ICP_UpdateParameters_PostMatch();
	T_ssm(matchPts, Tssm_matchPts);
	// compute sum of square distances of inliers
	sumSqrDist_Inliers = 0.0;
	for (unsigned int s = 0; s < nSamples; s++)
	{
		residuals_PostMatch.Element(s) = samplePtsXfmd.Element(s) - Tssm_matchPts.Element(s);
		sqrDist_PostMatch.Element(s) = residuals_PostMatch.Element(s).NormSquare();

		if (!outlierFlags[s])
		{
			sumSqrDist_Inliers += sqrDist_PostMatch.Element(s);
		}
	}
	//std::cout << "Here: " << samplePtsXfmd.Element(50) << " " << Tssm_matchPts.Element(50) << " "
		//<< residuals_PostMatch.Element(50) << " " << sqrDist_PostMatch.Element(50) << " " << sumSqrDist_Inliers << " " << std::endl;

	// update the match uncertainty factor
	sigma2 = sumSqrDist_Inliers / (nSamples - nOutliers);

	// apply max threshold
	if (sigma2 > sigma2Max)
	{
		sigma2 = sigma2Max;
	}

	// update noise models of the matches
	//std::cout << "Match Datums: " << matchDatums.Element(50) << " Sigma: " << sigma2 << std::endl;
	for (unsigned int s = 0; s < nSamples; s++)
	{
		//update target covariances
		Myi[s] = algICP::pTree->DatumCovPtr(matchDatums.Element(s));	// use pointer here for efficiency

		// target covariance with match uncertainty
		Myi_sigma2[s] = *Myi[s];
		Myi_sigma2[s].Element(0, 0) += sigma2;
		Myi_sigma2[s].Element(1, 1) += sigma2;
		Myi_sigma2[s].Element(2, 2) += sigma2;
	}
	//std::cout << "Myi: " << Myi_sigma2.Element(50) << std::endl;

	if (bFirstIter_Matches)
	{
		algICP_IMLP::UpdateNoiseModel_SamplesXfmd(FGuess);
	}

	vctRot3 R(FGuess.Rotation());
	for (unsigned int s = 0; s < nSamples; s++)
	{
		R_MsmtMxi_Rt[s] = R*MsmtMxi[s] * R.Transpose();
	}

#ifdef DEBUG_IMLP
	std::cout << "My0:" << std::endl << *Myi[0] << std::endl;
	std::cout << "My1:" << std::endl << *Myi[1] << std::endl;
#endif

	bFirstIter_Matches = false;
}

void algICP_DIMLP::ICP_UpdateParameters_PostRegister(vctFrm3 &Freg)
{
	// base class
	algICP_IMLP::ICP_UpdateParameters_PostRegister(Freg);
	//pTree = new PDTree_Mesh(*pMesh, 5, 5.0);
}

double algICP_DIMLP::ICP_EvaluateErrorFunction()
{

#ifdef COMPUTE_ERROR_FUNCTION
	//
	//	Compute negative log-likelihood of match probabilities for error function
	//
	//	likelihood function:	[Product_i( f_match(Xi, T_ssm(Yi), theta_i) )].f_shape(T_ssm(Y);s)
	//		where, in the case of Deformable IMLP,
	//			f_match	 = Product_i[ (1/sqrt((2*pi)^3*|Mi|) * exp( (-1/2)*(Yi-R*Xi-t)'*inv(Mi)*(Yi-R*Xi-t) ) ]
	//				(for more details on this equation, refer to algICP_IMLP.cpp)
	//			f_shape  = Product_i[ f_vertex(v_i; s) ]
	//			f_vertex = Product_i[ (1/sqrt(2*pi)^3) * exp( -(1/2)*||si||^2)
	//			si		 = wi^T(V - mean(V))
	//			wi		 = sqrt(lambda_i)*mi
	//			mi		 = modes of variation
	//			lambda_i = mode weights
	//			V		 = mean(V) + Sum_i[ bi*mi ]
	//			bi		 = mi^T(V - mean(V))
	//			theta_i	 = distribution parameters of the match likelihood function for IMLP
	//			T_ssm(Yi)= Sum_j[ mu_i^(j)*T_ssm(vi^(j)) ]
	//			mu_i	 = vertex weights
	//			v^(1), v^(2), v^(3)= vertices on a triangle face
	//			T_ssm( mean(v)_i ) = vi = mean(v)_i + Sum_j[ sj * wj^(i) ]
	//
	//	-loglik = Sum_i[ -loglik_IMLP ] + (3/2)*N*log(2*pi) + (n_data/2) * Sum_i [ ||s||^2 ]
	//
	//	Simplified Error = Sum_i[ -loglik_IMLP ] + n_data/2 * Sum_i [ ||s||^2 ]
	//

	vctDynamicVector<vct3x3>  Mi(nSamples);           // noise covariances of match (R*Mxi*Rt + Myi)
	vctDynamicVector<vct3x3>  inv_Mi(nSamples);       // inverse noise covariances of match (R*Mxi*Rt + Myi)^-1
	vctDynamicVector<double>  det_Mi(nSamples);       // determinant of noise covariances of match |R*Mxi*Rt + Myi|
	vctDynamicVector<double>  SqrMahalDist(nSamples); // square Mahalanobis distance of matches = (yi-Rxi-t)'*inv(Mi)*(yi-Rxi-t)

	// compute mahalanobis distances of the matches
	vct3 residual;
	T_ssm(matchPts, Tssm_matchPts);
	//std::cout << "R_Mxi_Rt: " << R_Mxi_Rt.Element(50) << " Myi_sigma2: " << Myi_sigma2.Element(50) << std::endl;
	for (unsigned int s = 0; s < nSamples; s++)
	{
		residual = samplePtsXfmd.Element(s) - Tssm_matchPts.Element(s);

		// match covariance
		Mi.Element(s) = R_Mxi_Rt.Element(s) + Myi_sigma2.Element(s);
		// match covariance decomposition
		ComputeCovDecomposition_NonIter(Mi.Element(s), inv_Mi.Element(s), det_Mi.Element(s));
		// match square Mahalanobis distance
		SqrMahalDist.Element(s) = residual*inv_Mi.Element(s)*residual;
	}

	// -- Here We Compute the Full Negative Log-Likelihood -- //
	static double nklog2PI = nSamples*3.0*log(2.0*cmnPI);
	double ssmCost = 0.0;
	double logCost = 0.0;
	double expCost = 0.0;
	for (unsigned int i = 0; i<nSamples; i++)
	{
		// Compute error contribution for this sample
		ssmCost += Si.Element(i).NormSquare();
		expCost += SqrMahalDist.Element(i);
		logCost += log(det_Mi.Element(i));
	}

	prevCostFuncValue = costFuncValue;
	//std::cout << costFuncValue << std::endl;
	//std::cout << nklog2PI << " " << logCost << " " << expCost << " " << ssmCost << std::endl;
	costFuncValue = nklog2PI + (logCost + expCost + nSamples*ssmCost) / 2.0;
	//std::cout << costFuncValue << std::endl;
	
	//-- Test for algorithm-specific termination --//

	// remove last iteration from monitoring variable
	//  by bit shifting one bit to the right
	costFuncIncBits >>= 1;
	if (costFuncValue > prevCostFuncValue)
	{
		// set 4th bit in monitoring variable
		//  (since we want to monitor up to 4 iterations)
		costFuncIncBits |= 0x08;

		// signal termination if cost function increased another time within 
		//  the past 3 trials and if the value has not decreased since that time
		//  TODO: better to test if the value is not the same value as before?
		if (costFuncIncBits > 0x08 && abs(prevIncCostFuncValue - costFuncValue) < 1e-10)
		{
			bTerminateAlgorithm = true;
		}
		prevIncCostFuncValue = costFuncValue;
	}
	else
	{ // record last time the cost function was non-increasing
		Fdec = Freg;
	}
	return costFuncValue;
#else
	return 0.0;
#endif
}


void algICP_DIMLP::T_ssm(vctDynamicVector<vct3> v, vctDynamicVector<vct3> &Tv, 
	vctDynamicVector<vct3> y, vctDynamicVector<vct3> &Ty)
{
	T_ssm(v, Tv);
	// temporary assignment of eta: equal weights for each vertex on triangle
	// TODO : For each y, find triangle v's closest to it, and compute eta
	// as the ratio of importance of the v's on y
	for (unsigned int s = 0; s < nSamples; s++)
	{
		Ty.Element(s) = eta.DotProduct(Tv.Element(s));
	}
}
void algICP_DIMLP::T_ssm(vctDynamicVector<vct3> &v, vctDynamicVector<vct3> &Tv)
{
	vctDynamicVector<vct3> tmpMean(nSamples);
	vctDynamicVector<vct3> tmpSi(nSamples);
	vctDynamicVector<vct3> tmpWi(nSamples);
	Tv.SetSize(nSamples);
	for (unsigned int s = 0; s < nSamples; s++)
	{
		tmpMean.Element(s) = meanShape.Element(pTree->MeshP->faces[matchDatums.Element(s)][0]); // Change this to average the vertices
		tmpWi.Element(s) = wi.Element(pTree->MeshP->faces[matchDatums.Element(s)][0]);
	}
	tmpSi = tmpWi.DotProduct(v - tmpMean); 
	for (unsigned int s = 0; s < nSamples; s++)
	{
		Tv.Element(s) = tmpMean.Element(s) + (tmpSi.DotProduct(tmpWi));
	}
}

bool algICP_DIMLP::ICP_Terminate(vctFrm3 &F)
{
	if (bTerminateAlgorithm)
	{
		// return registration for last iteration of non-increasing cost
		F = Fdec;
		return true;
	}
	else
	{
		return false;
	}
}

vctFrm3 algICP_DIMLP::ICP_RegisterMatches()
{
	vctFrm3 dF;
#if 1
	vct6 x0(0.0);
	vct6 x;

	// x_prev must be at a different value that x0
	x_prev.SetAll(std::numeric_limits<double>::max());

	//std::cout << "[1]" << std::endl;
	x = dlib.ComputeRegistration(x0);

	//std::cout << "[2]" << std::endl;
	//update transform
	vctFixedSizeVectorRef<double, 3, 1> alpha(x, 0);
	vctFixedSizeVectorRef<double, 3, 1> t(x, 3);
	dF.Rotation() = vctRot3(vctRodRot3(alpha));
	dF.Translation() = t;
	Freg = dF * Freg;
	//std::cout << "[3]" << std::endl;
#else
	RegisterP2P_TLS(samplePtsXfmd, Tssm_matchPts, //matchPts,
		R_Mxi_Rt, Myi_sigma2, dF);
	Freg = dF*Freg;
#endif
	return Freg;
}

void algICP_DIMLP::UpdateOptimizerCalculations(const vct6 &x)
{
	//std::cout << "[1.2.1]" << std::endl;
	a.Assign(x[0], x[1], x[2]);
	t.Assign(x[3], x[4], x[5]);

	// incremental rotation matrix
	Ra = vctRot3(vctRodRot3(a));

	vctDynamicVectorRef<vct3>   X_xfm(samplePtsXfmd); 
	vctDynamicVectorRef<vct3>   Tssm_Y(Tssm_matchPts);
	vctDynamicVector<vct3x3>  inv_Mxi(nSamples);       // inverse noise covariances of match Mxi^-1
	vctDynamicVector<double>  det_Mxi(nSamples);       // determinant of noise covariances of match |Mxi|

	//std::cout << "[1.2.2]" << std::endl;
	for (unsigned int i = 0; i < nSamples; i++)
	{
		Tssm_Y_t = Tssm_Y.Element(i) - t;
		//std::cout << "[1.2.3]" << std::endl;
		//std::cout << Tssm_Y.size() << ", " <</* samplePts.size() <<*/ std::endl;
		Rat_Tssm_Y_t_x.Element(i) = Ra.TransposeRef() * Tssm_Y_t.Element(i) - samplePts.Element(i);
		//std::cout << "[1.2.4]" << std::endl;
		ComputeCovDecomposition_NonIter(Mxi.Element(i), inv_Mxi.Element(i), det_Mxi.Element(i));
		//std::cout << "[1.2.5]" << std::endl;
		invMx_Rat_Tssm_Y_t_x.Element(i) = inv_Mxi.Element(i) * Rat_Tssm_Y_t_x.Element(i);
	}
	//std::cout << "[1.2.6]" << std::endl;
	x_prev = x;
}

double algICP_DIMLP::CostFunctionValue(const vct6 &x)
{
	// don't recompute these if already computed for gradient
	if (x.NotEqual(x_prev))
	{
		UpdateOptimizerCalculations(x);
	}
	//std::cout << "[1.2.1]" << std::endl;

	double f = 0.0;
	for (unsigned int i = 0; i < nSamples; i++)
	{
		f += Rat_Tssm_Y_t_x[i] * invMx_Rat_Tssm_Y_t_x[i]
			+ nSamples*Si[i].DotProduct(Si[i]);
	}

	//std::cout << "[1.2.2]" << std::endl;
	return f;
}

void algICP_DIMLP::CostFunctionGradient(const vct6 &x, vct9 &g)
{
	vctFixedSizeVector<vctRot3, 3> dRa;  // Rodrigues Jacobians of R(a) wrt ax,ay,az

	// don't recompute these if already computed for cost function value
	if (x.NotEqual(x_prev))
	{
		UpdateOptimizerCalculations(x);
	}

	//std::cout << "[1.2.3]" << std::endl;
	ComputeRodriguesJacobians(a, dRa);

	// form the cost function gradient
	g.SetAll(0.0);
	vctFixedSizeVectorRef<double, 3, 1> ga(g, 0);
	vctFixedSizeVectorRef<double, 3, 1> gt(g, 3);
	vctFixedSizeVectorRef<double, 3, 1> gs(g, 6);
	vct3x3 Jz_a;

	for (unsigned int s = 0; s < nSamples; s++)
	{
		for (unsigned int c = 0; c < 3; c++)
		{
			Jz_a.Column(c) = dRa[c].TransposeRef() * Tssm_Y_t[s];
		}

		ga += invMx_Rat_Tssm_Y_t_x[s].Multiply(2.0) * Jz_a;
		gt -= Ra * invMx_Rat_Tssm_Y_t_x[s].Multiply(2.0);
		gs += Ra * eta*wi[s] * invMx_Rat_Tssm_Y_t_x[s].Multiply(2.0);	// Cmatch component

		gs += 2 * (double)nSamples*Si[s];								// Cshape component
	}
	//std::cout << "[1.2.4]" << std::endl;
}

unsigned int algICP_DIMLP::ICP_FilterMatches()
{
	//
	// Filter Matches for Outliers
	//
	
	return algICP_IMLP::ICP_FilterMatches();
}


// PDTree Methods

// fast check if a datum might have smaller match error than error bound
int algICP_DIMLP::DatumMightBeCloser(const vct3 &v,
	int datum,
	double ErrorBound)
{
	return true;
}

double algICP_DIMLP::FindClosestPointOnDatum(const vct3 &v, vct3 &closest, int datum)
{
	// first iteration variables that don't change
	static const vct3x3 I_7071(vct3x3::Eye()*0.7071); // sqrt(1/2)*I
	static const vct3x3 I1_4142(vct3x3::Eye()*1.4142); // sqrt(2)*I
	static const vct3x3 I2(vct3x3::Eye()*2.0); // 2*I
	static const vct3x3 I_5(vct3x3::Eye()*0.5); // 0.5*I

	static vct3 d;
	static vct3x3 M, Minv, N, Ninv;
	double det_M;

	if (bFirstIter_Matches)
	{ // isotropic noise model for first iteration

		// M = Mx + NodeEigMax*I = 2*I = V*S*V'  =>  S = 2*I, V = I
		// Minv = N'*N = I*(1/2)*I  =>  N = sqrt(1/2)*I = 0.7071*I
		M = I2;
		Minv = I_5;
		N = I_7071;
		Ninv = I1_4142;
		det_M = 8;
	}
	else
	{ // Use noise model of node

		// compute noise model for this datum
		//  M = R*Mxi*R' + Myi
		M = sample_RMxRt_sigma2 + pMesh->TriangleCov[datum];
		ComputeCovDecomposition_NonIter(M, Minv, N, Ninv, det_M);
	}


#if 1
	// Find the closest point on this triangle in a Mahalanobis distance sense
	TCPS.FindMostLikelyPointOnTriangle(v, datum, N, Ninv, closest);
#else
	// Find the closest point on this triangle in a Euclidean distance sense
	//   Euclidean Distance:   ||x-v||
	TCPS.FindClosestPointOnTriangle(v, datum, closest);
#endif

	// return match error 
	//  log term plus square Mahalanobis distance
	d = (v - closest);

	//BoundingBox BB;
	//pTree->EnlargeBounds(Freg, datum, BB);

	return log(det_M) + vctDotProduct(d, Minv*d);
}

// fast check if a node might contain a datum having smaller match error
//  than the error bound
int algICP_DIMLP::NodeMightBeCloser(const vct3 &v,
	PDTreeNode *node,
	double ErrorBound)
{
	return algICP_IMLP::NodeMightBeCloser(v, node, ErrorBound);
}
