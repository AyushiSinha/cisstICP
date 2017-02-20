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
	double outlierChiSquareThreshold,
	double sigma2Max)
	: algICP_IMLP(pTree, samplePts, sampleCov, sampleMsmtCov, outlierChiSquareThreshold, sigma2Max),
	dlib(this),
	pTree(pTree),
	pMesh(pTree->MeshP),
	TCPS(*(pTree->MeshP))
{
	SetSamples(samplePts, sampleCov, sampleMsmtCov, meanShape);
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

	std::cout << "\nAverage Mahalanobis Distance = " << Avg << std::endl;
}

void algICP_DIMLP::SetSamples(
	vctDynamicVector<vct3> &argSamplePts,
	vctDynamicVector<vct3x3> &argMxi,
	vctDynamicVector<vct3x3> &argMsmtMxi,
	vctDynamicVector<vct3> &argMeanShape)
{
	//if (pMesh->mode.size() != nSamples)
	//{
	//	std::cout << "ERROR: number of sample mode matrices (" << pMesh->mode.size() << ") does not match number of samples (" << nSamples << ")" << std::endl;
	//}

	// base class
	algICP_IMLP::SetSamples(argSamplePts, argMxi, argMsmtMxi);

	meanShape = argMeanShape;
	nModes = pTree->MeshP->modeWeight.size();

	Si = pTree->MeshP->Si;
	wi = pTree->MeshP->wi;

	//Tssm_wi.resize(nSamples);
	Tssm_wi.resize(nModes);
	for (int i = 0; i < Tssm_wi.size(); i++)
		Tssm_wi[i].resize(nSamples);
	Tssm_matchPts.resize(nSamples);

	Tssm_Y_t.resize(nSamples);
	Rat_Tssm_Y_t_x.resize(nSamples);
	Rat_Tssm_Y_t_x_invMx.resize(nSamples);

	x_prev.SetSize(6+nModes); // 3 for rotation, 3 for translation, and n modes
	s.SetSize(nModes);
	
	eta = vct3(1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0);
}

void algICP_DIMLP::ICP_UpdateParameters_PostMatch()
{
	// base class
	algICP::ICP_UpdateParameters_PostMatch();
	T_ssm();

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

	// update the match uncertainty factor
	sigma2 = sumSqrDist_Inliers / (nSamples - nOutliers);

	// apply max threshold
	if (sigma2 > sigma2Max)
	{
		sigma2 = sigma2Max;
	}

	// update noise models of the matches
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
}

void algICP_DIMLP::ICP_ComputeMatches()
{
	///*** At the beginning of each correspondence phase, the positions 
	///*** of the representing the model shape must be recomputed based 
	///*** on the current values of the model-shape parameters, s 

	//static int itermesh = 0;
	//for (unsigned int i = 0; i < nModes; i++)
	//	if (bFirstIter_Matches)
	//		Si[i] = Si[i] / pTree->MeshP->modeWeight[i];

	// deformably transform each mesh vertex
	//pTree->MeshP->vertices = pMesh->meanShape;
	for (int s = 0; s < pMesh->NumVertices(); s++) 
		for (unsigned int i = 0; i < nModes; i++) 
			pTree->MeshP->vertices(s) += (Si[i] * wi[i].Element(s) / pTree->MeshP->modeWeight[i]);

		
	std::cout << "Si = " << Si[0] << " ";
	pMesh->Si = Si;
	//std::string imesh = std::to_string(itermesh) + ".ply";
	//cisstMesh currMesh;
	//currMesh.vertices = pTree->MeshP->vertices;
	//currMesh.faces = pTree->MeshP->faces;
	//currMesh.SavePLY(imesh);
	//itermesh++;

	// TODO: fix this so you're updating for all vertices,
	// TODO: and recursively for all parent bounding boxes
	for (unsigned int s = 0; s < pMesh->NumTriangles(); s++)
		pTree->EnlargeBounds(Freg, s, pTree->Bounds);
	
	// base class
	algICP_IMLP::ICP_ComputeMatches();
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
		expCost += SqrMahalDist.Element(i);
		logCost += log(det_Mi.Element(i));
	}
	ssmCost += Si.NormSquare();

	prevCostFuncValue = costFuncValue;
	costFuncValue = nklog2PI + (logCost + expCost + nSamples*ssmCost) / 2.0;
	//costFuncValue = (nklog2PI + logCost + expCost) / 2.0; // IMLP cost func
	
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

void algICP_DIMLP::T_ssm()
{
	vctDynamicVector<vct3> tmpMean(nSamples);
	vct3 tmpv0, tmpv1, tmpv2;
	vct3 tmpw0, tmpw1, tmpw2;
	vct3 tmpm0, tmpm1, tmpm2;
	vct3 tmpv0_matchPt, 
		tmpv1_matchPt, 
		tmpv2_matchPt;
	double tmpv0_norm, tmpv1_norm, tmpv2_norm;
	vctDynamicVector<double> tmp_si;
	tmp_si.resize(nModes);
	tmp_si.SetAll(0.0);

	//std::cout << "Before: " << matchPts[0] << std::endl;
	//std::cout << "Before: " << Si << std::endl;
	for (unsigned int s = 0; s < nSamples; s++)
	{
		// Find the 3 vertices of the triangle that the matchPoint lies on,
		// i.e., of the matchDatum on estimated shape 
		algICP_DIMLP::pTree->MeshP->FaceCoords(matchDatums[s], tmpv0, tmpv1, tmpv2);
		vctInt3 f = pTree->MeshP->faces[matchDatums.Element(s)];
		//tmpv0 = pMesh->vertices.Element(f[0]);
		//tmpv1 = pMesh->vertices.Element(f[1]);
		//tmpv2 = pMesh->vertices.Element(f[2]);

		// find distance between triangle vertices and matchPts
		tmpv0_matchPt = tmpv0 - matchPts.Element(s);
		tmpv1_matchPt = tmpv1 - matchPts.Element(s);
		tmpv2_matchPt = tmpv2 - matchPts.Element(s);

		double areaTri = vctCrossProduct(tmpv1 - tmpv0, tmpv2 - tmpv0).Norm();

		tmpv0_norm = vctCrossProduct(tmpv1_matchPt, tmpv2_matchPt).Norm() / areaTri;
		tmpv1_norm = vctCrossProduct(tmpv2_matchPt, tmpv0_matchPt).Norm() / areaTri;
		tmpv2_norm = vctCrossProduct(tmpv0_matchPt, tmpv1_matchPt).Norm() / areaTri;
		// Why aren't the barycentric coordinates summing to 1?
		//std::cout << eta[0] << "+" << eta[1] << "+" << eta[2] << "=" << eta[0] + eta[1] + eta[2] << std::endl;

		//// compute norm
		//tmpv0_norm = tmpv0_matchPt.Norm();
		//tmpv1_norm = tmpv1_matchPt.Norm();
		//tmpv2_norm = tmpv2_matchPt.Norm();

		// ensure sum of eta elements = 1
		eta.Element(0) = tmpv0_norm / (tmpv0_norm + tmpv1_norm + tmpv2_norm); // j = 1
		eta.Element(1) = tmpv1_norm / (tmpv0_norm + tmpv1_norm + tmpv2_norm); // j = 2
		eta.Element(2) = tmpv2_norm / (tmpv0_norm + tmpv1_norm + tmpv2_norm); // j = 3

		tmpm0 = meanShape.Element(f[0]); // j = 1
		tmpm1 = meanShape.Element(f[1]); // j = 2
		tmpm2 = meanShape.Element(f[2]); // j = 3
		tmpMean.Element(s) = eta.Element(0)*tmpm0 + eta.Element(1)*tmpm1 + eta.Element(2)*tmpm2; 
		
		tmpv0 = tmpm0;
		tmpv1 = tmpm1;
		tmpv2 = tmpm2;

		for (unsigned int i = 0; i < nModes; i++)
		{
			double curr_si = Si[i];
			double curr_modeweight = pTree->MeshP->modeWeight[i];

			tmpw0 = wi[i].Element(f[0]); // j = 1
			tmpw1 = wi[i].Element(f[1]); // j = 2
			tmpw2 = wi[i].Element(f[2]); // j = 3

			// TODO: Update the mean mesh not the current mesh
			tmpv0 += curr_si * tmpw0 /*/ curr_modeweight*/; // j = 1
			tmpv1 += curr_si * tmpw1 /*/ curr_modeweight*/; // j = 2
			tmpv2 += curr_si * tmpw2 /*/ curr_modeweight*/; // j = 3

			Tssm_wi[i][s] = eta.Element(0)*tmpw0 + eta.Element(1)*tmpw1 + eta.Element(2)*tmpw2;
		}

		Tssm_matchPts.Element(s) = eta.Element(0)*tmpv0 + eta.Element(1)*tmpv1 + eta.Element(2)*tmpv2;
		//matchPts[s] = Tssm_matchPts[s];

		for (unsigned int i = 0; i < nModes; i++)
			tmp_si[i] += Tssm_wi[i][s].DotProduct(Tssm_matchPts[s] - tmpMean[s]);
	}

	Si = tmp_si;
	//pTree->MeshP->Si = Si;
	//pTree->MeshP->wi = wi;

	//std::cout << "After: " << matchPts[0] << std::endl;
	//std::cout << "After: " << Si << std::endl;
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
	vctFrm3 F;
#if 1
	//vct7 x0(0.0);
	//vct7 x;
	vctDynamicVector<double> x0;
	vctDynamicVector<double> x;

	x0.SetSize(6+nModes);
	x.SetSize(6+nModes);
	
	x0.SetAll(0.0); 
	for (unsigned int i = 0; i < nModes; i++)
		x0[6+i] = Si[i];

	// x_prev must be at a different value that x0
	x_prev.SetAll(std::numeric_limits<double>::max());

	//printf("\nComputing registration... ");
	x = dlib.ComputeRegistration(x0);
	//printf("Done\n");

	//update transform
	vctFixedSizeVectorRef<double, 3, 1> alpha(x, 0);
	vctFixedSizeVectorRef<double, 3, 1> t(x, 3);
	vctDynamicVectorRef<double> s(x, 6, nModes);
	F.Rotation() = vctRot3(vctRodRot3(alpha));
	F.Translation() = t;
	Freg = F;
	Si = s;

	//std::cout << "Old = " << x0 << std::endl;
	//std::cout << "New = " << x << std::endl;
	//std::cout << "New Si = " << Si << std::endl;
#else
	RegisterP2P_TLS(samplePtsXfmd, Tssm_matchPts, //matchPts,
		R_Mxi_Rt, Myi_sigma2, F);
	Freg = F*Freg;
#endif
	return Freg;
}

//void algICP_DIMLP::UpdateOptimizerCalculations(const vct7 &x)
void algICP_DIMLP::UpdateOptimizerCalculations(const vctDynamicVector<double> &x)
{
	a.Assign(x[0], x[1], x[2]);
	t.Assign(x[3], x[4], x[5]);
	
	//s.Assign(x[6]);
	for (unsigned int i = 0; i < nModes; i++)
		s[i] = x[6 + i];

	// Rodrigues formulation
	Ra = vctRot3(vctRodRot3(a));

	vctDynamicVectorRef<vct3>   X_xfm(samplePtsXfmd); 
	vctDynamicVectorRef<vct3>   Tssm_Y(Tssm_matchPts);
	vctDynamicVector<vct3x3>  inv_Mxi(nSamples);       // inverse noise covariances of match Mxi^-1
	vctDynamicVector<double>  det_Mxi(nSamples);       // determinant of noise covariances of match |Mxi|

	for (unsigned int i = 0; i < nSamples; i++)
	{
		Tssm_Y_t = Tssm_Y.Element(i) - t;
		Rat_Tssm_Y_t_x.Element(i) = Ra.Transpose() * Tssm_Y_t.Element(i) - X_xfm.Element(i); 
		ComputeCovDecomposition_NonIter(Mxi.Element(i), inv_Mxi.Element(i), det_Mxi.Element(i));
		Rat_Tssm_Y_t_x_invMx.Element(i) = Rat_Tssm_Y_t_x.Element(i)*inv_Mxi.Element(i);
	}
	x_prev = x;
}

//double algICP_DIMLP::CostFunctionValue(const vct7 &x)
double algICP_DIMLP::CostFunctionValue(const vctDynamicVector<double> &x)
{
	// don't recompute these if already computed for gradient
	if (x.NotEqual(x_prev))
	{
		UpdateOptimizerCalculations(x);
	}

	double f = 0.0;
	for (unsigned int i = 0; i < nSamples; i++)
	{
		f += Rat_Tssm_Y_t_x_invMx.Element(i)*Rat_Tssm_Y_t_x.Element(i);
	}

	//f += nSamples*(Si.DotProduct(Si));
	f += nSamples*(s.DotProduct(s));

	return f;
}

//void algICP_DIMLP::CostFunctionGradient(const vct7 &x, vct7 &g)
void algICP_DIMLP::CostFunctionGradient(const vctDynamicVector<double> &x, vctDynamicVector<double> &g)
{
	vctFixedSizeVector<vctRot3, 3> dRa;  // Rodrigues Jacobians of R(a) wrt ax,ay,az

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
	vctDynamicVectorRef<double> gs(g, 6, nModes);

	vct3x3 Jz_a;

	for (unsigned int j = 0; j < nSamples; j++)
	{
		for (unsigned int c = 0; c < 3; c++)
			Jz_a.Column(c) = dRa[c].TransposeRef() * Tssm_Y_t[j];

		ga += Rat_Tssm_Y_t_x_invMx[j] * 2.0 * Jz_a;
		gt += Rat_Tssm_Y_t_x_invMx[j] * 2.0 * (-Ra.Transpose());

		for (unsigned int i = 0; i < nModes; i++)
			gs[i] += Rat_Tssm_Y_t_x_invMx[j] * 2.0 * (Ra.Transpose() * Tssm_wi[i][j]);	// Cmatch component
		
	}

	for (unsigned int i = 0; i < nModes; i++)
		gs[i] += 2.0 * (double)nSamples*s[i];	// Cshape component
		//gs[i] += 2.0 * (double)nSamples*Si[i];	// Cshape component
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

	return log(det_M) + vctDotProduct(d, Minv*d);
	//return vctDotProduct(d, Minv*d);   // rather than: log(det_M) + vctDotProduct(d,Minv*d);
}

// fast check if a node might contain a datum having smaller match error
//  than the error bound
//int algICP_DIMLP::NodeMightBeCloser(const vct3 &v,
//	PDTreeNode *node,
//	double ErrorBound)
//{
//	return algICP_IMLP::NodeMightBeCloser(v, node, ErrorBound);
//}
