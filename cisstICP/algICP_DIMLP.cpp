// ****************************************************************************
//
//    Copyright (c) 2017, Ayushi Sinha, Seth Billings, Russell Taylor, Johns Hopkins University. 
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
//#include <omp.h>

#define COMPUTE_ERROR_FUNCTION
#define ENABLE_PARALLELIZATION

algICP_DIMLP::algICP_DIMLP(
	PDTree_Mesh *pTree,
	vctDynamicVector<vct3> &samplePts,
	vctDynamicVector<vct3x3> &sampleCov,
	vctDynamicVector<vct3x3> &sampleMsmtCov,
	vctDynamicVector<vct3> &meanShape,
	double scale, bool bScale,
	double outlierChiSquareThreshold,
	double sigma2Max)
	: algICP_IMLP(pTree, samplePts, sampleCov, sampleMsmtCov, outlierChiSquareThreshold, sigma2Max),
	dlib(this),
	pTree(pTree),
	pMesh(pTree->MeshP),
	TCPS(*(pTree->MeshP))
{
	SetSamples(samplePts, sampleCov, sampleMsmtCov, meanShape, scale, bScale);
}

void algICP_DIMLP::SetConstraints(double argRotbounds,
									double argTransbounds,
									double argScalebounds,
									double argSPbounds)
{
	rb = argRotbounds;
	tb = argTransbounds;
	sb = argScalebounds;
	spb = argSPbounds;
}

void algICP_DIMLP::ComputeMatchStatistics(double &Avg, double &StdDev)
{
	// return the average mahalanobis distance of the matches
	//  based on point noise models only (measurement and surface model covariances)
	//  i.e. do not include sigma2
	// AS comment: maybe let this call default to the IMLP call?

	double sumSqrMatchDist = 0.0;
	double sumMatchDist = 0.0;
	double sqrMatchDist;

	totalSumSqrMahalDist = 0.0;
	sumSqrMahalDist = 0.0;
	double sumMahalDist = 0.0;
	double sqrMahalDist;
	nGoodSamples = 0;

	vct3x3 M, Minv;
	vct3 residual;

	for (unsigned int i = 0; i < nSamples; i++)
	{
		residual = Tssm_Y[i] - (Freg * samplePts[i]) * sc;
		M = Freg.Rotation() * Mxi[i] * Freg.Rotation().Transpose() + *Myi[i];
		ComputeCovInverse_NonIter(M, Minv);

		sqrMahalDist = residual*Minv*residual;
		totalSumSqrMahalDist += sqrMahalDist;

		if (outlierFlags[i])	continue;	// skip outliers

		sumSqrMahalDist += sqrMahalDist;
		sumMahalDist += sqrt(sqrMahalDist);

		sqrMatchDist = residual.NormSquare();
		sumSqrMatchDist += sqrMatchDist;
		sumMatchDist += sqrt(sqrMatchDist);
		nGoodSamples++;
	}

	Avg = sumMahalDist / nGoodSamples;
	StdDev = sqrt( (sumSqrMahalDist / nGoodSamples) + Avg*Avg );

	//std::cout << "\nFinal Scale = " << sc << std::endl;
	//std::cout << "\nAverage Match Distance = " << sumMatchDist / nGoodSamples << std::endl;
	//std::cout << "\nAverage Mahalanobis Distance = " << Avg << " (+/-" << StdDev << ")" << std::endl;
}

void algICP_DIMLP::PrintMatchStatistics(std::stringstream &tMsg)
{
	// For registration rejection purpose:
	tMsg << "\nSum square mahalanobis distance = " << totalSumSqrMahalDist << " over " << nSamples << " samples";
	tMsg << "\nSum square mahalanobis distance = " << sumSqrMahalDist << " over " << nGoodSamples << " inliers\n";
}

void algICP_DIMLP::SetSamples(
	vctDynamicVector<vct3> &argSamplePts,
	vctDynamicVector<vct3x3> &argMxi,
	vctDynamicVector<vct3x3> &argMsmtMxi,
	vctDynamicVector<vct3> &argMeanShape,
	double argScale, bool argbScale)
{
	// scale sample points (remove this from here when you move it to IMLP)
	sc = argScale;
	bScale = argbScale;
	if (bScale) {
		for (unsigned int i = 0; i < nSamples; i++)
			samplePts[i] = samplePts[i] * sc;
		nTrans = 7; // 3 for rotation, 3 for translation, 1 for scale
	}
	else
		nTrans = 6; // 3 for rotation, 3 for translation

	meanShape = argMeanShape;
	nModes = pTree->MeshP->modeWeight.size();

	Si = pTree->MeshP->Si;
	wi = pTree->MeshP->wi;

	Tssm_wi.resize(nModes);
	for (int i = 0; i < Tssm_wi.size(); i++)
		Tssm_wi[i].resize(nSamples);
	Tssm_matchPts.resize(nSamples);

	Tssm_Y.resize(nSamples);
	Tssm_Y_t.resize(nSamples);
	Rat_Tssm_Y_t_x.resize(nSamples);
	Rat_Tssm_Y_t_x_invMx.resize(nSamples);

	x_prev.SetSize(nTrans + nModes);  // transformation parameters, and n modes
	mu.SetSize(nSamples);
	f.SetSize(nSamples);
	s.SetSize(nModes);
}

void algICP_DIMLP::ICP_InitializeParameters(vctFrm3 &FGuess)
{
	// initialize base class
	algICP_IMLP::ICP_InitializeParameters(FGuess);
	this->FGuess = FGuess;

	// set x_prev to FGuess for Rotation and Translation 
	// components, and zero for shape components
	vct3 rot = vctRodRot3(FGuess.Rotation());
	vct3 trans = FGuess.Translation();
	double scale = sc;
	x_prev.SetAll(0.0); 
	for (int i = 0; i < 3; i++)
		x_prev[i] = rot[i];

	for (int i = 0; i < 3; i++)
		x_prev[3 + i] = trans[i];

	if (bScale)
		x_prev[6] = scale;

	for (unsigned int i = 0; i < nModes; i++)
		x_prev[nTrans + i] = Si[i];

	mu.SetAll(vct3(1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0));
}

void algICP_DIMLP::ICP_UpdateParameters_PostMatch()
{
	// base class
	algICP::ICP_UpdateParameters_PostMatch();

	// compute sum of square distances of inliers
	sumSqrDist_Inliers = 0.0;
	for (unsigned int s = 0; s < nSamples; s++)
	{
		residuals_PostMatch.Element(s) = samplePtsXfmd.Element(s) - Tssm_Y.Element(s); 
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

#ifdef DEBUG_DIMLP
	std::cout << "My0:" << std::endl << *Myi[0] << std::endl;
	std::cout << "My1:" << std::endl << *Myi[1] << std::endl;
#endif

	bFirstIter_Matches = false;
}

void algICP_DIMLP::ICP_UpdateParameters_PostRegister(vctFrm3 &Freg)
{
	// base class
	algICP_IMLP::ICP_UpdateParameters_PostRegister(Freg);
	if (bScale)
		for (unsigned int s = 0; s < nSamples; s++)
			samplePtsXfmd.Element(s) = sc * samplePtsXfmd.Element(s); // move this to IMLP also later

	UpdateShape(Si);
	UpdateTree();

	// Re-initialize to compute matches on updated mesh
	TCPS.init(pTree->MeshP->vertices, pTree->MeshP->faces);

#if 0
	static int count = 1;
	cisstMesh currentSamples;
	currentSamples.vertices.SetSize(nSamples);
	for (unsigned int i = 0; i < nSamples; i++)
		currentSamples.vertices[i] = samplePtsXfmd[i];
	currentSamples.SavePLY("currentSamples" + std::to_string(count) + ".ply");
	
	cisstMesh currentMesh;
	currentMesh.vertices = pMesh->vertices;
	currentMesh.faces = pMesh->faces;
	currentMesh.SavePLY("currentMesh" + std::to_string(count) + ".ply");

	count++;
#endif
}

void algICP_DIMLP::UpdateShape(vctDynamicVector<double>	&S)
{
	//static int itermesh = 0;

	// deformably transform each mesh vertex
	pTree->MeshP->vertices = meanShape;
	int s;
#ifdef ENABLE_PARALLELIZATION
#pragma omp parallel for
#endif
	for (s = 0; s < pMesh->NumVertices(); s++)
		for (unsigned int i = 0; i < nModes; i++)
			pTree->MeshP->vertices(s) += (S[i] * wi[i].Element(s));
}

void algICP_DIMLP::UpdateTree()
{
	vctFrm3 FId;
	FId.Assign(vctFrm3::Identity());

	pTree->EnlargeBounds(FId);
}

void algICP_DIMLP::ICP_ComputeMatches()
{
	//
	// At the beginning of each correspondence phase, the positions 
	// of the representing the model shape must be recomputed based 
	// on the current values of the model-shape parameters, s 
	//

	// base class
	algICP_IMLP::ICP_ComputeMatches();

	Tssm_Y = matchPts;
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
	//	-loglik = Sum_i[ -loglik_IMLP ] + (3/2)*N*log(2*pi) + (1/2) * Sum_i [ ||s||^2 ]
	//
	//	Simplified Error = Sum_i[ -loglik_IMLP ] + 1/2 * Sum_i [ ||s||^2 ]
	//

	vctDynamicVector<vct3x3>  Mi(nSamples);           // noise covariances of match (R*Mxi*Rt + Myi)
	vctDynamicVector<vct3x3>  inv_Mi(nSamples);       // inverse noise covariances of match (R*Mxi*Rt + Myi)^-1
	vctDynamicVector<double>  det_Mi(nSamples);       // determinant of noise covariances of match |R*Mxi*Rt + Myi|
	vctDynamicVector<double>  SqrMahalDist(nSamples); // square Mahalanobis distance of matches = (yi-Rxi-t)'*inv(Mi)*(yi-Rxi-t)

	vct3 residual;
	double nklog2PI = 0.0;
	double ssmCost	= 0.0;
	double logCost	= 0.0;
	double expCost	= 0.0;

	// compute mahalanobis distances of the matches
	for (unsigned int s = 0; s < nSamples; s++)
	{
		if (outlierFlags[s])	continue;	// skip outliers

		// Compute errors using the current match point on the 
		// deformed shape after computing updating Si
		residual = samplePtsXfmd.Element(s) - Tssm_Y.Element(s);  
		// match covariance
		Mi.Element(s) = R_Mxi_Rt.Element(s) + Myi_sigma2.Element(s);
		// match covariance decomposition
		ComputeCovDecomposition_NonIter(Mi.Element(s), inv_Mi.Element(s), det_Mi.Element(s));
		// match square Mahalanobis distance
		SqrMahalDist.Element(s) = residual*inv_Mi.Element(s)*residual;

		// -- Here We Compute the Full Negative Log-Likelihood -- //
		// Compute error contribution for this sample
		nklog2PI += 3.0*log(2.0*cmnPI);
		expCost += SqrMahalDist.Element(s);
		logCost += log(det_Mi.Element(s));
	}
	ssmCost += Si.NormSquare();

	prevCostFuncValue = costFuncValue;
	costFuncValue = nklog2PI + (logCost + expCost + ssmCost) / 2.0;
	
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

// Compute Mu for deformed shape after each optimization
void algICP_DIMLP::ComputeMu()
{
	vct3	v0, v1, v2;
	vct3	v0_Tssm_Y, 
			v1_Tssm_Y,	
			v2_Tssm_Y;

	for (unsigned int s = 0; s < nSamples; s++)
	{
		// Find the 3 vertices of the triangle that the matchPoint lies on,
		// i.e., of the matchDatum on estimated shape 
		f[s] = pTree->MeshP->faces[matchDatums.Element(s)];
		v0 = pTree->MeshP->vertices.Element(f[s][0]);
		v1 = pTree->MeshP->vertices.Element(f[s][1]);
		v2 = pTree->MeshP->vertices.Element(f[s][2]);

		// find distance between triangle vertices and matchPts
		v0_Tssm_Y = v0 - Tssm_Y.Element(s);
		v1_Tssm_Y = v1 - Tssm_Y.Element(s);
		v2_Tssm_Y = v2 - Tssm_Y.Element(s);

		double areaTri = vctCrossProduct(v1 - v0, v2 - v0).Norm(); 

		mu[s].Element(0) = vctCrossProduct(v1_Tssm_Y, v2_Tssm_Y).Norm() / areaTri;
		mu[s].Element(1) = vctCrossProduct(v2_Tssm_Y, v0_Tssm_Y).Norm() / areaTri;
		mu[s].Element(2) = vctCrossProduct(v0_Tssm_Y, v1_Tssm_Y).Norm() / areaTri;
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

void algICP_DIMLP::ReturnScale(double &scale)
{
	scale = sc;
}

void algICP_DIMLP::ReturnShapeParam(vctDynamicVector<double> &shapeParam)
{
	shapeParam = Si;
}

vctFrm3 algICP_DIMLP::ICP_RegisterMatches()
{
	vctFrm3 F;
	ComputeMu();

	vctDynamicVector<double> x0;
	vctDynamicVector<double> x;

	x0.SetSize(nTrans + nModes);
	x.SetSize(nTrans + nModes);

	x0 = x_prev; 

	// x_prev must be at a different value than x0
	x_prev.SetAll(std::numeric_limits<double>::max());

	for (unsigned int j = 0; j < nSamples; j++)
	{
		for (unsigned int i = 0; i < nModes; i++)
			Tssm_wi[i][j]	= mu[j][0] * wi[i][f[j][0]]
							+ mu[j][1] * wi[i][f[j][1]]
							+ mu[j][2] * wi[i][f[j][2]];
	}

	x = dlib.ComputeRegistration(x0);

	// update transform
	vctFixedSizeVectorRef<double, 3, 1> alpha(x, 0);
	vctFixedSizeVectorRef<double, 3, 1> t(x, 3); 
	double scale;
	if (bScale)
		scale = x[6];
	vctDynamicVectorRef<double> s(x, nTrans, nModes);
	F.Rotation() = vctRot3(vctRodRot3(alpha));
	F.Translation() = t;
	Freg = F;
	if (bScale)
		sc = scale;
	Si = s;

	pMesh->Si = Si;

	return Freg;
}

void algICP_DIMLP::UpdateOptimizerCalculations(const vctDynamicVector<double> &x)
{
	a.Assign(x[0], x[1], x[2]);
	t.Assign(x[3], x[4], x[5]);
	if (bScale)
		sc = x[6];
	
	for (unsigned int i = 0; i < nModes; i++)
		s[i] = x[nTrans + i];

	// Rodrigues formulation
	Ra = vctRot3(vctRodRot3(a));

	X = samplePts;
	vctDynamicVectorRef<vct3>   Mu(mu);

	vctDynamicVector<vct3x3>  inv_Mxi(nSamples);       // inverse noise covariances of match Mxi^-1
	vctDynamicVector<double>  det_Mxi(nSamples);       // determinant of noise covariances of match |Mxi|

	// Update shape based on current s (and wi and meanshape)
	// Compute Tssm_Y based on current Mu and shape
	UpdateShape(s);
	unsigned int j;
#ifdef ENABLE_PARALLELIZATION
#pragma omp parallel for
#endif
	for (j = 0; j < nSamples; j++)
	{
		Tssm_Y.Element(j)	= Mu[j][0] * pMesh->vertices[f[j][0]] 
							+ Mu[j][1] * pMesh->vertices[f[j][1]] 
							+ Mu[j][2] * pMesh->vertices[f[j][2]];

		Tssm_Y_t.Element(j) = Tssm_Y.Element(j) - t;
		Rat_Tssm_Y_t_x.Element(j) = Ra.Transpose() * Tssm_Y_t.Element(j) - sc * X.Element(j); 
		ComputeCovDecomposition_NonIter(Mxi.Element(j), inv_Mxi.Element(j), det_Mxi.Element(j));
		Rat_Tssm_Y_t_x_invMx.Element(j) = Rat_Tssm_Y_t_x.Element(j) * inv_Mxi.Element(j);
	}
	x_prev = x;
}

double algICP_DIMLP::CostFunctionValue(const vctDynamicVector<double> &x)
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
	for (i = 0; i < nSamples; i++)
	{
		if (outlierFlags[i])	continue;

		f += ( Rat_Tssm_Y_t_x_invMx.Element(i) * Rat_Tssm_Y_t_x.Element(i) ) / 2.0;
	}

#ifdef NOREGULARIZER
	f += 0;
#else
	f += (s.DotProduct(s)) / 2.0; 
#endif

	return f;
}

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
	vctFixedSizeVectorRef<double, 1, 1> gsc;
	if (bScale)
		gsc.SetRef(g, 6);
	vctDynamicVectorRef<double> gs(g, nTrans, nModes);

	vct3x3 Jz_a;

	unsigned int j;
#ifdef ENABLE_PARALLELIZATION
#pragma omp parallel for
#endif
	for (j = 0; j < nSamples; j++)
	{
		if (outlierFlags[j])	continue;

		for (unsigned int c = 0; c < 3; c++)
			Jz_a.Column(c) = dRa[c].TransposeRef() * Tssm_Y_t[j]; // check this computation

		ga += Rat_Tssm_Y_t_x_invMx[j] * Jz_a;
		gt += Rat_Tssm_Y_t_x_invMx[j] * (-Ra.Transpose());
		if (bScale)
			gsc += Rat_Tssm_Y_t_x_invMx[j] * (-X.Element(j));

		for (unsigned int i = 0; i < nModes; i++)
			gs[i] += Rat_Tssm_Y_t_x_invMx[j] * (Ra.Transpose() * Tssm_wi[i][j]);	// Cmatch component	
	}

#ifdef NOREGULARIZER
	gs += 0; 
#else
	gs +=  s;	// Cshape component 
#endif
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

	return vctDotProduct(d, Minv*d);   // rather than: log(det_M) + vctDotProduct(d,Minv*d);
}

