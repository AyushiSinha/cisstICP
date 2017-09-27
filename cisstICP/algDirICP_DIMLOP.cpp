// ****************************************************************************
//
//    Copyright (c) 2017, Ayushi Sinha, Russell Taylor, Johns Hopkins University
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

#include "algDirICP_DIMLOP.h"
#include "DirPDTreeNode.h"
#include "RegisterP2P.h"
#include "utilities.h"

#define EPS  1e-12

//void  algDirICP_DIMLOP::SetNoiseModel(
//  double initK, double initSigma2, double w_Rpos, bool dynamicallyEstParams)
//{
//  k_init = initK;
//  sigma2_init = initSigma2;
//  wRpos = w_Rpos;
//  dynamicParamEst = dynamicallyEstParams;
//}

algDirICP_DIMLOP::algDirICP_DIMLOP(
	//DirPDTreeBase *pDirTree,
	DirPDTree_Mesh *pDirTree,
	vctDynamicVector<vct3> &samplePts,
	vctDynamicVector<vct3> &sampleNorms,
	vctDynamicVector<vct3x3> &sampleCov,
	vctDynamicVector<vct3x3> &sampleMsmtCov,
	vctDynamicVector<vct3> &meanShape,
	double kinit, double sigma2init, double wRpos,
	double kfactor, double scale,
	bool dynamicallyEstParams,
	bool bScale)
   : algDirICP_IMLOP(pDirTree, samplePts, sampleNorms, kinit, sigma2init, wRpos, kfactor, dynamicallyEstParams),
   //algDirPDTree(pDirTree),
	dlib(this),
	pDirTree(pDirTree),
	pMesh(&pDirTree->mesh),
	TCPS(pDirTree->mesh)
{
   // Ensure SetSamples function of this derived class gets called
   SetSamples(samplePts, sampleNorms, sampleCov, sampleMsmtCov, meanShape, scale, bScale);
}

void algDirICP_DIMLOP::SetConstraints(double argSPbounds)
{
	spb = argSPbounds;
}

void algDirICP_DIMLOP::ComputeMatchStatistics(double &Avg, double &StdDev)
{
	double sumSqrMatchDist = 0.0;
	double sumMatchDist = 0.0;
	double sqrMatchDist;

	double sumSqrMahalDist = 0.0;
	double sumMahalDist = 0.0;
	double sqrMahalDist;
	int nGoodSamples = 0;

	vct3x3 M, Minv;
	vct3 residual;

	// NOTE: if using a method with outlier rejection, it may be desirable to
	//       compute statistics on only the inliers
	for (unsigned int i = 0; i < nSamples; i++)
	{
		if (outlierFlags[i]) continue;  // skip outliers

		residual = matchPts[i] - (Freg * samplePts[i]) * sc;
		M = Freg.Rotation() * Mxi[i] * Freg.Rotation().Transpose();// +*Myi[i];
		ComputeCovInverse_NonIter(M, Minv);

		sqrMahalDist = residual*Minv*residual;
		sumSqrMahalDist += sqrMahalDist;
		sumMahalDist += sqrt(sqrMahalDist);

		sqrMatchDist = residual.NormSquare();
		sumSqrMatchDist += sqrMatchDist;
		sumMatchDist += sqrt(sqrMatchDist);
		nGoodSamples++;
	}

	Avg = sumMahalDist / nGoodSamples;
	StdDev = sqrt( (sumSqrMahalDist / nGoodSamples) + Avg*Avg );

	std::cout << "\nFinal Scale = " << sc << std::endl;
	std::cout << "\n# good samples = " << nGoodSamples << std::endl;
	std::cout << "\nAverage Match Distance = " << sumMatchDist / nGoodSamples << std::endl;
	std::cout << "\nAverage Mahalanobis Distance = " << Avg << "(+/-" << StdDev << ")" << std::endl;

	//Avg = sumMatchDist / nSamples;
	//StdDev = (sumSqrMatchDist / nSamples) + Avg*Avg;
}

//void algDirICP_DIMLOP::ComputeMatchStatistics(
//	double &PosAvg, double &PosStdDev,
//	double &AngAvg, double &AngStdDev)
//{
//	double sumSqrMatchDist = 0.0;
//	double sumMatchDist = 0.0;
//	double sqrMatchDist;
//	double sumSqrMatchAngle = 0.0;
//	double sumMatchAngle = 0.0;
//	double matchAngle, sqrMatchAngle;
//
//	// return the average match distance of the inliers // AS comment: maybe return the mahalanobis distance instead?
//	for (unsigned int i = 0; i < nSamples; i++)
//	{
//		sqrMatchDist = (matchPts[i] - Freg * samplePts[i]).NormSquare();
//
//		sumSqrMatchDist += sqrMatchDist;
//		sumMatchDist += sqrt(sqrMatchDist);
//
//		matchAngle = acos(matchNorms[i] * (Freg.Rotation() * sampleNorms[i]));
//
//		sumMatchAngle += matchAngle;
//		sumSqrMatchAngle += matchAngle * matchAngle;
//	}
//
//	PosAvg = sumMatchDist / nSamples;
//	PosStdDev = (sumSqrMatchDist / nSamples) + PosAvg*PosAvg;
//	AngAvg = sumMatchAngle / nSamples;
//	AngStdDev = (sumSqrMatchAngle / nSamples) + AngAvg*AngAvg;
//}

void algDirICP_DIMLOP::SetSamples(
  const vctDynamicVector<vct3> &argSamplePts,
  const vctDynamicVector<vct3> &argSampleNorms,
  vctDynamicVector<vct3x3> &argMxi,
  vctDynamicVector<vct3x3> &argMsmtMxi,
  vctDynamicVector<vct3> &argMeanShape,
  double argScale, bool argbScale)
{
  // base class
  //algDirICP_IMLOP::SetSamples(argSamplePts, argSampleNorms);

  //std::cout << "algDirICP_IMLOP::SetSamples()" << std::endl;

  //unsigned int nSamples = samplePts.size();

  if (argMxi.size() != nSamples || argMsmtMxi.size() != nSamples)
  {
	  std::cout << "ERROR: number of covariances matrices does not match number of samples" << std::endl;
  }

  Mxi = argMxi;
  MsmtMxi = argMsmtMxi;
  meanShape = argMeanShape;

  //eigMxi.SetSize(nSamples);

  eigMxi.SetSize(nSamples);
  for (unsigned int i = 0; i < nSamples; i++)
  {
	  ComputeCovEigenValues_SVD(argMxi[i], eigMxi[i]);
  }

  R_Mxi_Rt.SetSize(nSamples);
  R_MsmtMxi_Rt.SetSize(nSamples);
  Myi_sigma2.SetSize(nSamples);
  Myi.SetSize(nSamples);

  outlierFlags.SetSize(nSamples);

  residuals_PostMatch.SetSize(nSamples);
  sqrDist_PostMatch.SetSize(nSamples);

  // scale sample points (remove this from here when you move it to IMLOP)
  sc = argScale;
  bScale = argbScale;
  if (bScale) {
	  for (unsigned int i = 0; i < nSamples; i++)
		  samplePts[i] = samplePts[i] * sc;
	  nTrans = 7; // 3 for rotation, 3 for translation, 1 for scale
  }
  else
	  nTrans = 6; // 3 for rotation, 3 for translation

  nModes = (unsigned int)pDirTree->mesh.modeWeight.size();

  Si = pDirTree->mesh.Si; 
  wi = pDirTree->mesh.wi;

  Tssm_wi.resize(nModes);
  for (int i = 0; i < Tssm_wi.size(); i++)
	  Tssm_wi[i].resize(nSamples);
  Tssm_matchPts.resize(nSamples);

  Tssm_Y.resize(nSamples);
  Tssm_Y_t.resize(nSamples);
  Rat_Tssm_Y_t_x.resize(nSamples);
  Rat_Tssm_Y_t_x_invMx.resize(nSamples);
  Yn_Rat_Xn.resize(nSamples); 

  x_prev.SetSize(nTrans + nModes); // transformation parameters, and n modes
  mu.SetSize(nSamples);
  f.SetSize(nSamples);
  s.SetSize(nModes);
}

void algDirICP_DIMLOP::UpdateShape(vctDynamicVector<double>	&S)
{
	// deformably transform each mesh vertex
	pDirTree->mesh.vertices = meanShape;
	int s;
#ifdef ENABLE_PARALLELIZATION
#pragma omp parallel for
#endif
	for (s = 0; s < pMesh->NumVertices(); s++)
		for (unsigned int i = 0; i < nModes; i++)
			pDirTree->mesh.vertices(s) += (S[i] * wi[i].Element(s));

}

void algDirICP_DIMLOP::UpdateTree()
{
	vctFrm3 FId;
	FId.Assign(vctFrm3::Identity());

	pDirTree->EnlargeBounds(FId);
}

void algDirICP_DIMLOP::ICP_ComputeMatches()
{
	//
	// At the beginning of each correspondence phase, the positions 
	// of the representing the model shape must be recomputed based 
	// on the current values of the model-shape parameters, s 
	//

	// base class
	algDirICP::ICP_ComputeMatches();

	Tssm_Y = matchPts;
}

double algDirICP_DIMLOP::ICP_EvaluateErrorFunction()
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

  vctDynamicVector<vct3x3>  Mi(nSamples);           // noise covariances of match (R*Mxi*Rt + Myi)
  vctDynamicVector<vct3x3>  inv_Mi(nSamples);       // inverse noise covariances of match (R*Mxi*Rt + Myi)^-1
  vctDynamicVector<double>  det_Mi(nSamples);       // determinant of noise covariances of match |R*Mxi*Rt + Myi|
  vctDynamicVector<double>  SqrMahalDist(nSamples); // square Mahalanobis distance of matches = (yi-Rxi-t)'*inv(Mi)*(yi-Rxi-t)

  vct3 residual;
  double sumSqrDist = 0.0;
  double nlogkappa	= 0.0;
  double nklog2PI	= 0.0;
  double ssmCost	= 0.0;
  double logCost	= 0.0;
  double expCost	= 0.0;
  double logExp		= 0.0;
  double sumNormProducts = 0.0;

  for (unsigned int s = 0; s < nSamples; s++)
  {
	//if (outlierFlags[s])	continue;	// skip outliers

	// TODO: Compute errors using the current match point on the 
	// deformed shape after computing updating Si
    residual = samplePtsXfmd.Element(s) - Tssm_Y.Element(s);	// matchPts.Element(s);
	// match covariance
	Mi[s] = R_Mxi_Rt[s] + Myi_sigma2[s];
	// match covariance decomposition
	ComputeCovDecomposition_NonIter(Mi[s], inv_Mi[s], det_Mi[s]);
	// match square distance
	sumSqrDist += residual.NormSquare();
	// match square Mahalanobis distance
	SqrMahalDist[s] = residual * inv_Mi[s] * residual;

	// -- Here We Compute the Full Negative Log-Likelihood -- //
	// Compute error contribution for this sample
	nlogkappa += log(k);
	nklog2PI += 5.0*log(2.0*cmnPI); // 1/2
	expCost += SqrMahalDist[s];		// 1/2
	logCost += log(det_Mi[s]);		// 1/2
	logExp += log(exp(k) - exp(-k));
    sumNormProducts += vctDotProduct(sampleNormsXfmd.Element(s), matchNorms.Element(s));
  }
  ssmCost += Si.NormSquare();		// 1/2

  // include match errors of good samples and of thresholded outliers
  //  to improve smoothness of cost function
  // NOTE: adding an extra k*nSamples to the cost produces a cost function >= 0

  // return B*(sumSqrDist)+k*(nSamples - sumNormProducts); // -logC*nSamples;	// IMLOP cost
  prevCostFuncValue = costFuncValue;
  costFuncValue = /*logExp +*/ k*(nSamples - sumNormProducts) /*- nlogkappa*/ + (nklog2PI + expCost + logCost + ssmCost) / 2.0;
  // ^ add logExp back if you can make sure it is stable - can be unstable if k is big => exp(k) blows up

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
}

// Compute Mu for deformed shape after each optimization
void algDirICP_DIMLOP::ComputeMu()
{
	vct3	v0, v1, v2;
	vct3	v0_Tssm_Y, 
			v1_Tssm_Y, 
			v2_Tssm_Y;

	for (unsigned int s = 0; s < nSamples; s++)
	{
		// Find the 3 vertices of the triangle that the matchPoint lies on,
		// i.e., of the matchDatum on estimated shape 
		f[s] = pDirTree->mesh.faces[matchDatums.Element(s)];
		v0 = pDirTree->mesh.vertices.Element(f[s][0]);
		v1 = pDirTree->mesh.vertices.Element(f[s][1]);
		v2 = pDirTree->mesh.vertices.Element(f[s][2]);

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

//bool algDirICP_DIMLOP::ICP_Terminate(vctFrm3 &F)
//{
//	if (bTerminateAlgorithm)
//	{
//		// return registration for last iteration of non-increasing cost
//		F = Fdec;
//		return true;
//	}
//	else
//	{
//		return false;
//	}
//}

void algDirICP_DIMLOP::ReturnScale(double &scale)
{
	scale = sc;
}

void algDirICP_DIMLOP::ReturnShapeParam(vctDynamicVector<double> &shapeParam)
{
	shapeParam = Si;
}

vctFrm3 algDirICP_DIMLOP::ICP_RegisterMatches()
{
	vctFrm3 F;
	ComputeMu();
	vctDynamicVector<double> x0;
	vctDynamicVector<double> x;

	x0.SetSize(nTrans + nModes);
	x.SetSize(nTrans + nModes);

	// initialize x_prev to FGuess where you initialize Freg
	x0 = x_prev;
	//for (int i = 6; i < x0.size(); i++)
	//{
	//	x0[i] = std::min(x0[i], 3.0);
	//	x0[i] = std::max(x0[i], -3.0);
	//}

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

void algDirICP_DIMLOP::UpdateOptimizerCalculations(const vctDynamicVector<double> &x)
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
	vctDynamicVectorRef<vct3>   X(samplePts);
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
		Tssm_Y.Element(j) = Mu[j][0] * pMesh->vertices[f[j][0]]
			+ Mu[j][1] * pMesh->vertices[f[j][1]]
			+ Mu[j][2] * pMesh->vertices[f[j][2]];

		Tssm_Y_t.Element(j) = Tssm_Y.Element(j) - t;
		Rat_Tssm_Y_t_x.Element(j) = Ra.Transpose() * Tssm_Y_t.Element(j) - sc * X.Element(j);
		ComputeCovDecomposition_NonIter(Mxi.Element(j), inv_Mxi.Element(j), det_Mxi.Element(j));
		Rat_Tssm_Y_t_x_invMx.Element(j) = Rat_Tssm_Y_t_x.Element(j) * inv_Mxi.Element(j); 
		Yn_Rat_Xn.Element(j) = vctDotProduct(Ra * sampleNorms/*Xfmd*/.Element(j), matchNorms.Element(j));
	}
	x_prev = x;
}

double algDirICP_DIMLOP::CostFunctionValue(const vctDynamicVector<double> &x)
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

		f += Rat_Tssm_Y_t_x_invMx.Element(i) * Rat_Tssm_Y_t_x.Element(i) + k*(1-Yn_Rat_Xn.Element(i));
	}

	f += s.DotProduct(s); // Comment out to test

	return f;
}

void algDirICP_DIMLOP::CostFunctionGradient(const vctDynamicVector<double> &x, vctDynamicVector<double> &g)
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
	//vctFixedSizeVectorRef<double, 1, 1> gs(g, 6);
	vctDynamicVectorRef<double> gs(g, nTrans, nModes);

	vct3x3 Jz_a;
	vct3 k_Yn_dRa_Xn;

	unsigned int j;
#ifdef ENABLE_PARALLELIZATION
#pragma omp parallel for
#endif
	for (j = 0; j < nSamples; j++)
	{
		if (outlierFlags[j])	continue;

		for (unsigned int c = 0; c < 3; c++)
		{
			Jz_a.Column(c) = dRa[c].TransposeRef() * Tssm_Y_t[j]; 
			k_Yn_dRa_Xn.Element(c) = -k*sampleNorms[j] * dRa[c] * matchNorms[j];
		}

		ga += Rat_Tssm_Y_t_x_invMx[j] * 2.0 * Jz_a + k_Yn_dRa_Xn;
		gt += Rat_Tssm_Y_t_x_invMx[j] * 2.0 * (-Ra.Transpose());
		if (bScale)
			gsc += Rat_Tssm_Y_t_x_invMx[j] * 2.0 * (-X.Element(j));

		for (unsigned int i = 0; i < nModes; i++)
			gs[i] += Rat_Tssm_Y_t_x_invMx[j] * 2.0 * (Ra.Transpose() * Tssm_wi[i][j]);	// Cmatch component	
	}

	gs += 2.0 * s;	// Cshape component
}

void algDirICP_DIMLOP::ICP_InitializeParameters(vctFrm3 &FGuess)
{
  // initialize base class
  algDirICP_IMLOP::ICP_InitializeParameters(FGuess);
  this->FGuess = FGuess;

  bFirstIter_Matches = true;
  nOutliers = 0;

  bTerminateAlgorithm = false;
  costFuncIncBits = 0;
  costFuncValue = std::numeric_limits<double>::max();

  sigma2 = 0.0;

  // begin with isotropic noise model for first match, 
  // since we don't yet have an approximation for sigma2

  R_Mxi_Rt.SetAll(vct3x3::Eye());
  R_MsmtMxi_Rt.SetAll(vct3x3(0.0));

  Myi_sigma2.SetAll(vct3x3::Eye());
  Myi.SetAll(NULL);

  outlierFlags.SetAll(0);

  if (nSamples != Mxi.size() || nSamples != eigMxi.size())
  {
	  std::cout << " ======> ERROR: noise model array for sample points does not match number of samples" << std::endl
		  << "  Did you forget to call algICP_IMLP::SetSampleCovariances() before starting ICP?" << std::endl;
	  assert(0);
  }

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

void algDirICP_DIMLOP::ICP_UpdateParameters_PostMatch() // CHECK IF YOU NEED NORMAL INFO HERE - should not for sigma because that's from positional data
{
  // base class
  algDirICP::ICP_UpdateParameters_PostMatch();

  // compute sum of square distance of inliers
  sumSqrDist_Inliers = 0.0;
  double sumNormProducts_Inliers = 0.0;
  for (unsigned int s = 0; s < nSamples; s++)
  {
	  residuals_PostMatch[s] = samplePtsXfmd[s] - Tssm_Y[s]; //matchPts[s];
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

  // update noise models of hte matches
  for (unsigned int s = 0; s < nSamples; s++)
  {
	  // update target covariances
	  Myi[s] = pDirTree->DatumCovPtr(matchDatums[s]);	// use pointer here for efficiency
	  //std::cout << matchDatums[s] << " = " << Myi[s] << std::endl;

	  // target covariance with match uncertainty
	  //Myi_sigma2[s] = *Myi[s];						// FIX THIS!!!
	  //Myi_sigma2[s].Element(0, 0) += sigma2;
	  //Myi_sigma2[s].Element(1, 1) += sigma2;
	  //Myi_sigma2[s].Element(2, 2) += sigma2;
  }

  if (bFirstIter_Matches)
  {
	// update noise model
	UpdateNoiseModel(sumSqrDist_Inliers, sumNormProducts_Inliers);

	UpdateNoiseModel_SamplesXfmd(FGuess);
  }

  vctRot3 R(FGuess.Rotation());
  for (unsigned int s = 0; s < nSamples; s++)
  {
	  R_MsmtMxi_Rt[s] = R*MsmtMxi[s] * R.Transpose();
  }

  bFirstIter_Matches = false;
}

void algDirICP_DIMLOP::ICP_UpdateParameters_PostRegister(vctFrm3 &Freg)
{
  // base class
	algDirICP_IMLOP::ICP_UpdateParameters_PostRegister(Freg);

	if (bScale)
		for (unsigned int s = 0; s < nSamples; s++)
			samplePtsXfmd.Element(s) = sc * samplePtsXfmd.Element(s); // move this to IMLOP also later

  UpdateNoiseModel_SamplesXfmd(Freg);

  UpdateShape(Si);
  UpdateTree();

  // Re-initialize to compute matches on updated mesh
  TCPS.init(pDirTree->mesh.vertices, pDirTree->mesh.faces);

#if 0
  static int count = 1; 
  cisstMesh currentSamples;
  currentSamples.vertices.SetSize(nSamples); 
  currentSamples.vertexNormals.SetSize(nSamples);
  currentSamples.vertices = samplePtsXfmd;
  currentSamples.vertexNormals = sampleNormsXfmd;
  currentSamples.SavePLY("currentSamples" + std::to_string(count) + ".ply");

  cisstMesh currentMesh;
  currentMesh.vertices = pMesh->vertices;
  currentMesh.faces = pMesh->faces;
  currentMesh.SavePLY("currentMesh" + std::to_string(count) + ".ply");

  count++;
#endif
}

void algDirICP_DIMLOP::UpdateNoiseModel_SamplesXfmd(vctFrm3 &Freg)
{
	// update noise models of the transformed sample points
	static vctRot3 R;
	R = Freg.Rotation();
	for (unsigned int s = 0; s < nSamples; s++)
	{
		R_Mxi_Rt[s] = R * Mxi[s] * R.Transpose();
	}
}

unsigned int algDirICP_DIMLOP::ICP_FilterMatches()
{
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
	ThetaThresh = ThetaThresh > cmnPI ? cmnPI : ThetaThresh;
	double NormProductThresh = cos(ThetaThresh);
	//std::cout << "(" << circSD << ", " << ThetaThresh << ", " << NormProductThresh << ") ";

	nOutliers = 0;
	//nPosOutliers = 0;
	//nNormOutliers = 0;
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
		Mo = R_Mxi_Rt.Element(s);
		Mo.Element(0, 0) += sigma2;
		Mo.Element(1, 1) += sigma2;
		Mo.Element(2, 2) += sigma2;
		//Mo = R_Mxi_Rt.Element(s) + Myi_sigma2.Element(s);
		//Mo.Element(0, 0) += outlier_alpha;
		//Mo.Element(1, 1) += outlier_alpha;
		//Mo.Element(2, 2) += outlier_alpha;

		// compute Mahalanobis distance
		ComputeCovInverse_NonIter(Mo, inv_Mo);
		sqrMahalDist = residuals_PostMatch.Element(s)*inv_Mo*residuals_PostMatch.Element(s); 
		normProduct = vctDotProduct(sampleNormsXfmd.Element(s), matchNorms.Element(s));

		// check if outlier
		if (sqrMahalDist > ChiSquareThresh)
		{ // an outlier
			nOutliers++;
			//nPosOutliers++;
			outlierFlags[s] = 1;
			// add isotropic outlier term to noise model for this match
			//  with magnitude of half the square match distance
			double outlierScale = 0.5 * sqrDist_PostMatch.Element(s) * varExpansionFactor;
			Myi_sigma2[s].Element(0, 0) += outlierScale;
			Myi_sigma2[s].Element(1, 1) += outlierScale;
			Myi_sigma2[s].Element(2, 2) += outlierScale;
			R_Mxi_Rt[s].Element(0, 0) += outlierScale;
			R_Mxi_Rt[s].Element(1, 1) += outlierScale;
			R_Mxi_Rt[s].Element(2, 2) += outlierScale;

			// This shouldn't be done here, because alpha term is not used
			//  in error function
			//// For Error Function Evaluation:
			//// update match covariance
			//Mi.Element(s) = R_Mxi_Rt.Element(s) + Myi.Element(s);
			//// match covariance decomposition
			//ComputeCovDecomposition(Mi.Element(s), inv_Mi.Element(s), det_Mi.Element(s));
			//// match Mahalanobis distance
			//SqrMahalDist.Element(s) = Residuals.Element(s)*inv_Mi.Element(s)*Residuals.Element(s);
		}
		else if (normProduct < NormProductThresh)
		{
			nOutliers++;
			//nNormOutliers++;
			outlierFlags[s] = 1;
			//std::cout << "\n(" << normProduct << " ? " << NormProductThresh << ") ";
		}
		else
		{
			outlierFlags[s] = 0;
		}
	}

	return nOutliers;
}

//unsigned int algDirICP_DIMLOP::ICP_FilterMatches()
//{
//
//  // Method 1: No Outlier Detection
//#if 1
//  nGoodSamples = 0;
//  nOutliers = 0;
//  nPosOutliers = 0;
//  nNormOutliers = 0;
//
//  // use all samples
//  for (unsigned int s = 0; s < nSamples; s++)
//  { // copy all samples to buffer
//    goodSamplePtsBuf.Element(s) = samplePts.Element(s);
//    goodMatchPtsBuf.Element(s) = matchPts.Element(s);
//
//    goodSampleNormsBuf.Element(s) = sampleNorms.Element(s);
//    goodMatchNormsBuf.Element(s) = matchNorms.Element(s);
//  }
//  nGoodSamples = nSamples;
//
//  // Non-destructively resize good sample reference vectors
//  goodSamplePts.SetRef(goodSamplePtsBuf, 0, nGoodSamples);
//  goodMatchPts.SetRef(goodMatchPtsBuf, 0, nGoodSamples);
//
//  goodSampleNorms.SetRef(goodSampleNormsBuf, 0, nGoodSamples);
//  goodMatchNorms.SetRef(goodMatchNormsBuf, 0, nGoodSamples);
//
//  //gSumSqrDist_PostMatch = sumSqrDist_PostMatch;
//  //gSumNormProducts_PostMatch = sumNormProducts_PostMatch;
//  //oSumSqrDist_PostMatch = 0.0;
//  //oSumNormProducts_PostMatch = 0.0;
//
//  return nOutliers;
//#endif 
//
//  // Method 2: Chi-Square Outlier Test
//#if 0
//  // Filer Matches for Outliers
//  //  
//  // Positions:
//  //  The Square Mahalanobis Distance of the matches follow a chi-square distribution
//  //  with 3 degrees of freedom (1 DOF for each spatial dimension) assuming that the
//  //  residuals of position are normally distributed about the mean. Since the
//  //  positional noise model is isotropic, the square Mahalanobis distance reduces
//  //  to the square Euclidean distance divided by the variance along one spatial dimension.
//  //
//  //  Detect position outliers as: ||Py - (R*Px+t)||^2 > ChiSquare(c) * Var(Pos)
//  //
//  //  Note:  ChiSquare(0.95) = 7.81     (1.96 Std Dev)
//  //         ChiSquare(0.975) = 9.35    (2.24 Std Dev)
//  //         ChiSquare(0.99) = 11.34    (2.56 Std Dev)
//  //         ChiSquare(0.9973) = 14.16  (3.0 Std Dev)     MATLAB: chi2inv(0.9973,3)
//  //
//  //  Note:  Var(Pos) = (1/3*N)*Sum_i(||Pyi - (R*Pxi+t)||^2)
//  //          - here we assume the residual error to have zero mean
//  //          - we use 3*N as the denominator rather than N because
//  //            the square distance is a sum of three independent Guassian
//  //            RVs (one for each axis) and it is the variance along a single
//  //            axis that appears in the Mahalanobis distance and Gaussian
//  //            probability equations.
//  //
//  // Orientations:
//  // Consider a point as outlier if orientation angle (acos(Na'Nb)) > 3 circular
//  //  standard deviations
//  //  
//  //  Detect Orientation outlier as:  acos(dot(Nx,Ny)) > 3*circStdDev
//  //
//  //  Copies only non outlier points to "GoodsamplePts" and 
//  //   nondestructively resizes "Good Sample" reference vectors 
//  //   to the number of points transferred.
//  //
//
//  //double ChiSquareThresh = 11.34;
//  double ChiSquareThresh = 14.16;
//  double S2 = pICP->SumSqrDist_PostMatch / (3*nSamples); // variance estimate along a single axis
//  double SquareDistThresh = S2 * ChiSquareThresh;
//
//  double ThetaThresh = 3.0*pICP->circSD_PostMatch;
//  ThetaThresh = ThetaThresh > cmnPI ? cmnPI : ThetaThresh;
//  double NormProductThresh = cos(ThetaThresh);
//
//  nGoodSamples = 0;
//  nOutliers = 0;
//  nPosOutliers = 0;
//  nNormOutliers = 0;
//
//  gSumSqrDist_PostMatch = 0.0;
//  gSumNormProducts_PostMatch = 0.0;
//  oSumSqrDist_PostMatch = 0.0;
//  oSumNormProducts_PostMatch = 0.0;
//  //gSumSqrNormDist_PostMatch = 0.0;
//  //oSumSqrNormDist_PostMatch = 0.0;
//
//  // for all model/sample sets
//  for (unsigned int set = 0; set < pICP->nSets; set++)
//  {
//    unsigned int nSamps_Set = pICP->nsamplePtsSets[set];
//    unsigned int nGood_Set = 0;
//
//    vctDynamicVectorRef<vct3>   samplePts( pICP->SampleSets[set] );
//    vctDynamicVectorRef<vct3>   sampleNorms( pICP->SampleNormSets[set] );
//    vctDynamicVectorRef<vct3>   matchPts( pICP->ClosestPointSets[set] );
//    vctDynamicVectorRef<vct3>   matchNorms( pICP->ClosestPointNormSets[set] );
//    vctDynamicVectorRef<double> SquareDistances( pICP->SquareDistanceSets_PostMatch[set] );
//    vctDynamicVectorRef<double> NormProducts( pICP->NormProductSets_PostMatch[set] );
//    //vctDynamicVectorRef<double> Distances( pICP->DistanceSets_PostMatch[set] );
//    //vctDynamicVectorRef<vct3>   TransformedsamplePts( pICP->TransformedSampleSets[set] );    
//    //vctDynamicVectorRef<vct3>   TransformedsampleNorms( pICP->TransformedSampleNormSets[set] );
//    //vctDynamicVectorRef<double> SquareDistances( pICP->SquareDistanceSets_PostMatch[set] );
//    //vctDynamicVectorRef<vct3>   Residuals( pICP->ResidualSets_PostMatch[set] );
//    //vctDynamicVectorRef<double> dThetas( pICP->dThetaSets_PostMatch[set] );
//
//    // Filter Outliers
//    for (unsigned int s = 0; s < nSamps_Set; s++)
//    {
//      // Positional Outlier Test
//      if ( SquareDistances.Element(s) > SquareDistThresh )
//      {	// a positional outlier
//        nPosOutliers += 1;
//
//        // Theshold the outlier error on position and (if applicable) orientation 
//        //  so that the cost function remains smooth, but the outlier does not
//        //  continue to increase the error.
//        //
//        // Dual error thresholds (one for position, one for orientation) are
//        //  not ideal, as the error of an outlier may then
//        //  change during its time as an outlier due to thresholding on
//        //  both positional and orientation error. The alternative is a single
//        //  threshold error value, but this would again create jumps in the error
//        //  when adding/removing outliers.
//        //
//        // Either way the optimization runs correctly, what this effects is the
//        //  error function value, which may be being monitored for the termination 
//        //  condition. We have two choices:
//        //   a) some change in outlier error affects the error function
//        //   b) jumpy cost function in between optimizations, when adding/removing
//        //      outliers, but cost of an outlier is always the same
//        // Since a) is smooth and b) may be jumpy, we choose situation a).
//        //
//        double normProduct;
//        //double sqrNormDist;
//        // check if orientation is outlier as well
//        //  if so, threshold orientation component of error as well
//        // Note: smaller norm products mean larger error
//        if ( NormProducts.Element(s) < NormProductThresh )
//        { // orientation is outlier as well => threshold its error
//          normProduct = NormProductThresh;
//          //sqrNormDist = sqrNormDistThresh;
//        }
//        else
//        {
//          normProduct = NormProducts.Element(s);
//          //sqrNormDist = SquareNormDistances.Element(s);
//        }
//        // apply match error to outlier error
//        oSumSqrDist_PostMatch += SquareDistThresh;
//        oSumNormProducts_PostMatch += normProduct;
//        //oSumSqrNormDist_PostMatch += sqrNormDist; 
//      }
//      // Orientation Outlier Test
//      else if ( NormProducts.Element(s) < NormProductThresh )
//      { // an orientation outlier
//        nNormOutliers += 1;
//
//        // Threshold the outlier error on orientation so that the cost function is
//        //  smooth, but the outlier does not continue to increase the error.
//        //  (Don't have to threshold position error, because the position was
//        //   already checked and cannot be an outlying value here)
//        // apply match error to outlier error
//        oSumSqrDist_PostMatch += SquareDistances.Element(s);
//        oSumNormProducts_PostMatch += NormProductThresh;
//        //oSumSqrNormDist_PostMatch += sqrNormDistThresh;
//      }
//      else
//      { // inlier
//        // copy this sample to the good sample buffers
//        goodSamplePtsBuf.Element(nGoodSamples) = samplePts.Element(s);        
//        goodSampleNormsBuf.Element(nGoodSamples) = sampleNorms.Element(s);
//        goodMatchPtsBuf.Element(nGoodSamples) = matchPts.Element(s);
//        goodMatchNormsBuf.Element(nGoodSamples) = matchNorms.Element(s);
//        //AllGoodTransformedsamplePtsBuf.Element(pICP->nGoodSamples) = TransformedsamplePts.Element(s);
//        //AllGoodTransformedsampleNormsBuf.Element(pICP->nGoodSamples) = TransformedsampleNorms.Element(s);
//        //GoodSampleResiduals_PostMatch.Row(pICP->nGoodSamples).Assign(Residuals.Element(s));
//
//        // apply match error to good sample error
//        gSumSqrDist_PostMatch += SquareDistances.Element(s);
//        gSumNormProducts_PostMatch += NormProducts.Element(s);
//        //gSumSqrNormDist_PostMatch += SquareNormDistances.Element(s);
//
//        nGoodSamples++;
//        nGood_Set++;
//      }
//    }
//
//    //// renormalize weights for good samples
//    //for (unsigned int k = nGoodSamples-nGood_Set; k < nGoodSamples; k++)
//    //{ // set sample weights
//    //  // must know # of good samples in the set before setting the sample weights
//    //  //  due to normalization factor
//    //  AllGoodSampleWeightsBuf.at(k) = pICP->Weights[set]/nGood_Set;
//    //}
//    //nOutliersSets[set] = nSamps_Set - nGood_Set;
//    nOutliers += nSamps_Set - nGood_Set;
//  }
//
//  goodSamplePts.SetRef( goodSamplePtsBuf,0,nGoodSamples );  
//  goodSampleNorms.SetRef( goodSampleNormsBuf,0,nGoodSamples );  
//  goodMatchPts.SetRef( goodMatchPtsBuf,0,nGoodSamples );
//  goodMatchNorms.SetRef( goodMatchNormsBuf,0,nGoodSamples );
//  //AllGoodSampleWeights.SetRef( AllGoodSampleWeightsBuf,0,nGoodSamples );
//  //AllGoodTransformedsamplePts.SetRef( AllGoodTransformedsamplePtsBuf,0,nGoodSamples );
//  //AllGoodTransformedsampleNorms.SetRef( AllGoodTransformedsampleNormsBuf,0,nGoodSamples );
//
//  pICP->iterData->nOutliersPos = nPosOutliers;
//  pICP->iterData->nOutliersNorm = nNormOutliers;
//
//  return nOutliers;
//
//  //if ( outlierDistFactor <= EPS && outlierNormFactor <= EPS )
//  //{	// outlier filtering is disabled => use all samples
//  //  for (unsigned int s = 0; s < nSamps_Set; s++)
//  //  { // copy all samples to buffer
//  //    goodSamplePtsBuf.Element(nGoodSamples) = samplePts.Element(s);
//  //    goodSampleNormsBuf.Element(nGoodSamples) = sampleNorms.Element(s);
//  //    goodMatchPtsBuf.Element(nGoodSamples) = matchPts.Element(s);        
//  //    goodMatchNormsBuf.Element(nGoodSamples) = matchNorms.Element(s);
//  //    AllGoodSampleWeightsBuf.Element(nGoodSamples) = pICP->Weights[set] / nSamps_Set;
//  //    nGoodSamples++;
//  //    //AllGoodTransformedsamplePtsBuf.Element(nGoodSamples) = TransformedsamplePts.Element(s);
//  //    //AllGoodTransformedsampleNormsBuf.Element(nGoodSamples) = TransformedsampleNorms.Element(s);
//  //  }
//  //  nGood_Set = nSamps_Set;
//  //  // apply all errors to good sample errors
//  //  gSumSqrDist_PostMatch += pICP->SumSquareDistanceSets_PostMatch[set];
//  //  gSumNormProducts_PostMatch += pICP->SumNormProductSets_PostMatch[set];
//  //  //gSumSqrNormDist_PostMatch += SumSquareNormDistanceSets_PostMatch[set];      
//  //}
//
//  //// compute outlier thresholds for position and orientation
//  //double sqrDistThresh = outlierDist2Factor * posVar;
//  //double thetaThresh = outlierNormFactor * pICP->circSD_PostMatch;    
//  //thetaThresh = thetaThresh > cmnPI ? cmnPI : thetaThresh;
//  //double normProductThresh = cos(thetaThresh);
//  ////// Use theta threshold to compute a threshold on the square geometric 
//  //////  distance between normal vectors; to compute this, compare a unit
//  //////  vector to itself rotated by theta threshold degrees
//  //////  (used for thresholding outlier square norm distance error, not
//  //////   used for detecting outlier)
//  ////vct3 t1(1.0,0.0,0.0);
//  ////vctAxAnRot3 AxAn(vct3(0.0,0.0,1.0),thetaThresh);
//  ////vct3 t2 = vctRot3(AxAn)*t1;   // *** why must convert to vctRot3 type?
//  ////double sqrNormDistThresh = (t1-t2).NormSquare();
//#endif
//}


// PD Tree Methods

int algDirICP_DIMLOP::NodeMightBeCloser(
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

// PD Tree Methods

double algDirICP_DIMLOP::FindClosestPointOnDatum(
	const vct3 &v, const vct3 &n,
	vct3 &closest, vct3 &closestNorm,
	int datum)
{
	// set closest point
	TCPS.FindClosestPointOnTriangle(v, datum, closest);

	// norm has same value everywhere on this datum
	closestNorm = pDirTree->mesh.faceNormals[datum];

	// return modified vMFG match error such that the match error
	//  is always >= 0
	return k*(1 - n.DotProduct(closestNorm)) + B*(v - closest).NormSquare();
}


int algDirICP_DIMLOP::DatumMightBeCloser(
	const vct3 &v, const vct3 &n,
	int datum,
	double ErrorBound)
{
	BoundingBox BB;
	vct3 v1, v2, v3;
	pDirTree->mesh.FaceCoords(datum, v1, v2, v3);
	BB.Include(v1);
	BB.Include(v2);
	BB.Include(v3);

	// We want to know if this point can produce a cost less than the error bound.
	//  Error bound is the best cost we have so far, and we know N for a triangle
	//  datum is the same everywhere => Nclosest is known for this datum.
	//  Some of the match error for this datum comes from Nclosest. Subtract this error
	//  from the error bound to get the remaining max error attributable to distance. Then
	//  check if the sample point lies w/in this distance of a bounding box
	//  around this datum.
	//
	//  error = k*(1-N'*Nclosest) + B*||v - closest||^2 
	//                                   (dist^2)
	//
	//    =>  maxDist = sqrt([ErrorBound - k*(1-N'*Nc)]/B)
	//

	double searchDist2 = (ErrorBound - k*(1 - n.DotProduct(pDirTree->mesh.faceNormals[datum]))) / B;
	// don't take square root of negative number
	if (searchDist2 > 0)
	{
		return BB.Includes(v, sqrt(searchDist2));
	}
	else
	{ // the difference in normal orientation alone creates an error
		//  greater than the error bound => this datum cannot be closer
		return 0;
	}
}


// Helper Methods:

void algDirICP_DIMLOP::UpdateNoiseModel(double sumSqrDist, double sumNormProducts)
{
  //// Position Parameters
  //// compute Gaussian variables
  ////  B = 1/(2*sigma2)
  //// divide sum of square distances by number of samples times 3 
  ////  to get the estimate of variance because square distance is 
  ////  a summation of three square Gaussian RVs, one for each coordinate axis
  //double S2 = sumSqrDist / (3.0*nSamples);
  //if (dynamicParamEst)
  //{
  //  B = 1.0 / (2.0*S2);
  //  B = B > threshB ? threshB : B;  // threshold in case of perfect matches
  //  sigma2 = 1.0 / (2.0*B);
  //}
  //
  //// Orientation Parameters
  //// compute Fisher variables
  ////
  ////  MLE estimate for k is an implicit ratio of
  ////   two Bessel functions => no analytical soln.
  ////
  ////  MLE for k:  set Ap(k) = R
  ////              
  ////    Ap(k) = Ip/2(k) / Ip/2-1(k)       <--- Bessel Functions
  ////    p =  spatial dimension (i.e. 3)
  ////    R =  Sum_i(dot(Ny,Nx))/N
  ////
  ////  Closed form approximation for k:  R(p-R^2)/(1-R^2)
  ////
  ////   NOTE: circular standard deviation of theta
  ////         circSD = sqrt(-2*ln(R))
  ////            where theta is angle between matched vectors
  ////            i.e. theta = acos(dot(Ny,Nx))
  ////

  // angular error of normal orientations
  ComputeCircErrorStatistics(sumNormProducts, Rnorm, circSD);

  //// angular error of positions (wrt ctr of rotation)
  //Rpos = ComputeRpos();
  //
  //// effective angular error
  //// only include positional error if it reduces the value of K
  //if (Rpos < Rnorm)
  //{
  //  R = wRpos*Rpos + (1.0 - wRpos)*Rnorm;
  //}
  //else
  //{
  //  R = Rnorm;
  //}
  //
  //double R2 = R*R;
  //if (dynamicParamEst)
  //{
  //  // protect from division by zero
  //  if (R2 >= 1.0)
  //  {
  //    k = threshK;  // set k to its max value
  //  }
  //  else
  //  {
  //    k = R*(3.0 - R2) / (1.0 - R2);  // approx for k
  //    k = k > threshK ? threshK : k;  // perfect match threshold
  //  }
  //
  //  // reduce k by a factor
  //  k = k * k_factor;
  //}

#ifdef TEST_STD_ICP
  k = 0.0;
#endif
}


//double algDirICP_DIMLOP::ComputeRpos()
//{
//
//#define NMLZ_METHOD 1   // Normalization method
//
//  // NOTE: could improve efficiency by computing this value as an 
//  //       optional output argument in the vMFG P2P Registration
//  //       function
//
//  // Compute angular match error of sample positions relative 
//  //  to the center of rotation; this is to include angular error
//  //  of position matches measured from center of rotation of
//  //  Procrustes problem as added indicator of rotational uncertainty
//
//  vctDynamicVectorRef<vct3> S(samplePtsXfmd);
//  vctDynamicVectorRef<vct3> C(/*matchPts*/Tssm_Y);
//  unsigned int N = S.size();
//  vct3 Smean = vctCentroid(S);
//  vct3 Cmean = vctCentroid(C);
//  vct3 Sp, Cp;
//  double dotProducts = 0.0;
//  double nmlzTerm = 0.0;
//  for (unsigned int i = 0; i < N; i++)
//  {
//    Sp.DifferenceOf(S[i], Smean);
//    Cp.DifferenceOf(C[i], Cmean);
//
//    //  Do not scale down to unit vectors, as points farther
//    //    away should have greater effect on rotation.
//    //  R calculation assumes unit vector normalization => need some
//    //    form of normalization in the end.
//    //  TODO: should normalization be a linear or square term?
//#if (NMLZ_METHOD == 1)
//    // use square normalization
//    //  this is the best solution, as it works with the data directly
//    //  i.e. the orientations do not have to be normalized prior to
//    //       taking their dot product
//    dotProducts += vctDotProduct(Sp, Cp);
//    nmlzTerm += Sp.Norm()*Cp.Norm();
//#elif (NMLZ_METHOD == 2)
//    // use linear normalization
//    double avgLen = (Sp.Norm() + Cp.Norm()) / 2.0;
//    dotProducts += avgLen*vctDotProduct(Sp.Normalized(), Cp.Normalized());
//    nmlzTerm += avgLen;
//#elif (NMLZ_METHOD == 3)
//    // use unit vectors
//    dotProducts += vctDotProduct(Sp.Normalized(), Cp.Normalized());
//#else
//    std::cout << "Error: normalization method unrecognized" << std::endl;
//#endif
//  }
//#if (NMLZ_METHOD == 3)
//  nmlzTerm = N;
//#endif
//
//  double Rpos = dotProducts / nmlzTerm;
//  // 0 <= Rpos <= 1
//  Rpos = Rpos < 0.0 ? 0.0 : Rpos;
//  Rpos = Rpos > 1.0 ? 1.0 : Rpos;
//  //double posCircSD = sqrt(-2*log(Rpos));
//  return Rpos;
//}


//// Below are various methods that were tried in an attempt to balance the k & B
////  weights s.t. k does not blow up too large
//
////switch (opt.errFunc)
////{
////case OptionsNorm::STDICP:
////  {
////    // orientation is not used in standard ICP => k=0
////    double sigma2 = this->S2;
////    k = 0.0;
//
////    // perfect match threshold
////    //  need to threshold both sigma2 & B, since sigma2 used for calculating
////    //  normalization constant.
////    //  Note:  if B > n, then sigma2 < 1/(2*n)
////    //B = B > 1e5 ? 1e5 : B;
////    sigma2 = sigma2 < 1.0/(2.0e5) ? 1.0/(2.0e5) : sigma2;
////    B = 1.0/(2.0*sigma2);
//
////    Compute_vMFGConstant( k,sigma2,logC );
////    err = ComputeErrorFunctionValue( B,k,logC,SumSquareDistances,SumNormProducts );
////    break;    
////  }
////case OptionsNorm::vMFG:
////  {
////    // update error function parameters
////    double sigma2;
////    Compute_vMFGParams( k,sigma2 );
//
////    // perfect match thresholds
////    //  need to threshold both sigma2 & B, since sigma2 used for calculating
////    //  normalization constant.
////    //  Note:  if B > n, then sigma2 < 1/(2*n)
////    //B = B > threshB ? threshB : B;
////    k = k > threshK ? threshK : k;
////    sigma2 = sigma2 < threshS2 ? threshS2 : sigma2;
////    B = 1.0/(2.0*sigma2);
//
////    // compute normalization constant for probability distribution
////    Compute_vMFGConstant( k,sigma2,logC );
//
////    // compute error function
////    err = ComputeErrorFunctionValue( B,k,logC,SumSquareDistances,SumNormProducts );
////    break;
////  }
////case OptionsNorm::kROTREG:
////  {
////    double sigma2 = this->S2;
//
////    // regularize k by dR
////    vctRodRot3 dR;
////    dR.From(iterData->dF.Rotation());   // convert to Rodrigues notation
////    double dAng = dR.Norm();
////    double cSD = this->circSD + dAng;
////    //  SD = sqrt(-2*log(R));
////    //   =>  R = exp(-(SD^2)/2)
////    double Rreg = exp(-(cSD*cSD)/2);
////    double Rreg2 = Rreg*Rreg;
////    k = Rreg*(3-Rreg2)/(1-Rreg2);
//
////    // perfect match thresholds
////    k = k > threshK ? threshK : k;
////    sigma2 = sigma2 < threshS2 ? threshS2 : sigma2;
////    B = 1.0/(2.0*sigma2);
//
////    Compute_vMFGConstant( k,sigma2,logC );
////    err = ComputeErrorFunctionValue( B,k,logC,SumSquareDistances,SumNormProducts );
////    break;
////  }
////case OptionsNorm::kREG:
////  {
////    double prevK = k; 
////    double sigma2;
////    Compute_vMFGParams( k,sigma2 );
//
////    // regularization on k
////    k = prevK + (k-prevK)*opt.kRegConst;
//
////    // perfect match thresholds
////    k = k > threshK ? threshK : k;
////    sigma2 = sigma2 < threshS2 ? threshS2 : sigma2;
////    B = 1.0/(2.0*sigma2);
//
////    Compute_vMFGConstant( k,sigma2,logC );
////    err = ComputeErrorFunctionValue( B,k,logC,SumSquareDistances,SumNormProducts );
////    break;
////  }
////case OptionsNorm::kbREG:
////  {
////    double prevK = k;
////    double prevB = B;
//
////    double sigma2;
////    Compute_vMFGParams( k,sigma2 );
//
////    // perfect match thresholds
////    k = k > threshK ? threshK : k;
////    sigma2 = sigma2 < threshS2 ? threshS2 : sigma2;
////    B = 1.0/(2.0*sigma2);
//
////    Compute_vMFGConstant( k,sigma2,logC );
////    err = ComputeErrorFunctionValue( B,k,logC,SumSquareDistances,SumNormProducts );
////    break;
////  }
////case OptionsNorm::kTHRESH:
////  {
////    double sigma2;
////    Compute_vMFGParams( k,sigma2 );
//
////    // Threshold K
////    //  Note: k = 600 equates to about 3.25 degrees circular standard deviation
////    k = k > opt.kThresh ? opt.kThresh : k;
//
////    // perfect match thresholds
////    sigma2 = sigma2 < threshS2 ? threshS2 : sigma2;
////    B = 1.0/(2.0*sigma2);
//
////    Compute_vMFGConstant( k,sigma2,logC );
////    err = ComputeErrorFunctionValue( B,k,logC,SumSquareDistances,SumNormProducts );
////    break;
////  }
////case OptionsNorm::kSCALE:
////  {
////    double sigma2;
////    Compute_vMFGParams( k,sigma2 );
//
////    // Scale K
////    k /= opt.kScale;
//
////    // perfect match thresholds
////    k = k > threshK ? threshK : k;
////    sigma2 = sigma2 < threshS2 ? threshS2 : sigma2;
////    B = 1.0/(2.0*sigma2);
//
////    Compute_vMFGConstant( k,sigma2,logC );
////    err = ComputeErrorFunctionValue( B,k,logC,SumSquareDistances,SumNormProducts );
////    break;
////  }
////case OptionsNorm::kFIXED:
////  {
////    k = opt.k;
////    double sigma2 = this->S2;
//
////    // perfect match thresholds
////    k = k > threshK ? threshK : k;
////    sigma2 = sigma2 < threshS2 ? threshS2 : sigma2;
////    B = 1.0/(2.0*sigma2);
//
////    Compute_vMFGConstant( k,sigma2,logC );      
////    err = ComputeErrorFunctionValue( B,k,logC,SumSquareDistances,SumNormProducts );
////    break;
////  }
////case OptionsNorm::kbFIXED:
////  {
////    k = opt.k;
////    B = opt.B;
////    double sigma2 = 1.0/(2.0*B);
////    Compute_vMFGConstant( k,sigma2,logC );
////    err = ComputeErrorFunctionValue( B,k,logC,SumSquareDistances,SumNormProducts );
////    break;
////  }
////case OptionsNorm::WRAPNORMAL:
////  {
////    // assign k according to the variance of a
////    //  wrapped normal distribution
////    double sigma2 = this->S2;
////    k = 1.0/(2.0*(this->circSD)*(this->circSD));  // ERROR: should be k = 1/circSD^2
//
////    // perfect match thresholds
////    k = k > threshK ? threshK : k;
////    sigma2 = sigma2 < threshS2 ? threshS2 : sigma2;
////    B = 1.0/(2.0*sigma2);
//
////    Compute_vMFGConstant( k,sigma2,logC );
////    err = ComputeErrorFunctionValue( B,k,logC,SumSquareDistances,SumNormProducts );
////    break;
////  }
////default:
////  {
////    std::cout << "No valid error function chosen" << std::endl;
////    assert(0);
////  }
////}