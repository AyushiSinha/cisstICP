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
#ifndef _algDirICP_DIMLOP_h
#define _algDirICP_DIMLOP_h

#include "algDirICP_IMLOP.h"
#include "DirPDTree_Mesh.h"
#include "TriangleClosestPointSolver.h"
#include "algDirICP_DIMLOP_dlibWrapper.h"


// Debug Modes
//#define TEST_STD_ICP

class algDirICP_DIMLOP : public algDirICP_IMLOP 
{
  //
  // This algorithm implements the von-Mises Fisher + Gaussian negative
  //   log likelihood cost function for position and orientation based
  //   registration in the presence of deformation making use of shape
  //   shape statistics.
  //


  //--- Algorithm Parameters ---//
public:
	algDirICP_DIMLOP_dlibWrapper dlib;

	vctDynamicVector<vct3>		Tssm_matchPts;
	vctDynamicVector<vct3>		mu;
	vctDynamicVector<vctInt3>	f;
	vctDynamicVector<vctDynamicVector<vct3>>	Tssm_wi;

	double rb, tb, sb, spb;		// rotation, translation, scale, and shape parameter bounds
	bool bScale;

	// -- Optimizer calculations common to both cost and gradient function
	vctRot3 Ra;
	vct3 a, t;
	double sc;
	vctDynamicVector<double> s;
	vctDynamicVector<double> x_prev;

	vctDynamicVector<vct3> X;
	vctDynamicVector<vct3> Tssm_Y;
	vctDynamicVector<vct3> Tssm_Y_t;
	vctDynamicVector<vct3> Rat_Tssm_Y_t_x;
	vctDynamicVector<vct3> Rat_Tssm_Y_t_x_invMx;
	vctDynamicVector<double> Yn_Rat_Xn;
	//vctDynamicVector<vct3> k_Yn_dRa_Xn;

protected:
	TriangleClosestPointSolver TCPS;
	cisstMesh *pMesh;
	DirPDTree_Mesh *pDirTree;

	vctFrm3 FGuess; // intiial guess for registration

	bool bFirstIter_Matches;  // flag that the first iteration is being run

	// dynamic noise model
	double sigma2;								// match uncertainty (added to My covariances as sigma2*I)
	double sigma2Max 
		= std::numeric_limits<double>::max();   // max threshold on match uncertainty
	vctDynamicVector<vct3>  residuals_PostMatch;    // Pclosest - Psample
	vctDoubleVec            sqrDist_PostMatch;      // ||Pclosest - Psample||^2

	// noise model
	// Note: Myi values are obtained from the mesh object
	vctDynamicVector<vct3x3>  Mxi;        // noise covariances of sample points
	vctDynamicVector<vct3>    eigMxi;     // eigenvalues of sample covariances
	vctDynamicVector<vct3x3>  R_Mxi_Rt;   // noise covariances of transformed sample points
	vctDynamicVector<vct3x3*> Myi;        // noise covariances of target correspondence points
	vctDynamicVector<vct3x3>  Myi_sigma2; // noise covariances of target correspondence points with match uncertainty added

	// measurement noise component of the sample noise model
	//  Note: if applying measurement noise to the target shape, then
	//        MsmtMyi should be used here as well
	vctDynamicVector<vct3x3> MsmtMxi;
	// covariance model for outlier tests
	//  (does not include planar noise model)
	vctDynamicVector<vct3x3> R_MsmtMxi_Rt;

	// algorithm-specific termination
	bool bTerminateAlgorithm;
	unsigned char costFuncIncBits;
	double costFuncValue, prevCostFuncValue, prevIncCostFuncValue;
	vctFrm3 Fdec;

	// outlier handling
	unsigned int nOutliers;
	double ChiSquareThresh = 7.81;
	double sumSqrDist_Inliers;
	vctDynamicVector<int>   outlierFlags;

	// Deformable Variables
	vctDynamicVector<vct3>	meanShape;
	unsigned int  nTrans; // # transformation parameters
	unsigned int  nModes;
	//vctDynamicVector<vct3>	sampleModes;
	//vctDynamicVector<double> sampleModeWts;	
	vctDynamicVector<vctDynamicVector<vct3>>	wi;
	vctDynamicVector<double>					Si;		// shape parameter


  //--- Algorithm Methods ---//
public:

  // constructor
	algDirICP_DIMLOP(
		//DirPDTreeBase *pDirTree,
		DirPDTree_Mesh *pDirTree,
		vctDynamicVector<vct3> &samplePts,
		vctDynamicVector<vct3> &sampleNorms,
		vctDynamicVector<vct3x3> &sampleCov,
		vctDynamicVector<vct3x3> &sampleMsmtCov,
		vctDynamicVector<vct3> &meanShape,
		double kinit = 0.0, double sigma2init = 1.0, double wRpos = 0.5,
		double kfactor = 1.0, double scale = 1.0,
		bool dynamicallyEstParams = true, 
		bool bScale = false);
 //   : algDirICP_IMLOP(pDirTree, samplePts, sampleNorms),
 //   //algDirPDTree(pDirTree),
	//pDirTree(pDirTree)
 // {
 //   // Ensure SetSamples function of this derived class gets called
 //   SetSamples(samplePts, sampleNorms);
 // }

  // destructor
  virtual ~algDirICP_DIMLOP() {}

  void    UpdateNoiseModel(double sumSqrDist, double sumNormProducts);
  //double  ComputeRpos();

  //void SetNoiseModel(
  //  double initK, double initSigma2, double wRpos, bool dynamicallyEstParams);

  void SetSamples(
    const vctDynamicVector<vct3> &argSamplePts,
	const vctDynamicVector<vct3> &argSampleNorms,
	vctDynamicVector<vct3x3> &sampleCov,
	vctDynamicVector<vct3x3> &sampleMsmtCov,
	vctDynamicVector<vct3> &meanShape,
	double argScale = 1.0,
	bool argbScale = false);

  void SetConstraints(double argRotbounds = DBL_MAX,
						double argTransbounds = DBL_MAX,
						double argScalebounds = 0.3,
						double argSPbounds = 3.0);

  virtual void algDirICP_DIMLOP::ComputeMatchStatistics(
	  double &Avg, double &StdDev);
  //virtual void algDirICP_DIMLOP::ComputeMatchStatistics(
	 // double &PosAvg, double &PosStdDev,
	 // double &AngAvg, double &AngStdDev);

  void	UpdateShape(vctDynamicVector<double> &si);
  void	UpdateTree();

  void    UpdateOptimizerCalculations(const vctDynamicVector<double> &x);
  void    CostFunctionGradient(const vctDynamicVector<double> &x, vctDynamicVector<double> &g);
  double  CostFunctionValue(const vctDynamicVector<double> &x);

  //--- ICP Interface Methods ---//

  virtual void          ICP_InitializeParameters(vctFrm3 &FGuess);
  virtual void          ICP_UpdateParameters_PostMatch();
  virtual void          ICP_UpdateParameters_PostRegister(vctFrm3 &Freg);

  virtual void			ICP_ComputeMatches();
  virtual vctFrm3       ICP_RegisterMatches();
  virtual double        ICP_EvaluateErrorFunction();
  virtual unsigned int  ICP_FilterMatches();
  //virtual bool			ICP_Terminate(vctFrm3 &Freg);

  virtual void ReturnScale(double &scale);
  virtual void ReturnShapeParam(vctDynamicVector<double> &shapeParam);

  //--- PD Tree Interface Methods ---//

  int  NodeMightBeCloser(
    const vct3 &v, const vct3 &n,
    DirPDTreeNode const *node,
    double ErrorBound);

  virtual double FindClosestPointOnDatum(const vct3 &v, const vct3 &n,
    vct3 &closest, vct3 &closestNorm,
    int datum);

  virtual int DatumMightBeCloser(const vct3 &v, const vct3 &n,
    int datum,
    double ErrorBound);

protected:
	// -- Deformable Methods -- //
	void ComputeMu();

	void algDirICP_DIMLOP::UpdateNoiseModel_SamplesXfmd(vctFrm3 &Freg);
};
#endif
