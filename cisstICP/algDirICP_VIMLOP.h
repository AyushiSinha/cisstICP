// ****************************************************************************
//
//    Copyright (c) 2016, Ayushi Sinha, Seth Billings, Russell Taylor, Johns Hopkins University
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
#ifndef _algDirICP_VIMLOP_h
#define _algDirICP_VIMLOP_h

#include "algICP_IMLP_Mesh.h"
#include "alg2D_DirICP.h"
#include "alg2D_DirPDTree_vonMises_Edges.h"
#include "Ellipsoid_OBB_Intersection_Solver.h"
#include "algDirICP_VIMLOP_dlibWrapper.h"
#include "utilities.h"
#include "camera.h"
#include <limits.h>

class algDirICP_VIMLOP : public algICP_IMLP_Mesh, public alg2D_DirPDTree_vonMises_Edges//, public camera
{
  //
  // This class implements algorithms for Iterative Most Likely Point
  //  registration, which is a modification of ICP to compute optimal registration
  //  in the presence of anisotropic noise.
  // Key differences in this algorithm compared to ICP are:
  //  1) total least squares framework for P2P registration (find registration that maximizes the match likelihood)
  //  2) assigned correspondences based on most likely match rather than closest match
  //
  //   likelihood function:   Product_i[ (1/sqrt((2*pi)^3*|Mi|) * exp( (-1/2)*(Yi-R*Xi-t)'*inv(Mi)*(Yi-R*Xi-t) ) ]
  //       where Mi  =  covariance of error in Yi-R*Xi-t  (Mi = R*Mxi*R' + Myi)
  //             Mxi =  covariance of measurement error for source point xi
  //             Myi =  covariance of measurement error for target point yi
  //
  //   Match Error function:  E = log|Mi| + (Yi-R*Xi-t)'*inv(Mi)*(Yi-R*Xi-t)
  //   (derived from negative log probability that Yi corresponds to R*Xi+t given the noise model Mi)
  //
  //   NOTE: negative values of the match error function 
  //         are possible (in case of very small |Mi|)
  //

public:

	algDirICP_VIMLOP_dlibWrapper dlib;

	// -- Optimizer calculations common to both cost and gradient function
	vct7 x_prev;
	vct3 a, t;
	double s;
	vctRot3 Ra;
	camera cam;

	vctDynamicVector<vct3>   Ysfm_t;
	vctDynamicVector<vct3>   Rat_Ysfm_t_x;
	vctDynamicVector<vct3>   R_Rat_Y3dp_t_st;
	vctDynamicVector<vct3>   invMx_Rat_Ysfm_t_x;

	vctDynamicVector<vct3>   Y3dp_t;
	vctDynamicVector<vct2>   F_R_Rat_Y3dp_t_x_st_2d_Xxfmd;
	vctDynamicVector<vct2>   invMx_F_R_Rat_Y3dp_t_x_st_2d_Xxfmd;

	vctDynamicVector<vct2>		S_R_Rat_Y3dn_2d;
	vctDynamicVector<double>	k_S_R_Rat_Y3dn_2d_Xxfmd;
  //-- Algorithm Parameters --//

protected:

  Ellipsoid_OBB_Intersection_Solver IntersectionSolver;

  vctFrm3 FGuess; // intiial guess for registration

  bool bFirstIter_Matches;  // flag that the first iteration is being run

  // dynamic noise model
  double sigma2;      // match uncertainty (added to My covariances as sigma2*I)
  double sigma2Max;   // max threshold on match uncertainty
  vctDynamicVector<vct2>  residuals_PostMatch;    // Pclosest - Psample
  vctDoubleVec            sqrDist_PostMatch;      // ||Pclosest - Psample||^2

  // noise model
  // Note: Myi values are obtained from the mesh object
  vctDynamicVector<vct3x3>  Mxi;        // noise covariances of sample points
  vctDynamicVector<vct3x3>  invMxi;
  vctDynamicVector<vct2>    eigMxi;     // eigenvalues of sample covariances
  vctDynamicVector<vct2x2>  R_Mxi_Rt;   // noise covariances of transformed sample points
  vctDynamicVector<vct3x3>  R_invMxi_Rt;
  vctDynamicVector<vct3x3*> Myi;        // noise covariances of target correspondence points
  vctDynamicVector<vct3x3>  Myi_sigma2; // noise covariances of target correspondence points with match uncertainty added
  //vctDynamicVector<vct3x3>  Mi;       // noise covariances of match (R*Mxi*Rt + Myi)
  //vctDynamicVector<vct3x3>  inv_Mi;   // inverse noise covariances of match (R*Mxi*Rt + Myi)^-1
  //vctDynamicVector<double>  det_Mi;   // determinant of noise covariances of match |R*Mxi*Rt + Myi|
  //vctDynamicVector<double>  SqrMahalDist;  // square Mahalanobis distance of matches = (yi-Rxi-t)'*inv(Mi)*(yi-Rxi-t)

  // outlier handling
  unsigned int nOutliers;
  double ChiSquareThresh; // Chi Square threshold for the outlier test  
  double sumSqrDist_Inliers;
  vctDynamicVector<int>   outlierFlags;
  // measurement noise component of the sample noise model
  //  Note: if applying measurement noise to the target shape, then
  //        MsmtMyi should be used here as well
  // 2D
  vctDynamicVector<vct2x2> MsmtMxi;
  vctDynamicVector<vct2x2> NormMxi;
  // covariance model for outlier tests
  //  (does not include planar noise model)
  vctDynamicVector<vct2x2> R_MsmtMxi_Rt;
  double concentration;

  // algorithm-specific termination
  bool bTerminateAlgorithm;
  unsigned char costFuncIncBits;
  double costFuncValue, prevCostFuncValue, prevIncCostFuncValue;
  vctFrm3 Fdec;
  
  // Temporary Match Variables
  //
  // These variables are recomputed each time a new sample point 
  //  is used in a PD tree search
  // Variables for current sample point undergoing a match search
  //  these are set once for each sample point by the pre-match function
  //  and used whenever a tree search routine is called thereafter.
  //  These do not change relative to the node being searched => precomputing
  //  these avoids multiple recomputations.
  vct2x2 sample_RMxRt_sigma2;      // covariance of transformed sample point plus registration uncertainty term
  vct2   sample_RMxRt_sigma2_Eig;  // eigenvalues of sample_RMxRt_sigma2 in order of decreasing magnitude

  vct2x2 M;       // effective measurement error covariance for a node & sample pair
  vct2x2 N;       // decomposition of inv(M) = N'N
  double Dmin;    // inverse sqrt of largest eigenvalue of M
                  //  (or the sqrt of the smallest eigenvalue of inv(M))
  double MinLogM; // lower bound on the log component of error for this node


  //-- Algorithm Methods --//

public:

  // constructor
  algDirICP_VIMLOP(
    PDTree_Mesh *pTree,
	DirPDTree2D_Edges *pDirTree,
    vctDynamicVector<vct3> &samplePts, 
    vctDynamicVector<vct3x3> &sampleCov,      // full noise model (measurement noise + surface model)
    vctDynamicVector<vct3x3> &sampleMsmtCov,  // partial noise model (measurement noise only)
	vctDynamicVector<vct2> &sampleEdgesV1,
	vctDynamicVector<vct2> &sampleEdgesV2,
	vctDynamicVector<vct2> &sampleEdgesNorm,
	vctDynamicVector<vct2x2> &sample2dMsmtCov,
	vctDynamicVector<vct2x2> &sample2dNormCov,
	camera cam,
    double outlierChiSquareThreshold = 7.81,
    double sigma2Max = std::numeric_limits<double>::max());

  virtual void  ComputeMatchStatistics(double &Avg, double &StdDev);

  // TODO: specify surface model independently and compute Mxi
  //       by adding surface and msmt covariances
  //  Mxi ~ full noise model (measurement noise + surface model)
  //  MsmtMxi ~ partial noise model (measurement noise only)
  virtual void  SetSamples(
    vctDynamicVector<vct3> &argSamplePts, 
    vctDynamicVector<vct3x3> &argMxi, 
	vctDynamicVector<vct3x3> &argMsmtMxi,
	vctDynamicVector<vct2> &argSampleEdgesV1,
	vctDynamicVector<vct2> &argSampleEdgesV2,
	vctDynamicVector<vct2> &argSampleEdgesNorm,
	vctDynamicVector<vct2x2> &arg2dMsmtMxi,
	vctDynamicVector<vct2x2> &arg2dNormMxi);

  //void SetSampleCovariances(vctDynamicVector<vct3x3> &Mi, vctDynamicVector<vct3x3> &MsmtMi);

  // Sets Chi Square threshold for the outlier test:
  // Note:    ChiSquare(0.6) = 2.95      
  //          ChiSquare(0.7) = 3.66      
  //          ChiSquare(0.8) = 4.64      
  //          ChiSquare(0.85) = 5.32     
  //          ChiSquare(0.9) = 6.25      
  //          ChiSquare(0.925) = 6.90    
  // default  ChiSquare(0.95) = 7.81     (1.96 Std Dev)
  //          ChiSquare(0.975) = 9.35    (2.24 Std Dev)
  //          ChiSquare(0.99) = 11.34    (2.56 Std Dev)
  //          ChiSquare(0.9973) = 14.16  (3.0 Std Dev)     MATLAB: chi2inv(0.9973,3)
  void SetChiSquareThreshold(double ChiSquareValue) { ChiSquareThresh = ChiSquareValue; }
  void SetSigma2Max(double sigma2MaxValue) { sigma2Max = sigma2MaxValue; }

  void    UpdateOptimizerCalculations(const vct7 &x); 
  void    CostFunctionGradient(const vct7 &x, vct7 &g);
  double  CostFunctionValue(const vct7 &x);

protected:

  void UpdateNoiseModel_SamplesXfmd(vctFrm3 &Freg);

  void ComputeNodeMatchCov(PDTreeNode *node, DirPDTree2DNode *dirNode);

  void ComputeCovDecomposition_NonIter(const vct3x3 &M, vct3x3 &Minv, double &det_M);
  void ComputeCovDecomposition_NonIter(const vct3x3 &M, vct3x3 &Minv, vct3x3 &N, vct3x3 &Ninv, double &det_M);

  void ComputeCovDecomposition_SVD(const vct3x3 &M, vct3x3 &Minv, double &det_M);
  void ComputeCovDecomposition_SVD(const vct3x3 &M, vct3x3 &Minv, vct3x3 &N, vct3x3 &Ninv, double &det_M);

  void algDirICP_VIMLOP::ComputeCovDecomposition_NonIter(const vct2x2 &M, vct2x2 &Minv, double &det_M);

  bool IntersectionSphereFace(const vct3 &n,
    const vct3 &v0, const vct3 &v1,
    const vct3 &v2, const vct3 &v3,
    double radius, double sqrRadius);

  int FindVisibleEdges(double q0, double q1,
    double v00, double v01,
    double v10, double v11,
    double v20, double v21,
    double v30, double v31,
    int *vsblEdges,
    bool ccwSequence);

  inline bool EdgeIsVisible(double q0, double q1,
    double v0, double v1,
    double n0, double n1);

  inline double SquareDistanceToEdge(const vct3 &p, const vct3 &r);

  // virtual standard routines for matching
  //void SamplePreMatch(unsigned int sampleIndex);


  //--- ICP Interface Methods ---//

public:

  virtual void    ICP_InitializeParameters(vctFrm3 &FGuess);
  virtual void    ICP_UpdateParameters_PostMatch();
  virtual void    ICP_UpdateParameters_PostRegister(vctFrm3 &Freg);

  virtual void ICP_ComputeMatches();
  virtual vctFrm3 ICP_RegisterMatches();
  virtual unsigned int ICP_FilterMatches();  

  virtual double  ICP_EvaluateErrorFunction();
  //virtual bool    ICP_Terminate(vctFrm3 &Freg);

  //virtual void  ICP_ComputeMatches();
  //virtual std::vector<cisstICP::Callback> ICP_GetIterationCallbacks();


  //--- PD Tree Interface Methods ---//

  int  NodeMightBeCloser(
    const vct3 &v,
    PDTreeNode *node,
	double ErrorBound);

  int  NodeMightBeCloser(
	  const vct2 &v,
	  DirPDTree2DNode *node,
	  double ErrorBound);

  double FindClosestPointOnDatum(
    const vct3 &v,
    vct3 &closest,
    int datum);

  int  DatumMightBeCloser(
    const vct3 &v,
    int datum,
    double ErrorBound);

  int  DatumMightBeCloser(
	  const vct2 &sample, const vct2 &sampleNorm,
	  int datum,
	  double ErrorBound);

  // finds the point on this datum with lowest match error
  //  and returns the match error and closest point
  double FindClosestPointOnDatum(
	  const vct2 &sample, const vct2 &sampleNorm,
	  vct2 &closest, vct2 &closestNorm,
	  int datum);

};
#endif
