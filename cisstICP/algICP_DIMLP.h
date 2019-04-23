// ****************************************************************************
//
//    Copyright (c) 2017, Ayushi Sinha, Seth Billings, Russell Taylor, Johns Hopkins University.
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
#ifndef _algICP_DIMLP_h
#define _algICP_DIMLP_h

#include "algICP_IMLP.h"
#include "PDTree_Mesh.h"
#include "TriangleClosestPointSolver.h"
#include "Ellipsoid_OBB_Intersection_Solver.h"
#include "algICP_DIMLP_dlibWrapper.h"
#include "utilities.h"
#include <limits.h>

class algICP_DIMLP : public algICP_IMLP
{
	//
	// This class implements algorithms for Deformable Iterative Most Likely Point
	// registration, which is a modification of ICP to compute optimal registration
	// in the presence of deformation.
	// 


  //--- Algorithm Parameters ---//
public:
	algICP_DIMLP_dlibWrapper dlib;

	vctDynamicVector<vct3>		Tssm_matchPts;
	vctDynamicVector<vct3>		mu;
	vctDynamicVector<vctInt3>	f;
	vctDynamicVector<vctDynamicVector<vct3>>	Tssm_wi;

	double rb, tb, sb, spb;		// rotation, translation, scale, and shape parameter bounds
	bool bScale;				// boolean for optional scale optimization

	// registration statistics
	double totalSumSqrMahalDist;
	double sumSqrMahalDist;

	int nGoodSamples;

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

protected:
	TriangleClosestPointSolver TCPS;
	PDTree_Mesh *pTree;
	cisstMesh *pMesh;

	// Deformable Variables
	vctDynamicVector<vct3>	meanShape;
	unsigned int nTrans; // # transformation parameters
	unsigned int nModes; // # shape parameters
	//vctDynamicVector<vct3>	sampleModes;
	//vctDynamicVector<double> sampleModeWts;	

	vctDynamicVector<vctDynamicVector<vct3>>		wi;
	vctDynamicVector<double>	Si;		// shape parameter


	// -- Algorithm Methods -- //
public:
	
	// constructor
	algICP_DIMLP(
		PDTree_Mesh *pTree,
		vctDynamicVector<vct3> &samplePts,
		vctDynamicVector<vct3x3> &sampleCov,      // full noise model (measurement noise + surface model)
		vctDynamicVector<vct3x3> &sampleMsmtCov,  // partial noise model (measurement noise only)
		vctDynamicVector<vct3> &meanShape,
		double scale = 1.0, bool bScale = false,
		double outlierChiSquareThreshold = 7.81,
		double sigma2Max = std::numeric_limits<double>::max());

	virtual ~algICP_DIMLP() {}

	virtual void SetSamples(
		vctDynamicVector<vct3> &argSamplePts,
		vctDynamicVector<vct3x3> &argMxi,
		vctDynamicVector<vct3x3> &argMsmtMxi,
		vctDynamicVector<vct3> &argmeanShape,
		double argScale = 1.0, 
		bool argbScale = false);

	void SetConstraints(double argRotbounds = DBL_MAX,
						double argTransbounds = DBL_MAX,
						double argScalebounds = 0.3,
						double argSPbounds = 3.0);

	virtual void ComputeMatchStatistics(double &Avg, double &stdDev);
	virtual void  PrintMatchStatistics(std::stringstream &tMsg);
	void	UpdateShape(vctDynamicVector<double> &si);
	void	UpdateTree();
	void    UpdateOptimizerCalculations(const vctDynamicVector<double> &x);
	void    CostFunctionGradient(const vctDynamicVector<double> &x, vctDynamicVector<double> &g);
	double	CostFunctionValue(const vctDynamicVector<double> &x);

protected:
	// -- Deformable Methods -- //
	void ComputeMu(); 

	// -- ICP Interface Methods -- //

public:
	virtual void ICP_InitializeParameters(vctFrm3 &FGuess);
	virtual void ICP_UpdateParameters_PostMatch();
	virtual void ICP_UpdateParameters_PostRegister(vctFrm3 &Freg);

	virtual void ICP_ComputeMatches();
	virtual vctFrm3 ICP_RegisterMatches();
	//virtual unsigned int ICP_FilterMatches();

	virtual double  ICP_EvaluateErrorFunction();
	virtual bool    ICP_Terminate(vctFrm3 &Freg);

	virtual void ReturnScale(double &scale);
	virtual void ReturnShapeParam(vctDynamicVector<double> &shapeParam);
	virtual void ReturnMatchPts(vctDynamicVector<vct3> &matchPts, vctDynamicVector<vct3> &matchNorms);

	// -- PD Tree Interface Methods -- //
	//int  NodeMightBeCloser(
	//	const vct3 &v,
	//	PDTreeNode *node,
	//	double ErrorBound);

	virtual double FindClosestPointOnDatum(
		const vct3 &v,
		vct3 &closest,
		int datum);

	virtual int  DatumMightBeCloser(
		const vct3 &v,
		int datum,
		double ErrorBound);
};
#endif