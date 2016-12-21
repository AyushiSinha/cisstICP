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
#ifndef _alg2D_DirPDTree_vonMises_Edges_h
#define _alg2D_DirPDTree_vonMises_Edges_h

//#include "alg2D_DirPDTree_vonMises.h"
#include "DirPDTree2D_Edges.h"

class alg2D_DirPDTree_vonMises_Edges : public DirPDTree2D_Edges
{
  //
  // Implements base class algorithms for a 2D edge shape
  //  (i.e. edge datum type)
  //

  //--- Algorithm Parameters ---//

protected:

  DirPDTree2D_Edges *pDirTree;

public:

  // This is no longer needed
  //// temporary buffer storage used to help determine matchLambdas
  ////  (since current architecture cannot store the matchLambdas directly)
  ////  matchPt = (edgeV2-edgeV1)*matchLambda + edgeV1
  ////  (stores matchLambda for each datum searched)
  //vctDoubleVec    searchLambdas;

	unsigned int nSamples;
	unsigned int nOutliers;

	vctDynamicVector<vct2>  sampleEdgesV1;
	vctDynamicVector<vct2>  sampleEdgesV2;
	vctDynamicVector<vct2>  sampleEdgesNorm;

	vctDynamicVector<vct2>  samplePtsXfmd;
	vctDynamicVector<vct2>  sampleNormsXfmd;

	double sigma2;  // variance of position
	double k;       // concentration of orientation

	vctDynamicVector<vct3>  matchPts;
	vctDynamicVector<vct3>  matchNorms;
	vctDynamicVector<int>  matchDatums;
	vctDoubleVec  matchErrors;

	vctDynamicVector<bool>	matchPermitted;

	int   minNodesSearched, maxNodesSearched, avgNodesSearched;

	double dThetaMax;               // maximum orientation error permitted for a match
	bool   bPermittedMatchFound;    // true if a permitted match found (i.e. w/in dThetaMax)
	//  this is used to search all matches until a permitted match is found
	//  if no permitted match is found, then the best match from the
	//  set of non-permitted matches is returned

#ifdef ValidatePDTreeSearch
	std::ofstream validFS;
	vct2 validPoint;
	vct2 validNorm;
	int  validDatum;
	double validDist;
	double validAng;
	double validError;
	double searchError;
	unsigned int numValidDatums;
	unsigned int numInvalidDatums;
	double validPercent;
	double doubleEps = 1e-16;
	int validIter = 0;
#endif

  //--- Algorithm Methods ---//

public:

  // constructor
  alg2D_DirPDTree_vonMises_Edges(
	  DirPDTree2D_Edges *pDirTree,
	  vctDynamicVector<vct2> &edgesV1,
	  vctDynamicVector<vct2> &edgesV2,
	  vctDynamicVector<vct2> &edgesNorm, 
	  int nThresh = 6, double diagThresh = 0.0, bool bUseOBB = false,
	  double k = 1.0, double sigma2 = 1.0, double thetaMax = cmnPI)
	  : DirPDTree2D_Edges(edgesV1, edgesV2, edgesNorm, nThresh, diagThresh, bUseOBB),
	  k(k), sigma2(sigma2), dThetaMax(thetaMax), bPermittedMatchFound(false)
  {
	  SetSamples(edgesV1, edgesV2, edgesNorm);
	  /*searchLambdas.SetSize(pDirTree->EdgeList.numEdges);*/ 
#ifdef ValidatePDTreeSearch
	  validFS.open("debugPDTreeSearchValidation.txt");
#endif
  }

  void SetSamples(vctDynamicVector<vct2> argEdgesV1, vctDynamicVector<vct2> argEdgesV2, vctDynamicVector<vct2> argEdgesNorm)
  {
	  if (argEdgesV1.size() != argEdgesV2.size())
	  {
		  CISST_THROW("ERROR: edge vertices should be the same size");
	  }

	  if (argEdgesV1.size() != argEdgesNorm.size())
	  {
		  CISST_THROW("ERROR: edges and normals are different sizes");
	  }

	  // copy sample points
	  sampleEdgesV1 = argEdgesV1;
	  sampleEdgesV2 = argEdgesV2;
	  sampleEdgesNorm = argEdgesNorm;

	  nSamples = sampleEdgesNorm.size();

	  // initialize variables dependent on sample size
	  samplePtsXfmd.SetSize(nSamples);
	  sampleNormsXfmd.SetSize(nSamples);

	  matchPts.SetSize(nSamples);
	  matchNorms.SetSize(nSamples);

	  matchDatums.SetSize(nSamples);
	  matchErrors.SetSize(nSamples);

	  matchPermitted.SetSize(nSamples);
  }

  // destructor
  virtual ~alg2D_DirPDTree_vonMises_Edges() {}

  virtual void ICP_InitializeParameters(vctFrm2 &FGuess);

  virtual void ICP_ComputeMatches();

  //--- PD Tree Interface Methods ---//

  // fast check if a datum might have smaller match error than error bound
  virtual int  DatumMightBeCloser(
    const vct2 &sample, const vct2 &sampleNorm,
    int datum, double ErrorBound);

  // finds the point on this datum with lowest match error
  //  and returns the match error and closest point
  virtual double FindClosestPointOnDatum(
    const vct2 &sample, const vct2 &sampleNorm,
    vct2 &closest, vct2 &closestNorm,
    int datum);
};

#endif
