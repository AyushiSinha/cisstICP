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
#ifndef _cisstAlgorithm2D_DirICP_StdICP_h
#define _cisstAlgorithm2D_DirICP_StdICP_h

#include "cisstAlgorithm2D_DirICP.h"
#include "cisstDirCovTree2DBase.h"
#include "cisstTriangleClosestPointSolver.h"

class cisstAlgorithm2D_DirICP_StdICP 
  : public cisstAlgorithm2D_DirICP  //, public cisstAlgorithm2D_DirCovTree
{
  //
  // This is a base class for the standard ICP algorithm
  //


  //--- Algorithm Parameters ---//

public:

  //
  //  Note: using buffers and buffer references for match filtering enables
  //        resizing the set of "good" points without altering memory allocation
  //

  // good samples (outliers removed)
  unsigned int  nGoodSamples;

  vctDynamicVectorRef<vct2>     goodSamplePts;  
  vctDynamicVectorRef<vct2>     goodSampleNorms;
  vctDynamicVectorRef<vct2>     goodMatchPts;
  vctDynamicVectorRef<vct2>     goodMatchNorms;

  vctDynamicVector<vct2>        goodSamplePtsBuf;  
  vctDynamicVector<vct2>        goodSampleNormsBuf;
  vctDynamicVector<vct2>        goodMatchPtsBuf;
  vctDynamicVector<vct2>        goodMatchNormsBuf;

  //// summary error statistics
  //double gSumSqrDist_PostMatch;   // good sample distances
  //double oSumSqrDist_PostMatch;   // thresholded distances on outlier samples (helpful for smoothing cost function)


  //--- Algorithm Methods ---//

public:

  // constructor
  cisstAlgorithm2D_DirICP_StdICP(
    cisstDirCovTree2DBase *pDirTree,
    vctDynamicVector<vct2> &samplePts,
    vctDynamicVector<vct2> &sampleNorms)
    : cisstAlgorithm2D_DirICP(pDirTree, samplePts, sampleNorms)
    //cisstAlgorithm2D_DirCovTree(pDirTree)
  {
    SetSamples(samplePts, sampleNorms);
  };

  void SetSamples(
    vctDynamicVector<vct2> &samplePts,
    vctDynamicVector<vct2> &sampleNorms)
  {
    cisstAlgorithm2D_DirICP::SetSamples(samplePts, sampleNorms);

    // size sample buffers
    goodSamplePtsBuf.SetSize(nSamples);
    goodSampleNormsBuf.SetSize(nSamples);
    goodMatchPtsBuf.SetSize(nSamples);
    goodMatchNormsBuf.SetSize(nSamples);
  }


  //--- ICP Interface Methods ---//

  void          ICP_InitializeParameters(vctFrm2 &FGuess);
  vctFrm2       ICP_RegisterMatches();
  double        ICP_EvaluateErrorFunction();
  unsigned int  ICP_FilterMatches();

  //void  ICP_UpdateParameters_PostMatch();
  //void  ICP_UpdateParameters_PostRegister(vctFrm2 &Freg);
  //void  ICP_ComputeMatches();
  //std::vector<cisstICP::Callback> ICP_GetIterationCallbacks();


  ////--- Covariance Tree Interface Methods ---//

  //int  NodeMightBeCloser(
  //  const vct2 &v, const vct2 &n,
  //  cisstDirCovTree2DNode const *node,
  //  double ErrorBound);

  //virtual double FindClosestPointOnDatum(
  //  const vct2 &v, const vct2 &n,
  //  vct2 &closest, vct2 &closestNorm,
  //  int datum) = 0;

  //virtual int DatumMightBeCloser(
  //  const vct2 &v, const vct2 &n,
  //  int datum,
  //  double ErrorBound) = 0;
};
#endif
