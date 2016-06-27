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
#ifndef _alg2D_DirPDTree_vonMises_h
#define _alg2D_DirPDTree_vonMises_h

#include "alg2D_DirPDTree.h"
#include "DirPDTree2DBase.h"

class alg2D_DirPDTree_vonMises : public alg2D_DirPDTree
{
  //
  // Implements the isotropic oriented point algorithm for 2D PD tree search
  //   using the isotropic vonMises and isotropic Gaussian distributions
  //

public:

  // Noise Parameters
  double k;       // concentration of orientation
  double sigma2;  // variance of position
  
  double dThetaMax;               // maximum orientation error permitted for a match
  bool   bPermittedMatchFound;    // true if a permitted match found (i.e. w/in dThetaMax)
                                  //  this is used to search all matches until a permitted match is found
                                  //  if no permitted match is found, then the best match from the
                                  //  set of non-permitted matches is returned
  //bool  bPermittedMatchOverride;  // override to return errors of non-permitted matches
  
protected:

  // constructor
  alg2D_DirPDTree_vonMises(
    DirPDTree2DBase *pDirTree, 
    double k = 1.0, double sigma2 = 1.0, double thetaMax = cmnPI)
    : alg2D_DirPDTree(pDirTree),
    k(k), sigma2(sigma2), dThetaMax(thetaMax), bPermittedMatchFound(false)
  {}

  // destructor
  virtual ~alg2D_DirPDTree_vonMises() {}

  //// Helper Methods
  //int FastInitializeProximalDatum(
  //  const vct2 &v, const vct2 &n,
  //  vct2 &proxPoint, vct2 &proxNorm);


  //--- PD Tree Interface Methods ---//

  // fast check if a node might contain a datum having smaller match error
  //  than the error bound
  virtual int  NodeMightBeCloser(
    const vct2 &sample, const vct2 &sampleNorm,
    DirPDTree2DNode const *node,
    double ErrorBound);

  // fast check if a datum might have smaller match error than error bound
  virtual int  DatumMightBeCloser(
    const vct2 &sample, const vct2 &sampleNorm,
    int datum,
    double ErrorBound) = 0;

  // finds the point on this datum with lowest match error
  //  and returns the match error and closest point
  virtual double FindClosestPointOnDatum(
    const vct2 &sample, const vct2 &sampleNorm,
    vct2 &closest, vct2 &closestNorm,
    int datum) = 0;
};

#endif
