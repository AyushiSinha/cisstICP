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
#ifndef _cisstAlgorithmCovTree_h
#define _cisstAlgorithmCovTree_h

#include <cisstVector.h>
#include "cisstCovTreeNode.h"

class cisstCovTreeBase;     // forward decleration for mutual dependency


class cisstAlgorithmCovTree
{
  //
  // This is the base class for a family of Covariance Tree search algorithms
  //


  //--- Algorithm Parameters ---//

public:

  cisstCovTreeBase  *pTree;   // the covariance tree


  //--- Algorithm Methods ---//

public:

  // constructor
  cisstAlgorithmCovTree(cisstCovTreeBase *pTree);


  //--- Covariance Tree Interface Methods ---//

  // finds the point on this datum with lowest match error
  //  and returns the match error and closest point
  virtual double FindClosestPointOnDatum(
    const vct3 &sample,
    vct3 &closest,
    int datum) = 0;

  // fast check if a datum might have smaller match error than the error bound
  virtual int  DatumMightBeCloser(
    const vct3 &sample,
    int datum,
    double ErrorBound) = 0;

  // fast check if a node might contain a datum having smaller match error
  //  than the error bound
  virtual int  NodeMightBeCloser(
    const vct3 &sample,
    cisstCovTreeNode *node,
    double ErrorBound) = 0;
};

#endif