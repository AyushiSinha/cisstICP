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

#ifndef _DirPDTreeNode_h
#define _DirPDTreeNode_h

#include <stdio.h>
#include <assert.h>

#include <cisstVector.h>

#include "BoundingBox.h"

class DirPDTreeBase;      // forward declerations for mutual dependency
class algPDTree;    //  ''


class DirPDTreeNode
{
  //
  // This class implements a node of a directed PD tree;
  //  there are no derived versions of this class.
  //


  //--- Variables ---//

public:

  DirPDTreeBase*	pMyTree;
  DirPDTreeNode*  pParent;
  DirPDTreeNode   *pLEq, *pMore;  // child nodes

  int* pDataIndices;
  int NData;
  int myDepth;

  vctFrm3 F;                // transforms world -> local node coords
  BoundingBox Bounds;  // bounding box for this node

  vct3    Navg;         // avg orientation of the datum in this node
  double  dThetaMax;    // max deviation (in radians) from the avg orientation



  //--- Methods ---//

public:

  // constructor
  DirPDTreeNode(
    int *pDataIndexArray, int numIndexes,
    DirPDTreeBase* pTree, DirPDTreeNode* pParent);

  // debug constructor
  DirPDTreeNode(double dummy) :
    pMyTree(NULL),
    pParent(NULL),
    pLEq(NULL),
    pMore(NULL),
    pDataIndices(NULL),
    NData(0)
  {};

  // destructor
  ~DirPDTreeNode();

  DirPDTreeNode* GetChildSplitNode(const vct3 &datumPos);

  // Check if a datum in this node has a lower match error than the error bound
  //  If a lower match error is found, set the new closest point, update error
  //  bound, and return the global datum index of the closest datum.
  //  Otherwise, return -1.
  int FindClosestDatum(
    const vct3 &v, const vct3 &n,
    vct3 &closestPoint, vct3 &closestPointNorm,
    double &ErrorBound,
    unsigned int &numNodesVisited, 
    unsigned int &numNodesSearched);

  int   NumData() const { return NData; };
  int   IsTerminalNode() const { return pLEq == NULL; };

  int     SortNodeForSplit();
  vctFrm3 ComputeCovFrame(int i0, int i1);
  int     ConstructSubtree(int CountThresh, double DiagThresh);

  void  AccumulateCentroid(int datum, vct3 &sum) const;
  void  AccumulateVariances(int datum, const vct3 &mean, vctDouble3x3 &C) const;

  // Return the global datum index of the ith datum in this node
  int Datum(int i) const
  {
    assert(i >= 0 && i < NData);
    return pDataIndices[i];
  };
  int& Datum(int i)
  {
    assert(i >= 0 && i < NData);
    return pDataIndices[i];
  };

protected:

  // computes orientation parameters for this node and
  // returns the running sum of orientations in the node
  //   (return value is for speedup when computing parent node params)
  vct3 ComputeOrientationParams();

public:

  // debug routines
  int   FindTerminalNode(int datum, DirPDTreeNode **termNode);
  void  PrintTerminalNodes(std::ofstream &fs);
  bool  NodeContainsDatum(int datum)
  {
    for (int i = 0; i < NData; i++)
    {
      if (datum == pDataIndices[i]) return true;
    }
    return false;
  }

  //void  RecomputeBoundingBoxesUsingExistingCovFrames();
  //void  ReSortUsingExistingCovFrames();
  //void  RecomputeCovFrames();

  //virtual void Print(FILE* chan, int indent);

};

#endif

//// Check if a datum in this node has a lower match error than the error bound,
////  and an, in addition, subject to the Euclidean distance constraint from the sample.
////  If a lower match error is found, set the new closest point, update error
////  bound, and return the global datum index of the closest datum.
////  Otherwise, return -1.
////  Returns datum closer than "ErrorBound" to v subject to a geometric distance
////  constraint, else return -1
////  Sets the new closest point and updates error bound
//int FindClosestDatum( const vct3 &v, const vct3 &n,
//                       vct3 &closestPoint, vct3 &closestPointNorm,
//                       double posWeight, double normWeight,
//                       double &ErrorBound,
//                       unsigned int &numNodesVisited,
//                       unsigned int &numNodesSearched,
//                       double geomDistBound );