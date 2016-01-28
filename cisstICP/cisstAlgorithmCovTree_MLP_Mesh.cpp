// ****************************************************************************
//
//    Copyright (c) 2014, Seth Billings, Russell Taylor, Johns Hopkins University
//    All rights reserved.
//
//    Redistribution and use in source and binary forms, with or without
//    modification, are NOT permitted without specific prior written permission
//    of the copyright holders.
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
#include "cisstAlgorithmCovTree_MLP_Mesh.h"
#include "utilities.h"

// finds the point on this datum with lowest match error
double cisstAlgorithmCovTree_MLP_Mesh::FindClosestPointOnDatum(const vct3 &v,
  vct3 &closest,
  int datum)
{

  static vct3 d, v0, v1, v2;
  static vct3x3 M, Minv, N, Ninv;
  double det_M;

  // compute noise model for this datum
  M = sampleXfm_M + pTree->MeshP->TriangleCov[datum];
  ComputeCovDecomposition_NonIter(M, Minv, N, Ninv, det_M);

  // Find the closest point on this triangle in a Mahalanobis distance sense
  static vct3 p0, p1, p2, c;
  pTree->TriangleVerticesCoords(datum, v0, v1, v2);
  p0 = N*(v0 - v);
  p1 = N*(v1 - v);
  p2 = N*(v2 - v);
  TCPS.FindClosestPointOnTriangle(vct3(0.0), p0, p1, p2, -1, c);
  closest = Ninv*c + v;

  d = (v - closest);
  return log(det_M) + vctDotProduct(d, Minv*d);
}


int cisstAlgorithmCovTree_MLP_Mesh::DatumMightBeCloser(const vct3 &v,
  int datum,
  double ErrorBound)
{
  return true;
}