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
#include "alg2D_DirPDTree_CP_Edges.h"

double alg2D_DirPDTree_CP_Edges::FindClosestPointOnDatum(
  const vct2 &v, const vct2 &n,
  vct2 &closest, vct2 &closestNorm,
  int datum)
{
  // Project point onto the edge
  //  (storing lambda value in temp buffer)
  closest = pDirTree->GetEdge(datum).ProjectOnEdge(v, &searchLambdas.Element(datum));
  closestNorm = pDirTree->GetEdge(datum).Norm;

  // return distance as match error
  //  NOTE: distance is more convenient than square distance
  //        since the distance is required for doing the
  //        bounding box proximity tests.
  return (v - closest).Norm();
}

int alg2D_DirPDTree_CP_Edges::DatumMightBeCloser(
  const vct2 &v, const vct2 &n,
  int datum,
  double ErrorBound)
{
  // doing a decent proximity check is complicated enough that it is
  //  better to just compute the full error directly
  return 1;
}
