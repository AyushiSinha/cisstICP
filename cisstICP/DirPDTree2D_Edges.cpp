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

#include "DirPDTree2D_Edges.h"


DirPDTree2D_Edges::DirPDTree2D_Edges( 
                            const vctDynamicVector<vct2> &edgesV1,
                            const vctDynamicVector<vct2> &edgesV2,
                            const vctDynamicVector<vct2> &edgesNorm,
                            int countThresh, double diagThresh, bool bUseOBB )
{ 
  EdgeList.SetEdges( edgesV1, edgesV2, edgesNorm );
  NData = EdgeList.numEdges;

  DataIndices = new int[NData];
  for (int i=0;i<NData;i++) 
  { 
    DataIndices[i]=i;
  }

  Top = new DirPDTree2DNode(DataIndices, NData, this, NULL, bUseOBB, 0);
  NNodes = 1;
  treeDepth = Top->ConstructTree(countThresh,diagThresh);

#ifdef DebugDirPDTree2D
  fprintf(debugFile, "Directional Mesh Cov Tree built: NNodes=%d  NData=%d  TreeDepth=%d\n", NumNodes(), NumData(), TreeDepth());
#endif
}

DirPDTree2D_Edges::~DirPDTree2D_Edges()
{
  if (Top) delete Top;
  if (DataIndices) delete DataIndices;
}

vct2 DirPDTree2D_Edges::DatumSortPoint(int datum) const
{ // use first vertex as sort point
  return GetEdge(datum).V1;
}

vct2 DirPDTree2D_Edges::DatumNorm(int datum) const
{
  return GetEdge(datum).Norm;
}

void DirPDTree2D_Edges::EnlargeBounds(const vctFrm2& F, int datum, BoundingBox2D& BB) const
{ 
  BB.Include( F*GetEdge(datum).V1 );
  BB.Include( F*GetEdge(datum).V2 );
}

void DirPDTree2D_Edges::EnlargeBounds(int datum, BoundingBox2D& BB) const
{ 
  BB.Include( GetEdge(datum).V1 );
  BB.Include( GetEdge(datum).V2 );
}
