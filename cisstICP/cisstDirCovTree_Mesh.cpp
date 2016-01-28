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
#include "cisstDirCovTree_Mesh.h"


cisstDirCovTree_Mesh::cisstDirCovTree_Mesh( cisstMesh &mesh, int countThresh, double diagThresh )
 : MeshP(&mesh)
{ 
  NData = MeshP->NumTriangles();
  DataIndices = new int[NData];
  for (int i=0;i<NData;i++) 
  { 
    DataIndices[i]=i;
  }
  Top = new cisstDirCovTreeNode(DataIndices,NData,this,NULL);
  NNodes = 0; NNodes++;
  treeDepth = Top->ConstructSubtree(countThresh,diagThresh);

#ifdef DebugDirCovTree
  fprintf(debugFile, "Directional Mesh Cov Tree built: NNodes=%d  NData=%d  TreeDepth=%d\n", NumNodes(), NumData(), TreeDepth());
#endif
}

cisstDirCovTree_Mesh::~cisstDirCovTree_Mesh()
{
  if (Top) delete Top;
  if (DataIndices) delete DataIndices;
}

vct3 cisstDirCovTree_Mesh::DatumSortPoint(int datum) const
{ // use vertex 0 as the sort point
  return TriangleVertexCoord(datum,0);
}

vct3 cisstDirCovTree_Mesh::DatumNorm(int datum) const
{
  return TriangleNorm(datum);
}

void cisstDirCovTree_Mesh::EnlargeBounds(const vctFrm3& F, int datum, cisstBoundingBox& BB) const
{ for (int vx=0;vx<3;vx++)
	{ BB.Include(F*TriangleVertexCoord(datum,vx));
	}
}

void cisstDirCovTree_Mesh::PrintDatum(FILE* chan,int level,int datum)
{   int i;
	for (i=0;i<level;i++)
		fprintf(chan," ");
	int v0=MeshP->TriangleVertexIndex(datum,0); 
	int v1=MeshP->TriangleVertexIndex(datum,1);
	int v2=MeshP->TriangleVertexIndex(datum,2);
	fprintf(chan,"%5d (%4d,%4d,%4d):",datum,v0,v1,v2);
	// to std output
  //for (i=0;i<3;i++)
	//{ fprintf(stdout," ["); fprintfVct3(stdout,TriangleVertexCoords(datum,i)); fprintf(stdout,"] ");};
	//fprintf(stdout,"\n");
}


// These functions were creating an inbalanced covariance tree
//  since it was adding additional data to the centroid calculation
//  in addition to the sort point position.
//void cisstDirCovTree_Mesh::AccumulateCentroid(int datum, vct3 &sum) const
//{ 
//  vct3 s=TriangleVertexCoords(datum,0);
//  s+=TriangleVertexCoords(datum,1);
//  s+=TriangleVertexCoords(datum,2);
//  s*=(1./3.0);
//  sum +=s;
//}
//void cisstDirCovTree_Mesh::AccumulateVariances(int datum, const vct3 &mean, vctDoubleMat &C, vctDoubleMat &M) const
//{ 
//  for (int vx=0; vx<3; vx++)
//	{ 
//    vct3 p = TriangleVertexCoords(datum,vx);
//	  vct3 d = p - mean;
//    vctOuter(d,d,M); 
//	  C+=M;
//	}
//}