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

#include "cisstDirCovTree_PointCloud.h"

// Build tree from point cloud input
//  this constructor does not define the noise model for the point cloud;
//  the user must do this manually by calling the noise model methods of this class
cisstDirCovTree_PointCloud::cisstDirCovTree_PointCloud( vctDynamicVector<vct3> &pointCloud,
                                                      vctDynamicVector<vct3> &pointCloudNorms,
                                                      int nThresh, double diagThresh )
{ 
  points = pointCloud;
  pointNorms = pointCloudNorms;
  NData = pointCloud.size();
  DataIndices = new int[NData];
  for (int i=0;i<NData;i++) 
  { 
    DataIndices[i]=i;
  }
  Top = new cisstDirCovTreeNode(DataIndices,NData,this,NULL);
  NNodes = 0; NNodes++;
  treeDepth = Top->ConstructSubtree(nThresh,diagThresh);

#ifdef DebugDirCovTree
  fprintf(debugFile, "Directional Point Cloud Cov Tree built: NNodes=%d  NData=%d  TreeDepth=%d\n", NumNodes(), NumData(), TreeDepth());
#endif
}


// build tree from mesh, converting mesh to point cloud by choosing
//  the centerpoint from each triangle as a point in the cloud
cisstDirCovTree_PointCloud::cisstDirCovTree_PointCloud( cisstMesh &mesh, 
                                                      int nThresh, double diagThresh )
{ 
  vct3 c,v0,v1,v2;

  // build point cloud from triangle centers
  NData = mesh.NumTriangles();
  points.SetSize(NData);
  pointNorms.SetSize(NData);
  DataIndices = new int[NData];
  for (int i=0;i<NData;i++) 
  { 
    DataIndices[i]=i;
  }
  for (int i=0; i<NData; i++)
  {
    points.at(i) = mesh.Triangles.at(i).Midpoint();
    pointNorms.at(i) = mesh.Triangles.at(i).norm;
  }

  // Build tree nodes
  Top = new cisstDirCovTreeNode(DataIndices,NData,this,NULL);
  NNodes = 0; NNodes++;
  treeDepth = Top->ConstructSubtree(nThresh,diagThresh);

#ifdef DebugCovTree
  fprintf(debugFile, "Point Cloud Cov Tree built: NNodes=%d  NData=%d  TreeDepth=%d\n", NumNodes(), NumData(), TreeDepth());
#endif
}


cisstDirCovTree_PointCloud::~cisstDirCovTree_PointCloud()
{
  if (Top) delete Top;
  if (DataIndices) delete DataIndices;
}

vct3 cisstDirCovTree_PointCloud::DatumSortPoint(int datum) const
{ // the sort point is the point itself
  return points.Element(datum);
}

vct3 cisstDirCovTree_PointCloud::DatumNorm(int datum) const
{
  return pointNorms.Element(datum);
}

void cisstDirCovTree_PointCloud::EnlargeBounds(const vctFrm3& F, int datum, cisstBoundingBox& BB) const
{ 
  BB.Include(F*points.Element(datum));
}

void cisstDirCovTree_PointCloud::PrintDatum(FILE* chan,int level,int datum)
{ 
  for (int i=0;i<level;i++)
		fprintf(chan," ");
	fprintf(chan,"%5d (%f,%f,%f):",datum,points.at(datum).X(),points.at(datum).Y(),points.at(datum).Z());
}


// These functions could create an inbalanced covariance tree
//  since it was adding additional data to the centroid calculation
//  in addition to the sort point position.
//void cisstDirCovTree_PointCloud::AccumulateCentroid(int datum, vct3 &sum) const
//{ 
//  vct3 s=TriangleVertexCoords(datum,0);
//  s+=TriangleVertexCoords(datum,1);
//  s+=TriangleVertexCoords(datum,2);
//  s*=(1./3.0);
//  sum +=s;
//}
//void cisstDirCovTree_PointCloud::AccumulateVariances(int datum, const vct3 &mean, vctDoubleMat &C, vctDoubleMat &M) const
//{ 
//  for (int vx=0; vx<3; vx++)
//	{ 
//    vct3 p = TriangleVertexCoords(datum,vx);
//	  vct3 d = p - mean;
//    vctOuter(d,d,M); 
//	  C+=M;
//	}
//}