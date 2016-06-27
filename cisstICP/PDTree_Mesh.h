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
#ifndef _PDTree_Mesh_h
#define _PDTree_Mesh_h

#include "PDTreeBase.h"
#include "cisstMesh.h"
#include <limits>

class PDTree_Mesh : public PDTreeBase
{ 
  //
  // This class implements a PD tree for a mesh shape
  //


  //--- Variables ---//

public:

  cisstMesh *MeshP;


  //--- Methods ---//

public:

  // constructor
  //  mesh       - target shape from which to construct the tree
  //               (each triangle of the mesh becomes a datum in the tree)
  //  nThresh    - min number of datums to subdivide a node
  //  diagThresh - min physical size to subdivide a node
	PDTree_Mesh(cisstMesh &mesh, int nThresh, double diagThresh);

  // destructor
  virtual ~PDTree_Mesh();


  //--- Base Class Virtual Methods ---//

  virtual vct3 DatumSortPoint(int datum) const;

  virtual void EnlargeBounds(const vctFrm3& F, int datum, BoundingBox& BB) const;


  //--- Noise Model Methods ---//

  vct3x3& DatumCov(int datum)       // return measurement noise model for this datum
  {
    return MeshP->TriangleCov[datum];
  }
  vct3x3* DatumCovPtr(int datum)    // return measurement noise model for this datum
  {
    return &(MeshP->TriangleCov[datum]);
  }

  vct3& DatumCovEig(int datum)         // return measurement noise model for this datum
  {
    return MeshP->TriangleCovEig[datum];
  }
  vct3* DatumCovEigPtr(int datum)      // return measurement noise model for this datum
  {
    return &(MeshP->TriangleCovEig[datum]);
  }

};

#endif
