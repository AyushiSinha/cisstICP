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
#ifndef _algICP_RobustICP_PointCloud_h
#define _algICP_RobustICP_PointCloud_h

#include "algICP_RobustICP.h"
#include "PDTree_PointCloud.h"
#include "algPDTree_CP_PointCloud.h"

class algICP_RobustICP_PointCloud : public algICP_RobustICP, public algPDTree_CP_PointCloud
{ 
  //
  // This class implements the standard ICP algorithm for 
  //  a point cloud target shape.
  //

  //--- Algorithm Parameters ---//


  //--- Algorithm Methods ---//

public:

  // constructor
  algICP_RobustICP_PointCloud(
    PDTree_PointCloud *pTree, 
    vctDynamicVector<vct3> &samplePts,
    double D, double D0max)
    : algICP_RobustICP(pTree, samplePts, D, D0max),
    algPDTree_CP_PointCloud(pTree)
  {}

  // destructor
  ~algICP_RobustICP_PointCloud() {}

};

#endif
