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
#ifndef _cisstPointCloud2D_h_
#define _cisstPointCloud2D_h_

#include <stdio.h>
#include <cisstVector.h>
#include <cisstCommon.h>

class cisstPointCloud2D
{

public:

  vctDynamicVector<vct2> points;

  // optional point cloud properties
  vctDynamicVector<vct2> pointOrientations;

  // point cloud noise model
  vctDynamicVector<vct2x2>  pointCov;       // covariance of measurement noise
  vctDynamicVector<vct2>    pointCovEig;    // eigenvalues of covariance
  
  //--- Methods ---//

  // constructors
  cisstPointCloud2D() {};

  cisstPointCloud2D( vctDynamicVector<vct2> &points) :
    points(points)
  {};

  cisstPointCloud2D(
	  vctDynamicVector<vct2> &points,
	  vctDynamicVector<vct2> &pointOrientations);


  // Point Set I/O
  static int WritePointCloudToFile(vctDynamicVector<vct2> &points, std::string &filePath);
  static int ReadPointCloudFromFile(vctDynamicVector<vct2> &points, std::string &filePath);
  static int AppendPointCloudFromFile(vctDynamicVector<vct2> &points, std::string &filePath);

  static int WritePointCloudToFile(vctDynamicVector<vct2> &points, vctDynamicVector<vct2> &orientaitons, std::string &filePath);
  static int ReadPointCloudFromFile(vctDynamicVector<vct2> &points, vctDynamicVector<vct2> &orientaitons, std::string &filePath);
  static int AppendPointCloudFromFile(vctDynamicVector<vct2> &points, vctDynamicVector<vct2> &orientaitons, std::string &filePath);

  int WritePointCloudToFile(std::string &filePath)
  {
    WritePointCloudToFile(points, filePath);
  }

  int ReadPointCloudFromFile(std::string &filePath)
  {
    ReadPointCloudFromFile(points, filePath);
  }

  int AppendPointCloudFromFile(std::string &filePath)
  {
    AppendPointCloudFromFile(points, filePath);
  }

};

#endif // _cisstPointCloud2D_h_
