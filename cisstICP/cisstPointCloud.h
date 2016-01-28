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
#ifndef _cisstPointCloud_h_
#define _cisstPointCloud_h_

#include <cisstVector.h>

#include "cisstMesh.h"

class cisstPointCloud
{

public:

  //--- Variables ---//

  vctDynamicVector<vct3> points;

  // optional point cloud properties
  vctDynamicVector<vct3> pointOrientations;

  // point cloud noise model
  vctDynamicVector<vct3x3>  pointCov;       // covariance of measurement noise
  vctDynamicVector<vct3>    pointCovEig;    // eigenvalues of covariance


  //--- Methods ---//
  
  // constructors
  cisstPointCloud() {};

  cisstPointCloud(vctDynamicVector<vct3> &points);

  cisstPointCloud(
    vctDynamicVector<vct3> &points,
    vctDynamicVector<vct3> &pointOrientations);

  // build point cloud from a mesh by choosing the centerpoint 
  // from each mesh triangle as a point in the cloud; also
  // the noise model of each mesh triangle is copied as the
  // noise model for each point in the point cloud
  cisstPointCloud( cisstMesh &mesh );

  // build point cloud from a mesh by choosing the centerpoint 
  // from each mesh triangle as a point in the cloud; the
  // noise model of the point cloud is initialized such that
  // the variance along the point orientation is as specified
  // and the variance perpendicular to the point orientation
  // is equal to the variance of the distance from the
  // center of each triangle to its 3 vertices
  cisstPointCloud( cisstMesh &mesh, double noisePerpPlaneSD);


  // initializes all point cloud properties to empty (default initializer);
  //  this is a useful routine to use while building a point cloud,
  //  since some point cloud properties are optional and may not be
  //  initialized by the data used to build the point cloud; calling this 
  //  ensures that unused properties are emptied rather than left with
  //  possibly invalid values
  void ResetPointCloud();

  // initializes point cloud noise model to zero (default initializer)
  void InitializeNoiseModel();


  // I/O

  // save point cloud noise model covariance matrices to file
  void SavePointCloudCov(std::string &file);

  // methods for writing / reading point cloud objects to / from a file
  int WritePointCloudToFile(std::string &filePath);
  int ReadPointCloudFromFile(std::string &filePath);
  int AppendPointCloudFromFile(std::string &filePath);

  // static methods for writing / reading point cloud points and 
  // optionall point orientations to / from file without 
  // referencing a point cloud object

  static int WritePointCloudToFile( 
    std::string &filePath, 
    vctDynamicVector<vct3> &points, 
    vctDynamicVector<vct3> &orientaitons);
  
  static int ReadPointCloudFromFile( 
    std::string &filePath, 
    vctDynamicVector<vct3> &points, 
    vctDynamicVector<vct3> &orientaitons);
  
  static int AppendPointCloudFromFile( 
    std::string &filePath, 
    vctDynamicVector<vct3> &points, 
    vctDynamicVector<vct3> &orientaitons);


  static int WritePointCloudToFile(
    std::string &filePath,
    vctDynamicVector<vct3> &points);

  static int ReadPointCloudFromFile(
    std::string &filePath,
    vctDynamicVector<vct3> &points);

  static int AppendPointCloudFromFile(
    std::string &filePath,
    vctDynamicVector<vct3> &points);

};

#endif // _cisstPointCloud_h_
