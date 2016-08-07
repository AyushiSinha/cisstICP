// ****************************************************************************
//
//    Copyright (c) 2016, Seth Billings, Johns Hopkins University
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

#ifndef _ply_io_h
#define _ply_io_h

#include <cisstVector.h>

class ply_io
{
public:
  
  // read ply file storing the data in data arrays
  // returns  1: success  0: error
  int read_ply( const std::string &input_ply,
    vctDynamicVector<vct3>    &vertices,
    vctDynamicVector<vctInt3> &faces,
    vctDynamicVector<vct3>    &face_normals,
    vctDynamicVector<vctInt3> &face_neighbors,
    vctDynamicVector<vct3>    &vertex_normals
    );

  // write array data to ply file
  // returns  1: success  0: error
  int write_ply( const std::string &output_ply,
    const vctDynamicVector<vct3>    &vertices,
    const vctDynamicVector<vctInt3> &faces,
    const vctDynamicVector<vct3>    &face_normals,
    const vctDynamicVector<vctInt3> &face_neighbors,
    const vctDynamicVector<vct3>    &vertex_normals
    );
};


#endif // _ply_io_h
