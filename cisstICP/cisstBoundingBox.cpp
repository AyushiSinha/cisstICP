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
#include "cisstBoundingBox.h"

cisstBoundingBox& cisstBoundingBox::Include(const vct3& V)
{ if (MinCorner[0] > V[0]) MinCorner[0]=V[0];
  if (MinCorner[1] > V[1]) MinCorner[1]=V[1];
  if (MinCorner[2] > V[2]) MinCorner[2]=V[2];
  if (MaxCorner[0] < V[0]) MaxCorner[0]=V[0];
  if (MaxCorner[1] < V[1]) MaxCorner[1]=V[1];
  if (MaxCorner[2] < V[2]) MaxCorner[2]=V[2];
  // TODO: would be faster to call this once after calling all includes
  //       but more complicated --> putting this here for now.
  //ComputeHalfExtents();
  return *this;
};

cisstBoundingBox& cisstBoundingBox::Include(const cisstBoundingBox& him)
{ if (MinCorner[0] > him.MinCorner[0]) MinCorner[0]=him.MinCorner[0];
  if (MinCorner[1] > him.MinCorner[1]) MinCorner[1]=him.MinCorner[1];
  if (MinCorner[2] > him.MinCorner[2]) MinCorner[2]=him.MinCorner[2];
  if (MaxCorner[0] < him.MaxCorner[0]) MaxCorner[0]=him.MaxCorner[0];
  if (MaxCorner[1] < him.MaxCorner[1]) MaxCorner[1]=him.MaxCorner[1];
  if (MaxCorner[2] < him.MaxCorner[2]) MaxCorner[2]=him.MaxCorner[2];
  // TODO: would be faster to call this once after calling all includes
  //       but more complicated --> putting this here for now.
  //ComputeHalfExtents();
  return *this;
};

//void cisstBoundingBox::ComputeHalfExtents()
//{
//  HalfExtents = (MaxCorner-MinCorner)/2.0;
//};