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

#ifndef _camera_h
#define _camera_h

#include <cisstVector.h>

class camera
{
  //
  // This is the base class for storing data pertaining to a calibrated camera
  //

  //--- Standard Algorithm Parameters ---//

public:

  // Intrinsic Parameters (K) Matrix
	//| fx  0  Ox |   | f*Sx  0   Ox |
	//| 0  fy  Oy | = |   0  f*Sy Oy |
	//| 0   0   1 |   |   0   0    1 |

	//Note: For calibrated parameters -
	//	  f(mm), Sx(pix / mm), Sy(pix/mm)
	//then -
	//	  fx = f*Sx, fy = f*Sy

	//Note: Focal length matrix
	//	  FL = | fx   0 |
	//		   |  0  fy |

	int imWidth;		// width of image in pixels
	int imHeight;		// height of image in pixels
	int fx;				// focal length in X in pixels
	int fy;				// focal length in Y in pixels
	int Ox;				// offset for optical center in X in pixels
	int Oy;				// offset for optical center in Y in pixels

	vctDynamicVector<vct2> pImPix;

private:

	vctFrm3 F;
	vct3 cpos;
	vct3 cup;
	vct3 ctar;

  //--- Standard Algorithm Methods ---//

public:

  // constructor
  camera(int imWidth, int imHeight, int fx, int fy, int Ox, int Oy);

  virtual void  SetSamples(int imWidth, int imHeight, int fx, int fy, int Ox, int Oy);

  virtual void SetPose_Xfm(vctFrm3 F);
  virtual void camGetPose_Xfm2Cam(vctFrm3 F);
  virtual void SetPose_Cam(vct3 cpos, vct3 ctar, vct3 cup);
  virtual void camGetPose_Cam2Xfm(vct3 cpos, vct3 ctar, vct3 cup);
  virtual void GetCoord_Pix2imPix(vctDynamicVector<vct2> pPixel);

protected:


  //--- ICP Interface Methods ---//

public:

  
};

#endif
