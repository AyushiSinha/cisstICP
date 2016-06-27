// ****************************************************************************
//
//    Copyright (c) 2016, Ayushi Sinha, Seth Billings, Russell Taylor, Johns Hopkins University
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

#include "camera.h"

// constructor
camera::camera(int imWidth, int imHeight, int fx, int fy, int Ox, int Oy)
{
	SetSamples(imWidth, imHeight, fx, fy, Ox, Oy);
	SetPose_Xfm(vctFrm3::Identity());
}


void camera::SetSamples(int argImWidth, int argImHeight, int argFx, int argFy, int argOx, int argOy)
{
	imWidth = argImWidth;
	imHeight = argImHeight;
	fx = argFx;
	fy = argFy;
	Ox = argOx;
	Oy = argOy;
}

void camera::SetPose_Xfm(vctFrm3 argF)
{
	F = argF;
	camGetPose_Xfm2Cam(argF);
}

void camera::camGetPose_Xfm2Cam(vctFrm3 argF)
{
	vct3 tmpCup, tmpCtar;
	tmpCup.Zeros();
	tmpCup(1) = -1; // tmpCup = [0, -1,  0]
	tmpCtar.Zeros();
	tmpCtar(2) = 1; //tmpCtar = [0,  0,  1]

	vctRot3 Rinv(argF.Rotation().Transpose());
	vct3 tinv = -Rinv*argF.Translation();
	cpos = tinv;
	cup = tmpCup*Rinv.Transpose();
	ctar = cpos + tmpCtar*Rinv.Transpose();
}

void camera::SetPose_Cam(vct3 argCpos, vct3 argCtar, vct3 argCup)
{
	cpos = argCpos;
	ctar = argCtar;
	cup = argCup;
	camGetPose_Cam2Xfm(argCpos, argCtar, argCup);
}
void camera::camGetPose_Cam2Xfm(vct3 argCpos, vct3 argCtar, vct3 argCup)
{
	vct3 zdir = argCtar - argCpos;
	vct3 xdir;
	xdir.CrossProductOf(zdir, argCup);
	vct3 ydir;
	ydir.CrossProductOf(zdir, xdir);
	vct3 x = xdir / xdir.Norm();
	vct3 y = ydir / ydir.Norm();
	vct3 z = zdir / zdir.Norm();

	vctRot3 R;
	R.at(0, 0) = x(0); R.at(0, 1) = x(1); R.at(0, 2) = x(2);
	R.at(1, 0) = y(0); R.at(1, 1) = y(1); R.at(1, 2) = y(2);
	R.at(2, 0) = z(0); R.at(2, 1) = z(1); R.at(2, 2) = z(2);
	vct3 t = -R*argCpos;
	F.Assign(R, t);
}

void camera::GetCoord_Pix2imPix(vctDynamicVector<vct2> pPixel)
{
	vct2 offset;
	offset(0) = Ox + 1; offset(1) = Oy + 1;
	for (unsigned int s = 0; s < pPixel.size(); s++)
		pImPix(s) = pPixel(s) - offset;
}
