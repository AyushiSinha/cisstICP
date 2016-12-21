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
 // Compute  camera pose parameters from a camera's extrinsic 
 // parameter matrix.

 //  Input
 //   F = [R, t]   ~rotation and translation(world->cam)

 //  Output
 //    cpos   ~camera position
 //    ctar   ~some point along the camera's optical axis
 //    cup    ~camera up direction within image

 // define camera frame : x ~right
 //                      y ~down
 //                      z ~outwards along optical axis

 // cpos = inv(Fcam)*[0 0 0]'
 // cup = inv(Rcam)*[0 - 1 0]
 // ctar = cpos + inv(Rcam)*[0 0 1]'

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
 // Compute extrinsic parameter matrix for a camera given camera pose
 // parameters from Matlab.

 //  Input
 //    cpos   ~camera position
 //    ctar   ~camera view direction
 //    cup    ~camera up direction within image

 //  Output
 //   F = [R, t]   ~rotation and translation

 // define camera frame : x ~right
 //                      y ~down
 //                      z ~outwards along optical axis

 // compute orientation of each camera axis in world coords

	vct3 zdir = argCtar - argCpos;
	vct3 xdir;
	xdir.CrossProductOf(zdir, argCup);
	vct3 ydir;
	ydir.CrossProductOf(zdir, xdir);
	vct3 x = xdir / xdir.Norm();
	vct3 y = ydir / ydir.Norm();
	vct3 z = zdir / zdir.Norm();

// compute tranformation to take camera axis to standard camera frame
	vctRot3 R;
	R.at(0, 0) = x(0); R.at(0, 1) = x(1); R.at(0, 2) = x(2);
	R.at(1, 0) = y(0); R.at(1, 1) = y(1); R.at(1, 2) = y(2);
	R.at(2, 0) = z(0); R.at(2, 1) = z(1); R.at(2, 2) = z(2);
	vct3 t = -R*argCpos;
	F.Assign(R, t);
}

void camera::GetPose_Xfm(vctFrm3 &argF)
{
	argF = F;
}

void camera::GetPose_Cam(vct3 &argCpos, vct3 &argCtar, vct3 &argCup)
{
	argCpos = cpos;
	argCtar = ctar;
	argCup = cup;
}

void camera::MoveScene_Xfm(vctFrm3 argdF)
{
	//   move entire scene(camera and scene elements) by incremental[R, t]
	//   the camera's view matrix changes as the inverse of the scene motion
	//   initial view matrix : Tcam = [Rcam, tcam] = >   Xcam = Tcam*Xlocal
	//   scene moves by dT = [dR, dt]	= >   Xworld = dT*Xlocal
	//   new camera view matrix : Xcam = Tcam2*Xworld = Tcam2*dT*Xlocal
	//									= Tcam*Xlocal = >  Tcam2*dT = Tcam
	//												  = >  Tcam2 = Tcam*inv(dT)

	GetXfmScenePose_Xfm(argdF);
	camGetPose_Xfm2Cam(F);
}

void camera::GetXfmScenePose_Xfm(vctFrm3 argdF)
{
	// return camera view matrix for a scene transformed by[R, t]
	// (see documentation for move scene method of this class)

	F = F*argdF.Inverse();
}

void camera::GetScaledExtrinsics(float scale, vctFrm3 &argFs)
{
	// Return the scaled extrinsic parameters of this camera

	vctRot3 Rcam;
	vct3 tcam;
	Rcam = F.Rotation();
	tcam = F.Translation();
	argFs.Assign(Rcam, tcam.Multiply(scale));
}

void camera::GetFocalLengthMatrix(vct2x2 &argFL)
{
	argFL.Zeros();
	argFL.at(0, 0) = fx;
	argFL.at(1, 1) = fy;
}

void camera::GetCoord_Pix2World(vctDynamicVector<vct2> pPixel, vctDynamicVector<vct1> Z, vctDynamicVector<vct3> &pWorld)
{
	 // convert 2d pixel coordinate(s) to 3d "world" coordinate(s)
	 // for an image pixel at depth Z
	 // inputs:
	 //  pPixel	    ~point in pixel coords with (1, 1) referencing the pixel origin
	 //				 located in the center of the pixel at upper - left corner of image (nPts x 2)
	 //	Z			~depth in world units (nPts x 1)
	 // outputs :
	 //  pWorld		~point in world coords (nPts x 3)

	vctDoubleRot3 K;
	K.Zeros();
	K.at(0, 0) = fx; K.at(0, 2) = Ox;
	K.at(1, 1) = fy; K.at(1, 2) = Oy;
	K.at(2, 2) = 1;

	// pCamera = K^-1 * pPixel
	int nPts = Z.size();
	pWorld.resize(nPts);
	pWorld.Zeros();

	// convert from base-1 to base-0 indexing
	vctDynamicVector<vct2> pPixel_base0;
	vct3 pPixel3D;
	vct3 pCamera;
	for (int i = 0; i < nPts; i++)
	{
		pPixel_base0[i][0] = pPixel[i][0] - 1;
		pPixel_base0[i][1] = pPixel[i][1] - 1;
		pPixel3D.at(0) = pPixel_base0[i][0];
		pPixel3D.at(1) = pPixel_base0[i][1];
		pPixel3D.at(2) = Z[i][0];
		pCamera = K.Inverse()*pPixel3D;

		// Apply transform F^-1 to pCamera
		pWorld.at(i) = pCamera*F.Inverse().Rotation() + F.Inverse().Translation();
	}
}

void camera::GetCoord_Pix2imPix(vctDynamicVector<vct2> pPixel, vctDynamicVector<vct2> &pImPix)
{
 // convert 2d pixel coordinate(s) to 2d "image pixel" coordinate(s)

 // NOTE: "image pixel" coords are "image coordinates" with "pixel" units,
 //      i.e, the origin is located at the optical center(the "image coordinate" part)
 // but the units are in pixels(the "Pixel" part) rather than metric units

 // inputs :
 //  pPixel  ~point in pixel coords with(1, 1) referencing the pixel origin
 //            located in the center of the pixel at upper - left corner of image(nPts x 2)
 // outputs :
 //  pImPix  ~point(s) in image pixel coords[nPts x 2]

	// re-index pixel units from base-1 to base-0 indexing
	vct2 offset;
	offset(0) = Ox + 1; offset(1) = Oy + 1;

	for (unsigned int s = 0; s < pPixel.size(); s++)
		pImPix(s) = pPixel(s) - offset;
}

void camera::GetCoord_imPix2Pix(vctDynamicVector<vct2> pImPixel, vctDynamicVector<vct2> &pPix)
{
 // convert 2d "image pixel" coordinate(s) to 2d pixel coordinate(s)

 // NOTE: "image pixel" coords are "image coordinates" with "pixel" units,
 //      i.e, the origin is located at the optical center(the "image coordinate" part)
 // but the units are in pixels(the "Pixel" part) rather than metric units

 // inputs :
 //  pPixel  ~point in pixel coords with(1, 1) referencing the pixel origin
 //            located in the center of the pixel at upper - left corner of image(nPts x 2)
 // outputs :
 //  pImPix  ~point(s) in image pixel coords[nPts x 2]

	vct2 offset;
	offset(0) = Ox + 1; offset(1) = Oy + 1;

	// re-index pixel units from base-0 to base-1 indexing
	for (unsigned int s = 0; s < pImPixel.size(); s++)
		pPix(s) = pImPixel(s) + offset;
}

void camera::camGetProjectionMatrix(vct4x4& xfmProj, float zNear = 0.5, float zFar = 500)
{
 // Compute DirectX projection matrix based on the intrinsic camera
 // parameters

 // cam     ~instance of objCamera defining parameters of a physical camera
 // zNear   ~near cutting plane
 // zFar    ~far cutting plane

 // NOTE: Returns projection matrix assuming a matrix - vector multiplication
 //       convention of placing vector to right of matrix.Note that the
 //       DirectX library assumes the opposite convention.

	float m00 = 2 * fx / imWidth; // matlab m11
	float m11 = 2 * fy / imHeight; // m22
	float m22 = (zFar + zNear) / (zFar - zNear); // m33
	float m23 = -2.0 * zFar * zNear / (zFar - zNear); //m34

	xfmProj.SetAll(0.0);
	xfmProj.at(0, 0) = m00;
	xfmProj.at(1, 1) = m11;
	xfmProj.at(2, 2) = m22;
	xfmProj.at(2, 3) = m23;
	xfmProj.at(3, 2) = 1;
}

void camera::camGetPerspectiveProjection(vctDynamicVector<vct2>& point2d, vctDynamicVector<vct3> point3d)
{
	point2d.SetSize(point3d.size());
	for (int i = 0; i < point3d.size(); i++) {
		point2d[i].Element(0) = point3d[i].Element(0) / point3d[i].Element(2);
		point2d[i].Element(1) = point3d[i].Element(1) / point3d[i].Element(2);
	}
}

void camera::camGetPerspectiveProjectionJacobian(vctDynamicVector<vct2x3>& J_PP, vctDynamicVector<vct3> point3d)
{
	J_PP.SetSize(point3d.size());
	for (int i = 0; i < point3d.size(); i++) {
		J_PP[i].SetAll (0.0);
		J_PP[i].Element(0, 0) = 1 / point3d[i].Element(2);
		J_PP[i].Element(0, 2) = -point3d[i].Element(0) / (point3d[i].Element(2)*point3d[i].Element(2));
		J_PP[i].Element(1, 1) = 1 / point3d[i].Element(2);
		J_PP[i].Element(1, 2) = -point3d[i].Element(1) / (point3d[i].Element(2)*point3d[i].Element(2));
	}
}

void camera::camGetOrthographicProjection(vctDynamicVector<vct2>& point2d, vctDynamicVector<vct3> point3d)
{
	point2d.SetSize(point3d.size());
	for (int i = 0; i < point3d.size(); i++) {
		point2d[i].Element(0) = point3d[i].Element(0) ; 
		point2d[i].Element(1) = point3d[i].Element(1) ; 
	}
}