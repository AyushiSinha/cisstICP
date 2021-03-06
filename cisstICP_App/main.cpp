// ****************************************************************************
//
//    Copyright (c) 2014, Ayushi Sinha, Seth Billings, Russell Taylor, Johns Hopkins University.
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
#include <stdio.h>
#include <iostream>
#include <vector>

//#include <cisstCommon.h>
#include <cisstVector.h>
#include "utility.h"

// Registration Tests
#include "testICP.h"
#include "testICPNormals.h"

// Command Line Options
#include "CmdLineParser.h"

cmdLineString	Alg("alg"), Target("target"), In("in"),
				Out("out"), Xfm("xfm"), SSM("ssm"),
				Cov("cov"), Axes("axes"),
				ModeWeights("modewts"), 
				WorkingDir("workdir");
cmdLineInt 		targetType("targettype"), 
				nModes("modes"), nSamples("samples"),
				nThresh("nthresh"), 					// PD-tree variables
				nIters("iters");
cmdLineFloat	Scale("scale"),							// transformation offsets
				MinPos("minpos"), MaxPos("maxpos"),
				MinAng("minang"), MaxAng("maxang"),
				pOutliers("poutliers"),					// outlier variables
				outMinPos("outminpos"), outMaxPos("outmaxpos"),
				outMinAng("outminang"), outMaxAng("outmaxang"),
				NoiseInPlane("noiseinplane"),			// noise variables
				NoisePerpPlane("noiseperpplane"),
				NoiseDeg("noisedeg"), NoiseEcc("noiseecc"),
				diagThresh("diagthresh"),				// PD-tree variables
				RotationBounds("rbounds"),				// optimization bounds
				TranslationBounds("tbounds"),
				ScaleBounds("sbounds"),
				ShapeParamBounds("spbounds");
cmdLineReadable bScale("bscale"),
				h("h"), help("help");

cmdLineReadable* params[] =
{
	&Alg, &Target, &In,		// strings
	&Out, &Xfm, &SSM,
	&Cov, &Axes,			/* Positional and angular noise */
	&ModeWeights,
	&WorkingDir,
	&targetType,			// ints
	&nModes, &nSamples,
	&nIters,
	&Scale,					// floats
	&MinPos, &MaxPos,
	&MinAng, &MaxAng,
	&pOutliers,
	&outMinPos, &outMaxPos,
	&outMinAng, &outMaxAng,
	&NoiseInPlane, 
	&NoisePerpPlane,
	&NoiseDeg, &NoiseEcc,
	&nThresh, &diagThresh,
	&RotationBounds,
	&TranslationBounds,
	&ScaleBounds,
	&ShapeParamBounds,
	&bScale,				// readable 
	&h, &help,				// help
	NULL
};

void SetParams() 
{
	int i = 0;
	// Algorithm
	params[i]->description = strdup("Select algorithm from the following:\n"
									"\t\t\tStdICP: Implements the standard ICP algorithm\n"
									"\t\t\tIMLP: Implements the iterative most likely point algorithm (default)\n"
									"\t\t\tDIMLP: Implements the deformable IMLP algorithm\n"
									"\t\t\tDirICP: Implements ICP with orientation algorithm\n"
									"\t\t\tIMLOP: Implements the iterative most likely oriented point algorithm\n"
									"\t\t\tDIMLOP: Implements the deformable IMLOP algorithm\n"
									"\t\t\tGIMLOP: Implements the generalized IMLOP algorithm\n"
									"\t\t\tGDIMLOP: Implements the deformable generalized IMLOP algorithm\n"
									"\t\t\tPIMLOP: Implements the projected IMLOP algorithm\n\n"
									/*"\t\t\tVIMLOP: Implements the video IMLOP algorithm\n\n"*/);
	i++;
	// Target location
	params[i]->description = strdup("Enter the location of the target mesh or point cloud\n\n");
	i++;
	// Input location
	params[i]->description = strdup("Enter the location of the moving mesh or point cloud to be registered to the target\n"
									"\t\t(default = sample points from the target mesh\n"
									"\t\tor deformed target for deformable algorithms)\n\n");
	i++;
	// Directory extention for output
	params[i]->description = strdup("Enter the directory extention that will be appended to the working directory. The output will be saved here\n\n");
	i++;
	// Guess transform
	params[i]->description = strdup("Enter the location of the initial transformation to be applied to the moving data\n\n");
	i++;
	// SSM location
	params[i]->description = strdup("Enter the location of the statistical shape model (SSM) that is homologous to the target\n"
									"\t\tNOTE:\tAn SSM is required for the deformable algorithms. If an SSM is not provided,\n"
									"\t\t\tthe corresponding rigid registration will be computed\n\n");
	i++;
	// File containing covariance matrices
	params[i]->description = strdup("Enter the location of the covariance matrices for positional noise\n\n");
	i++;
	// File containing major/minor axes
	params[i]->description = strdup("Enter the location of the major/minor axes for angular noise\n\n");
	i++;
	// File containing mode weights
	params[i]->description = strdup("Enter the location of the mode weights\n\n");
	i++;
	// Working directory
	params[i]->description = strdup("Enter the new working directory (default = \"..\\..\\..\\test_data\\LastRun_<algorithm name>\"\n\n");
	i++;
	// Target type (mesh or point cloud)
	params[i]->description = strdup("Specify the target type:\n"
									"\t\tMesh: 1 (default)\n"
									"\t\tPoint cloud: 0\n\n");
	i++;
	// Number of modes
	params[i]->description = strdup("Enter the number of modes you want to use (default = 3)\n\n");
	i++;
	// Number of samples
	params[i]->description = strdup("When no input point cloud is provided, number of points to sample from model shape (default = 300), or\n"
									"\t\tWhen an input point cloud is provided, number of points to subsample from it (default = point cloud size)\n\n");
	i++;
	// Maximum number of iterations
	params[i]->description = strdup("Enter the maximum number of iterations to be performed (default = 100)\n\n");
	i++;
	// Scaling factor for input
	params[i]->description = strdup("Enter the initial scaling factor for the input (default = 1.0)\n\n");
	i++;
	// Minimum positional offset
	params[i]->description = strdup("Enter the lower bound for positional offset (default = 10)\n\n");
	i++;
	// Maximum positional offset
	params[i]->description = strdup("Enter the upper bound for positional offset (default = 20)\n\n");
	i++;
	// Minimum orientation offset
	params[i]->description = strdup("Enter the lower bound for orientation offset (default = 6)\n\n");
	i++;
	// Maximum orientation offset
	params[i]->description = strdup("Enter the upper bound for orientation offset (default = 12)\n\n");
	i++;
	// Percent outliers
	params[i]->description = strdup("Enter the percentage of outliers to add to sampled points (default = 0%)\n\n");
	i++;
	// Minimum positional offset for outliers
	params[i]->description = strdup("Enter the lower bound for outlier positional offset (default = 5)\n\n");
	i++;
	// Maximum positional offset for outliers
	params[i]->description = strdup("Enter the upper bound for outlier positional offset (default = 10)\n\n");
	i++;
	// Minimum orientation offset for outliers
	params[i]->description = strdup("Enter the lower bound for outlier orientation offset (default = 5)\n\n");
	i++;
	// Maximum orientation offset for outliers
	params[i]->description = strdup("Enter the upper bound for outlier orientation offset (default = 10)\n\n");
	i++;
	// In-plane positional noise SD
	params[i]->description = strdup("Standard deviation of positional noise in plane (default = 1)\n\n");
	i++;
	// Out-of-plane positional noise SD
	params[i]->description = strdup("Standard deviation of positional noise out of plane (default = 1)\n\n");
	i++;
	// Angular noise SD
	params[i]->description = strdup("Standard deviation of angular noise in degrees (default = 2)\n\n");
	i++;
	// Angular noise eccentricity
	params[i]->description = strdup("Eccentricity of angular noise (default = 0.5)\n\n");
	i++;
	// Min number of datums in each PD-tree node
	params[i]->description = strdup("Min number of datums in each PD-tree node (default = 5 (without orientaion) or 15 (with orientation))\n\n");
	i++;
	// Min size of PD-tree node
	params[i]->description = strdup("Min size of PD-tree node (default = 5 (without orientaion) or 15 (with orientation))\n\n");
	i++;
	// Rotation constraint
	params[i]->description = strdup("Constrain rotation component search between [prev-n, prev+n] (default n = DBL_MAX)\n\n");
	i++;
	// Translation constraint
	params[i]->description = strdup("Constrain translation component search between [prev-n, prev+n] (default n = DBL_MAX)\n\n");
	i++;
	// Scale constraint
	params[i]->description = strdup("Constrain scale component search between [1-n, 1+n] (default n = 0.3)\n\n");
	i++;
	// Shape parameter constraint
	params[i]->description = strdup("Constrain shape parameter search between [-n, n] (default n = 3.0)\n\n");
	i++;
	// Scale optimization
	params[i]->description = strdup("Optimize over scale in addition to [R,t] and shape parameters (default = false)\n"
									"\t\tOnly available for D-IMLP, D-IMLOP, G-IMLOP, and GD-IMLOP algorithms\n\n");
	i++;
	// Brief usage directions
	params[i]->description = strdup("Prints short usage directions\n\n");
	i++;
	// Detailed usage directions
	params[i]->description = strdup("Prints detailed usage directions\n\n");
}

void Usage(const char* exec)
{
	printf("Usage %s: \n", exec);
	printf("\t--%s <algorithm type>\n", Alg.name);
	printf("\t--%s <target>\n", Target.name);
	printf("\t--%s <target type>\n", targetType.name); 
	printf("\t--%s <input>\n", In.name);
	printf("\t--%s <output extention>\n", Out.name);
	printf("\t--%s <initial guess transform>\n", Xfm.name);
	printf("\t--%s <ssm>\n", SSM.name);
	printf("\t--%s <cov (positional noise)>\n", Cov.name);
	printf("\t--%s <axes (angular noise)>\n", Axes.name);
	printf("\t--%s <mode weights>\n", ModeWeights.name);
	printf("\t--%s <working directory>\n", WorkingDir.name);
	printf("\t--%s <number of modes>\n", nModes.name);
	printf("\t--%s <number of samples>\n", nSamples.name);
	printf("\t--%s <max iterations>\n", nIters.name);
	printf("\t--%s <scale>\n", Scale.name);
	printf("\t--%s \n", bScale.name);
	printf("\t--%s <min pos offset>\n", MinPos.name);
	printf("\t--%s <max pos offset>\n", MaxPos.name);
	printf("\t--%s <min ang offset>\n", MinAng.name);
	printf("\t--%s <min ang offset>\n", MaxAng.name);
	printf("\t--%s <percent outliers>\n", pOutliers.name);
	printf("\t--%s <outlier min pos offset>\n", outMinPos.name);
	printf("\t--%s <outlier max pos offset>\n", outMaxPos.name);
	printf("\t--%s <outlier min ang offset>\n", outMinAng.name);
	printf("\t--%s <outlier min ang offset>\n", outMaxAng.name);
	printf("\t--%s <in plane noise>\n", NoiseInPlane.name);
	printf("\t--%s <out of plane noise>\n", NoisePerpPlane.name);
	printf("\t--%s <angular noise (deg)>\n", NoiseDeg.name);
	printf("\t--%s <angular noise eccentricity>\n", NoiseEcc.name);
	printf("\t--%s <min datums in PD-tree node>\n", nThresh.name);
	printf("\t--%s <min size of PD-tree node>\n", diagThresh.name);
	printf("\t--%s <rotation constraints>\n", RotationBounds.name);
	printf("\t--%s <translation constraints>\n", TranslationBounds.name);
	printf("\t--%s <scale constraints>\n", ScaleBounds.name);
	printf("\t--%s <shape parameter constraints>\n", ShapeParamBounds.name);
	printf("\t--%s (Prints usage directions)\n", h.name);
	printf("\t--%s (Prints detailed usage directions)\n", help.name);
}

void Help(const char* exec)
{
	printf("Usage %s: \n", exec);
	printf("\t--%s <algorithm type>\n\t\t%s", Alg.name, Alg.description);
	printf("\t--%s <target>\n\t\t%s", Target.name, Target.description);
	printf("\t--%s <target type>\n\t\t%s", targetType.name, targetType.description);
	printf("\t--%s <input>\n\t\t%s", In.name, In.description);
	printf("\t--%s <output directory>\n\t\t%s", Out.name, Out.description);
	printf("\t--%s <initial guess transform>\n\t\t%s", Xfm.name, Xfm.description);
	printf("\t--%s <ssm>\n\t\t%s", SSM.name, SSM.description);
	printf("\t--%s <cov (positional noise)>\n\t\t%s", Cov.name, Cov.description);
	printf("\t--%s <axes (angular noise)>\n\t\t%s", Axes.name, Axes.description);
	printf("\t--%s <mode weights>\n\t\t%s", ModeWeights.name, ModeWeights.description);
	printf("\t--%s <working directory>\n\t\t%s", WorkingDir.name, WorkingDir.description);
	printf("\t--%s <number of modes>\n\t\t%s", nModes.name, nModes.description);
	printf("\t--%s <number of samples>\n\t\t%s", nSamples.name, nSamples.description);
	printf("\t--%s <max iterations>\n\t\t%s", nIters.name, nIters.description);
	printf("\t--%s <scale>\n\t\t%s", Scale.name, Scale.description);
	printf("\t--%s \n\t\t%s", bScale.name, bScale.description);
	printf("\t--%s <min pos offset>\n\t\t%s", MinPos.name, MinPos.description);
	printf("\t--%s <max pos offset>\n\t\t%s", MaxPos.name, MaxPos.description);
	printf("\t--%s <min ang offset>\n\t\t%s", MinAng.name, MinAng.description);
	printf("\t--%s <min ang offset>\n\t\t%s", MaxAng.name, MaxAng.description);
	printf("\t--%s <percent outliers>\n\t\t%s", pOutliers.name, pOutliers.description);
	printf("\t--%s <outlier min pos offset>\n\t\t%s", outMinPos.name, outMinPos.description);
	printf("\t--%s <outlier max pos offset>\n\t\t%s", outMaxPos.name, outMaxPos.description);
	printf("\t--%s <outlier min ang offset>\n\t\t%s", outMinAng.name, outMinAng.description);
	printf("\t--%s <outlier min ang offset>\n\t\t%s", outMaxAng.name, outMaxAng.description);
	printf("\t--%s <in plane noise>\n\t\t%s", NoiseInPlane.name, NoiseInPlane.description);
	printf("\t--%s <out of plane noise>\n\t\t%s", NoisePerpPlane.name, NoisePerpPlane.description);
	printf("\t--%s <angular noise (deg)>\n\t\t%s", NoiseDeg.name, NoiseDeg.description);
	printf("\t--%s <angular noise eccentricity>\n\t\t%s", NoiseEcc.name, NoiseEcc.description);
	printf("\t--%s <min datums in PD-tree node>\n\t\t%s", nThresh.name, nThresh.description);
	printf("\t--%s <min size of PD-tree node>\n\t\t%s", diagThresh.name, diagThresh.description);
	printf("\t--%s <rotation constraints>\n\t\t%s", RotationBounds.name, RotationBounds.description);
	printf("\t--%s <translation constraints>\n\t\t%s", TranslationBounds.name, TranslationBounds.description);
	printf("\t--%s <scale constraints>\n\t\t%s", ScaleBounds.name, ScaleBounds.description);
	printf("\t--%s <shape parameter constraints>\n\t\t%s", ShapeParamBounds.name, ShapeParamBounds.description);
	printf("\t--%s \t%s", h.name, h.description);
	printf("\t--%s \t%s", help.name, help.description);
}

int main( int argc, char* argv[] )
{
	// initialize variables 
	ICPAlgType algType;
	ICPDirAlgType dirAlgType;
	cisstICP::CmdLineOptions cmdLineOpts;

	// set defaults
	bool TargetShapeAsMesh = true;

	// read in command line options
	std::vector< std::string > nonoptArgs;
	cmdLineParse(argc, argv, params, nonoptArgs);

	if (h.set)
	{
		Usage(argv[0]);
		return 0;
	}

	if (help.set)
	{
		SetParams();
		Help(argv[0]);
		return 0;
	}

	if (!Alg.set)
	{
		 // By default, not providing any input will run the IMLP algorithm with default settings
		 // But, you may uncomment below to manually pick an algorithm to run with default settings
		 // It is recommended NOT to do this - you can set the algorithm to run via command line

		//-- Registration Test Runs --//
		algType = AlgType_IMLP;		//Alg.value = "IMLP";
		testICP(TargetShapeAsMesh, algType, cmdLineOpts);
		
		return 0;
	}
	else
	{
		//std::cout << Alg.value << "\n";
		if (!strcmp(Alg.value, "StdICP"))
			algType = AlgType_StdICP;
		else if (!strcmp(Alg.value, "IMLP"))
			algType = AlgType_IMLP;
		else if (!strcmp(Alg.value, "DIMLP"))
		{
			algType = AlgType_DIMLP;
			cmdLineOpts.deformable = true;
		}
		else if (!strcmp(Alg.value, "VIMLOP"))
			algType = AlgType_VIMLOP;
		else if (!strcmp(Alg.value, "DirICP"))
			dirAlgType = DirAlgType_StdICP;
		else if (!strcmp(Alg.value, "IMLOP"))
			dirAlgType = DirAlgType_IMLOP;
		else if (!strcmp(Alg.value, "DIMLOP"))
		{
			dirAlgType = DirAlgType_DIMLOP;
			cmdLineOpts.deformable = true;
		}
		else if (!strcmp(Alg.value, "GIMLOP")) {		// Comment these two lines to remove G-IMLOP as an option in case 
			dirAlgType = DirAlgType_GIMLOP;			// you do not have access to Numerical Recipes, which is a dependency
		}
		else if (!strcmp(Alg.value, "GDIMLOP"))
		{
			dirAlgType = DirAlgType_GDIMLOP;
			cmdLineOpts.deformable = true;
		}
		else if (!strcmp(Alg.value, "PIMLOP"))
			dirAlgType = DirAlgType_PIMLOP;
	}

	if (Target.set) {
		cmdLineOpts.target = Target.value;
		cmdLineOpts.useDefaultTarget = false;
	}

	if (In.set) {
		cmdLineOpts.input = In.value;
		cmdLineOpts.useDefaultInput = false;
	}

	if (ModeWeights.set) {
		cmdLineOpts.modeweights = ModeWeights.value;
		cmdLineOpts.readModeWeights = true;
		cmdLineOpts.useDefaultNumModes = false;
	}

	if (Out.set) {
		cmdLineOpts.output = Out.value;
		cmdLineOpts.useDefaultOutput = false;
	}

	if (Xfm.set) {
		cmdLineOpts.xfm = Xfm.value;
		cmdLineOpts.useDefaultXfm = false;
	}

	if (SSM.set) {
		cmdLineOpts.ssm = SSM.value;
		cmdLineOpts.useDefaultSSM = false;
	}

	if (Cov.set) {
		cmdLineOpts.cov = Cov.value;
		cmdLineOpts.useDefaultCov = false;
	}

	if (Axes.set) {
		cmdLineOpts.axes = Axes.value;
		cmdLineOpts.useDefaultAxes = false;
	}

	if (WorkingDir.set) {
		cmdLineOpts.workingdir = WorkingDir.value;
		cmdLineOpts.useDefaultWorkingDir = false;
	}

	if (targetType.set) {
		if (targetType.value == 0)
			TargetShapeAsMesh = false;
		else if (targetType.value == 1)
			TargetShapeAsMesh = true;
		else {
			std::cerr << "Invalid option for target type: " << targetType.value << "\nExiting...\n\n";
			return 0;
		}
	}
	
	if (nModes.set) {
		cmdLineOpts.modes = nModes.value; 
		cmdLineOpts.useDefaultNumModes = false;
	}

	if (nSamples.set) {
		cmdLineOpts.samples = nSamples.value;
		cmdLineOpts.useDefaultNumSamples = false;
	}

	if (nIters.set) {
		cmdLineOpts.niters = nIters.value;
		cmdLineOpts.useDefaultNumIters = false;
	}

	if (Scale.set) {
		cmdLineOpts.scale = Scale.value;
		cmdLineOpts.useDefaultScale = false;
	}

	if (bScale.set)
	{
		cmdLineOpts.bScale = true;
	}

	if (MinPos.set) {
		cmdLineOpts.minpos = MinPos.value;
		cmdLineOpts.useDefaultMinPos = false;
	}

	if (MaxPos.set) {
		cmdLineOpts.maxpos = MaxPos.value;
		cmdLineOpts.useDefaultMaxPos = false;
	}

	if (MinAng.set) {
		cmdLineOpts.minang = MinAng.value;
		cmdLineOpts.useDefaultMinAng = false;
	}

	if (MaxAng.set) {
		cmdLineOpts.maxang = MaxAng.value;
		cmdLineOpts.useDefaultMaxAng = false;
	}

	if (pOutliers.set) {
		cmdLineOpts.poutliers = pOutliers.value;
		cmdLineOpts.useDefaultNumOutliers = false;
	}

	if (outMinPos.set) {
		cmdLineOpts.outminpos = outMinPos.value;
		cmdLineOpts.useDefaultOutMinPos = false;
	}

	if (outMaxPos.set) {
		cmdLineOpts.outmaxpos = outMaxPos.value;
		cmdLineOpts.useDefaultOutMaxPos = false;
	}

	if (outMinAng.set) {
		cmdLineOpts.outminang = outMinAng.value;
		cmdLineOpts.useDefaultOutMinAng = false;
	}

	if (outMaxAng.set) {
		cmdLineOpts.outmaxang = outMaxAng.value;
		cmdLineOpts.useDefaultOutMaxAng = false;
	}

	if (NoiseInPlane.set) {
		cmdLineOpts.noiseinplane = NoiseInPlane.value;
		cmdLineOpts.useDefaultNoiseInPlane = false;
	}

	if (NoisePerpPlane.set) {
		cmdLineOpts.noiseperpplane = NoisePerpPlane.value;
		cmdLineOpts.useDefaultNoisePerpPlane = false;
	}

	if (NoiseDeg.set) {
		cmdLineOpts.noisedeg = NoiseDeg.value;
		cmdLineOpts.useDefaultNoiseDeg = false;
	}

	if (NoiseEcc.set) {
		cmdLineOpts.noiseecc = NoiseEcc.value;
		cmdLineOpts.useDefaultNoiseEcc = false;
	}

	if (nThresh.set) {
		cmdLineOpts.nthresh = nThresh.value;
		cmdLineOpts.useDefaultNThresh = false;
	}

	if (diagThresh.set) {
		cmdLineOpts.diagthresh = diagThresh.value;
		cmdLineOpts.useDefaultDiagThresh = false;
	}

	if (RotationBounds.set) {
		cmdLineOpts.rbounds = RotationBounds.value;
		cmdLineOpts.useDefaultRotationBounds = false;
	}

	if (TranslationBounds.set) {
		cmdLineOpts.tbounds = TranslationBounds.value;
		cmdLineOpts.useDefaultTranslationBounds = false;
	}

	if (ScaleBounds.set) {
		cmdLineOpts.sbounds = ScaleBounds.value;
		cmdLineOpts.useDefaultScaleBounds = false;
	}

	if (ShapeParamBounds.set) {
		cmdLineOpts.spbounds = ShapeParamBounds.value;
		cmdLineOpts.useDefaultShapeParamBounds = false;
	}

	if (!strcmp(Alg.value, "StdICP") || !strcmp(Alg.value, "IMLP") 
		|| !strcmp(Alg.value, "DIMLP") || !strcmp(Alg.value, "VIMLOP"))
		testICP(TargetShapeAsMesh, algType, cmdLineOpts);
	else if (!strcmp(Alg.value, "DirICP") || !strcmp(Alg.value, "IMLOP")
		|| !strcmp(Alg.value, "DIMLOP")  || !strcmp(Alg.value, "GIMLOP") 
		|| !strcmp(Alg.value, "GDIMLOP") || !strcmp(Alg.value, "PIMLOP"))
		testICPNormals(TargetShapeAsMesh, dirAlgType, cmdLineOpts);

	return 0;
}

#if 0
int main(void)
{

	//-- Registration Test Runs --//

	bool TargetShapeAsMesh = true;
	ICPAlgType algType;
	//algType = AlgType_StdICP;
	algType = AlgType_IMLP;
	//algType = AlgType_DIMLP;
	//algType = AlgType_VIMLOP;
	testICP(TargetShapeAsMesh, algType);

	TargetShapeAsMesh = true;
	ICPDirAlgType dirAlgType;
	//dirAlgType = DirAlgType_StdICP;
	//dirAlgType = DirAlgType_IMLOP;
	//dirAlgType = DirAlgType_GIMLOP;
	//dirAlgType = DirAlgType_PIMLOP;
	//testICPNormals(TargetShapeAsMesh, dirAlgType);

	return 0;
}
#endif