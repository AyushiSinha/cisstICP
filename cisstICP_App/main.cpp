// ****************************************************************************
//
//    Copyright (c) 2014, Seth Billings, Ayushi Sinha, Russell Taylor, 
//	  Johns Hopkins University.
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
				Out("out"), Xfm("xfm"), SSM("ssm");
cmdLineInt 		nModes("modes"), nSamples("samples");
cmdLineFloat	Scale("scale");
cmdLineReadable h("h"), help("help");

cmdLineReadable* params[] =
{
	&Alg, &Target, &In,		// strings
	&Out, &Xfm, &SSM,
	&nModes, &nSamples,		// ints
	&Scale,					// floats
	&h, &help,
	NULL
};

void SetParams() 
{
	int i = 0;
	params[i]->description = strdup("Select algorithm from the following:\n"
									"\t\t\tStdICP: Implements the standard ICP algorithm\n"
									"\t\t\tIMLP: Implements the iterative most likely point algorithm (default)\n"
									"\t\t\tDIMLP: Implements the deformable IMLP algorithm\n"
									"\t\t\tDirICP: Implements ICP with orientation algorithm\n"
									"\t\t\tIMLOP: Implements the iterative most likely oriented point algorithm\n"
									"\t\t\tDIMLOP: Implements the deformable IMLOP algorithm\n"
//									"\t\t\tGIMLOP: Implements the generalized IMLOP algorithm\n"
									"\t\t\tPIMLOP: Implements the projected IMLOP algorithm\n"
									"\t\t\tVIMLOP: Implements the video IMLOP algorithm\n\n");
	i++;
	params[i]->description = strdup("Enter the location of the target mesh or point cloud\n\n");
	i++;
	params[i]->description = strdup("Enter the location of the moving mesh or point cloud to be registered to the target\n"
									"\t\t(default = sample points from the target mesh\n"
									"\t\tor deformed target for deformable algorithms)\n\n");
	i++;
	params[i]->description = strdup("Enter the location where the output will be saved\n\n");
	i++;
	params[i]->description = strdup("Enter the location of the initial transformation to be applied to the moving data\n\n");
	i++;
	params[i]->description = strdup("Enter the location of the statistical shape model (SSM) that is homologous to the target\n"
									"\t\tNOTE:\tAn SSM is required for the deformable algorithms. If an SSM is not provided,\n"
									"\t\t\tthe corresponding rigid registration will be computed\n\n");
	i++;
	params[i]->description = strdup("Enter the number of modes you want to use (default = 3)\n\n");
	i++;
	params[i]->description = strdup("Enter the sample size to use (default = 300), or\n"
									"\t\tthe # vertices to subsample from input\n"
									"\t\t(default = input mesh or point cloud size)\n\n");
	i++;
	params[i]->description = strdup("Enter the initial scaling factor for the input (default = 1)\n\n");
	i++;
	params[i]->description = strdup("Prints short usage directions\n\n");
	i++;
	params[i]->description = strdup("Prints detailed usage directions\n\n");
	i++;
}

void Usage(const char* exec)
{
	printf("Usage %s: \n", exec);
	printf("\t--%s <algorithm type>\n", Alg.name);
	printf("\t--%s <target>\n", Target.name);
	printf("\t--%s <input>\n", In.name);
	printf("\t--%s <output directory>\n", Out.name);
	printf("\t--%s <initial guess transform>\n", Xfm.name);
	printf("\t--%s <ssm>\n", SSM.name);
	printf("\t--%s <number of modes>\n", nModes.name);
	printf("\t--%s <number of subsamples>\n", nSamples.name);
	printf("\t--%s <scale>\n", Scale.name);
	printf("\t--%s (Prints usage directions)\n", h.name);
	printf("\t--%s (Prints detailed usage directions)\n", help.name);
}

void Help(const char* exec)
{
	printf("Usage %s: \n", exec);
	printf("\t--%s <algorithm type>\n\t\t%s", Alg.name, Alg.description);
	printf("\t--%s <target>\n\t\t%s", Target.name, Target.description);
	printf("\t--%s <input>\n\t\t%s", In.name, In.description);
	printf("\t--%s <output directory>\n\t\t%s", Out.name, Out.description);
	printf("\t--%s <initial guess transform>\n\t\t%s", Xfm.name, Xfm.description);
	printf("\t--%s <ssm>\n\t\t%s", SSM.name, SSM.description);
	printf("\t--%s <number of modes>\n\t\t%s", nModes.name, nModes.description);
	printf("\t--%s <number of subsamples>\n\t\t%s", nSamples.name, nSamples.description);
	printf("\t--%s <scale>\n\t\t%s", Scale.name, Scale.description);
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
		//algType = AlgType_StdICP;	Alg.value = "StdICP";
		algType = AlgType_IMLP;		//Alg.value = "IMLP";
		//algType = AlgType_DIMLP;	Alg.value = "DIMLP";
		//algType = AlgType_VIMLOP; Alg.value = "VIMLOP";

		testICP(TargetShapeAsMesh, algType, cmdLineOpts);

		// If you switch to one of the oriented point algorithms,
		// don't forget to uncomment testICPNormals, and comment testICP

		//dirAlgType = DirAlgType_StdICP;	Alg.value = "DirICP";
		//dirAlgType = DirAlgType_IMLOP;	Alg.value = "IMLOP";
		////dirAlgType = DirAlgType_GIMLOP;	Alg.value = "GIMLOP";
		//dirAlgType = DirAlgType_PIMLOP;	Alg.value = "PIMLOP";
		//testICPNormals(TargetShapeAsMesh, dirAlgType);
		
		return 0;
	}
	else
	{
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
		//else if (strcmp(Alg.value, "GIMLOP"))		// We did not build G-IMLOP in this instance of the code build because
		//	dirAlgType = DirAlgType_GIMLOP;			// we did not have access to Numerical Recipes, which is a dependency
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

	
	if (nModes.set) {
		cmdLineOpts.modes = nModes.value; 
		cmdLineOpts.useDefaultNumModes = false;
	}

	if (nSamples.set) {
		cmdLineOpts.samples = nSamples.value;
		cmdLineOpts.useDefaultNumSamples = false;
	}

	if (Scale.set)
	{
		cmdLineOpts.scale = Scale.value;
		cmdLineOpts.useDefaultScale = false;
	}

	if (!strcmp(Alg.value, "StdICP") || !strcmp(Alg.value, "IMLP") 
		|| !strcmp(Alg.value, "DIMLP") || !strcmp(Alg.value, "VIMLOP"))
		testICP(TargetShapeAsMesh, algType, cmdLineOpts);
	else if (!strcmp(Alg.value, "DirICP") || !strcmp(Alg.value, "IMLOP")
		|| !strcmp(Alg.value, "DIMLOP") /*|| !strcmp(Alg.value, "GIMLOP")*/ 
		|| !strcmp(Alg.value, "PIMLOP"))
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