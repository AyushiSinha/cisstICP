// ****************************************************************************
//
//    Copyright (c) 2014, Ayushi Sinha, Seth Billings, Russell Taylor, Johns Hopkins University. 
//	  All rights reserved.
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
#ifndef _testICP_H
#define _testICP_H

#include <stdio.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <limits.h>

#include "utility.h"
#include "camera.h"
#include "cisstICP.h"
#include "cisstMesh.h"
#include "cisstPointCloud.h"
#include "PDTree_Mesh.h"
#include "PDTree_PointCloud.h"

#include "algICP_StdICP_Mesh.h"
#include "algICP_IMLP_Mesh.h"
#include "algDirICP_VIMLOP.h"
#include "algICP_DIMLP.h"

#include "algICP_StdICP_PointCloud.h"
#include "algICP_IMLP_PointCloud.h"


enum ICPAlgType { AlgType_StdICP, AlgType_IMLP, AlgType_DIMLP, AlgType_VIMLOP };

void Callback_TrackRegPath_testICP(cisstICP::CallbackArg &arg, void *userData)
{
	// Save to file:
	//  - error function
	//  - incremental transform

	// output format:
	//  iter error r00 r01 r02 r10 r11 r12 r20 r21 r22 tx ty tz s sp00 sp01 ... spn
	std::ofstream *fs = (std::ofstream *)(userData);
	(*fs) << arg.iter << ", , " << arg.E << ", , "
		<< arg.Freg.Rotation().Row(0)(0) << ", " << arg.Freg.Rotation().Row(0)(1) << ", " << arg.Freg.Rotation().Row(0)(2) << ", , "
		<< arg.Freg.Rotation().Row(1)(0) << ", " << arg.Freg.Rotation().Row(1)(1) << ", " << arg.Freg.Rotation().Row(1)(2) << ", , "
		<< arg.Freg.Rotation().Row(2)(0) << ", " << arg.Freg.Rotation().Row(2)(1) << ", " << arg.Freg.Rotation().Row(2)(2) << ", , "
		<< arg.Freg.Translation()(0) << ", " << arg.Freg.Translation()(1) << ", " << arg.Freg.Translation()(2) << ", ,"
		<< arg.scale << ", ,";
	for (int m = 0; m < arg.S.size(); m++)
		(*fs) << arg.S(m) << ",";
	(*fs) << std::endl;
}

void Callback_SaveIterationsToFile_testICP(cisstICP::CallbackArg &arg, void *userData)
{
	std::ofstream *fs = (std::ofstream *)(userData);

	vctRodRot3 dR(arg.dF.Rotation());
	std::stringstream ss;
	//ss << cmnPrintf("iter=%u  E=%.3f  tolE=%.4f  (dAng/dPos)= %.2f/%.2f  t=%.3f  NNodes=%u/%u/%u  NOut=%u")
	ss << cmnPrintf("iter=%u  E=%.3f  tolE=%.4f  (dAng/dPos)= %.2f/%.2f  t=%.3f  NOut=%u")
		<< arg.iter
		<< arg.E
		<< arg.tolE
		<< dR.Norm() * 180 / cmnPI << arg.dF.Translation().Norm()
		<< arg.time
		//<< arg.maxNodesSearched << arg.avgNodesSearched << arg.minNodesSearched
		<< arg.nOutliers;

	(*fs) << ss.str() << std::endl;
}

void Callback_SaveIterationsToCSV_testICP(cisstICP::CallbackArg &arg, void *userData)
{
	std::ofstream *fs = (std::ofstream *)(userData);

	vctRodRot3 dR(arg.dF.Rotation());
	std::stringstream ss;
	ss << cmnPrintf("%u ,, %.3f ,, %.4f ,, %.2f,%.2f ,, %.3f ,, %u")
		<< arg.iter
		<< arg.E
		<< arg.tolE
		<< dR.Norm() * 180 / cmnPI << arg.dF.Translation().Norm()
		<< arg.time
		//<< arg.maxNodesSearched << arg.avgNodesSearched << arg.minNodesSearched
		<< arg.nOutliers;

	(*fs) << ss.str() << std::endl;
}

// TargetShapeAsMesh    true - uses mesh to represent target shape
//                      false - uses point cloud (taken from mesh) to represent target shape
void testICP(bool TargetShapeAsMesh, ICPAlgType algType, cisstICP::CmdLineOptions cmdOpts)
{
	// Set output directories
	std::string workingDir, outputDir,
				algDir,		loadMeshPath;

	// Set directories per algorithm
	switch (algType)
	{
	case AlgType_StdICP:
	{		std::cout << "\nRunning standard ICP" << std::endl;		algDir = "LastRun_StdICP/";		break;	}
	case AlgType_IMLP:
	{		std::cout << "\nRunning IMLP" << std::endl;				algDir = "LastRun_IMLP/";		break;	}
	case AlgType_DIMLP:
	{		std::cout << "\nRunning D-IMLP" << std::endl;			algDir = "LastRun_DIMLP/";		break;	}
	case AlgType_VIMLOP:
	{		std::cout << "\nRunning V-IMLOP" << std::endl;			algDir = "LastRun_VIMLOP/";		break;	}
	default:
	{
		std::cout << "ERROR: unknown algorithm type" << std::endl;
		assert(0);
	}
	}

	// Create directories
	workingDir	= cmdOpts.workingdir;
	outputDir	= workingDir + algDir + cmdOpts.output + "/";
	if (CreateDirectory(outputDir.c_str(), NULL))
		std::cout << "Directory created... \n" << std::endl;
	else if (ERROR_ALREADY_EXISTS == GetLastError())
		std::cout << "Directory already exists... \n" << std::endl;
	else
		std::cout << "Directory failed to be created... \n" << std::endl;


	// input files
	std::string normRVFile			 = workingDir + "GaussianValues.txt";	// Random Numbers: Normal RV's

	// output files
	std::string saveMeshPath		 = outputDir + "SaveMesh";
	std::string saveSamplesPath		 = outputDir + "SaveSamples.pts";
	std::string saveSubSamplesPath	 = outputDir + "SaveSubSamples.pts";
	std::string saveNoisySamplesPath = outputDir + "SaveNoisySamples.pts";
	std::string saveCovPath			 = outputDir + "SaveSampleCov.txt";
	std::string saveOffsetXfmPath	 = outputDir + "SaveOffsetXfm.txt";
	std::string saveRegXfmPath		 = outputDir + "SaveRegXfm.txt";
	std::string saveModeWeightsPath	 = outputDir + "saveModeWeights.txt";
	std::string meanMesh			 = outputDir + "meanMesh.ply";
	std::string meshDir				 = outputDir + "estimatedMesh.ply";
	//std::string saveLPath			 = outputDir + "SaveSampleL";

	// Declare variables
	cisstMesh						mesh,	mesh_ssm,	pts;
	PDTreeBase*						pTree;

	vctDynamicVector<unsigned int>  sampleDatums;
	vctDynamicVector<double>		weight;

	vctDynamicVector<vct3>			samples, sampleNorms,
		noisySamples, noisySampleNorms;
	vctDynamicVector<vct3x3>		sampleNoiseCov, sampleNoiseInvCov;
	vctDynamicVector<vct3x2>		sampleNoiseL;


	// Declare and initialize variables (Default values set in cisstICP.h) 
	int		nSamples	= cmdOpts.samples;	// 300	(Default number of samples for default input)
	int		maxIters	= cmdOpts.niters;	// 100
	double	scale		= cmdOpts.scale; 	// 1.0;	
	bool	bScale		= cmdOpts.bScale;	// false;

	int		modes		= 1;				// +1 (for mean)
	int		nThresh		= 5;				// Cov Tree Params
	if (!cmdOpts.useDefaultNThresh)
		nThresh			= cmdOpts.nthresh;
	double	diagThresh	= 5.0;				//  ''
	if (!cmdOpts.useDefaultDiagThresh)
		diagThresh		= cmdOpts.diagthresh;

	double minOffsetPos = (double)cmdOpts.minpos;	// 10.0
	double maxOffsetPos = (double)cmdOpts.maxpos;	// 20.0
	double minOffsetAng = (double)cmdOpts.minang;	//  6.0
	double maxOffsetAng = (double)cmdOpts.maxang;	// 12.0 

	double percentOutliers		= cmdOpts.poutliers / 100.0;	//  0.0
	double minPosOffsetOutlier	= (double)cmdOpts.outminpos;	//  5.0
	double maxPosOffsetOutlier	= (double)cmdOpts.outmaxpos;	// 10.0
	double minAngOffsetOutlier	= (double)cmdOpts.outminang;	//  5.0
	double maxAngOffsetOutlier	= (double)cmdOpts.outmaxang;	// 10.0

	// generate random seeds
	std::srand(time(NULL)); unsigned int randSeed1		= std::rand();	// generates samples
	std::srand(time(NULL)); unsigned int randSeqPos1	= std::rand();
	std::srand(time(NULL)); unsigned int randSeed2		= std::rand();	// generates offsets
	std::srand(time(NULL)); unsigned int randSeqPos2	= std::rand();
	std::srand(time(NULL)); unsigned int randSeed3		= std::rand();	// generates shape parameters
	std::srand(time(NULL)); unsigned int randSeqPos3	= std::rand();
	// use specific seeds if testing for release
	if (cmdOpts.output == "testingforrelease")
	{
		randSeed1	= 0;			// generates samples
		randSeqPos1	= 0;		
		randSeed2	= 17;			// generates offsets
		randSeqPos2	= 28;		
		randSeed3	= 28;			// generates shape parameters
		randSeqPos3	= 8;		
	}

	// Samples Noise Model
	//  NOTE: this is a generative noise model 
	//		  (i.e. noise is generated according
	//        to the noise properties defined here)
	//double noiseSampsSD[3] = {1.0, 1.0, 1.0};		// noise model for samples (std dev along each axis)
	double sampleNoiseInPlane	= (double)cmdOpts.noiseinplane;		// 1.0	(Standard deviation of noise in and out of plane)
	double sampleNoisePerpPlane = (double)cmdOpts.noiseperpplane;	// 1.0					   ''

	// Target Noise Model (for point cloud target only)
	//  NOTE: this is a descriptive model, not a generative one
	//        i.e. no noise is added to the point cloud, the noise model is merely
	//        allow for errors at intermediate locations between the points and penalize
	//        errors offset from the surface
	double PointCloudNoisePerpPlane = 1.0;				// 1.0	(Noise model for point cloud using mesh constructor)
	//  Note: in-plane noise set automatically relative to triangle size

	double rotbounds		= cmdOpts.rbounds;			// DBL_MAX;
	double transbounds		= cmdOpts.tbounds;			// DBL_MAX;
	double scalebounds		= cmdOpts.sbounds;			// 0.3;
	double shapeparambounds = cmdOpts.spbounds;			// 3.0;	

	// load mesh
	loadMeshPath = cmdOpts.target;
	CreateMesh(mesh, loadMeshPath, &saveMeshPath);
	mesh_ssm	 = mesh;								// Initializing mesh_ssm to mesh for non-deformable algorithms 

	if (cmdOpts.deformable)
	{
		// read mode weights
		if (cmdOpts.readModeWeights) {
			shapeparam_read(weight, cmdOpts.modeweights);
			modes += weight.size();
		}
		else
			modes += cmdOpts.modes;
		weight.SetSize(modes - 1);
		std::cout << "\nNumber of modes = " << modes - 1;

		// read ssm
		std::string loadModelPath = cmdOpts.ssm;
		int check = ReadShapeModel(mesh, loadModelPath, modes);		// makes sure mesh is mean mesh
		if (check != 1) {
			std::cout << "\nWARNING: Could not read model data, switching to IMLP..." << std::endl;
			algType = AlgType_IMLP;
		}
		mesh_ssm.vertices = mesh.meanShape;
		mesh_ssm.SavePLY("currentMesh0.ply");

		// if using default input (i.e., mean shape), 
		// create an instance of the shape model to estimate
		if (cmdOpts.useDefaultInput)
		{
			GenerateRandomShapeParams(randSeed3, randSeqPos3, modes - 1, weight);
			// add deformation to mean shape
			for (int j = 0; j < modes - 1; j++)
				for (int i = 0; i < mesh.NumVertices(); i++)
					mesh_ssm.vertices[i] += weight[j] * mesh.wi[j][i];
			mesh_ssm.SavePLY(meshDir);
			shapeparam_write(weight, saveModeWeightsPath);
		}
		// TODO: if reading in an instance of the shape model,
		// compute weights by projecting instance onto shape model
		else if (!cmdOpts.readModeWeights)
			weight.SetAll(0.0);	
	}

	// Create target shape from mesh (as a PD tree)
	if (TargetShapeAsMesh)
	{
		// build PD tree on the mesh directly
		// Note: defines measurement noise to be zero
		printf("\nBuilding mesh PD tree with nThresh: %d and diagThresh: %.2f... \n", nThresh, diagThresh);
		pTree = new PDTree_Mesh(mesh, nThresh, diagThresh);
		//tree.RecomputeBoundingBoxesUsingExistingCovFrames();      // *** is this ever needed?
		printf("Tree built: NNodes=%d  NData=%d  TreeDepth=%d\n\n", pTree->NumNodes(), pTree->NumData(), pTree->TreeDepth());
	}
	else
	{
		printf("\nBuilding point cloud PD tree... \n");
		// build Point Cloud PD tree from mesh
		//  Note: the mesh constructor for the point cloud PD tree
		//        assumes the specified noise in direction perpendicular to surface
		//        and sets the in-plane noise based on the triangle size in order
		//        to allow for greater freedom of match anywhere along the triangle surface
		//        even though each triangle surface is represented by only one point.

		// Uncomment to use triangle center points as the point cloud
		cisstPointCloud pointCloud(mesh, PointCloudNoisePerpPlane);

		// Uncomment to use vertices as the point cloud
		/*for (int i = 0; i < mesh.NumVertices(); i++)
		{
			mesh.vertices[i][0] = mesh.vertices[i][0] * 10;
			mesh.vertices[i][1] = mesh.vertices[i][1] * 10;
			mesh.vertices[i][2] = mesh.vertices[i][2] * 10;
		}
		cisstPointCloud pointCloud(mesh.vertices); */

		PDTree_PointCloud *pPointCloudTree;
		pPointCloudTree = new PDTree_PointCloud(pointCloud, nThresh, diagThresh);
		pTree = pPointCloudTree;
		//tree.RecomputeBoundingBoxesUsingExistingCovFrames();      //*** is this ever needed?
		printf("Tree built: NNodes=%d  NData=%d  TreeDepth=%d\n\n", pTree->NumNodes(), pTree->NumData(), pTree->TreeDepth());
		//printf(" Point Cloud Noise Model:\n  perp-plane variance = %f\n  in-plane variance = %f (avg)\n\n", 
		//  PointCloudNoisePerpPlane, pPointCloudTree->avgVarInPlane);
	}

	// subsample input if need be
	if (!cmdOpts.useDefaultInput)
	{
		// load source mesh
		std::string loadPtsPath;
		loadPtsPath = cmdOpts.input;
		CreateMesh(pts, loadPtsPath, &saveSamplesPath);

		if (!cmdOpts.useDefaultNumSamples)
		{
			nSamples = cmdOpts.samples;
			if (nSamples > pts.NumVertices())
			{
				printf("WARNING: Cannot sample more points than number of vertices in input point cloud/mesh\n");
				nSamples	= pts.NumVertices();
				samples		= pts.vertices;
				sampleNorms = pts.vertexNormals;
			}
			else
				GenerateSamples(pts, randSeed1, randSeqPos1, nSamples,
				samples, sampleNorms,
				&saveSubSamplesPath);
		}
		else
		{
			nSamples	= pts.NumVertices(); 
			samples		= pts.vertices;
			sampleNorms = pts.vertexNormals;
		}
	}
	// Generate random samples from target mesh_ssm 
	else if (TargetShapeAsMesh)
		GenerateSamples(mesh_ssm, randSeed1, randSeqPos1, nSamples,		// mesh_ssm = mesh, if registration algorithm is not deformable
		samples, sampleNorms, sampleDatums,
		&saveSamplesPath);
	else
		// Generate random samples from target point cloud 
		GenerateSamples(mesh_ssm, randSeed1, randSeqPos1, nSamples,		// mesh_ssm = mesh, if registration algorithm is not deformable
		samples, sampleNorms,
		&saveSamplesPath);

	// scale input if need be
	if (!cmdOpts.useDefaultScale)
	{
		for (int i = 0; i < nSamples; i++)
			samples[i] = scale * samples[i];
	}

	std::ifstream randnStream(normRVFile.c_str());  // streams N(0,1) RV's

	// Add noise to samples
	if (!cmdOpts.useDefaultCov || !cmdOpts.useDefaultAxes)
	{
		// Read noise if noise model available
		ReadSampleSurfaceNoise(cmdOpts.useDefaultCov, cmdOpts.useDefaultAxes,
			randSeed1, randSeqPos1, randnStream,
			sampleNoiseInPlane, sampleNoisePerpPlane, 0.0, 0.0,
			samples, sampleNorms,
			noisySamples, noisySampleNorms,
			sampleNoiseCov, sampleNoiseInvCov, sampleNoiseL,
			percentOutliers,
			minPosOffsetOutlier, maxPosOffsetOutlier,
			minAngOffsetOutlier, maxAngOffsetOutlier,
			&saveNoisySamplesPath,
			cmdOpts.cov, &saveCovPath);
	}
	else
	{
		// Generate noise given standard deviations
		GenerateSampleSurfaceNoise(randSeed1, randSeqPos1, randnStream,
			sampleNoiseInPlane, sampleNoisePerpPlane, 0.0, 0.0,
			samples, sampleNorms,
			noisySamples, noisySampleNorms,
			sampleNoiseCov, sampleNoiseInvCov, sampleNoiseL,
			percentOutliers,
			minPosOffsetOutlier, maxPosOffsetOutlier,
			minAngOffsetOutlier, maxAngOffsetOutlier,
			&saveNoisySamplesPath,
			&saveCovPath);
	}

	// If reading in points, assume that noise is already added
	if (!cmdOpts.useDefaultInput) {
		noisySamples	 = samples; 
		noisySampleNorms = sampleNorms; 
	}
	std::cout << "Using " << nSamples << " noisy sample points with " << percentOutliers * 100 << "% outliers...\n";


	// Generate random initial offset
	vctFrm3 Fi;
	if (cmdOpts.useDefaultXfm) {
		GenerateRandomTransform(randSeed2, randSeqPos2,
			minOffsetPos, maxOffsetPos,
			minOffsetAng, maxOffsetAng,
			Fi);
		// save initial offset
		transform_write(Fi, saveOffsetXfmPath);
	}

	// creating ICP solver
	cisstICP ICP;

	// user callbacks for ICP
	std::vector<cisstICP::Callback> userCallbacks;

	//  callback: iteration file (txt)
	cisstICP::Callback iterCallback;
	iterCallback.cbFunc = Callback_SaveIterationsToFile_testICP;
	std::stringstream iterFile;
	iterFile << outputDir << "/SaveIterations.txt";
	std::ofstream iterFileStream(iterFile.str().c_str());
	iterCallback.userData = (void*)(&iterFileStream);
	userCallbacks.push_back(iterCallback);

	//  callback: iteration file (csv)
	cisstICP::Callback iterCallbackCSV;
	iterCallbackCSV.cbFunc = Callback_SaveIterationsToCSV_testICP;
	std::stringstream iterCSV;
	iterCSV << outputDir << "/SaveIterations.csv";
	std::ofstream iterCsvStream(iterCSV.str().c_str());
	iterCallbackCSV.userData = (void*)(&iterCsvStream);
	userCallbacks.push_back(iterCallbackCSV);

	//  callback: track path file
	cisstICP::Callback xfmCallback;
	xfmCallback.cbFunc = Callback_TrackRegPath_testICP;
	std::stringstream trackPathFile;
	trackPathFile << outputDir << "/SaveTrackRegPath.csv"; //"SaveTrackRegPath.txt";
	std::ofstream xfmFileStream(trackPathFile.str().c_str());
	xfmCallback.userData = (void*)(&xfmFileStream);
	userCallbacks.push_back(xfmCallback);

	// Set initial transform
	vctFrm3 FGuess = vctFrm3::Identity();
	if (!cmdOpts.useDefaultXfm)
	{
		transform_read(Fi, cmdOpts.xfm);
		std::cout << std::endl << "Setting initial transform guess: " << std::endl << Fi << std::endl;
		FGuess = Fi;
	}
	else
	{
		std::cout << std::endl << "Applying Sample Offset Fi: " << std::endl << Fi << std::endl;
		FGuess = Fi;
	}
	std::cout << "scale:\n " << scale << std::endl;

	// ICP Algorithm
	algICP *pICPAlg = NULL;
	switch (algType)
	{
	case AlgType_StdICP:
	{
		if (TargetShapeAsMesh)
		{ // target shape is a mesh
			PDTree_Mesh *pTreeMesh = dynamic_cast<PDTree_Mesh*>(pTree);
			pICPAlg = new algICP_StdICP_Mesh(pTreeMesh, noisySamples);
		}
		else
		{ // target shape is a point cloud
			PDTree_PointCloud *pTreePointCloud = dynamic_cast<PDTree_PointCloud*>(pTree);
			pICPAlg = new algICP_StdICP_PointCloud(pTreePointCloud, noisySamples);
		}
		break;
	}
	case AlgType_IMLP:
	{
		if (TargetShapeAsMesh)
		{ // target shape is a mesh
			PDTree_Mesh *pTreeMesh = dynamic_cast<PDTree_Mesh*>(pTree);
			algICP_IMLP_Mesh *pAlg;
			pAlg = new algICP_IMLP_Mesh(pTreeMesh, noisySamples, sampleNoiseCov, sampleNoiseCov);
			// set sample noise model
			//pAlg->SetSampleCovariances( sampleNoiseCov );
			// set mesh noise model to zero noise
			mesh.TriangleCov.SetSize(mesh.NumTriangles());
			mesh.TriangleCovEig.SetSize(mesh.NumTriangles());
			mesh.TriangleCov.SetAll(vct3x3(0.0));
			mesh.TriangleCovEig.SetAll(vct3(0.0));
			pTreeMesh->ComputeNodeNoiseModels();
			//double noiseSDInPlane = 0.5;
			//double noiseSDPerpPlane = 1.0;
			//SetMeshTriangleCovariances( mesh, noiseSDInPlane, noiseSDPerpPlane );
			pICPAlg = pAlg;
		}
		else
		{ // target shape is a point cloud
			PDTree_PointCloud *pTreePointCloud = dynamic_cast<PDTree_PointCloud*>(pTree);
			algICP_IMLP_PointCloud *pAlg;
			pAlg = new algICP_IMLP_PointCloud(pTreePointCloud, noisySamples, sampleNoiseCov, sampleNoiseCov);
			// set sample noise model
			//pAlg->SetSampleCovariances( sampleNoiseCov );
			pICPAlg = pAlg;
			// NOTE: (PD tree noise model was already defined
			//       by the point cloud cov tree constructor)
		}
		break;
	}
	case AlgType_DIMLP:
	{
		std::cout << "shape parameters:\n" << weight << std::endl;
		if (!TargetShapeAsMesh)
		{
			std::cout << "ERROR: Currently only mesh target supported for DIMLP" << std::endl;
			assert(0);
		}
		PDTree_Mesh *pTreeMesh = dynamic_cast<PDTree_Mesh*>(pTree);
		algICP_DIMLP *pAlg;
		pAlg = new algICP_DIMLP(pTreeMesh, noisySamples, sampleNoiseCov, sampleNoiseCov, mesh.meanShape, 1, bScale);

		pAlg->SetConstraints(rotbounds, transbounds, scalebounds, shapeparambounds);

		mesh.TriangleCov.SetSize(mesh.NumTriangles());
		mesh.TriangleCovEig.SetSize(mesh.NumTriangles());
		mesh.TriangleCov.SetAll(vct3x3(0.0));
		pTreeMesh->ComputeNodeNoiseModels();
		pICPAlg = pAlg;

		break;
	}
	default:
	{
		std::cout << "ERROR: unknown algorithm type" << std::endl;
		assert(0);
	}
	}

	cisstMesh samplePts;
	samplePts.vertices.SetSize(noisySamples.size());
	samplePts.vertexNormals.SetSize(noisySampleNorms.size());
	samplePts.vertices		= noisySamples;
	samplePts.vertexNormals = noisySampleNorms;
	samplePts.SavePLY(outputDir + "/Pts.ply");

	// ICP Options
	cisstICP::Options opt;							// TODO: Make these command line inputs
	opt.auxOutputDir	= outputDir;
	opt.maxIter			= maxIters;
	opt.termHoldIter	= 2;
	opt.numShapeParams	= modes - 1;
	opt.minE			= -std::numeric_limits<double>::max();
	opt.tolE			= 0.0;
	opt.dPosThresh		= 0.1;
	opt.dAngThresh		= 0.1*(cmnPI / 180);
	opt.dShapeThresh	= 0.1;
	opt.dPosTerm		= 0.001;					// 0.01;
	opt.dAngTerm		= 0.001*(cmnPI / 180);		// 0.01*(cmnPI / 180);
	opt.dShapeTerm		= 0.001;					// 0.01; 
	opt.deformable		= cmdOpts.deformable;

	// Run ICP
	vctFrm3 Freg;
	int numRuns		= 1;
	double runtime	= 0.0;
	cisstICP::ReturnType rv;
	for (int i = 0; i < numRuns; i++)
	{
		rv = ICP.RunICP(pICPAlg, opt, FGuess, &userCallbacks);
		//rv = ICP.RunICP_Rigid( samples, *pTree, opt, FGuess, Freg );
		std::cout << rv.termMsg;
		iterFileStream << rv.termMsg;
		runtime += rv.runTime;
		Freg = rv.Freg;
	}
	std::cout << std::endl << " ===> Avg RunTime: " << runtime / numRuns << std::endl;

	// save registration result
	transform_write(Freg, saveRegXfmPath);

	// Freg now includes Fi as FGuess => Freg should be identity for perfect registration
	vctFrm3 Ferr	= Freg;
	//vctFrm3 Ferr	= Fi * Freg;
	vctRodRot3 Rerr(Ferr.Rotation());
	double terr		= Ferr.Translation().Norm();
	double rerr		= Rerr.Norm();
	vctRodRot3 Rinit(Fi.Rotation());
	double tinit	= Fi.Translation().Norm();
	double rinit	= Rinit.Norm();

	std::stringstream resultStream;
	resultStream << std::endl;
	if (cmdOpts.deformable)
	{
		resultStream << "Starting Offset:   \tdAng: " << rinit * 180 / cmnPI << "\tdPos: " << tinit << "\tdShp: " << weight;
		if (cmdOpts.bScale)
			resultStream << "\tdScl: " << scale << std::endl;
		else
			resultStream << std::endl;
		pICPAlg->ReturnScale(scale);
		resultStream << "Registration Error:\tdAng: " << rerr * 180 / cmnPI << "\tdPos: " << terr << "\tdShp: " << weight - mesh.Si;
		if (cmdOpts.bScale)
			resultStream << "\tdScl: " << scale << std::endl;
		else
			resultStream << std::endl;
	}
	else
	{
		resultStream << "Starting Offset:   \tdAng: " << rinit * 180 / cmnPI << "\tdPos: " << tinit << std::endl;
		resultStream << "Registration Error:\tdAng: " << rerr * 180 / cmnPI << "\tdPos: " << terr << std::endl;
	}
	resultStream << std::endl << "Average Match Error (Mahalanobis):\t" << rv.MatchPosErrAvg << " (+/-" << rv.MatchPosErrSD << ")\n" << std::endl;
	std::cout << resultStream.str();
	iterFileStream << resultStream.str();

	std::cout << "=============================================================\n" << std::endl;

	if (cmdOpts.deformable)
	{
		mesh.vertices = mesh.meanShape;
		for (int s = 0; s < mesh.NumVertices(); s++)
			for (unsigned int i = 0; i < (unsigned int)(modes - 1); i++)
				mesh.vertices(s) += (mesh.Si[i] * mesh.wi[i].Element(s));
	}
	mesh.SavePLY(outputDir + "/finalMesh.ply");

	for (int i = 0; i < mesh.NumVertices(); i++)
		mesh.vertices(i) = Freg * mesh.vertices(i);
	mesh.SavePLY(outputDir + "/finalRegMesh.ply");

	for (int i = 0; i < noisySamples.size(); i++)
		samplePts.vertices[i] = Fi * noisySamples[i];
	samplePts.SavePLY("currentSamples0.ply");
	samplePts.SavePLY(outputDir + "/initPts.ply");

	for (int i = 0; i < noisySamples.size(); i++) 
		samplePts.vertices[i] = Freg * noisySamples[i] * scale;
	samplePts.SavePLY(outputDir + "/finalPts.ply");

	if (pICPAlg) delete pICPAlg;
}

#endif // _testICP_H
