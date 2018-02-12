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


enum ICPAlgType { AlgType_StdICP, AlgType_IMLP, AlgType_DIMLP , AlgType_VIMLOP};

void Callback_TrackRegPath_testICP( cisstICP::CallbackArg &arg, void *userData )
{
  // Save to file:
  //  - error function
  //  - incremental transform
  // output format:
  //  error r00 r01 r02 r10 r11 r12 r20 r21 r22 tx ty tz

  //std::ofstream *fs = (std::ofstream *)(userData);
  //(*fs) << arg.E << " " << arg.dF.Rotation().Row(0) << " " << arg.dF.Rotation().Row(1) << " " 
  //  << " " << arg.dF.Rotation().Row(2) << " " << arg.dF.Translation() << std::endl;

	
  // output format:
  //  iter error r00 r01 r02 r10 r11 r12 r20 r21 r22 tx ty tz s1 ... sn
  std::ofstream *fs = (std::ofstream *)(userData);
  (*fs) << arg.iter << ", , " << arg.E << ", , "
	  << arg.Freg.Rotation().Row(0)(0) << ", " << arg.Freg.Rotation().Row(0)(1) << ", " << arg.Freg.Rotation().Row(0)(2) << ", , "
	  << arg.Freg.Rotation().Row(1)(0) << ", " << arg.Freg.Rotation().Row(1)(1) << ", " << arg.Freg.Rotation().Row(1)(2) << ", , "
	  << arg.Freg.Rotation().Row(2)(0) << ", " << arg.Freg.Rotation().Row(2)(1) << ", " << arg.Freg.Rotation().Row(2)(2) << ", , "
	  << arg.Freg.Translation()(0) << ", " << arg.Freg.Translation()(1) << ", " << arg.Freg.Translation()(2) << ", ,"
	  << arg.scale << ", ," ;
  for (int m = 0; m < arg.S.size(); m++)
	  (*fs) << arg.S(m) << ",";
  (*fs) << std::endl;
}
void Callback_SaveIterationsToFile_testICP( cisstICP::CallbackArg &arg, void *userData )
{
  std::ofstream *fs = (std::ofstream *)(userData);

  vctRodRot3 dR(arg.dF.Rotation());
  std::stringstream ss;
  //ss << cmnPrintf("iter=%u  E=%.3f  tolE=%.4f  (dAng/dPos)= %.2f/%.2f  t=%.3f  NNodes=%u/%u/%u  NOut=%u")
  ss << cmnPrintf("iter=%u  E=%.3f  tolE=%.4f  (dAng/dPos)= %.2f/%.2f  t=%.3f  NOut=%u")
    << arg.iter 
    << arg.E 
    << arg.tolE
    << dR.Norm()*180/cmnPI << arg.dF.Translation().Norm()
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
	//char k;
	//std::cout << "press key then enter to start:" << std::endl;
	//std::cin >> k;

	// Set output directories
	std::string workingDir;
	std::string outputDir;
	std::string algDir;

	//if (cmdOpts.useDefaultOutput)
		workingDir = "../../../test_data/";
	//else
	//	workingDir = cmdOpts.output;

	switch (algType)
	{
	case AlgType_StdICP:
	{
		std::cout << "\nRunning standard ICP\n" << std::endl;
		algDir = "LastRun_StdICP/";
		break;
	}
	case AlgType_IMLP:
	{
		std::cout << "\nRunning IMLP\n" << std::endl;
		algDir = "LastRun_IMLP/";
		break;
	}
	case AlgType_DIMLP:
	{
		std::cout << "\nRunning D-IMLP\n" << std::endl;
		algDir = "LastRun_DIMLP/";
		break;
	}
	case AlgType_VIMLOP:
	{
		std::cout << "\nRunning V-IMLOP\n" << std::endl;
		algDir = "LastRun_VIMLOP/";
		break;
	}
	default:
	{
		std::cout << "ERROR: unknown algorithm type" << std::endl;
		assert(0);
	}
	}

	if (cmdOpts.useDefaultOutput)
		outputDir = workingDir + algDir;
	else
		outputDir = workingDir + algDir + cmdOpts.output + "/";
	CreateDirectory(outputDir.c_str(), NULL);

	// input files
	std::string normRVFile				= workingDir + "GaussianValues.txt";	// Random Numbers: Normal RV's

	// output files
	std::string saveMeshPath			= outputDir + "SaveMesh";
	std::string saveSamplesPath			= outputDir + "SaveSamples.pts";
	std::string saveSubSamplesPath		= outputDir + "SaveSubSamples.pts";
	std::string saveNoisySamplesPath	= outputDir + "SaveNoisySamples.pts";
	std::string saveCovPath				= outputDir + "SaveSampleCov.txt";
	std::string saveOffsetXfmPath		= outputDir + "SaveOffsetXfm.txt";
	std::string saveRegXfmPath			= outputDir + "SaveRegXfm.txt";
	std::string saveModeWeightsPath		= outputDir + "saveModeWeights.txt";
	std::string meanMesh				= outputDir + "meanMesh.ply";
	std::string meshDir					= outputDir + "estimatedMesh.ply";
	//std::string saveLPath				= outputDir + "SaveSampleL";

	// Declare variables
	cisstMesh						mesh, 
									mesh_ssm, 
									pts;
	PDTreeBase*						pTree;
	vctDynamicVector<vct3>			samples;
	vctDynamicVector<vct3>			sampleNorms;
	vctDynamicVector<vct3>			noisySamples;
	vctDynamicVector<vct3>			noisySampleNorms;
	vctDynamicVector<unsigned int>  sampleDatums;
	vctDynamicVector<vct3x3>		sampleNoiseCov;
	vctDynamicVector<vct3x3>		sampleNoiseInvCov;
	vctDynamicVector<vct3x2>		sampleNoiseL;

	vctDynamicVector<vct2>			noisyEdgesV1;
	vctDynamicVector<vct2>			noisyEdgesV2;
	vctDynamicVector<vct2>			noisyNorms2d;
	vctDynamicVector<vct2x2>		sampleNoiseCov2d;
	vctDynamicVector<vct2x2>		normNoiseCov2d;

	vctDynamicVector<double>		weight;

	// Set default values/initialize variables
	int		nSamples	= 300;	// default number of samples (for default input)
	int		maxIters	= 100;
	double	scale		= 1.0;
	bool	bScale		= false;

	int		nThresh		= 5;	// 5;	// Cov Tree Params
	double	diagThresh = 5.0; // 5.0;	// 5.0;	//  ''

	double minOffsetPos = 10.0;	// 50.0; 
	double maxOffsetPos = 20.0;	// 100.0; 
	double minOffsetAng = 6.0;	// 30.0; 
	double maxOffsetAng = 12.0;	// 60.0;  

	double percentOutliers		= 0.0;	// 0.0; // 0.05;
	double minPosOffsetOutlier	= 2.0;	// 5.0;
	double maxPosOffsetOutlier	= 5.0;	// 10.0;
	double minAngOffsetOutlier	= 2.0;	// 0.0;
	double maxAngOffsetOutlier	= 5.0;	// 0.0;

	std::srand(time(NULL)); unsigned int randSeed1	 = /*0;		*/std::rand();	// generates samples
	std::srand(time(NULL)); unsigned int randSeqPos1 = /*0;		*/std::rand();
	std::srand(time(NULL)); unsigned int randSeed2	 = /*17;	*/std::rand();	// generates offsets
	std::srand(time(NULL)); unsigned int randSeqPos2 = /*28;	*/std::rand();
	std::srand(time(NULL)); unsigned int randSeed3	 = /*28;	*/std::rand();	// generates shape parameters
	std::srand(time(NULL)); unsigned int randSeqPos3 = /*8;		*/std::rand();

	// Samples Noise Model
	//  NOTE: this is a generative noise model 
	//		  (i.e. noise is generated according
	//        to the noise properties defined here)
	//double noiseSampsSD[3] = {1.0, 1.0, 1.0};		// noise model for samples (std dev along each axis)
	double sampleNoiseInPlane = 1.0;				// standard deviation of noise in and out of plane
	double sampleNoisePerpPlane = 1.0;				//   ''

	// Target Noise Model (for point cloud target only)
	//  NOTE: this is a descriptive model, not a generative one
	//        i.e. no noise is added to the point cloud, the noise model is merely
	//        allow for errors at intermediate locations between the points and penalize
	//        errors offset from the surface
	double PointCloudNoisePerpPlane = 1.0; // 0.5;			// noise model for point cloud using mesh constructor
	//  Note: in-plane noise set automatically relative to triangle size

	double rotbounds = DBL_MAX;
	double transbounds = DBL_MAX;
	double scalebounds = 0.3;
	double shapeparambounds = 3.0;
	
	// Modify default values
	if (!cmdOpts.useDefaultNumSamples)
		nSamples = cmdOpts.samples;
	if (!cmdOpts.useDefaultNumIters)
		maxIters = cmdOpts.niters;

	if (!cmdOpts.useDefaultMinPos)
		minOffsetPos = (double)cmdOpts.minpos;
	if (!cmdOpts.useDefaultMaxPos)
		maxOffsetPos = (double)cmdOpts.maxpos;
	if (!cmdOpts.useDefaultMinAng)
		minOffsetAng = (double)cmdOpts.minang;
	if (!cmdOpts.useDefaultMaxAng)
		maxOffsetAng = (double)cmdOpts.maxang;

	if (!cmdOpts.useDefaultNumOutliers)
		percentOutliers = cmdOpts.poutliers/100.0;
	if (!cmdOpts.useDefaultOutMinPos)
		minPosOffsetOutlier = (double)cmdOpts.outminpos;
	if (!cmdOpts.useDefaultOutMaxPos)
		maxPosOffsetOutlier = (double)cmdOpts.outmaxpos;
	if (!cmdOpts.useDefaultOutMinAng)
		minAngOffsetOutlier = (double)cmdOpts.outminang;
	if (!cmdOpts.useDefaultOutMaxAng)
		maxAngOffsetOutlier = (double)cmdOpts.outmaxang;

	if (!cmdOpts.useDefaultNoiseInPlane)
		sampleNoiseInPlane = (double)cmdOpts.noiseinplane;
	if (!cmdOpts.useDefaultNoisePerpPlane)
		sampleNoisePerpPlane = (double)cmdOpts.noiseperpplane;

#if 1
	// load mesh
	std::string loadMeshPath;
	if (cmdOpts.useDefaultTarget)
		loadMeshPath = workingDir + "MT.ply";
		//std::string loadMeshPath = workingDir + "ProximalFemur.ply";
		//std::string loadMeshPath = workingDir + "RIGHTHEMIPELVIS_centered.mesh";
		//std::string loadMeshPath = workingDir + "RIGHTHEMIPELVIS.mesh";
		//std::string loadMeshPath = workingDir + "CTBreastImage_Dec20000_Shell.mesh";  
	else
		loadMeshPath = cmdOpts.target;
	CreateMesh(mesh, loadMeshPath, &saveMeshPath);

	mesh_ssm = mesh;	// Initializing mesh_ssm to mesh for non-deformable algorithms 
	int modes = 1;		// +1 (for mean)
	
	if (algType == AlgType_DIMLP)
	{
		// load ssm
		std::string loadModelPath;
		if (cmdOpts.useDefaultSSM)
			loadModelPath = workingDir + "atlas_mt.txt";
		else
			loadModelPath = cmdOpts.ssm;

		if (cmdOpts.useDefaultNumModes)
			modes += 3;
		else	  
		{
			if (cmdOpts.readModeWeights) {
				shapeparam_read(weight, cmdOpts.modeweights);
				modes += weight.size();
			}
			else
				modes += cmdOpts.modes;
		}
		std::cout << "Number of modes = " << modes << std::endl;
		
		int check = ReadShapeModel(mesh, loadModelPath, modes);		// makes sure mesh is mean mesh
		if (check != 1) {
			std::cout << "Unsuccessful in reading model data, switching to IMLP..." << std::endl;
			algType = AlgType_IMLP;
		}
		mesh_ssm.vertices = mesh.meanShape;
		mesh_ssm.SavePLY("currentMesh0.ply");

		if (cmdOpts.useDefaultInput)
		{
			weight.SetSize(modes - 1);
			GenerateRandomShapeParams(randSeed3, randSeqPos3, modes - 1, weight);
			// add deformation to mean shape
			for (int j = 0; j < modes - 1; j++)
				for (int i = 0; i < mesh.NumVertices(); i++)
					mesh_ssm.vertices[i] += weight[j] * mesh.wi[j][i];
			mesh_ssm.SavePLY(meshDir);
			shapeparam_write(weight, saveModeWeightsPath);
		}
		else
		{
			weight.SetSize(modes - 1);
			if (!cmdOpts.readModeWeights) {
				weight.SetAll(0.0);
			}
		}

		if (!cmdOpts.useDefaultShapeParamBounds)
			shapeparambounds = cmdOpts.spbounds;
	}

	if (!cmdOpts.useDefaultRotationBounds)
		rotbounds = cmdOpts.rbounds;

	if (!cmdOpts.useDefaultTranslationBounds)
		transbounds = cmdOpts.tbounds;

	if (!cmdOpts.useDefaultScaleBounds)
		scalebounds = cmdOpts.sbounds;

	// Create target shape from mesh (as a PD tree)
	if (TargetShapeAsMesh)
	{
		// build PD tree on the mesh directly
		//  Note: defines measurement noise to be zero
		printf("Building mesh PD tree... ");
		pTree = new PDTree_Mesh(mesh, nThresh, diagThresh);
		//tree.RecomputeBoundingBoxesUsingExistingCovFrames();      //*** is this ever needed?
		printf("Tree built: NNodes=%d  NData=%d  TreeDepth=%d\n", pTree->NumNodes(), pTree->NumData(), pTree->TreeDepth());
	}
	else
	{
		// build Point Cloud PD tree from mesh
		// (uses triangle center points as the point cloud)
		//  Note: the mesh constructor for the point cloud PD tree
		//        assumes the specified noise in direction perpendicular to surface
		//        and sets the in-plane noise based on the triangle size in order
		//        to allow for greater freedom of match anywhere along the triangle surface
		//        even though each triangle surface is represented by only one point.
		printf("Building point cloud PD tree... ");
		cisstPointCloud pointCloud(mesh, PointCloudNoisePerpPlane);
		PDTree_PointCloud *pPointCloudTree;
		pPointCloudTree = new PDTree_PointCloud(pointCloud, nThresh, diagThresh);
		pTree = pPointCloudTree;
		//tree.RecomputeBoundingBoxesUsingExistingCovFrames();      //*** is this ever needed?
		printf("Tree built: NNodes=%d  NData=%d  TreeDepth=%d\n", pTree->NumNodes(), pTree->NumData(), pTree->TreeDepth());
		//printf(" Point Cloud Noise Model:\n  perp-plane variance = %f\n  in-plane variance = %f (avg)\n\n", 
		//  PointCloudNoisePerpPlane, pPointCloudTree->avgVarInPlane);
	}

	// sub-sample and scale input if need be
	if (!cmdOpts.useDefaultInput)
	{
		std::string loadPtsPath;
		loadPtsPath = cmdOpts.input;
		CreateMesh(pts, loadPtsPath, &saveSamplesPath);

		if (!cmdOpts.useDefaultNumSamples)
		{
			int nSubSamples = cmdOpts.samples;
			vctDynamicVector<vct3> subsampledPts;
			vctDynamicVector<bool> selectedPts;

			subsampledPts.SetSize(nSubSamples);
			selectedPts.SetSize(pts.NumVertices());

			selectedPts.SetAll(false);
			subsampledPts.Zeros();

			GenerateSubSamples(pts, selectedPts, subsampledPts, nSubSamples, &saveSubSamplesPath); // change this to be more similar to generatesamples
			//GenerateSamples(pts, randSeed1, randSeqPos1, nSubSamples,		
			//	samples, sampleNorms, sampleDatums,
			//	&saveSamplesPath);
		}
		else
			nSamples = pts.NumVertices(); //samples.size();

		if (!cmdOpts.useDefaultScale)
		{
			scale = cmdOpts.scale; // (((double)rand() / (double)RAND_MAX) /** (1.5 - 0.5)*/) + 0.5; 
			for (int i = 0; i < nSamples; i++)
				//samples[i] = scale * samples[i];
				pts.vertices[i] = scale * pts.vertices[i];
		}
		if (cmdOpts.bScale)
			bScale = true;
	}
	std::ifstream randnStream(normRVFile.c_str());  // streams N(0,1) RV's

	// Generate random samples from target mesh_ssm 
	GenerateSamples(mesh_ssm, randSeed1, randSeqPos1, nSamples,		// mesh_ssm = mesh, if registration algorithm is not deformable
		samples, sampleNorms, sampleDatums,
		&saveSamplesPath);

	if (!cmdOpts.useDefaultInput) {
		for (int i = 0; i < nSamples; i++)
		{
			samples[i] = pts.vertices[i];
			//sampleNorms[i] = pts.vertexNormals[i];
		}
	}

	// Add noise to samples
	if (!cmdOpts.useDefaultCov || !cmdOpts.useDefaultAxes)
	{
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

	if (!cmdOpts.useDefaultInput && cmdOpts.useDefaultNumSamples) {
		for (int i = 0; i < nSamples; i++)
		{
			noisySamples[i] = pts.vertices[i];
			//noisySampleNorms[i] = pts.vertexNormals[i];
		}
	}
	std::cout << "Using " << nSamples << " noisy sample points";
	if (!cmdOpts.useDefaultNumOutliers)
		std::cout << " with " << percentOutliers << " outliers...\n";
	else
		std::cout << "...\n";

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

	//noisySamples2d.SetSize(1);
	noisyEdgesV1.SetSize(1);
	noisyEdgesV2.SetSize(1);
	noisyNorms2d.SetSize(1);
	sampleNoiseCov2d.SetSize(1);
	normNoiseCov2d.SetSize(1);

#else
	// Replay Randomized Trial
	vctFrm3 Fi;
	std::string baseFolder = "..\\ICP_TestData\\RandomTests\\";
	std::string noiseDir = "NoiseIP0.5PP0.5_OffsetR0-10P0-10\\";
	std::string trialName = "IMLP";     // this can be any of the trials (since Finit and samples are same for all trials)
	std::string trialNum = "1";
	std::string commonDir = "CommonFiles\\";
	std::string loadMeshFile = baseFolder + commonDir + "SaveMesh.mesh";
	//std::string loadSamplesFile = baseFolder + noiseLevel + commonDir + "SaveSamples_" + trialNum + ".pts";
	std::string loadNoisySamplesFile = baseFolder + noiseDir + commonDir + "SaveNoisySamples_" + trialNum + ".pts";
	// need this to get FGuess
	std::string loadTrackPathFile = baseFolder + noiseDir + trialName + "\\SaveTrackRegPath_" + trialName + "_" + trialNum + ".txt";

	//std::string workingFolder = "..\\ICP_TestData\\LastRun_CovEst\\";
	//std::string loadMeshFile = workingFolder + "SaveMesh.mesh";
	//std::string loadSamplesFile = workingFolder + "SaveSamples" + ".pts";
	//std::string loadNoisySamplesFile = workingFolder + "SaveNoisySamples" + ".pts";
	//std::string loadTrackPathFile = workingFolder + "SaveTrackRegPath" + ".txt";

	mesh.LoadMeshFile( loadMeshFile );
	cisstICP::MeshSave(mesh, saveMeshFile);
	pTree = new PDTree_PointCloud(mesh, nThresh, diagThresh);
	printf("Tree built: NNodes=%d  NData=%d  TreeDepth=%d\n", pTree->NumNodes(), pTree->NumData(), pTree->TreeDepth());
	cisstICP::SamplesLoad( noisySamples, loadNoisySamplesFile );
	cisstICP::SamplesSave( noisySamples, saveNoisySamplesFile );
	std::ifstream fs(loadTrackPathFile.c_str());
	double temp;
	double r00,r01,r02,r10,r11,r12,r20,r21,r22,p1,p2,p3;
	fs >> temp >> r00 >> r01 >> r02 >> r10 >> r11 >> r12 >> r20 >> r21 >> r22 >> p1 >> p2 >> p3;
	vctRot3 Rm( 
		r00, r01, r02,
		r10, r11, r12,
		r20, r21, r22,
		VCT_NORMALIZE );
	vct3 pm( p1,p2,p3 );
	Fi.Assign(Rm,pm);

	std::cout << Fi << std::endl;

	double meshstdDevPerpPlanePlane = 0.0;         // noise model for mesh
	double meshStdDevInPlane = 0.0;
	double noisePosSD[3] = {0.5, 0.5, 0.5};   // noise model for samples (std dev along each axis)
#endif

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
	std::cout << "scale:\n   " << scale << std::endl;

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
		if (bScale)
			pAlg = new algICP_DIMLP(pTreeMesh, noisySamples, sampleNoiseCov, sampleNoiseCov, mesh.meanShape, scale, bScale);
		else
			pAlg = new algICP_DIMLP(pTreeMesh, noisySamples, sampleNoiseCov, sampleNoiseCov, mesh.meanShape, scale, bScale); 
																		// ^ for cases when scale is specified, but must not be optimized over

		pAlg->SetConstraints(rotbounds, transbounds, scalebounds, shapeparambounds);

		mesh.TriangleCov.SetSize(mesh.NumTriangles());
		mesh.TriangleCovEig.SetSize(mesh.NumTriangles());
		mesh.TriangleCov.SetAll(vct3x3(0.0));
		pTreeMesh->ComputeNodeNoiseModels();
		pICPAlg = pAlg;

		break;
	}
	//case AlgType_VIMLOP:
	//{
	//	PDTree_Mesh *pTreeMesh = dynamic_cast<PDTree_Mesh*>(pTree);
	//	DirPDTree2D_Edges *pDirTree = dynamic_cast<DirPDTree2D_Edges*>(pDirTree);
	//	camera cam(500, 500, 1, 1, 0, 0);
	//	algDirICP_VIMLOP *pAlg;
	//	pAlg = new algDirICP_VIMLOP(pTreeMesh, pDirTree, noisySamples, sampleNoiseCov, sampleNoiseCov, noisyEdgesV1, noisyEdgesV2, noisyNorms2d, sampleNoiseCov2d, normNoiseCov2d, cam);
	//
	//	mesh.TriangleCov.SetSize(mesh.NumTriangles());
	//	mesh.TriangleCovEig.SetSize(mesh.NumTriangles());
	//	mesh.TriangleCov.SetAll(vct3x3(0.0));
	//	pTreeMesh->ComputeNodeNoiseModels();
	//	pICPAlg = pAlg;
	//
	//	break;
	//}
	default:
	{
		std::cout << "ERROR: unknown algorithm type" << std::endl;
		assert(0);
	}
	}

	cisstMesh samplePts;
	samplePts.vertices.SetSize(noisySamples.size());
	samplePts.vertices = noisySamples;
	samplePts.SavePLY(outputDir + "/Pts.ply");

	// ICP Options
	cisstICP::Options opt;
	opt.auxOutputDir = outputDir;
	opt.maxIter = maxIters;
	opt.termHoldIter = 2;
	opt.numShapeParams = modes - 1;
	opt.minE = -std::numeric_limits<double>::max();
	opt.tolE = 0.0;
	opt.dPosThresh = 0.1; 
	opt.dAngThresh = 0.1*(cmnPI / 180); 
	opt.dShapeThresh = 0.1;  
	opt.dPosTerm = 0.001; // 0.01;
	opt.dAngTerm = 0.001*(cmnPI / 180); // 0.01*(cmnPI / 180);
	opt.dShapeTerm = 0.001; // 0.01; 
	if (cmdOpts.deformable)
		opt.deformable = true;

	// Run ICP
	int numRuns = 1;
	vctFrm3 Freg;
	double runtime = 0.0;
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
	vctFrm3 Ferr = Freg;
	//vctFrm3 Ferr = Fi * Freg;
	vctRodRot3 Rerr(Ferr.Rotation());
	double terr = Ferr.Translation().Norm();
	double rerr = Rerr.Norm();
	vctRodRot3 Rinit(Fi.Rotation());
	double tinit = Fi.Translation().Norm();
	double rinit = Rinit.Norm();

	std::stringstream resultStream;
	resultStream << std::endl;
	if (algType == AlgType_DIMLP)
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
	
	if (algType == AlgType_DIMLP)
	{
		mesh.vertices = mesh.meanShape;
		for (int s = 0; s < mesh.NumVertices(); s++)
		{
			for (unsigned int i = 0; i < (unsigned int)(modes - 1); i++)
			{
				mesh.vertices(s) += (mesh.Si[i] * mesh.wi[i].Element(s));
			}
		}
	}
	mesh.SavePLY(outputDir + "/finalMesh.ply");

	for (int i = 0; i < mesh.NumVertices(); i++)
		mesh.vertices(i) = Freg * mesh.vertices(i);
	mesh.SavePLY(outputDir + "/finalRegMesh.ply");

	for (int i = 0; i < noisySamples.size(); i++)
		samplePts.vertices[i] = Fi * noisySamples[i];
	samplePts.SavePLY("currentSamples0.ply");
	samplePts.SavePLY(outputDir + "/initPts.ply");

	for (int i = 0; i < noisySamples.size(); i++) {
		if (bScale)
			samplePts.vertices[i] = Freg * noisySamples[i] * scale;
		else
			samplePts.vertices[i] = Freg * noisySamples[i];
	}
	samplePts.SavePLY(outputDir + "/finalPts.ply");

	if (pICPAlg) delete pICPAlg;

  //std::cout << "press key then enter to quit:" << std::endl;
  //std::cin >> k;
}

#endif // _testICP_H
