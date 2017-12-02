#ifndef _testICPNormals_H
#define _testICPNormals_H

#include <stdio.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <limits.h>

#include "utility.h"
#include "cisstICP.h"
#include "cisstMesh.h"
#include "DirPDTree_Mesh.h"
#include "DirPDTree_PointCloud.h"

#include "algDirICP_StdICP_PointCloud.h"
#include "algDirICP_StdICP_Mesh.h"
#include "algDirICP_IMLOP_Mesh.h"
#include "algDirICP_GIMLOP_Mesh.h"
#include "algDirICP_PIMLOP_Mesh.h"
#include "algDirICP_DIMLOP.h"
#include "algDirICP_GDIMLOP.h"

// disable for run-time tests
#define ENABLE_CALLBACKS


enum ICPDirAlgType { DirAlgType_StdICP, DirAlgType_IMLOP, DirAlgType_DIMLOP, DirAlgType_GIMLOP, DirAlgType_GDIMLOP, DirAlgType_PIMLOP };

void Callback_TrackRegPath_testICPNormals(cisstICP::CallbackArg &arg, void *userData)
{
  // cast to norm callback arg type
  //cisstICP::CallbackArg *argp = &arg;
  //argp = dynamic_cast<cisstICP::CallbackArgNormals*>(&arg);
  //if (argp==0)
  //{
  //  std::cerr << "ERROR: cannot cast callback argument to cisstICP type" << std::endl;
  //  assert(argp);   // terminate program
  //}

  // Save to file:
  //  - error function (-loglik)
  //  - incremental transform
  //  - IMLOP params
  //  - residual match errors
  // output format:
  //  err r00 r01 r02 r10 r11 r12 r20 r21 r22 tx ty tz normWeight posWeight avgNormError avgPosError 

  //std::ofstream *fs = (std::ofstream *)(userData);
  //(*fs) << argp->E << " " << argp->dF.Rotation().Row(0) << " " << argp->dF.Rotation().Row(1) << " "
  //  << argp->dF.Rotation().Row(2) << " " << argp->dF.Translation() << " "
  //  //<< argp->normWeight << " " << argp->posWeight << " "
  //  //<< argp->MatchPosErrorAvg << " " << argp->MatchPosErrorSD << " "
  //  //<< argp->MatchNormErrorAvg << " " << argp->MatchNormErrorSD 
  //  << std::endl;

  // output format:
  //  iter err r00 r01 r02 r10 r11 r12 r20 r21 r22 tx ty tz s1 ... sn
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
void Callback_SaveIterationsToFile_testICPNormals(cisstICP::CallbackArg &arg, void *userData)
{
  std::ofstream *fs = (std::ofstream *)(userData);

  // cast to norm callback arg type
  //cisstICP::CallbackArg *argp;
  //argp = dynamic_cast<cisstICP::CallbackArgNormals*>(&arg);
  //if (argp==0)
  //{
  //  std::cerr << "ERROR: cannot cast callback argument to cisstICP type" << std::endl;
  //  assert(argp);   // terminate program
  //}

  vctRodRot3 dR(arg.dF.Rotation());
  std::stringstream ss;
  //ss << cmnPrintf("iter=%u  E=%.3f  tolE=%.4f  nW/pW=%.4f/%.4f  (cSD/dSD)=%.2f/%.2f  dTheta=%.2f/%.2f/%.2f  (dAng/dPos)= %.2f/%.2f  t=%.4f NNodes=%u/%u/%u NOut=%u/%u")
  ss << cmnPrintf("iter=%u  E=%.3f  tolE=%.4f  (dAng/dPos)= %.2f/%.2f  t=%.4f  NOut=%u/%u")
	  // (RMS/Res)=%.4f/%.4f
	  //  dTheta=%.2f/%.2f/%.2f  
	  << arg.iter //<< argp->iter
	  << arg.E //<< argp->E
	  << arg.tolE //<< argp->tolE
	  //<< argp->normWeight
	  //<< argp->posWeight
	  //<< argp->circSD*180.0 / cmnPI << sqrt(argp->posVar)
	  //<< argp->dThetaMin*180.0 / cmnPI << argp->dThetaMax*180.0 / cmnPI << argp->dThetaAvg*180.0 / cmnPI
	  << dR.Norm() * 180 / cmnPI << arg.dF.Translation().Norm()
	  << arg.time //<< argp->time
	  << arg.nOutliers; //<< argp->nOutliers;
  //<< argp->maxNodesSearched << argp->avgNodesSearched << argp->minNodesSearched
  //<< argp->nOutliersPos
  //<< argp->nOutliersNorm;    

  (*fs) << ss.str() << std::endl;
}
void Callback_SaveIterationsToCSV_testICPNormals(cisstICP::CallbackArg &arg, void *userData)
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
void testICPNormals(bool TargetShapeAsMesh, ICPDirAlgType algType, cisstICP::CmdLineOptions cmdOpts)
{
  //char k;
  //std::cout << "press key then enter to start:" << std::endl;
  //std::cin >> k;

  // Set output directories
  std::string workingDir;
  std::string outputDir;
  std::string algDir; // = "LastRun//";

  //if (cmdOpts.useDefaultOutput)
  workingDir = "F:/Research/SinusProject/Seth_code/cisstICP/cissticp/test_data/";
  //else
	 // workingDir = cmdOpts.output;

  switch (algType)
  {
  case DirAlgType_StdICP:
  {
	  std::cout << "\nRunning directional ICP\n" << std::endl;
	  algDir = "LastRun_DirICP/";
	  break;
  }
  case DirAlgType_IMLOP:
  {
	  std::cout << "\nRunning IMLOP\n" << std::endl;
	  algDir = "LastRun_IMLOP/";
	  break;
  }
  case DirAlgType_DIMLOP:
  {
	  std::cout << "\nRunning D-IMLOP\n" << std::endl;
	  algDir = "LastRun_DIMLOP/";
	  break;
  }  
  case DirAlgType_GIMLOP:
  {
	  std::cout << "\nRunning GIMLOP\n" << std::endl;
	  algDir = "LastRun_GIMLOP/";
	  break;
  }
  case DirAlgType_GDIMLOP:
  {
	  std::cout << "\nRunning GDIMLOP\n" << std::endl;
	  algDir = "LastRun_GDIMLOP/";
	  break;
  }
  case DirAlgType_PIMLOP:
  {
	  std::cout << "\nRunning PIMLOP\n" << std::endl;
	  algDir = "LastRun_PIMLOP/";
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
  std::string normRVFile = workingDir + "GaussianValues.txt";	// Random Numbers: Normal RV's

  // output files
  std::string saveSourceMeshPath	= outputDir + "SaveMeshSource";
  std::string saveTargetMeshPath	= outputDir + "SaveMeshTarget";
  std::string saveSamplesPath		= outputDir + "SaveSamples";
  std::string saveSubSamplesPath	= outputDir + "SaveSubSamples.pts";
  std::string saveNoisySamplesPath	= outputDir + "SaveNoisySamples";
  std::string saveNoisySamplesPath2	= outputDir + "SaveEccNoisySamples";
  std::string saveOffsetXfmPath		= outputDir + "SaveOffsetXfm.txt";
  std::string saveRegXfmPath		= outputDir + "SaveRegXfm.txt";
  std::string saveModeWeightsPath	= outputDir + "saveModeWeights.txt";
  std::string savePath_L			= outputDir + "SaveL.pts";
  std::string savePath_L2			= outputDir + "SaveEccL.pts";
  std::string savePath_Cov			= outputDir + "SaveCov.pts";
  std::string meshDir				= outputDir + "estimatedMesh.ply";
  std::string initmeshDir			= outputDir + "initMesh.ply";

  // declare variables
  cisstMesh         mesh_source;
  cisstMesh         mesh_target;
  cisstMesh			mesh_ssm_target;
  DirPDTreeBase*      pTree;
  vctDynamicVector<vct3>    samples;
  vctDynamicVector<vct3>    sampleNorms;
  vctDynamicVector<vct3>    noisySamples;
  vctDynamicVector<vct3>    noisySampleNorms;
  vctDynamicVector<vct3>    noisySampleNorms2;
  vctDynamicVector<unsigned int>  sampleDatums;
  vctDynamicVector<vct3x3>  sampleNoiseCov;
  vctDynamicVector<vct3x3>  sampleNoiseInvCov;
  vctDynamicVector<vct3x2>  sampleNoiseL;
  vctDynamicVector<vct3x2>  sampleNoiseL2;

  vctDynamicVector<double> weight;

#if 1
  // Set default values/initialize variables
  int nSamples = 300;		// default number of samples (for default input)
  int maxIters = 100;
  double scale = 1.0;
  bool	bScale = false;

  int    nThresh = 15;			// Cov Tree Params
  double diagThresh = 15.0;		//  ''
  //int    nThresh = 5;			// Cov Tree Params
  //double diagThresh = 5.0;	//  ''

  double minOffsetPos = 10.0; // 5.0;
  double maxOffsetPos = 20.0; // 10.0;
  double minOffsetAng = 6.0;  // 20.0;
  double maxOffsetAng = 12.0; // 60.0;

  double percentOutliers	 = 0.0; // 0.05;
  double minPosOffsetOutlier = 5.0;  // 15.0;
  double maxPosOffsetOutlier = 10.0; // 20.0;
  double minAngOffsetOutlier = 0.0;  // 15.0;
  double maxAngOffsetOutlier = 0.0;  // 20.0;

  std::srand(time(NULL)); unsigned int randSeed1	= /*0;			*/std::rand();  // generates samples
  std::srand(time(NULL)); unsigned int randSeqPos1	= /*0;			*/std::rand();
  std::srand(time(NULL)); unsigned int randSeed2	= /*17; // 1; 	*/std::rand();  // generates offsets
  std::srand(time(NULL)); unsigned int randSeqPos2	= /*28; // 0;	*/std::rand();
  std::srand(time(NULL)); unsigned int randSeed3	= /*28;			*/std::rand();	// generates shape parameters
  std::srand(time(NULL)); unsigned int randSeqPos3	= /*8;			*/std::rand();

  // Samples Noise Model
  //  NOTE: this is a generative noise model (i.e. noise is generated according
  //        to the noise properties defined here)
  //double noiseSampsSD[3] = {1.0, 1.0, 1.0};   // noise model for samples (std dev along each axis)
  double sampleNoiseInPlane = 1.0;      // standard deviation of noise in and out of plane
  double sampleNoisePerpPlane = 1.0;    //   ''
  double sampleNoiseCircSDDeg = 2.0;    // noise to apply to sample orientations 
  double sampleNoiseEccentricity = 0.5; // eccentricity of orientation noise

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

  if (!cmdOpts.useDefaultNoiseInPlane)
	  sampleNoiseInPlane = (double)cmdOpts.noiseinplane;
  if (!cmdOpts.useDefaultNoisePerpPlane)
	  sampleNoisePerpPlane = (double)cmdOpts.noiseperpplane;
  if (!cmdOpts.useDefaultNoiseDeg)
	  sampleNoiseCircSDDeg = (double)cmdOpts.noisedeg;
  if (!cmdOpts.useDefaultNoiseEcc)
	  sampleNoiseEccentricity = (double)cmdOpts.noiseecc;

  double sampleNoiseCircSD = sampleNoiseCircSDDeg*cmnPI / 180.0;
  
  // load target mesh
  std::string loadTargetMeshPath;
  if (cmdOpts.useDefaultTarget)
	  loadTargetMeshPath = workingDir + "MT.ply";
  else
	  loadTargetMeshPath = cmdOpts.target;
	  //std::string loadTargetMeshPath = workingDir + "ProximalFemur.ply";
	  //std::string loadSourceMeshPath = workingDir + "RIGHTHEMIPELVIS_centered.mesh";
	  //std::string loadTargetMeshPath = workingDir + "RIGHTHEMIPELVIS_centered.mesh";
  //std::string loadMeshPath = workingDir + "RIGHTHEMIPELVIS_centered.mesh";
  //std::string loadMeshPath = workingDir + "RIGHTHEMIPELVIS.mesh";
  //std::string loadMeshPath = workingDir + "CTBreastImage_Dec20000_Shell.mesh";  
  CreateMesh(mesh_target, loadTargetMeshPath, &saveTargetMeshPath);

  mesh_ssm_target = mesh_target;	// Initializing mesh_ssm to mesh for non-deformable algorithms 
  int modes = 1;					// +1 (for mean)

  if (cmdOpts.deformable/*algType == DirAlgType_DIMLOP*/)
  {
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

	  int check = ReadShapeModel(mesh_target, loadModelPath, modes);
	  if (check != 1) {
		  std::cout << "Unsuccessful in reading model data, switching to IMLOP..." << std::endl;
		  algType = DirAlgType_IMLOP;
	  }
	  mesh_ssm_target.vertices = mesh_target.meanShape;
	  //mesh_ssm_target.SavePLY("currentMesh0.ply");

	  if (cmdOpts.useDefaultInput)
	  {
		  weight.SetSize(modes - 1);
		  GenerateRandomShapeParams(randSeed3, randSeqPos3, modes - 1, weight);
		  // add deformation to mean shape
		  for (int j = 0; j < modes - 1; j++)
			  for (int i = 0; i < mesh_target.NumVertices(); i++)
				  mesh_ssm_target.vertices[i] += weight[j] * mesh_target.wi[j][i];
		  mesh_ssm_target.SavePLY(meshDir);
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

  //mesh_target.SavePLY(initmeshDir);

  if (TargetShapeAsMesh)
  {
    // build PD tree on the mesh directly
    //  Note: defines measurement noise to be zero
    printf("Building mesh PD tree .... ");
    pTree = new DirPDTree_Mesh(mesh_target, nThresh, diagThresh); 
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
    printf("Building point cloud PD tree .... ");
    cisstPointCloud pointCloud(mesh_target);
    DirPDTree_PointCloud *pPointCloudTree;
    pPointCloudTree = new DirPDTree_PointCloud(pointCloud, nThresh, diagThresh);
    pTree = pPointCloudTree;
    //tree.RecomputeBoundingBoxesUsingExistingCovFrames();      //*** is this ever needed?
	printf("Tree built: NNodes=%d  NData=%d  TreeDepth=%d\n", pTree->NumNodes(), pTree->NumData(), pTree->TreeDepth());
  }

  // sub-sample and scale input if need be
  if (!cmdOpts.useDefaultInput)
  {
	  // load source mesh
	  std::string loadSourceMeshPath;
	  loadSourceMeshPath = cmdOpts.input;
	  CreateMesh(mesh_source, loadSourceMeshPath, &saveSourceMeshPath);

	  if (!cmdOpts.useDefaultNumSamples)
	  {
		  int nSubSamples = cmdOpts.samples;
		  vctDynamicVector<vct3> subsampledPts;
		  vctDynamicVector<bool> selectedPts;

		  subsampledPts.SetSize(nSubSamples);
		  selectedPts.SetSize(mesh_source.NumVertices());

		  selectedPts.SetAll(false);
		  subsampledPts.Zeros();

		  GenerateSubSamples(mesh_source, selectedPts, subsampledPts, nSubSamples, &saveSubSamplesPath); // change this to be more similar to generatesamples
		  //GenerateSamples(pts, randSeed1, randSeqPos1, nSubSamples,		
		  //	samples, sampleNorms, sampleDatums,
		  //	&saveSamplesPath);
	  }
	  else
		  nSamples = mesh_source.NumVertices(); //samples.size();

	  if (!cmdOpts.useDefaultScale)
	  {
		  scale = (double)cmdOpts.scale;
		  for (int i = 0; i < nSamples; i++)
			  //samples[i] = scale * samples[i];
			  mesh_source.vertices[i] = scale * mesh_source.vertices[i];
	  }
	  if (cmdOpts.bScale)
		  bScale = true;
  }
  std::ifstream randnStream(normRVFile.c_str());  // streams N(0,1) RV's

  // Generate random samples from mesh_ssm_target
  GenerateSamples(mesh_ssm_target, randSeed1, randSeqPos1, nSamples,		// mesh_ssm_target = mesh_target, if registration algorithm is not deformable
    samples, sampleNorms, sampleDatums,
	&saveSamplesPath);

  if (!cmdOpts.useDefaultInput) {
	  for (int i = 0; i < nSamples; i++)
	  {
		  samples[i] = mesh_source.vertices[i];
		  sampleNorms[i] = mesh_source.vertexNormals[i];
	  }
  }

  //noisySampleNorms.SetSize(nSamples);
  //noisySampleNorms2.SetSize(nSamples);
  // Add noise to samples
  if (!cmdOpts.useDefaultCov || !cmdOpts.useDefaultAxes)
  {
	  if (algType == DirAlgType_GIMLOP || algType == DirAlgType_GDIMLOP)
		  ReadSampleSurfaceNoise(cmdOpts.useDefaultCov, cmdOpts.useDefaultAxes,
			  randSeed1, randSeqPos1, randnStream,
			  sampleNoiseInPlane, sampleNoisePerpPlane,
			  sampleNoiseCircSDDeg*cmnPI / 180.0, sampleNoiseEccentricity,
			  samples, sampleNorms,
			  noisySamples, noisySampleNorms2,
			  sampleNoiseCov, sampleNoiseInvCov, sampleNoiseL,
			  percentOutliers,
			  minPosOffsetOutlier, maxPosOffsetOutlier,
			  minAngOffsetOutlier, maxAngOffsetOutlier,
			  &saveNoisySamplesPath,
			  cmdOpts.cov, &savePath_Cov,
			  cmdOpts.axes, &savePath_L2); 
	  else
		  ReadSampleSurfaceNoise(cmdOpts.useDefaultCov, cmdOpts.useDefaultAxes,
		  randSeed1, randSeqPos1, randnStream,
		  sampleNoiseInPlane, sampleNoisePerpPlane,
		  sampleNoiseCircSDDeg*cmnPI / 180.0, sampleNoiseEccentricity,
		  samples, sampleNorms,
		  noisySamples, noisySampleNorms,
		  sampleNoiseCov, sampleNoiseInvCov, sampleNoiseL,
		  percentOutliers,
		  minPosOffsetOutlier, maxPosOffsetOutlier,
		  minAngOffsetOutlier, maxAngOffsetOutlier,
		  &saveNoisySamplesPath,
		  cmdOpts.cov, &savePath_Cov,
		  cmdOpts.axes, &savePath_L);
  }
  else
  {
	  if (algType == DirAlgType_GIMLOP || algType == DirAlgType_GDIMLOP)
		  GenerateSampleSurfaceNoise2(randSeed1, randSeqPos1, randnStream, 
			  sampleNoiseInPlane, sampleNoisePerpPlane,
			  sampleNoiseCircSDDeg*cmnPI / 180.0, sampleNoiseEccentricity,
			  samples, sampleNorms,
			  noisySamples, noisySampleNorms2, noisySampleNorms,
			  sampleNoiseCov, sampleNoiseInvCov,
			  sampleNoiseL, sampleNoiseL2,
			  percentOutliers,
			  minPosOffsetOutlier, maxPosOffsetOutlier,
			  minAngOffsetOutlier, maxAngOffsetOutlier,
			  &saveNoisySamplesPath2, &saveNoisySamplesPath,
			  &savePath_Cov, &savePath_L2, &savePath_L);
	  else 
		  GenerateSampleSurfaceNoise(randSeed1, randSeqPos1, randnStream,
		  sampleNoiseInPlane, sampleNoisePerpPlane,
		  sampleNoiseCircSDDeg*cmnPI / 180.0, sampleNoiseEccentricity,
		  samples, sampleNorms,
		  noisySamples, noisySampleNorms,
		  sampleNoiseCov, sampleNoiseInvCov,
		  sampleNoiseL,
		  percentOutliers,
		  minPosOffsetOutlier, maxPosOffsetOutlier,
		  minAngOffsetOutlier, maxAngOffsetOutlier,
		  &saveNoisySamplesPath,
		  &savePath_Cov, &savePath_L);
  }
  
  if (!cmdOpts.useDefaultInput && cmdOpts.useDefaultNumSamples) {
	  for (int i = 0; i < nSamples; i++)
	  {
		  noisySamples[i] = mesh_source.vertices[i];
		  if (algType == DirAlgType_GIMLOP || algType == DirAlgType_GDIMLOP)
			  noisySampleNorms2[i] = mesh_source.vertexNormals[i];
		  else
			  noisySampleNorms[i] = mesh_source.vertexNormals[i];
	  }
  }

  std::cout << "Using " << nSamples << " noisy sample points...\n";

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

#else
  // Replay Randomized Trial
  vctFrm3 Fi;
  //std::string baseFolder = "..\\ICP_TestData\\RandomTests\\";
  std::string baseFolder = "//jhdfs/data/lcsr$/CIIS/sbillin3/Research/ICOP/Data/RandomTests_FemurPatch/TestD_BivariateNormalTheta_OffsetR10-20P10-20_NoOutlierDetection/";  
  std::string testDir = "Noise_Plane1Norm1_Offset_R10-20T10-20_NSamps50/";
  std::string trialName = "IMLOP";     // Note: Finit and samples are same for all trials
  std::string trialNum = "19";
  std::string commonDir = "CommonFiles/";  
  std::string loadTargetMeshFile = baseFolder + commonDir + "SaveMeshTarget.mesh";
  std::string loadSourceMeshFile = baseFolder + commonDir + "SaveMeshSource.mesh";
  std::string loadNoisySamplesFile = baseFolder + testDir + commonDir + "SaveNoisySamples_" + trialNum + ".pts";  
  // need this to get FGuess
  std::string loadTrackPathFile = baseFolder + testDir + trialName + "/SaveTrackRegPath_" + trialName + "_" + trialNum + ".txt";

  //std::string workingFolder = "..\\ICP_TestData\\LastRun_CovEst\\";
  //std::string loadMeshFile = workingFolder + "SaveMesh.mesh";
  //std::string loadSamplesFile = workingFolder + "SaveSamples" + ".pts";
  //std::string loadNoisySamplesFile = workingFolder + "SaveNoisySamples" + ".pts";
  //std::string loadTrackPathFile = workingFolder + "SaveTrackRegPath" + ".txt";

  mesh_target.LoadMeshFile( loadTargetMeshFile );
  cisstICP::MeshSave(mesh_target, saveTargetMeshPath);
  if (TargetShapeAsMesh)
  { // target shape is a mesh
    pTree = new DirPDTree_Mesh(mesh_target, nThresh, diagThresh);
  }
  else
  { // target shape is a point cloud
    pTree = new DirPDTree_PointCloud(mesh_target, nThresh, diagThresh);
  }
  printf("Tree built: NNodes=%d  NData=%d  TreeDepth=%d\n", pTree->NumNodes(), pTree->NumData(), pTree->TreeDepth());
  cisstICP::SamplesLoad( noisySamples, noisySampleNorms, loadNoisySamplesFile );
  cisstICP::SamplesSave( noisySamples, noisySampleNorms, saveNoisySamplesPath );
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

  //double meshstdDevPerpPlanePlane = 0.0;         // noise model for mesh
  //double meshStdDevInPlane = 0.0;
  //double noisePosSD[3] = {0.5, 0.5, 0.5};   // noise model for samples (std dev along each axis)
#endif

  // creating ICP solver
  cisstICP ICP;

#ifdef ENABLE_CALLBACKS
  // adding ICP callbacks
  std::vector<cisstICP::Callback> callbacks;

  //  callback: iteration file (txt)
  cisstICP::Callback iterCallback;
  iterCallback.cbFunc = Callback_SaveIterationsToFile_testICPNormals;
  std::stringstream iterFile;
  iterFile << outputDir << "SaveIterations.txt";
  std::ofstream iterFileStream(iterFile.str().c_str());
  iterCallback.userData = (void*)(&iterFileStream);
  callbacks.push_back(iterCallback);

  //  callback: iteration file (csv)
  cisstICP::Callback iterCallbackCSV;
  iterCallbackCSV.cbFunc = Callback_SaveIterationsToCSV_testICPNormals;
  std::stringstream iterCSV;
  iterCSV << outputDir << "SaveIterations.csv";
  std::ofstream iterCsvStream(iterCSV.str().c_str());
  iterCallbackCSV.userData = (void*)(&iterCsvStream);
  callbacks.push_back(iterCallbackCSV);

  //  callback: track path file
  cisstICP::Callback xfmCallback;
  xfmCallback.cbFunc = Callback_TrackRegPath_testICPNormals;
  std::stringstream trackPathFile;
  trackPathFile << outputDir << "SaveTrackRegPath.csv"; //"SaveTrackRegPath.txt";
  std::ofstream xfmFileStream(trackPathFile.str().c_str());
  xfmCallback.userData = (void*)(&xfmFileStream);
  callbacks.push_back(xfmCallback);
#endif

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
  algDirICP *pICPAlg = NULL;
  switch (algType)
  {
  case DirAlgType_StdICP:
  {
    if (TargetShapeAsMesh)
    { // target shape is a mesh
      DirPDTree_Mesh *pTreeMesh = dynamic_cast<DirPDTree_Mesh*>(pTree);
      pICPAlg = new algDirICP_StdICP_Mesh(pTreeMesh, noisySamples, noisySampleNorms);
    }
    else
    { // target shape is a point cloud
      DirPDTree_PointCloud *pTreePointCloud = dynamic_cast<DirPDTree_PointCloud*>(pTree);
      pICPAlg = new algDirICP_StdICP_PointCloud(pTreePointCloud, noisySamples, noisySampleNorms);
    }
    break;
  }
  case DirAlgType_IMLOP:
  {
    if (!TargetShapeAsMesh)
    {
      std::cout << "ERROR: Currently only mesh target supported for IMLOP" << std::endl;
      assert(0);
    }
    DirPDTree_Mesh *pTreeMesh = dynamic_cast<DirPDTree_Mesh*>(pTree);
    //algDirICP_IMLOP_Mesh *pAlg = new algDirICP_IMLOP_Mesh( pTreeMesh, &ICP );
    //pAlg->k_init = 0.0;
    //pAlg->sigma2_init = 1.0;
    //pAlg->wRpos = 0.5;
    double k = 1.0 / (sampleNoiseCircSD*sampleNoiseCircSD);
    double sigma2 = sampleNoiseInPlane*sampleNoiseInPlane;
    double wRpos = 0.5;
    double kfactor = 1.0;
    bool dynamicParamEst = true;
    algDirICP_IMLOP_Mesh *pAlg = new algDirICP_IMLOP_Mesh(
      pTreeMesh, noisySamples, noisySampleNorms,
      k, sigma2, wRpos, kfactor, dynamicParamEst);
    pICPAlg = pAlg;
    break;
  }
  case DirAlgType_DIMLOP:
  {
	  std::cout << "shape parameters:\n" << weight << std::endl;
	  if (!TargetShapeAsMesh)
	  {
		  std::cout << "ERROR: Currently only mesh target supported for DIMLOP" << std::endl;
		  assert(0);
	  }
	  DirPDTree_Mesh *pTreeMesh = dynamic_cast<DirPDTree_Mesh*>(pTree);
	  double k = 1.0 / (sampleNoiseCircSD*sampleNoiseCircSD);
	  double sigma2 = sampleNoiseInPlane*sampleNoiseInPlane;
	  double wRpos = 0.5;
	  double kfactor = 1.0;
	  bool dynamicParamEst = true;
	  algDirICP_DIMLOP *pAlg;
	  if (bScale)
		  pAlg = new algDirICP_DIMLOP(
			  pTreeMesh, noisySamples, noisySampleNorms,
			  sampleNoiseCov, sampleNoiseCov, mesh_target.meanShape,
			  k, sigma2, wRpos, kfactor, scale,
			  dynamicParamEst, bScale);
	  else
		  pAlg = new algDirICP_DIMLOP(
			  pTreeMesh, noisySamples, noisySampleNorms,
			  sampleNoiseCov, sampleNoiseCov, mesh_target.meanShape,
			  k, sigma2, wRpos, kfactor, scale,
			  dynamicParamEst);											// for cases when scale is specified, but must not be optimized over
	  pAlg->SetConstraints(shapeparambounds);

	  pICPAlg = pAlg;
	  break;
  }
  case DirAlgType_GIMLOP:
  {
    if (!TargetShapeAsMesh)
    {
      std::cout << "ERROR: Currently only mesh target supported for GIMLOP" << std::endl;
      assert(0);
    }
    DirPDTree_Mesh *pTreeMesh = dynamic_cast<DirPDTree_Mesh*>(pTree);
    // define GIMLOP parameters
    double k = 1.0 / (sampleNoiseCircSD*sampleNoiseCircSD);
    double B = sampleNoiseEccentricity*k / 2.0;
    std::cout << "k: " << k << " B: " << B << std::endl;
    vctDynamicVector<double> argK(nSamples, k);
	vctDynamicVector<double> argB(nSamples, sampleNoiseEccentricity); // B);
    //vctDynamicVector<double> argB( nSamples,0.0 );

    //double sigma2 = sampleNoiseInPlane*sampleNoiseInPlane;
    //vctDynamicVector<double> argSigma2( nSamples,sigma2 );
    //vctDynamicVector<vctFixedSizeMatrix<double,3,2>> argL( nSamples );
    //vct3 xProd, L1,L2;
    //vct3 xAxis( 1.0,0.0,0.0 );
    //vct3 yAxis( 0.0,1.0,0.0 );
    //for (unsigned int i=0; i<nSamples; i++)
    //{ // set argL as isotropic
    //  xProd = vctCrossProduct( noisySampleNorms(i),xAxis );
    //  if (xProd.Norm() < 0.01)
    //  {
    //    xProd = vctCrossProduct( noisySampleNorms(i),yAxis );
    //  }
    //  L1 = xProd.Normalized();
    //  L2 = vctCrossProduct(noisySampleNorms(i),L1).Normalized();
    //  argL(i).Column(0) = L1;
    //  argL(i).Column(1) = L2;
    //}

    // create algorithm
    algDirICP_GIMLOP_Mesh *pAlg = new algDirICP_GIMLOP_Mesh(
      pTreeMesh, noisySamples, noisySampleNorms2,
	  argK, argB, sampleNoiseL, sampleNoiseCov /*sampleNoiseInvCov*/ );
    //pAlg->SetNoiseModel(argK, argB, sampleNoiseL, sampleNoiseInvCov);
    pICPAlg = pAlg;
    break;
  }
  case DirAlgType_GDIMLOP:
  {
	  //enum  algDirICP_GDIMLOP::PARAM_EST_TYPE{ PARAMS_FIXED, PARAMS_DYN_THRESH, PARAMS_FIXED_1stITER_POS };
	  //PARAM_EST_TYPE paramEstMethod;
	  if (!TargetShapeAsMesh)
	  {
		  std::cout << "ERROR: Currently only mesh target supported for GIMLOP" << std::endl;
		  assert(0);
	  }
	  DirPDTree_Mesh *pTreeMesh = dynamic_cast<DirPDTree_Mesh*>(pTree);
	  // define GIMLOP parameters
	  double k = 1.0 / (sampleNoiseCircSD*sampleNoiseCircSD);
	  double B = sampleNoiseEccentricity*k / 2.0;
	  std::cout << "k: " << k << " B: " << B << std::endl;
	  vctDynamicVector<double> argK(nSamples, k);
	  vctDynamicVector<double> argB(nSamples, sampleNoiseEccentricity); // B);
	  //vctDynamicVector<double> argB( nSamples,0.0 );

	  //double sigma2 = sampleNoiseInPlane*sampleNoiseInPlane;
	  //vctDynamicVector<double> argSigma2( nSamples,sigma2 );
	  //vctDynamicVector<vctFixedSizeMatrix<double,3,2>> argL( nSamples );
	  //vct3 xProd, L1,L2;
	  //vct3 xAxis( 1.0,0.0,0.0 );
	  //vct3 yAxis( 0.0,1.0,0.0 );
	  //for (unsigned int i=0; i<nSamples; i++)
	  //{ // set argL as isotropic
	  //  xProd = vctCrossProduct( noisySampleNorms(i),xAxis );
	  //  if (xProd.Norm() < 0.01)
	  //  {
	  //    xProd = vctCrossProduct( noisySampleNorms(i),yAxis );
	  //  }
	  //  L1 = xProd.Normalized();
	  //  L2 = vctCrossProduct(noisySampleNorms(i),L1).Normalized();
	  //  argL(i).Column(0) = L1;
	  //  argL(i).Column(1) = L2;
	  //}

	  // create algorithm
	  algDirICP_GDIMLOP *pAlg;
	  if (bScale)
		  pAlg = new algDirICP_GDIMLOP(
			  pTreeMesh, noisySamples, noisySampleNorms2,
			  argK, argB, sampleNoiseL, sampleNoiseCov,
			  sampleNoiseCov, mesh_target.meanShape, scale, bScale);
	  else
		  pAlg = new algDirICP_GDIMLOP(
			  pTreeMesh, noisySamples, noisySampleNorms2,
			  argK, argB, sampleNoiseL, sampleNoiseCov,
			  sampleNoiseCov, mesh_target.meanShape, scale, bScale);	// for cases when scale is specified, but must not be optimized over
	  pAlg->SetConstraints(shapeparambounds);
	  //pAlg->SetNoiseModel(argK, argB, sampleNoiseL, sampleNoiseInvCov);
	  pICPAlg = pAlg;
	  break;
  }
  case DirAlgType_PIMLOP:
  {
    if (!TargetShapeAsMesh)
    {
      std::cout << "ERROR: Currently only mesh target supported for PIMLOP" << std::endl;
      assert(0);
    }
    DirPDTree_Mesh *pTreeMesh = dynamic_cast<DirPDTree_Mesh*>(pTree);
    // define noise model parameters
    vctDynamicVector<vct2> Xpln(nSamples, vct2(0.0));
    double k = 0.0;
    vct3x3 M = vct3x3::Eye();
    vctRot3 Rx_pln;
    vctDynamicVector<double> argK(nSamples, k);
    vctDynamicVector<vct3x3> argM(nSamples, M);
    vctDynamicVector<vctRot3> argRx_pln(nSamples, Rx_pln);

    // create algorithm
    algDirICP_PIMLOP_Mesh *pAlg = new algDirICP_PIMLOP_Mesh(
      pTreeMesh, 
      noisySamples, Xpln, argRx_pln, argK, argM );
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
  samplePts.vertexNormals.SetSize(noisySamples.size());
  samplePts.vertices = noisySamples;
  if (algType == DirAlgType_GIMLOP || algType == DirAlgType_GDIMLOP)
  {
	  samplePts.vertexNormals = noisySampleNorms2;
	  samplePts.SavePLY(outputDir + "/EccPts.ply");
  }
  if (!noisySampleNorms.empty())
  {
	  samplePts.vertexNormals = noisySampleNorms;
	  samplePts.SavePLY(outputDir + "/Pts.ply");
  }


  // ICP Options
  //cisstICP::OptionsNormals opt;
  cisstICP::Options opt;
  opt.auxOutputDir = outputDir;
  opt.maxIter = maxIters;
  opt.termHoldIter = 2;
  opt.numShapeParams = modes - 1;
  opt.minE = -std::numeric_limits<double>::max();
  opt.tolE = 0.0;
  opt.dPosThresh = 0.1; // 0.25;
  opt.dAngThresh = 0.1*(cmnPI / 180); // 0.25*(cmnPI / 180);
  opt.dShapeThresh = 0.1;
  opt.dPosTerm = 0.001;
  opt.dAngTerm = 0.001*(cmnPI / 180);
  opt.dShapeTerm = 0.001;
  //opt.alg = pICPAlg;
  if (cmdOpts.deformable)
	  opt.deformable = true;

  // Run ICP
  int numRuns = 1;
  vctFrm3 Freg;
  double runtime = 0.0;
  cisstICP::ReturnType rv;
  for (int i = 0; i < numRuns; i++)
  {
    rv = ICP.RunICP(pICPAlg, opt, FGuess, &callbacks);
    //rv = ICP.RunICP_Rigid( samples, sampleNorms, *pTree, opt, FGuess, Freg );
    std::cout << rv.termMsg;
#ifdef ENABLE_CALLBACKS
    iterFileStream << rv.termMsg;
#endif
    runtime += rv.runTime;
    Freg = rv.Freg;
  }
  std::cout << std::endl << " ===> Avg RunTime: " << runtime / numRuns << std::endl;

  //// save registration result
  transform_write(Freg, saveRegXfmPath);

  // Freg now includes Fi as FGuess => Freg should be identity for perfect registration
  vctFrm3 Ferr = Freg;
  vctRodRot3 Rerr(Ferr.Rotation());
  double terr = Ferr.Translation().Norm();
  double rerr = Rerr.Norm();
  vctRodRot3 Rinit(Fi.Rotation());
  double tinit = Fi.Translation().Norm();
  double rinit = Rinit.Norm();

  std::stringstream resultStream;
  resultStream << std::endl;
  if (cmdOpts.deformable/*algType == DirAlgType_DIMLOP*/)
  {
	  resultStream << "Starting Offset:   \tdAng: " << rinit * 180 / cmnPI << "\tdPos: " << tinit << "\tdShp: " << weight;
	  if (cmdOpts.bScale)
		  resultStream << "\tdScl: " << scale << std::endl;
	  else
		  resultStream << std::endl;
	  pICPAlg->ReturnScale(scale);
	  pICPAlg->ReturnShapeParam(mesh_target.Si);
	  resultStream << "Registration Error:\tdAng: " << rerr * 180 / cmnPI << "\tdPos: " << terr << "\tdShp: " << weight - mesh_target.Si;
	  if (cmdOpts.bScale)
		  resultStream << "\tdScl: " << scale << std::endl;
	  else
		  resultStream << std::endl;
  }
  else
  {
	  resultStream << "Starting Offset:   \tdAng: " << rinit * 180 / cmnPI << "\tdPos: " << tinit << std::endl;
	  resultStream << "Registration Error:\tdAng: " << rerr * 180 / cmnPI << "\tdPos: " << terr << std::endl << std::endl;
  }
  resultStream << std::endl << "Average Match Error (Mahalanobis):\t" << rv.MatchPosErrAvg << " (+/-" << rv.MatchPosErrSD << ")\n" << std::endl;
  std::cout << resultStream.str();
#ifdef ENABLE_CALLBACKS
  iterFileStream << resultStream.str();
#endif

  std::cout << "=============================================================\n" << std::endl;

  if (cmdOpts.deformable/*algType == DirAlgType_DIMLOP*/)
  {
	  mesh_target.vertices = mesh_target.meanShape;
	  for (int s = 0; s < mesh_target.NumVertices(); s++)
	  {
		  for (unsigned int i = 0; i < (unsigned int)(modes - 1); i++)
		  {
			  mesh_target.vertices(s) += (mesh_target.Si[i] * mesh_target.wi[i].Element(s));
		  }
	  }
  }
  mesh_target.SavePLY(outputDir + "/finalMesh.ply");

  for (int i = 0; i < mesh_target.NumVertices(); i++)
	  mesh_target.vertices(i) = Freg * mesh_target.vertices(i);
  mesh_target.SavePLY(outputDir + "/finalRegMesh.ply");

  //for (int i = 0; i < mesh_source.NumVertices(); i++) {
	 // mesh_source.vertices(i) = pICPAlg->matchPts[i];
	 // mesh_source.vertexNormals(i) = pICPAlg->matchNorms[i];
  //}
  //mesh_source.SavePLY(outputDir + "/matchedMesh.ply");

  samplePts.vertices.SetSize(noisySamples.size());
  samplePts.vertexNormals.SetSize(noisySamples.size());
  if (algType == DirAlgType_GIMLOP || algType == DirAlgType_GDIMLOP)
  {
	  for (int i = 0; i < noisySamples.size(); i++) {
		  samplePts.vertices[i] = Fi * noisySamples[i];
		  samplePts.vertexNormals[i] = Fi.Rotation() * noisySampleNorms2[i];
	  }
	  //samplePts.SavePLY("currentSamples0.ply");
	  samplePts.SavePLY(outputDir + "/initPts.ply");

	  for (int i = 0; i < noisySamples.size(); i++) {
		  if (bScale)
			  samplePts.vertices[i] = Freg * noisySamples[i] * scale;
		  else
			  samplePts.vertices[i] = Freg * noisySamples[i];
		  samplePts.vertexNormals[i] = Freg.Rotation() * noisySampleNorms2[i];
	  }
	  samplePts.SavePLY(outputDir + "/finalPts.ply");
  }
  else
  {
	  for (int i = 0; i < noisySamples.size(); i++) {
		  samplePts.vertices[i] = Fi * noisySamples[i];
		  samplePts.vertexNormals[i] = Fi.Rotation() * noisySampleNorms[i];
	  }
	  //samplePts.SavePLY("currentSamples0.ply");
	  samplePts.SavePLY(outputDir + "/initPts.ply");

	  for (int i = 0; i < noisySamples.size(); i++) {
		  if (bScale)
			  samplePts.vertices[i] = Freg * noisySamples[i] * scale;
		  else
			  samplePts.vertices[i] = Freg * noisySamples[i];
		  samplePts.vertexNormals[i] = Freg.Rotation() * noisySampleNorms[i];
	  }
	  samplePts.SavePLY(outputDir + "/finalPts.ply");
  }

  if (pICPAlg) delete pICPAlg;

  //std::cout << "press key then enter to quit:" << std::endl;
  //std::cin >> k;
}

#endif // _testICPNormals_H
