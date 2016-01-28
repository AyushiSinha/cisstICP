#ifndef _cisstAlgorithm2D_DirCovTree_CP_Edges_MexInterface_h
#define _cisstAlgorithm2D_DirCovTree_CP_Edges_MexInterface_h

#include "matlabClassHandle.h"
#include "mex.h"
#include "matlabParser.h"
#include "matlabExtras.h"

#include "cisstAlgorithm2D_DirCovTree_CP_Edges.h"


// debug
//#define DEBUG_MEX

#ifdef DEBUG_MEX
#define MEX_DEBUG(x) MEX_PRINT((std::string("MEX_DEBUG: ") + x + "\n").c_str())
#else
#define MEX_DEBUG(x)
#endif;


// Matlab gateway function
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);

// Command Interface
void CommandNew(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
void CommandDelete(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
void CommandInitialize(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
void CommandComputeMatches(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);

class cisstAlgorithm2D_DirCovTree_CP_Edges_MexInterface
{

  //--- Parameters ---//

public:

  cisstDirCovTree2D_Edges *pDirTree;
  cisstAlgorithm2D_DirCovTree_CP_Edges *pAlg;

  // covariance tree inputs
  vctDynamicVector<vct2> edgesV1;
  vctDynamicVector<vct2> edgesV2;
  vctDynamicVector<vct2> edgesNorm;

  vctDynamicVector<vct3>   edgesV1_3D;
  vctDynamicVector<vct3>   edgesV2_3D;
  vctDynamicVector<vct3>   edgesNorm_3D;


  // per sample-set I/O
  unsigned int nSamples;

  vctDynamicVector<vct2> samplePtsXfmd;
  vctDynamicVector<vct2> sampleNormsXfmd;

  vctDynamicVector<vct2> matchPts;
  vctDynamicVector<vct2> matchNorms;
  vctDynamicVector<int>  matchDatums;
  vctDoubleVec  matchErrors;
  vctDoubleVec  matchLambdas;
  
  vctDynamicVector<vct3>   MatchPts_3D;
  //vctDynamicVector<vct3>   MatchNorms_3D;

  double minNodesSearched, maxNodesSearched, avgNodesSearched;



  //--- Methods ---//

public:

  // constructor
	cisstAlgorithm2D_DirCovTree_CP_Edges_MexInterface() 
    : pAlg(NULL)
  {}

	~cisstAlgorithm2D_DirCovTree_CP_Edges_MexInterface()
  {
    if (pAlg)
    {
      delete pAlg;
      pAlg = NULL;
    }
  }


  void ConvertMatchesTo3D()
  {
    vct3 V1, V2;
    for (unsigned int i = 0; i < nSamples; i++)
    {
      V1 = edgesV1_3D.Element(matchDatums[i]);
      V2 = edgesV2_3D.Element(matchDatums[i]);
      MatchPts_3D[i] = (V2 - V1)*matchLambdas[i] + V1;
      //MatchNorms_3D[i] = edgesNorm_3D.Element(matchDatums[i]);
    }
  }

  void ComputeMatches();

};


#endif