#ifndef _cisstAlgorithmDirICP_vMFG_Mesh_MexInterface_h
#define _cisstAlgorithmDirICP_vMFG_Mesh_MexInterface_h

#include "matlabClassHandle.h"
#include "mex.h"
#include "matlabParser.h"
#include "matlabExtras.h"

#include "cisstAlgorithmDirICP_vMFG_Mesh.h"


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
void CommandSetSamples(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
void CommandICP_InitializeParameters(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
void CommandICP_ComputeMatches(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
void CommandICP_RegisterMatches(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
void CommandICP_EvaluateErrorFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
void CommandICP_Terminate(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
//void CommandTest(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
//void CommandTest2(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);

class cisstAlgorithmDirICP_vMFG_Mesh_MexInterface
{

  //--- Parameters ---//

public:

  cisstAlgorithmDirICP_vMFG_Mesh *pAlg;


  //--- Methods ---//

public:

  // constructor
	cisstAlgorithmDirICP_vMFG_Mesh_MexInterface() 
    : pAlg(NULL)
  {}

	~cisstAlgorithmDirICP_vMFG_Mesh_MexInterface()
  {
    if (pAlg)
    {
      delete pAlg;
      pAlg = NULL;
    }
  }

};


#endif