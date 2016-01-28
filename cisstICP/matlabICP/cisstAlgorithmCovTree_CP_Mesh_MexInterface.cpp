// ****************************************************************************
//
//  $Id: cisstAlgorithmCovTree_CP_Mesh_MexInterface.cpp, v0.1 Exp $
//
//  Description:
//	        
//
//  Usage:
//	        
//
//  See Also:
//	        
//
//  Author(s):  Seth Billings
//
//  Created on:
//
//
//              Developed by the Engineering Research Center for
//          Computer-Integrated Surgical Systems & Technology (cisstST)
//
//               Copyright(c) 2001, The Johns Hopkins University
//                          All Rights Reserved.
//
//  Experimental Software - Not certified for clinical use.
//  For use only by cisstST sponsored research projects.  For use outside of
//  cisstST, please contact Dr. Russell Taylor, cisstST Director, at rht@cs.jhu.edu
//  
// ****************************************************************************

#include "cisstAlgorithmCovTree_CP_Mesh_MexInterface.h"


// NOTES:
//
// mexGet: returns value of the specified property in the specified graphics
//         object on success. Returns NULL on failure. 
//         Do not modify the return argument from mexGet. Changing the data in 
//         a const (read-only) mxArray can produce undesired side effects.
//

// Matlab gateway function
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  //mexPrintf("I'm Running!\n");

  // Get the command string
  char cmd[64];
  if (nrhs < 1 || mxGetString(prhs[0], cmd, sizeof(cmd)))
  {
    MEX_ERROR("First input should be a command string less than 64 characters long.");
  }

  //--- Standard Functions ---//

  // Command: New
  if (!strcmp("new", cmd))
  {
    return CommandNew(nlhs, plhs, nrhs, prhs);
  }

  // All Other Comamnds:
  //  2nd input must be class instance handle
  if (nrhs < 2)
  {
    MEX_ERROR("ERROR: Second input should be a class instance handle.");
  }

  // Command: Delete
  if (!strcmp("delete", cmd))
  {
    return CommandDelete(nlhs, plhs, nrhs, prhs);
  }


  //--- ICP Methods ---//

  // Command: Initialize 
  if (!strcmp("Initialize", cmd))
  {
    return CommandInitialize(nlhs, plhs, nrhs, prhs);
  }

  // Command: Match
  if (!strcmp("ComputeMatches", cmd))
  {
    return CommandComputeMatches(nlhs, plhs, nrhs, prhs);
  }

  MEX_ERROR("ERROR: Command not recognized.");


  //// Command: Template 
  //if (!strcmp("template", cmd))
  //{
  //  return CommandInitialize(  nlhs, plhs, nrhs, prhs );
  //}
}


void CommandNew(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{ // create new class instance

  // Expected Input:
  //  cmd

  // Check parameters
  if (nlhs != 1 || nrhs != 1)
  {
    MEX_ERROR("New: requires 1 outputs, 1 inputs.");
  }
  // Return a handle to a new C++ instance
  plhs[0] = convertPtr2Mat < cisstAlgorithmCovTree_CP_Mesh_MexInterface >
    (new cisstAlgorithmCovTree_CP_Mesh_MexInterface);
}


void CommandDelete(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{ // delete class instance

  // Expected Input:
  //  cmd
  //  class handle

  // Warn if other commands were ignored
  if (nlhs != 0 || nrhs != 2)
  {
    MEX_WARNING("Delete: requires 0 outputs, 2 inputs.");
  }

  // Destroy the C++ object
  destroyObject<cisstAlgorithmCovTree_CP_Mesh_MexInterface>(prhs[1]);
}

void CommandInitialize(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
  // Expected Input:
  //  cmd
  //  class handle
  //  mesh
  //    V ~ 3D vertex positions     (3 x Nv double)
  //    T ~ triangle vertex indices (3 x Nt integer)
  //    N ~ triangle normals        (3 x Nt double)
  //  

  std::stringstream ss;

  // Get the class instance from the second input
  cisstAlgorithmCovTree_CP_Mesh_MexInterface& obj =
    *convertMat2Ptr<cisstAlgorithmCovTree_CP_Mesh_MexInterface>(prhs[1]);

  // Check parameters
  if (nlhs != 0 || nrhs != 5)
  {
    MEX_ERROR("Cmd Initialize: requires 0 outputs, 5 inputs.");
  }


  //--- Inputs ---//

  MEX_DEBUG("Parsing Inputs...\n");

  // Parse Mesh  
  //vctDynamicVector<vct3> V;
  //vctDynamicVector<vctInt3> T;
  //vctDynamicVector<vct3> Tn;
  MEX_DEBUG(" ...V\n");
  Parse3DVectorArray_Double(prhs[2], obj.V);
  MEX_DEBUG(" ...T\n");
  Parse3DVectorArray_Int(prhs[3], obj.T);
  MEX_DEBUG(" ...N\n");
  Parse3DVectorArray_Double(prhs[4], obj.Tn);

  MEX_DEBUG(" ...done\n");


  //--- Processing ---//

  // build mesh
  MEX_DEBUG("Building mesh...\n");
  cisstMesh *pMesh = new cisstMesh();
  pMesh->LoadMesh(obj.V, obj.T, obj.Tn);
  if (pMesh->NumVertices() == 0)
  {
    MEX_ERROR("ERROR: Build mesh failed\n");
  }
  // build covariance tree
  int    nThresh = 5;       // Cov Tree Params
  double diagThresh = 5.0;  //  ''
  MEX_DEBUG("Building mesh covariance tree...\n");
  obj.pTree = new cisstCovTree_Mesh(*pMesh, nThresh, diagThresh);
  //cisstCovTree_Mesh *pCovTree = new cisstCovTree_Mesh(*pMesh, nThresh, diagThresh);  
  ss.str("");
  ss << "Tree built:  NNodes = " << obj.pTree->NumNodes() << "  NDatums = "
    << obj.pTree->NumData() << "  TreeDepth = " << obj.pTree->TreeDepth() << "\n";
  MEX_DEBUG(ss.str());

  // create & initialize algorithm
  if (obj.pAlg)
  {
    delete obj.pAlg;
    obj.pAlg = NULL;
  }
  obj.pAlg = new cisstAlgorithmCovTree_CP_Mesh(obj.pTree);

}

void CommandComputeMatches(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
  // Expected Input:
  //  cmd
  //  class handle
  //  samplePtsXfmd       ~ sample points          (3 x nSamps double)
  //  matchDatumsInit     ~ initial match datums   (1 x nSamps int32)   [optional]
  // 
  // Output:
  //  MatchPts3D      ~ 3D match points (in camera coords)    (3 x nSamps double)
  //  MatchNorms3D    ~ 3D match norms  (in camera coords)    (3 x nSamps double)  [optional]
  //  matchDatums     ~ final match datums                    (1 x nSamps int32)   [optional]
  // 

  // Check parameters
  if ((nlhs < 1 || nlhs > 3) || (nrhs < 3 || nrhs > 4))
  {
    MEX_ERROR("ComputeMatches: requires 1 or 3 outputs, 3 or 4 inputs.");
  }

  // Get the class instance from the second input
  cisstAlgorithmCovTree_CP_Mesh_MexInterface& obj =
    *convertMat2Ptr<cisstAlgorithmCovTree_CP_Mesh_MexInterface>(prhs[1]);

  cisstAlgorithmCovTree_CP_Mesh& alg = *obj.pAlg;


  //--- Inputs ---//

  MEX_DEBUG("Parsing Inputs...\n");

  // Parse Samples
  MEX_DEBUG(" ...samplePts\n");
  Parse3DVectorArray_Double(prhs[2], obj.samplePtsXfmd);

  // Parse Match Datum Initializers
  bool bDatumsInitialized = false;
  if (nrhs > 3)
  {    
    MEX_DEBUG(" ...matchDatumsInit\n");
    ParseVector_Int(prhs[3], obj.matchDatums);
    bDatumsInitialized = true;
  }

  MEX_DEBUG(" ...done\n");


  //--- Computation ---//

  obj.nSamples = obj.samplePtsXfmd.size();
  obj.matchPts.SetSize(obj.nSamples);
  obj.matchNorms.SetSize(obj.nSamples);
  obj.matchDatums.SetSize(obj.nSamples);
  obj.matchErrors.SetSize(obj.nSamples);

  if (!bDatumsInitialized)
  {
    // initialize match datums with accelerated approximate search
    for (unsigned int i = 0; i < obj.nSamples; i++)
    {
      obj.matchDatums[i] = obj.pTree->FastInitializeProximalDatum(
        obj.samplePtsXfmd[i], obj.matchPts[i]);
    }
  }

  // compute matches
  obj.ComputeMatches();


  //--- Outputs ---//

  unsigned int nSamps = obj.nSamples;

  // outputs 1 & 2: match points and normals
  plhs[0] = mxCreateDoubleMatrix(3, nSamps, mxREAL);
  plhs[1] = mxCreateDoubleMatrix(3, nSamps, mxREAL);
  double *matchPts = mxGetPr(plhs[0]);
  double *matchNorms = mxGetPr(plhs[1]);
  for (unsigned int i = 0; i < nSamps; i++)
  {
    (*matchPts) = obj.matchPts[i].Element(0);
    matchPts++;
    (*matchPts) = obj.matchPts[i].Element(1);
    matchPts++;
    (*matchPts) = obj.matchPts[i].Element(2);
    matchPts++;

    (*matchNorms) = obj.matchNorms[i].Element(0);
    matchNorms++;
    (*matchNorms) = obj.matchNorms[i].Element(1);
    matchNorms++;
    (*matchNorms) = obj.matchNorms[i].Element(2);
    matchNorms++;
  }

  // output 3: match datums
  if (nlhs >= 3)
  {
    plhs[2] = mxCreateNumericMatrix(nSamps, 1, mxINT32_CLASS, mxREAL);
    int *matchDatums = static_cast<int*> (mxGetData(plhs[2]));
    for (unsigned int i = 0; i < nSamps; i++)
    {
      (*matchDatums) = obj.matchDatums[i];
      matchDatums++;
    }
  }

}

//void CommandTemplate( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
//{
//  // Expected Input:
//  //  cmd
//  //  class handle
//
//  // Get the class instance pointer from the second input
//  cisstAlgorithmCovTree_CP_Mesh* pAlg = 
//    convertMat2Ptr<cisstAlgorithmCovTree_CP_Mesh>( prhs[1] );
//}