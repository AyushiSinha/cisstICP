#ifndef _mexInterface_AlgPDTree_MLP_Mesh_h
#define _mexInterface_AlgPDTree_MLP_Mesh_h

#include "matlabClassHandle.h"
#include "mex.h"
#include "matlabParser.h"
#include "matlabExtras.h"
#include "utilities.h"

#include "algPDTree_MLP_Mesh.h"


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

class mexInterface_AlgPDTree_MLP_Mesh
{

  //--- Parameters ---//

public:

  PDTree_Mesh             *pTree;
  algPDTree_MLP_Mesh *pAlg;

  // covariance tree inputs
  vctDynamicVector<vct3>      V;    // mesh vertices
  vctDynamicVector<vctInt3>   T;    // mesh triangles
  vctDynamicVector<vct3>      Tn;   // mesh triangle normals

  vctDynamicVector<vct3x3>  TriangleCov;        // mesh triangle covariances
  vctDynamicVector<vct3>    TriangleCovEig;     // mesh triangle covariance eigenvalues (in decreasing order)


  // per sample-set I/O
  unsigned int nSamples;

  vctDynamicVector<vct3>    samplePtsXfmd;
  vctDynamicVector<vct3x3>  sampleCovXfmd;

  vctDynamicVector<vct3>  matchPts;
  vctDynamicVector<vct3>  matchNorms;
  vctDynamicVector<int>   matchDatums;
  vctDoubleVec            matchErrors;

  double minNodesSearched, maxNodesSearched, avgNodesSearched;



  //--- Methods ---//

public:

  // constructor
	mexInterface_AlgPDTree_MLP_Mesh() 
    : pAlg(NULL)
  {}

	~mexInterface_AlgPDTree_MLP_Mesh()
  {
    if (pAlg)
    {
      delete pAlg;
      pAlg = NULL;
    }
  }


  void ComputeMatches()
  {
    // Find the point on the model having lowest match error
    //  for each sample point

#ifdef ValidatePDTreeSearch  
    numInvalidDatums = 0;
    numValidDatums = 0;
#endif

    unsigned int nSamples = samplePtsXfmd.size();

    unsigned int nodesSearched = 0;
    minNodesSearched = std::numeric_limits<unsigned int>::max();
    maxNodesSearched = std::numeric_limits<unsigned int>::min();
    avgNodesSearched = 0;

    vct3x3 dummyMat;
    vct3 sampleCovEig;
    for (unsigned int s = 0; s < nSamples; s++)
    {
      // inform algorithm beginning new match
      //SamplePreMatch(s);

      ComputeCovEigenDecomposition_NonIter(sampleCovXfmd[s], sampleCovEig, dummyMat);
      pAlg->InitializeSampleSearch(sampleCovXfmd[s], sampleCovEig);

      // Find best match for this sample
      matchDatums.Element(s) = pTree->FindClosestDatum(
        samplePtsXfmd.Element(s), matchPts.Element(s),
        matchDatums.Element(s),
        matchErrors.Element(s),
        nodesSearched);

      avgNodesSearched += nodesSearched;
      minNodesSearched = (nodesSearched < minNodesSearched) ? nodesSearched : minNodesSearched;
      maxNodesSearched = (nodesSearched > maxNodesSearched) ? nodesSearched : maxNodesSearched;

#ifdef ValidatePDTreeSearch        
      validDatum = pDirTree->ValidateClosestDatum(
        samplePtsXfmd.Element(s), sampleNormsXfmd.Element(s),
        validPoint, validNorm);
      if (validDatum != matchDatums.Element(s))
      {
        // It is possible to have different datums for same point if the match
        //  lies on a datum edge; if this is the case, then the search didn't
        //  actually fail
        // Note: cannot compare two double values for exact equality due to
        //       inexact representation of decimal values in binary arithmetic
        searchError = (validPoint - matchPts.Element(s)).NormSquare();
        if (searchError > doubleEps)
        {
          numInvalidDatums++;
          //matchDist = (matchPts.Element(s)-samplePtsXfmd.Element(s)).Norm();
          //matchAng  = acos( vctDotProduct(matchNorms.Element(s),sampleNormsXfmd.Element(s)) );
          validDist = (validPoint - samplePtsXfmd.Element(s)).Norm();
          validAng = acos(vctDotProduct(validNorm, sampleNormsXfmd.Element(s)));
          vct2 tmp1, tmp2;
          //searchError = algorithm->FindClosestPointOnDatum(
          //                    samplePtsXfmd.Element(s), sampleNormsXfmd.Element(s),
          //                    tmp1, tmp2, matchDatums.Element(s));
          validError = pDirTree->pAlgorithm->FindClosestPointOnDatum(
            samplePtsXfmd.Element(s), sampleNormsXfmd.Element(s),
            tmp1, tmp2, validDatum);
          validFS << "Correspondence Search Data:  (foundMatch / validMatch)" << std::endl
            << " MatchError = " << matchErrors.Element(s) << "/" << validError << std::endl
            << " dPos = " << dist_PostMatch.Element(s) << "/" << validDist << std::endl
            << " dAng = " << dTheta*180.0 / cmnPI << "/" << validAng*180.0 / cmnPI << std::endl
            << " XfmSamplePoint = [" << samplePtsXfmd.Element(s) << "]" << std::endl
            << " XfmSampleNorm  = [" << sampleNormsXfmd.Element(s) << "]" << std::endl
            << " MatchPoint =     [" << matchPts.Element(s) << "]" << std::endl
            << " MatchNorm =      [" << matchNorms.Element(s) << "]" << std::endl
            << " MatchDatum = " << matchDatums.Element(s) << std::endl
            << " ValidPoint =     [" << validPoint << "]" << std::endl
            << " ValidNorm =      [" << validNorm << "]" << std::endl
            << " ValidDatum = " << validDatum << std::endl
            << " SampleIndex = " << s << std::endl;

          //validFS << "Fact = [" << std::endl << Fact << "]" << std::endl;

          DirPDTreeNode *termNode = 0;
          pDirTree->FindTerminalNode(validDatum, &termNode);
          if (!termNode)
          {
            std::cout << "ERROR: did not find terminal node for datum: " << validDatum << std::endl;
            assert(0);
          }
          validFS << " Valid Terminal Node:" << std::endl;
          validFS << "   MinCorner: " << termNode->Bounds.MinCorner << std::endl;
          validFS << "   MaxCorner: " << termNode->Bounds.MaxCorner << std::endl;
          validFS << "   NData: " << termNode->NData << std::endl;
          validFS << "Fnode_valid = [" << std::endl << termNode->F << "]" << std::endl;
        }
        else
        {
          numValidDatums++;
        }
      }
      else
      {
        numValidDatums++;
      }
#endif
      
      // save lambda value of match for computing 3D match equivalent
      matchNorms.Element(s) = pTree->Triangle(matchDatums.Element(s)).norm;
      //SamplePostMatch(s);
    }

    avgNodesSearched /= nSamples;

#ifdef ValidatePDTreeSearch  
    validPercent = (double)numValidDatums / (double)nSamples;
    validFS << "iter " << validIter << ":  NumMatches(valid/invalid): "
      << numValidDatums << "/" << numInvalidDatums << "  valid% = "
      << validPercent << std::endl;
    validIter++;
#endif

  }

};


#endif