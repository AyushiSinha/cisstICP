// ****************************************************************************
//
//    Copyright (c) 2014, Seth Billings, Russell Taylor, Johns Hopkins University
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
#include "alg2D_DirPDTree_vonMises_Edges.h"

#include <limits>

void alg2D_DirPDTree_vonMises_Edges::ICP_InitializeParameters(vctFrm2 &FGuess)
{
	// set starting sample positions
	//UpdateSampleXfmPositions(FGuess);

	// initialize matches with accelerated approximate search
	for (unsigned int i = 0; i < nSamples; i++)
	{
		//FIX THIS!!!
		//matchDatums[i] = pDirTree->FastInitializeProximalDatum(
		//	samplePtsXfmd[i], sampleNormsXfmd[i],
		//	matchPts[i], matchNorms[i]);
	}

	nOutliers = 0;
	matchPermitted.SetAll(false);
}

void alg2D_DirPDTree_vonMises_Edges::ICP_ComputeMatches()
{
	// Find the point on the model having lowest match error
	//  for each sample point

#ifdef ValidatePDTreeSearch  
	numInvalidDatums = 0;
	numValidDatums = 0;
#endif

	unsigned int nodesSearched = 0;
	minNodesSearched = std::numeric_limits<unsigned int>::max();
	maxNodesSearched = std::numeric_limits<unsigned int>::min();
	avgNodesSearched = 0;

	for (unsigned int s = 0; s < nSamples; s++)
	{
		bPermittedMatchFound = false;

		// set sigma and k per point?
		
		/// FIX THIS!
		//matchDatums.Element(s) = pDirTree->FindClosestDatum(
		//	samplePtsXfmd.Element(s), sampleNormsXfmd.Element(s),
		//	matchPts.Element(s), matchNorms.Element(s),
		//	matchDatums.Element(s),
		//	matchErrors.Element(s),
		//	nodesSearched);

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
				validError = FindClosestPointOnDatum(
					samplePtsXfmd.Element(s), sampleNormsXfmd.Element(s),
					tmp1, tmp2, validDatum);
				validFS << "Correspondence Search Data:  (foundMatch / validMatch)" << std::endl
					<< " MatchError = " << matchErrors.Element(s) << "/" << validError << std::endl
					//<< " dPos = " << dist_PostMatch.Element(s) << "/" << validDist << std::endl
					//<< " dAng = " << dTheta*180.0 / cmnPI << "/" << validAng*180.0 / cmnPI << std::endl
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

				DirPDTree2DNode *termNode = 0;
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

		matchPermitted.Element(s) = bPermittedMatchFound;
	}

	avgNodesSearched /= nSamples;

#ifdef ValidatePDTreeSearch  
	validPercent = (double)numValidDatums / (double)nSamples * 100.0;
	validFS << "iter " << validIter << ":  NumMatches(valid/invalid): "
		<< numValidDatums << "/" << numInvalidDatums << "  valid% = "
		<< validPercent << std::endl;
	validIter++;
#endif
}

double alg2D_DirPDTree_vonMises_Edges::FindClosestPointOnDatum(
  const vct2 &v, const vct2 &n,
  vct2 &closest, vct2 &closestNorm,
  int datum)
{
  //  cost:  k*(1-N'*Nclosest) + ||v - closest||^2 / (2*sigma2)

  // Project point onto the edge
  //  (storing lambda value in temp buffer)
  closest = pDirTree->GetEdge(datum).ProjectOnEdge(v); //, &searchLambdas.Element(datum));
  closestNorm = pDirTree->GetEdge(datum).Norm;

  // Is this a permitted match?
  double dTheta = acos(vctDotProduct(n, closestNorm));
  if (dTheta > dThetaMax)
  {
    return std::numeric_limits<double>::max();
  }
  else
  {
    bPermittedMatchFound = true;
    return k*(1.0 - n*closestNorm) + (v - closest).NormSquare() / (2.0*sigma2);
  }

  // This doesn't work, i.e. we can't keep returning the best match among the non-permitted
  //  matches because then when a permitted match finally comes along, then it may have
  //  higher error, which would prevent it from being chosen. The only way to accomplish the
  //  above is to modify the core PD tree search routine, which I don't want to do.
  //  => only return match errors for permitted matches.
  //// is this a permitted match?
  //double matchError = k*(1.0 - n*closestNorm) + (v - closest).NormSquare() / (2.0*sigma2);
  //double dTheta = acos(vctDotProduct(n, closestNorm));
  //if (dTheta > dThetaMax)
  //{
  //  if (bPermittedMatchFound)
  //  { // skip this match as long as some other permitted match has been already been found
  //    // do this by returning an astronomical match error
  //    matchError = std::numeric_limits<double>::max();
  //  }
  //}
  //else
  //{
  //  bPermittedMatchFound = true;
  //}
  //return matchError;
}

int alg2D_DirPDTree_vonMises_Edges::DatumMightBeCloser(
  const vct2 &v, const vct2 &n,
  int datum,
  double ErrorBound)
{
  // doing a decent proximity check is complicated enough that it is
  //  better to just compute the full error directly
  return 1;
}
