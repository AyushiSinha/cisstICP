// ****************************************************************************
//
//    Copyright (c) 2014, Seth Billings, Ayushi Sinha, Russell Taylor, 
//    Johns Hopkins University. 
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


#include <stdio.h>
#include <sstream>
#include <limits>

#include <cisstOSAbstraction.h>

#include "cisstMesh.h"
#include "cisstICP.h"
#include "algICP.h"

// debug
//#define ENABLE_CODE_TRACE
//#define ENABLE_CODE_PROFILER

cisstICP::ReturnType cisstICP::RunICP(
  algICP *pAlg,
  const Options &opt,
  const vctFrm3 &FGuess,
  std::vector<Callback> *pUserCallbacks,
  bool bEnableAlgorithmCallbacks)
{
  if (pAlg == NULL)
  {
    std::cout << "ERROR: no registration algorithm specified" << std::endl;
    assert(0);
    return ReturnType();
  }
  this->pAlgorithm = pAlg;
  this->opt = opt;
  this->FGuess = FGuess;

  // setup iteration callbacks
  ClearIterationCallbacks();
  if (pUserCallbacks)
  {
    // add user callbacks
    AddIterationCallbacks(*pUserCallbacks);
  }
  if (bEnableAlgorithmCallbacks)
  {
    // add algorithm callbacks
    // NOTE: for some reason, linux requires callbacks to be stored
    //       to a local variable before use as a function argument
    std::vector<Callback> callbacks = pAlg->ICP_GetIterationCallbacks();
    AddIterationCallbacks(callbacks);
  }

  // begin registration
  return IterateICP();
}


cisstICP::ReturnType cisstICP::IterateICP()
{
  //bool JustDidAccelStep = false;
  std::stringstream termMsg;
  double dAng, dPos;
  double dAng01, dAng12 = 0.0;
  double dPos01, dPos12 = 0.0;
  vctFrm3 dF;
  double dS;
  vctFrm3 F01, F12, F02;
  vctFrm3 Freg0, Freg1, Freg2;
  vctRodRot3 dR;
  int numModes;
  double E0, E1, E2;
  double tolE;
  ReturnType rt;
  unsigned int terminateIter = 0;  // consecutive iterations satisfying termination
  vctDynamicVector<double> sp(opt.numShapeParams);
  double scale;

#ifdef ENABLE_CODE_PROFILER
  osaStopwatch codeProfiler;
  double time_Callbacks = 0.0;
  double time_Extras = 0.0;
  double time_Match = 0.0;
  double time_UpdateParams_PostMatch = 0.0;
  double time_FilterMatches = 0.0;
  double time_EvalErrorFunc = 0.0;
  double time_Register = 0.0;
  double time_UpdateParams_PostRegister = 0.0;
  codeProfiler.Reset();
  codeProfiler.Start();
#endif 

  if (opt.printOutput)
  {
    std::cout << "\n===================== Beginning Registration ==================\n";
  }

  osaStopwatch totalTimer;
  osaStopwatch iterTimer;
  totalTimer.Reset();
  totalTimer.Start();
  iterTimer.Reset();
  iterTimer.Start();

  //--- ICP Initialize ---//

  // initialize algorithm
  Freg0 = Freg1 = vctFrm3::Identity();
  dF = Freg2 = Freg = FGuess;
  sp.SetAll(0.0);
  ShapeNorm = 0.0;
  prevShapeNorm = 0.0;
#ifdef ENABLE_CODE_TRACE
  std::cout << "InitializeParameters()" << std::endl;
#endif
  pAlgorithm->ICP_InitializeParameters(FGuess);

#ifdef ENABLE_CODE_PROFILER
  time_Extras = codeProfiler.GetElapsedTime();
  codeProfiler.Reset();
  codeProfiler.Start();
#endif

  //------------ ICP Iterate ----------------//

  // compute matches
  unsigned int iter;
  for (iter = 1; iter <= opt.maxIter; iter++)
  {

#ifdef ENABLE_CODE_TRACE
    std::cout << "ComputeMatches()" << std::endl;
#endif
	pAlgorithm->ICP_ComputeMatches();

#ifdef ENABLE_CODE_PROFILER
    time_Match = codeProfiler.GetElapsedTime();
    codeProfiler.Reset();
    codeProfiler.Start();
#endif

#ifdef ENABLE_CODE_TRACE
    std::cout << "UpdateParameters_PostMatch()" << std::endl;
#endif
	pAlgorithm->ICP_UpdateParameters_PostMatch();

#ifdef ENABLE_CODE_PROFILER
    time_UpdateParams_PostMatch = codeProfiler.GetElapsedTime();
    codeProfiler.Reset();
    codeProfiler.Start();
#endif

#ifdef ENABLE_CODE_TRACE
    std::cout << "FilterMatches()" << std::endl;
#endif
	nOutliers = pAlgorithm->ICP_FilterMatches();

#ifdef ENABLE_CODE_PROFILER
    time_FilterMatches = codeProfiler.GetElapsedTime();
    codeProfiler.Reset();
    codeProfiler.Start();
#endif

    //--- First Iteration: report initial match statistics ---//
    if (iter == 1)
    {
#ifdef ENABLE_CODE_TRACE
      std::cout << "EvaluateErrorFunction()" << std::endl;
#endif
	  E = pAlgorithm->ICP_EvaluateErrorFunction();
      E0 = E1 = std::numeric_limits<double>::max();
      E2 = E;
      tolE = 0.0;
      Ebest = E;
      iterBest = 0;
      Fbest = FGuess;

#ifdef ENABLE_CODE_PROFILER
      time_EvalErrorFunc = codeProfiler.GetElapsedTime();
      codeProfiler.Reset();
      codeProfiler.Start();
#endif

#ifdef ENABLE_CODE_TRACE
      std::cout << "Callbacks" << std::endl;
#endif

	  pAlgorithm->ReturnScale(scale);
	  pAlgorithm->ReturnShapeParam(sp);

      // initial callback
      iterData.iter = 0;
      iterData.E = E;
      iterData.tolE = tolE;
      iterData.Freg.Assign(FGuess);
	  iterData.dF.Assign(FGuess);
	  iterData.scale = scale;
	  iterData.S = sp;
      iterData.time = iterTimer.GetElapsedTime();
      iterData.nOutliers = nOutliers;
      //iterData.isAccelStep = false;
      std::vector<Callback>::iterator cbIter;
      for (cbIter = this->iterationCallbacks.begin(); cbIter != this->iterationCallbacks.end(); cbIter++)
      {
        cbIter->cbFunc(iterData, cbIter->userData);
      }      
      iterTimer.Reset();
      iterTimer.Start();

      rt.runTimeFirstMatch = iterData.time;

#ifdef ENABLE_CODE_PROFILER
      time_Callbacks = codeProfiler.GetElapsedTime();
      codeProfiler.Reset();
      codeProfiler.Start();
#endif

#ifdef ENABLE_CODE_PROFILER
      std::cout
        << " time_Match:              " << time_Match << std::endl
        << " time_FilterMatches:      " << time_FilterMatches << std::endl
        << " time_UpdateParams_PostMatch: " << time_UpdateParams_PostMatch << std::endl
        << " time_EvalErrorFunc:      " << time_EvalErrorFunc << std::endl
        << " time_Callbacks:          " << time_Callbacks << std::endl
        << " time_Extras:             " << time_Extras << std::endl;
      codeProfiler.Reset();
      codeProfiler.Start();
#endif
    }

#ifdef ENABLE_CODE_TRACE
    std::cout << "RegisterMatches()" << std::endl;
#endif
	Freg = pAlgorithm->ICP_RegisterMatches();
    Freg0 = Freg1;
    Freg1 = Freg2;
    Freg2 = Freg;

	pAlgorithm->ReturnScale(scale);
	pAlgorithm->ReturnShapeParam(sp); 

	prevShapeNorm = ShapeNorm;
	ShapeNorm = sp.Norm();

#ifdef ENABLE_CODE_TRACE
    std::cout << "F:" << std::endl << Freg << std::endl;
#endif

    // dF = xfm from Freg1 to Freg2
    // first go back along Freg1 then go forward along Freg2
	dF = Freg2 * Freg1.Inverse();
	dS = abs(prevShapeNorm - ShapeNorm);

#ifdef ENABLE_CODE_PROFILER
    time_Register = codeProfiler.GetElapsedTime();
    codeProfiler.Reset();
    codeProfiler.Start();
#endif

#ifdef ENABLE_CODE_TRACE
    std::cout << "UpdateParameters_PostRegister()" << std::endl;
#endif

    // update algorithm's post-registration step parameters
    pAlgorithm->ICP_UpdateParameters_PostRegister(Freg);

#ifdef ENABLE_CODE_PROFILER
    time_UpdateParams_PostRegister = codeProfiler.GetElapsedTime();
    codeProfiler.Reset();
    codeProfiler.Start();
#endif

#ifdef ENABLE_CODE_TRACE
    std::cout << "EvaluateErrorFunction()" << std::endl;
#endif

    // On a rare iteration (typically one that involves adding a new outlier) the cost function
    //  has been seen to slightly increase, despite using the outlier threshold error correction.
    // Note: One idea is this may occur when the avg error over the entire set increases,
    //       thereby increasing the outlier threshold (and hence the outlier error contribution)
    //       for all outliers in the set. This has not been confirmed, and the true source is
    //       yet unknown.
    // Note: On the rare iterations when the rms error increases, the weighted point match
    //       registration function still reduces the rms error for the points being matched.
    //       So the increased error has something to do with how the outliers are handled.
    //if (E2 > E1)
    //{ std::cout << "  ---> Runtime Warning: cost function increased!" << std::endl; }

    // compute error function value
    E = pAlgorithm->ICP_EvaluateErrorFunction();
    E0 = E1;
    E1 = E2;
    E2 = E;
    tolE = fabs((E2 - E1) / E1);
    if (E <= Ebest)
    {
      Ebest = E;
      Fbest = Freg;
      iterBest = iter;
    }

#ifdef ENABLE_CODE_PROFILER
    time_EvalErrorFunc = codeProfiler.GetElapsedTime();
    codeProfiler.Reset();
    codeProfiler.Start();
#endif

#ifdef ENABLE_CODE_TRACE
    std::cout << "Callbacks" << std::endl;
#endif

    //-- Callbacks --//
    iterData.iter = iter;
    iterData.E = E;
    iterData.tolE = tolE;
    iterData.Freg.Assign(Freg);
    iterData.dF.Assign(dF);
	iterData.scale = scale;
	iterData.S = sp;
    iterData.time = iterTimer.GetElapsedTime();
    iterData.nOutliers = nOutliers;
    //iterData.isAccelStep = JustDidAccelStep;
    std::vector<Callback>::iterator cbIter;
    for (cbIter = this->iterationCallbacks.begin(); cbIter != this->iterationCallbacks.end(); cbIter++)
    {
      cbIter->cbFunc(iterData, cbIter->userData);
    }
    iterTimer.Reset();
    iterTimer.Start();

#ifdef ENABLE_CODE_PROFILER
    time_Callbacks = codeProfiler.GetElapsedTime();
    codeProfiler.Reset();
    codeProfiler.Start();
#endif

#ifdef ENABLE_CODE_PROFILER
    std::cout
      << " time_Register:                  " << time_Register << std::endl
      << " time_UpdateParams_PostRegister: " << time_UpdateParams_PostRegister << std::endl
      << " time_Match:                     " << time_Match << std::endl
      << " time_UpdateParams_PostMatch:    " << time_UpdateParams_PostMatch << std::endl
      << " time_FilterMatches:             " << time_FilterMatches << std::endl
      << " time_EvalErrorFunc:             " << time_EvalErrorFunc << std::endl
      << " time_Callbacks:                 " << time_Callbacks << std::endl
      << " time_Extras:                    " << time_Extras << std::endl;
    codeProfiler.Reset();
    codeProfiler.Start();
#endif


#ifdef ENABLE_CODE_TRACE
    std::cout << "Termination Test" << std::endl;
#endif

#if 0
	vctDynamicVector<vct3> matchPts;
	vctDynamicVector<vct3> matchNorms;
	pAlgorithm->ReturnMatchPts(matchPts, matchNorms);

	cisstMesh samplePts;
	samplePts.vertices.SetSize(matchPts.size());
	samplePts.vertexNormals.SetSize(matchNorms.size());
	samplePts.vertices = matchPts;
	samplePts.vertexNormals = matchNorms;
	samplePts.SavePLY(std::to_string(iter) + ".ply");
#endif

    //-- Termination Test --//

    dR.From(dF.Rotation());   // convert rotation to Rodrigues form
    dAng01 = dAng12;
    dAng12 = dAng;
    dAng = dR.Norm();
    dPos01 = dPos12;
    dPos12 = dPos;
    dPos = dF.Translation().Norm();

    // Algorithm specific termination
    //  also enables algorithm to update the registration
    //  to a different iteration if desired
    if (pAlgorithm->ICP_Terminate(Freg))
    {
      totalTimer.Stop();
      termMsg << std::endl << "Termination Condition:  Termination Requested by Algorithm" << std::endl;
      break;  // exit iteration loop
    }

    // Consider termination
    if (dAng < opt.dAngThresh && dPos < opt.dPosThresh && dS < opt.dShapeThresh)
    {
      // Termination Test
      //  Note: max iterations is enforced by for loop
      if ((dAng < opt.dAngTerm && dPos < opt.dPosTerm && dS < opt.dShapeTerm) // TODO: modify the termination condition
        || E < opt.minE
        || tolE < opt.tolE)
      {
        // termination condition must be satisfied for min number of consecutive iterations
        terminateIter++;
        if (terminateIter >= opt.termHoldIter)
        {
          // prepare termination message
          totalTimer.Stop();
          termMsg << std::endl << "Termination Condition satisfied for " << opt.termHoldIter << " iterations: " << std::endl;
          if (E < opt.minE) termMsg << "reached minE (" << opt.minE << ")" << std::endl;
          else if (tolE < opt.tolE) termMsg << "reached min dE/E (" << opt.tolE << ")" << std::endl;
		  else termMsg << "reached min dAngle & dTrans & dShape (" << opt.dAngTerm * 180 / cmnPI << "/" << opt.dPosTerm << "/" << opt.dShapeTerm << ")" << std::endl;
          break;  // exit iteration loop
        }
      }
	  else if ( (E == E0) && /*(dAng == dAng01) &&*/ (dPos == dPos01) )
	  {
		  terminateIter++;
		  if (terminateIter > opt.termHoldIter)
		  {
			  // prepare termination message
			  totalTimer.Stop();
			  termMsg << std::endl << "Termination Condition satisfied for " << opt.termHoldIter << " iterations: " << std::endl;
			  termMsg << "oscilating between " << E << "/" << E1 << ";" << dAng << "/" << dAng12 << ";" << dPos << "/" << dPos12 << std::endl;
			  break;  // exit iteration loop
		  }
	  }
      else
      {
        terminateIter = 0;
      }
    }
    else
    {
      terminateIter = 0;
    }

    if (iter == opt.maxIter)
    {
      // prepare termination message
      totalTimer.Stop();
      termMsg << std::endl << "Termination Condition: reached max iteration (" << opt.maxIter << ")" << std::endl;
    }

#ifdef ENABLE_CODE_PROFILER
    time_Extras = codeProfiler.GetElapsedTime();
    codeProfiler.Reset();
    codeProfiler.Start();
#endif

  }
  iterTimer.Stop();

  // complete termination message
  termMsg << " E: " << E << std::endl;
  termMsg << " dE/E: " << tolE << std::endl;
  termMsg << " dAng: " << dAng12 * 180 / cmnPI << " " << dAng01*180.0 / cmnPI << " (deg)" << std::endl;
  termMsg << " dPos: " << dPos12 << " " << dPos01 << std::endl;
  termMsg << " dShp: " << prevShapeNorm << " " << ShapeNorm << std::endl;
  termMsg << " iter: " << iter << std::endl;
  termMsg << " runtime: " << totalTimer.GetElapsedTime() << std::endl;
  termMsg << std::endl << Freg << std::endl;
  termMsg << "scale:\n   " << scale << std::endl;
  if (opt.deformable)
	termMsg << "shape parameters:\n" << sp << std::endl;
  if (iterData.iter != iterBest)
  {
    termMsg << std::endl << "WARNING: best iteration (" << iterBest << ") is not the final iteration" << std::endl;
  }
  //std::cout << termMsg.str().c_str();

  // compute final match distance
  pAlgorithm->ComputeMatchStatistics(rt.MatchPosErrAvg, rt.MatchPosErrSD);
  pAlgorithm->PrintMatchStatistics(termMsg);

  rt.termMsg = termMsg.str();
  rt.Freg = Freg;    
  rt.runTime = totalTimer.GetElapsedTime();
  rt.numIter = iter;  
  rt.nOutliers = nOutliers;

  return rt;
}


void cisstICP::AddIterationCallback(Callback &callback)
{
  this->iterationCallbacks.push_back(callback);
}

void cisstICP::AddIterationCallbacks(std::vector<Callback> &callbacks)
{
  if (callbacks.size() > 0)
  {
    Callback *callbackArray = callbacks.data();
    iterationCallbacks.insert(iterationCallbacks.end(), callbackArray, callbackArray + callbacks.size());
  }
}

void cisstICP::ClearIterationCallbacks()
{
  this->iterationCallbacks.clear();
}

