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
#ifndef _cisstICP_h
#define _cisstICP_h

#include <limits>

#include "PDTreeBase.h"
//#include "PDTree_Mesh.h"
//#include "cisstPointCloud.h"

class algICP;      // forward declerations for mutual dependency


class cisstICP
{
  //
  // This class implements the high-level structure of an ICP-style
  //  registration algorithm. Low-level implementations of key routines
  //  are implemented in a set of algorithm classes, whose methods are 
  //  called from here.
  //  


  //--- TYPES ---//

public:

	// ICP command line options
	struct CmdLineOptions
	{
		// command line options
		std::string target;			// file name and location of target
		std::string input;			// file name and location of input
		std::string output;			// directory extention added to working dir where output will be stored
		std::string modeweights;	// file name and location of mode weights
		std::string workingdir;		// replaces the default working direction, where input and output data is stored
		std::string xfm;			// file name and location of initial guess transform
		std::string ssm;			// file name and location of statistical shape model
		std::string cov;			// file name and location of covariance matrices of positional noise model
		std::string axes;			// file name and location of major/minor axes of angular noise model

		int modes;
		int samples;
		int niters;

		float scale;
		float minpos, maxpos;
		float minang, maxang;
		float poutliers;
		float outminpos, outmaxpos;
		float outminang, outmaxang;
		float noiseinplane, noiseperpplane;
		float noisedeg, noiseecc;
		float nthresh, diagthresh;
		float rbounds, tbounds, sbounds, spbounds;

		bool deformable;	// is algorithm deformable?
		bool bScale;
		bool useDefaultTarget;
		bool useDefaultInput;
		bool readModeWeights;
		bool useDefaultOutput;
		bool useDefaultWorkingDir;
		bool useDefaultXfm;
		bool useDefaultSSM;
		bool useDefaultCov;
		bool useDefaultAxes;

		bool useDefaultNumModes;
		bool useDefaultNumSamples;
		bool useDefaultNumIters;

		bool useDefaultScale;
		bool useDefaultMinPos;
		bool useDefaultMaxPos;
		bool useDefaultMinAng;
		bool useDefaultMaxAng;
		bool useDefaultNumOutliers;
		bool useDefaultOutMinPos;
		bool useDefaultOutMaxPos;
		bool useDefaultOutMinAng;
		bool useDefaultOutMaxAng;
		bool useDefaultNoiseInPlane;
		bool useDefaultNoisePerpPlane;
		bool useDefaultNoiseDeg;
		bool useDefaultNoiseEcc;
		bool useDefaultNThresh;
		bool useDefaultDiagThresh;
		bool useDefaultRotationBounds;
		bool useDefaultTranslationBounds;
		bool useDefaultScaleBounds;
		bool useDefaultShapeParamBounds;

		// default constructor
		CmdLineOptions() :
			target("../../../test_data/MT.ply"),
			input(""), 
			modeweights(""),
			output(""),
			workingdir("../../../test_data/"),
			xfm(""),
			ssm("../../../test_data/atlas_mt.txt"),
			cov(""),
			axes(""),
			modes(3),
			samples(300),
			niters(100),
			scale(1.0),
			minpos(10.0),
			maxpos(20.0),
			minang(6.0),
			maxang(12.0),
			poutliers(0.0),
			outminpos(5.0),
			outmaxpos(10.0),
			outminang(5.0),
			outmaxang(10.0),
			noiseinplane(1.0),
			noiseperpplane(1.0),
			noisedeg(2.0),
			noiseecc(0.5),
			rbounds(DBL_MAX),
			tbounds(DBL_MAX),
			sbounds(0.3),
			spbounds(3.0),
			bScale(false),
			deformable(false),
			useDefaultTarget(true),
			useDefaultInput(true),
			readModeWeights(false),
			useDefaultOutput(true),
			useDefaultWorkingDir(true),
			useDefaultXfm(true),
			useDefaultSSM(true),
			useDefaultCov(true),
			useDefaultAxes(true),
			useDefaultNumModes(true),
			useDefaultNumSamples(true),
			useDefaultNumIters(true),
			useDefaultScale(true),
			useDefaultMinPos(true),
			useDefaultMaxPos(true),
			useDefaultMinAng(true),
			useDefaultMaxAng(true),
			useDefaultNumOutliers(true),
			useDefaultOutMinPos(true),
			useDefaultOutMaxPos(true),
			useDefaultOutMinAng(true),
			useDefaultOutMaxAng(true),
			useDefaultNoiseInPlane(true),
			useDefaultNoisePerpPlane(true),
			useDefaultNoiseDeg(true),
			useDefaultNoiseEcc(true),
			useDefaultRotationBounds(true),
			useDefaultTranslationBounds(true),
			useDefaultScaleBounds(true),
			useDefaultShapeParamBounds(true)
		{};
	};

  // ICP run-time options
  struct Options 
  {
    // run-time options
    std::string auxOutputDir; // directory for saving run-time logs
    bool    printOutput;      // print runtime output
	bool	deformable;
    // termination conditions
    unsigned int  maxIter;        // max iterations
    unsigned int  termHoldIter;   // min iterations for which termination condition must be satisfied
	unsigned int  numShapeParams; // number of shape parameters
    double  minE;             // min error value (E)                                                          (RHT: 0.00000000001)
    double  tolE;             // min % change in error dE/E
    double  dPosThresh;       // min change in position to consider termination (sample/model distance units) (RHT: 0.00005)
    double  dAngThresh;       // min change in angle to consider termination (radians)                        (RHT: 0.00005)
	double	dShapeThresh;
    double  dPosTerm;         // terminate if position and angle are less than
    double  dAngTerm;         //  these termination values
	double	dShapeTerm;		  // terminate if shape parameter doesn't change
    //double  errorRatioThresh; // error ratio constraint (E/Eprev > ratio && Eprev/E < ratio);               (RHT: 0.999) 
                            
    // default constructor
    Options() :
      auxOutputDir(""),
      printOutput(true),
	  deformable(false),
      maxIter(100),
      termHoldIter(2),
	  numShapeParams(0),
      minE(-std::numeric_limits<double>::max()),
      tolE(0.0),  // 0.005 good number if using this
      dPosThresh(0.1), dAngThresh(0.1*(cmnPI/180)),
	  dShapeThresh(0.1),
      dPosTerm(0.1), dAngTerm(0.1*(cmnPI/180)),
	  dShapeTerm(0.1)
      {};

    virtual std::string toString()
    {
      std::stringstream ss;
      ss
        << "ICP Options:" << std::endl
        << " maxIter:\t" << maxIter << std::endl
        << " termHoldIter:\t" << termHoldIter << std::endl
        << " minE:\t" << minE << std::endl
        << " tolE:\t" << tolE << std::endl
        << " dPosThresh:\t" << dPosThresh << std::endl
        << " dAngThresh (deg):\t" << dAngThresh * 180.0/cmnPI << std::endl
        << " dPosTerm:\t" << dPosTerm << std::endl
        << " dAngTerm (deg):\t" << dAngTerm * 180.0 / cmnPI << std::endl
        << " auxOutputDir:\t" << auxOutputDir << std::endl
        << " printOutput:\t" << printOutput << std::endl
        ;
      return ss.str();
    };
  };

  // ICP Return Type
  struct ReturnType
  {
    vctFrm3 Freg;         // final transformation
    std::string termMsg;  // termination message
    double runTime;
    double runTimeFirstMatch;
    unsigned int numIter;
    double  MatchPosErrAvg;
    double  MatchNormErrAvg;
    double  MatchPosErrSD;
    double  MatchNormErrSD;
    unsigned int nOutliers;

    ReturnType()
      : termMsg(""), runTime(0.0), runTimeFirstMatch(0.0), numIter(0),
        MatchPosErrAvg(0.0), MatchNormErrAvg(0.0), 
        MatchPosErrSD(0.0), MatchNormErrSD(0.0),
        nOutliers(0)
        {};
    ReturnType( std::string msg, double time, double timeFirstIter, unsigned int numIterations,
                double PosErrAvg, double NormErrAvg,
                double PosErrorSD, double NormErrSD,
                double nOutliers)
      : termMsg(msg), runTime(time), runTimeFirstMatch(timeFirstIter), numIter(numIterations),
        MatchPosErrAvg(PosErrAvg), MatchNormErrAvg(NormErrAvg),
        MatchPosErrSD(PosErrorSD), MatchNormErrSD(NormErrSD),
        nOutliers((unsigned int)nOutliers)
        {};
  };

  // callback function standard payload
  struct CallbackArg
  {
    unsigned int  iter;               // current iteration
    vctFrm3       Freg;               // current rigid transform estimate
    vctFrm3       dF;                 // incremental change in rigid transform on this iteration
	double		  scale;			  // current scale
	vctDynamicVector<double> S;		  // current shape params
    double        E;                  // error function value (RMS for standard ICP)
    double        tolE;               // percent change in error function value
    double        time;               // time transpired for this iteration
    unsigned int  nOutliers;          // number of outliers this iteration
    
    virtual void ThisMakesMePolymorphic(){};

    // Constructor
    CallbackArg() :
      // defaults
      iter(0),
      Freg(vctFrm3::Identity()),
      dF(vctFrm3::Identity()),
	  scale(0.0),
	  S(0.0),
      E(0.0),
      tolE(0.0),
      time(0.0),
      nOutliers(0)
      {};
  };

  // callback function prototype
  typedef void(*cbFuncType)(CallbackArg &, void *userData);

  // container for callback function and user-defined payload
  struct Callback
  {
    cbFuncType cbFunc;
    void *userData;

    // constructors
    Callback()
      : cbFunc(0),
        userData(0) {};
    Callback( cbFuncType cbFunc, void *userData )
      : cbFunc(cbFunc),
        userData(userData) {};
  };


  //--- VARIABLES ---//

protected:

  Options opt;
  algICP *pAlgorithm;

  vctFrm3 FGuess;
  vctFrm3 Freg, Fbest;        // rigid body transforms (model = F * sample)
  double ShapeNorm, prevShapeNorm;
  double  E, Ebest;           // registration error
  unsigned int iterBest;      // iteration with best registration error
  unsigned int nOutliers;     // number of outliers in last iteration

  std::vector< Callback > iterationCallbacks;
  CallbackArg iterData;


  //--- METHODS ---//

public:

  // constructor
  cisstICP()
    : pAlgorithm(NULL)
  {}

  // destructor
  ~cisstICP() {}


  ReturnType RunICP( 
    algICP *pAlg, const Options &opt, const vctFrm3 &FGuess,
    std::vector<Callback> *pUserCallbacks = NULL,
    bool bEnableAlgorithmCallbacks = true);

protected:

  ReturnType IterateICP();

  void AddIterationCallback(Callback &callback);
  void AddIterationCallbacks(std::vector<Callback> &callbacks);
  void AddDefaultIterationCallback();
  void ClearIterationCallbacks();

};

#endif // _cisstICP_h
