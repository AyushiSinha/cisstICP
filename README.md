# Most-Likely Point Paradigm

This repository provides data and source code for algorithms developed within the Most-Likely Point and Deformable Most-Likely Point paradigms described in the following papers:

**IMLP**: Seth D. Billings, Emad M. Boctor, and Russell H. Taylor, "Iterative Most-Likely Point Registration (IMLP): A Robust Algorithm for Computing Optimal Shape Alignment", PLOS One 10 (3): e0117688 (2015)

**IMLOP**: Seth D. Billings, Russell H. Taylor, "Iterative most likely oriented point registration", Medical Image Computing and Computer-Assisted Intervention, Lecture Notes in Computer Science 8673: 178--185 (2014)

**G-IMLOP**: Seth D. Billings, Russell H. Taylor, "Generalized iterative most likely oriented-point (G-IMLOP) registration", International Journal of Computer Assisted Radiology and Surgery 10 (8): 1213--1226 (2015)

**P-IMLOP**: Seth D. Billings, Hyun J. Kang, Alexis Cheng, Emad M. Boctor, Peter Kazanzides, Russell H. Taylor, "Minimally invasive registration for computer-assisted orthopedic surgery: combining tracked ultrasound and bone surface points via the P-IMLOP algorithm", International Journal of Computer Assisted Radiology and Surgery 10 (6): 761--771 (2015)

<!-- **V-IMLOP**: Seth D. Billings, Ayushi Sinha, Austin Reiter, Simon Leonard, Masaru Ishii, Gregory D. Hager, Russell H. Taylor, "Anatomically Constrained Video-CT Registration via the V-IMLOP Algorithm", Medical Image Computing and Computer-Assisted Intervention, Lecture Notes in Computer Science 9902: 133--141 (2016) -->

**D-IMLP**: Ayushi Sinha, Seth D. Billings, Austin Reiter, Masaru Ishii, Gregory D. Hager, Russell H. Taylor, "The deformable most likely point paradigm", TBP

**D-IMLOP**: Ayushi Sinha, Seth D. Billings, Austin Reiter, Masaru Ishii, Gregory D. Hager, Russell H. Taylor, "The deformable most likely point paradigm", TBP

**GD-IMLOP**: Ayushi Sinha, Xington Liu, Austin Reiter, Masaru Ishii, Gregory D. Hager, Russell H. Taylor, "Endoscopic navigation in the absence of CT imaging", Medical Image Computing and Computer-Assisted Intervention, Lecture Notes in Computer Science (2018)

NOTE: This repository also contains an implementation of ICP, directional ICP and robust ICP. Other algorithms evaluated in these papers, for which source code is not included in this repository, include:

 - GICP

    > downloaded from http://www.robots.ox.ac.uk/~avsegal/generalized_icp.html

    > see IMLP paper for minor modifications required for termination condition and to set covariance matrices

 - CPD

    > downloaded from https://sites.google.com/site/myronenko/research/cpd

    > see IMLP paper for minor modifications required for termination condition

NOTE: To run experiments described in the IMLP paper, go to https://github.com/sbillin/IMLP. Compile the source code following instructions there (or here) and run the experiments in the `PLOSONE` folder.

## Required Software Dependencies

 - CMake (http://www.cmake.org)

 - Visual Studio for Windows (http://www.visualstudio.com/downloads/download-visual-studio-vs)
    > application was developed in Visual Studio 2013

 - CISST Libraries (https://github.com/jhu-cisst/cisst/)
    > must compile this from source (instructions below)

 - WildMagic5 Libraries (http://www.geometrictools.com/Downloads/Downloads.html)
    > requires: Wm5Core.lib, Wm5Mathematics.lib

    > must compile these from source (instructions below)

    > included in `dependencies` folder

 - DLib (http://dlib.net/)
    > for iterative solutions to our minimization problems

    > included in `dependencies` folder

 - Numerical Recipes (http://numerical.recipes/)
    > needed only for G-IMLOP and GD-IMLOP algorithm; other algorithms don't use it and the code can be built without it

 - RPly Library
    > for reading/writing PLY files

    > included in `dependencies` folder

 - PLY_IO Library (if using Matlab)
    > for reading/writing PLY files in Matlab

    > included in `dependencies/Matlab Dependencies` folder

NOTE: Visual Studio 12 2013 Win64 was used successfully when creating this `README`. Visual Stuido 14 2015 Win64 did not work due to errors such as `cisstNetlib_f2c.lib(endfile.obj) : error LNK2001: unresolved external symbol __imp_sprintf`. This issue has not been solved and it is recommended to build the `cisstICP` application using Visual Studio 2013.

## Instructions for Compiling Dependencies
 
### CISST Libraries: 
https://github.com/jhu-cisst/cisst/
 - Clone `cisst` repository from github

    > git clone https://github.com/jhu-cisst/cisst/ 
    
    > only jhu-cisst/cisst required (not cisst-saw, etc.)
 - Create a build directory `cisst_build` alongside (not inside) the `cisst` directory that was cloned from git
 - Create Visual Studio solution for CISST Libraries using CMake <br />
    > run CMake on `cisst` source directory with build directory set to `cisst_build`

       -- check boxes for libraries: CISST_cisstCommon, CISST_cisstVector, CISST_cisstNumerical, CISST_cisstOSAbstraction

       -- uncheck boxes for cisstMultitask, cisstMutliTask_, cisstParameterTypes, cisstRobot_ 

       -- check box for CISST_HAS_CISSTNETLIB 

    > configure with `Visual Studio 12 2013 Win64` as the compiler 
       -- check box: CISSTNETLIB_DOWNLOAD_NOW 

    > configure and generate project 
    
    > open the `cisst_build/cisst.sln` file and build in Visual Studio with build mode set to `Release` and `x64`

### WildMagic5 Libraries: 
http://www.geometrictools.com/Downloads/Downloads.html 
These libraries are used solely for the closed-form implementation of computing the decomposition of a 3x3 covariance matrix. The libraries are no longer available online from their author, but they are included in the dependencies folder of the cisstICP repo.
 - Build WildMagic5 <br />
    > extract source from the `WildMagic5p13.zip` file located in the `dependencies` folder
    
    > open the `GeometricTools\WildMagic5\WildMagic5Wgl_VC120.sln` file in Visual Studio and set build configurations to `Release` and `x64` 
    
    > compile Wm5Core.lib and Wm5Mathematics.lib by right clicking `LibCore_VC120` and `LibMathematics_VC120` and selecting build 
 - See website of [WildMagic5](http://www.geometrictools.com) for more instructions on this

### Dlib Libraries:
 - Extract source from the `dlib-a8.6.7z` file located in the `dependencies` folder

### Numerical Recipes:
 - License to use the source code for this must be bought (usually available on amazon along with the book)
 - Contact sinha(at)jhu(dot)edu for modifications in numerical recipes code once license is obtained

### RPly Libraries:
 - Extract source from the `rply-1.1.4.zip` file located in the `dependencies` folder

### PLY_IO Libraries: 
If using Matlab
 - Extract source from the `PLY_IO` file located in the `dependencies/Matlab Dependencies` folder

## Instructions for Compiling Source Code
These instructions are based on the Windows platform with VisualStudio 2013. The output is a static library (cisstICP.lib).
 
### Compile cisstICP Library:
 - Clone `cisstICP` repository from github

    > git clone https://github.com/AyushiSinha/cisstICP
 - Create a build directory `cisstICP_build` within the `cisstICP` directory that was cloned from git
 - run CMake on source directory `cisstICP`

    > specify paths to cisstICP source code and build directories (`cisstICP` and `cisstICP_build`) and configure
    
    > set cisst_DIR to location of CISST Library build `cisst/cisst_build` and configure again. If the compiled CISST libraries were found correctly by CMake, then `CISST_USE_FILE` should now automatically contain a path to the file `Usecisst.cmake` 
    
    > specify path to WildMagic5 base directory by setting `WM5_BASE_DIR` to `cisstICP/dependencies/WildMagic5p13/GeometricTools/WildMagic5` (other WildMagic5 fields should auto-detect) 
    
    > set `DLIB_INCLUDE` to the path to the extracted dlib folder 
    
    > set `RPLY_DIR` to the path to `cisstICP/dependencies/rply-1.1.4` 
    
    > optionally, check `USE_EXTRA_ALGORITHMS` to include IMLP_CP (closest point), IMLP_MD (mahalanobis distance) and RobustICP in the build 
    
    > configure and generate project
 - Build cisstICP library in Visual Studio <br />
    > open `cisstICP_build/cisstICP.sln` 
    
    > build with build mode set to `Release` and `x64`

### Compile cisstICP_App:
 - Create a build directory `cisstICP_App_build` within the `cisstICP_App` directory

    > specify paths to cisstICP source code and build directories (`cisstICP_App` and `cisstICP_App_build`) 
 - run CMake on source directory `cisstICP_App` and configure <br />
    > set `cisst_DIR` to the path to the `cisst/cisst_build` folder and configure again. If the compiled CISST libraries were found correctly by CMake, then `CISST_USE_FILE` should now automatically contain a path to the file `Usecisst.cmake`
    
    > set `DLIB_INCLUDE` to the path to the extracted dlib folder. 
    
    > set `WM5_BASE_DIR` to the path to `cisstICP/dependencies/WildMagic5p13/GeometricTools/WildMagic5`  
    
    > set `MATLAB_ENGINE_INCLUDE_DIR` and `MATLAB_ENGINE_LIB_DIR` to your system paths for Matlab, such as `C:/Program Files/MATLAB/R2015a/extern/include` and `C:/Program Files/MATLAB/R2015a/extern/lib/win64/microsoft` (these paths will change depending on the Matlab version). 
    
    > set `cisstICP_LIB` to the path to `/cisstICP/cisstICP_build/Release/cisstICP.lib` 
    
    > set `cisstICP_LIB_INCLUDE` to the path to `cisstICP/cisstICP` 
    
    > configure again and Generate
 - Build cisstICP_App library in Visual Studio <br />
    > open the `cisstICP_App_build/cisstICP_App.sln` file in Visual Studio, set the build options to 'Release' and 'x64', and build the solution.

### Compile matlabICP Library:
This library provides a Matlab interface to the C++ code of the cisstICP Library. This library is not maintained by me.

 - Create a `matlabICP_build` folder alongside the `matlabICP` folder in the cloned repo

    > set CMake's source and build locations to these folders (`matlabICP` and `matlabICP_build`)
 - run CMake on source directory `matlabICP` and configure <br />
    > set `cisst_DIR` to the path to the `cisst/cisst_build` folder and configure again. If the compiled CISST libraries were found correctly by CMake, then `CISST_USE_FILE` should now automatically contain a path to the file `Usecisst.cmake`. 
    
    > set `DLIB_INCLUDE` to the path to the extracted dlib folder. 
    
    > set `WM5_BASE_DIR` to the path to `cisstICP/dependencies/WildMagic5p13/GeometricTools/WildMagic5` 
    
    > set `MATLAB_ENGINE_INCLUDE_DIR` and `MATLAB_ENGINE_LIB_DIR` to your system paths for Matlab, such as  `C:/Program Files/MATLAB/R2015a/extern/include` and `C:/Program Files/MATLAB/R2015a/extern/lib/win64/microsoft` (these paths will change depending on the Matlab version). 
    
    > set `cisstICP_LIB` to the path to `/cisstICP/cisstICP_build/Release/cisstICP.lib` 
    
    > set `cisstICP_LIB_INCLUDE` to the path to `cisstICP/cisstICP` 
    
    > configure again and Generate
 - Build matlabICP library in Visual Studio <br />
    > open the `matlabICP_build/matlabICP.sln` file in Visual Studio, set the build options to 'Release' and 'x64', and build the solution.
 - Add to Matlab
    > extract the `/cisstICP/dependencies/MatlabDependencies/mtimesx.zip` and add it to your Matlab path. 
    
    > add `/cisstICP/dependencies/MatlabDependencies/PLY_IO` to your Matlab path.  Note that the [PLY_IO library code](http://people.sc.fsu.edu/~jburkardt/m_src/ply_io/ply_io.html) included here containes a bug fix for `ply_read.m`, changing line 571 from
```matlab
        if ( ( nargin > 1 & strcmpi(Str,'Tri') ) || nargout > 2 )
```
to
```matlab
        if ( ( nargin > 1 && strcmpi(Str,'Tri') ) || nargout > 2 )
```

NOTE: If you do not have the Numerical Recipes code, then the Matlab interface for GIMLOP will fail to compile (since GIMLOP will not have been compiled in the C++ library), giving an error. This is not a problem, as the Matlab interfaces for the other algorithms won't be affected and can still be used. 

## Test Run
 - `cd` to `cisstICP_App/cisstICP_App_build/Release/` and run `ICP_App` via command line

    > this should run the executable with default settings
 - Test whether everything is running correctly <br />
    > run the following commands: 
```cmd
    ICP_App --alg IMLP --out testingforrelease
    ICP_App --alg DIMLP --out testingforrelease
    ICP_App --alg IMLOP --out testingforrelease
    ICP_App --alg DIMLOP --out testingforrelease
    ICP_App --alg GIMLOP --out testingforrelease
    ICP_App --alg GDIMLOP --out testingforrelease
    ICP_App --alg PIMLOP --out testingforrelease 
```
 - `cd` into `cisstICP_App/tests` and run `compare_outputs_script.m` in Matlab. If all algorithms `pass`, then everything is running correctly. This is a good test to run when developing to make sure changes to the code have not broken anything

 - Test if the matlabICP libraries are working (if using Matlab)
    > `cd` to `matlabICP/TestApps/` and run `App_Test_IMLOP.m` in Matlab (not maintained)
