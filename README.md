Most-Likely Point Paradigm
==========================

This repository provides data and source code for algorithms developed within the Most-Likely Point and Deformable Most-Likely Point paradigms described in the following papers:

**IMLP**: Seth Billings, Emad Boctor, and Russell Taylor, "Iterative Most-Likely Point Registration (IMLP): A Robust Algorithm for Computing Optimal Shape Alignment", PLOS One, 10 (3), 2015

**IMLOP**: Seth Billings, Russell Taylor, "Iterative most likely oriented point registration", MICCAI 2014

**G-IMLOP**: Seth Billings, Russell Taylor, "Generalized iterative most likely oriented-point (G-IMLOP) registration", IJCARS 2015

**P-IMLOP**: Seth Billings, Hyun Kang, Alexis Cheng, Emad Boctor, Peter Kazanzides, Russell Taylor, "Minimally invasive registration for computer-assisted orthopedic surgery: combining tracked ultrasound and bone surface points via the P-IMLOP algorithm", IJCARS 2015

**V-IMLOP**: Seth Billings, Ayushi Sinha, Austin Reiter, Simon Leonard, Masaru Ishii, Gregory Hager, Russell Taylor, "Anatomically Constrained Video-CT Registration via the V-IMLOP Algorithm", MICCAI 2016

**D-IMLP**: Ayushi Sinha, Seth Billings, Austin Reiter, Masaru Ishii, Gregory Hager, Russell Taylor, "The deformable most likely point paradigm", TBP

**D-IMLOP**: Ayushi Sinha, Seth Billings, Austin Reiter, Masaru Ishii, Gregory Hager, Russell Taylor, "The deformable most likely point paradigm", TBP

**GD-IMLOP**: Ayushi Sinha, Xington Liu, Austin Reiter, Masaru Ishii, Gregory Hager, Russell Taylor, "Endoscopic navigation in the absence of CT imaging", MICCAI 2018

NOTE: This repository also contains an implementation of ICP, directional ICP, robust ICP.

## Required Software Dependencies

 - CMake (http://www.cmake.org)
 - Visual Studio for Windows: http://www.visualstudio.com/downloads/download-visual-studio-vs
    > application was developed in Visual Studio 2013
 - CISST Libraries (https://trac.lcsr.jhu.edu/cisst)
    > must compile this from source (instructions below)
 - WildMagic5 Libraries: http://www.geometrictools.com/Downloads/Downloads.html
    > requires: Wm5Core.lib, Wm5Mathematics.lib<br />
    > must compile these from source
 - DLib: 
http://dlib.net/
    > for iterative solutions to our minimization problems
 - Numerical Recipes: 
http://numerical.recipes/
    > needed only for G-IMLOP and GD-IMLOP algorithm; other algorithms don't use it and the code can be built without it
 - RPly Library
    > for reading/writing PLY files
 - PLY_IO Library (if using Matlab)
    > for reading/writing PLY files in Matlab


## Instructions for Compiling Dependencies
 
### CISST Libraries: 
https://github.com/jhu-cisst/cisst/
 - Clone `cisst` repository from github <br />
    > git clone https://github.com/jhu-cisst/cisst/ <br />
    > only jhu-cisst/cisst required (not cisst-saw, etc.)
 - Create a build directory `cisst_build` alongside (not inside) the `cisst` directory that was cloned from git
 - Create Visual Studio solution for CISST Libraries using CMake <br />
    > run CMake on `cisst` source directory with build directory set to `cisst_build`  <br />
       -- check boxes for libraries: CISST_cisstCommon, CISST_cisstVector, CISST_cisstNumerical, CISST_cisstOSAbstraction <br />
       -- uncheck boxes for cisstMultitask, cisstMutliTask_, cisstParameterTypes, cisstRobot_ <br />
       -- check box for CISST_HAS_CISSTNETLIB <br />
    > configure with `Visual Studio 12 2013 Win64` as the compiler <br />
       -- check box: CISSTNETLIB_DOWNLOAD_NOW <br />
    > configure and generate project <br />
    > open the `cisst_build/cisst.sln` file and build in Visual Studio with build mode set to `Release` and `x64`

### WildMagic5 Libraries: 
http://www.geometrictools.com/Downloads/Downloads.html <br />
These libraries are used solely for the closed-form implementation of computing the decomposition of a 3x3 covariance matrix. The libraries are no longer available online from their author, but they are included in the dependencies folder of the cisstICP repo.
 - Build WildMagic5 <br />
    > extract source from the `WildMagic5p13.zip` file located in the `dependencies` folder <br />
    > open the `GeometricTools\WildMagic5\WildMagic5Wgl_VC120.sln` file in Visual Studio and set build configurations to `Release` and `x64` <br />
    > compile Wm5Core.lib and Wm5Mathematics.lib by right clicking `LibCore_VC120` and `LibMathematics_VC120` and selecting build <br />
 - See website of WildMagic5 for more instructions on this (http://www.geometrictools.com)

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
 - Clone `cisstICP` repository from github <br />
    > git clone https://github.com/AyushiSinha/cisstICP
 - Create a build directory `cisstICP_build` within the `cisstICP` directory that was cloned from git
 - run CMake on source directory `cisstICP` <br />
    > specify paths to cisstICP source code and build directories (`cisstICP` and `cisstICP_build`) and configure <br />
    > set cisst_DIR to location of CISST Library build `cisst/cisst_build` and configure again. If the compiled CISST libraries were found correctly by CMake, then `CISST_USE_FILE` should now automatically contain a path to the file `Usecisst.cmake` <br />
    > specify path to WildMagic5 base directory by setting `WM5_BASE_DIR` to `cisstICP/dependencies/WildMagic5p13/GeometricTools/WildMagic5` (other WildMagic5 fields should auto-detect) <br />
    > set `DLIB_INCLUDE` to the path to the extracted dlib folder <br />
    > set `RPLY_DIR` to the path to `cisstICP/dependencies/rply-1.1.4` <br />
    > optionally, check `USE_EXTRA_ALGORITHMS` to include IMLP_CP (closest point), IMLP_MD (mahalanobis distance) and RobustICP in the build <br />
    > configure and generate project
 - Build cisstICP library in Visual Studio <br />
    > open `cisstICP_build/cisstICP.sln` <br />
    > build with build mode set to `Release` and `x64`

### Compile cisstICP_App:
 - Create a build directory `cisstICP_App_build` within the `cisstICP_App` directory <br />
    > specify paths to cisstICP source code and build directories (`cisstICP_App` and `cisstICP_App_build`) 
 - run CMake on source directory `cisstICP_App` and configure <br />
    > set `cisst_DIR` to the path to the `cisst/cisst_build` folder and configure again. If the compiled CISST libraries were found correctly by CMake, then `CISST_USE_FILE` should now automatically contain a path to the file `Usecisst.cmake`. <br />
    > set `DLIB_INCLUDE` to the path to the extracted dlib folder. <br />
    > set `WM5_BASE_DIR` to the path to `cisstICP/dependencies/WildMagic5p13/GeometricTools/WildMagic5`  <br />
    > set `MATLAB_ENGINE_INCLUDE_DIR` and `MATLAB_ENGINE_LIB_DIR` to your system paths for Matlab, such as `C:/Program Files/MATLAB/R2015a/extern/include` and `C:/Program Files/MATLAB/R2015a/extern/lib/win64/microsoft` (these paths will change depending on the Matlab version). <br />
    > set `cisstICP_LIB` to the path to `/cisstICP/cisstICP_build/Release/cisstICP.lib` <br />
    > set `cisstICP_LIB_INCLUDE` to the path to `cisstICP/cisstICP` <br />
    > configure again and Generate
 - Build cisstICP_App library in Visual Studio <br />
    > open the `cisstICP_App_build/cisstICP_App.sln` file in Visual Studio, set the build options to 'Release' and 'x64', and build the solution.

### Compile matlabICP Library:
This library provides a Matlab interface to the C++ code of the cisstICP Library. This library is not maintained by me.

 - Create a `matlabICP_build` folder alongside the `matlabICP` folder in the cloned repo <br />
    > set CMake's source and build locations to these folders 
 - run CMake on source directory `matlabICP` and configure <br />
    > set `cisst_DIR` to the path to the `cisst/cisst_build` folder and configure again. If the compiled CISST libraries were found correctly by CMake, then `CISST_USE_FILE` should now automatically contain a path to the file `Usecisst.cmake`. <br />
    > set `DLIB_INCLUDE` to the path to the extracted dlib folder. <br />
    > set `WM5_BASE_DIR` to the path to `cisstICP/dependencies/WildMagic5p13/GeometricTools/WildMagic5` <br />
    > set `MATLAB_ENGINE_INCLUDE_DIR` and `MATLAB_ENGINE_LIB_DIR` to your system paths for Matlab, such as  <br />`C:/Program Files/MATLAB/R2015a/extern/include` and `C:/Program Files/MATLAB/R2015a/extern/lib/win64/microsoft` (these paths will change depending on the Matlab version). <br />
    > set `cisstICP_LIB` to the path to `/cisstICP/cisstICP_build/Release/cisstICP.lib` <br />
    > set `cisstICP_LIB_INCLUDE` to the path to `cisstICP/cisstICP` <br />
    > configure again and Generate
 - Build matlabICP library in Visual Studio <br />
    > open the `matlabICP_build/matlabICP.sln` file in Visual Studio, set the build options to 'Release' and 'x64', and build the solution.
 - Add to Matlab
    > extract the `/cisstICP/dependencies/MatlabDependencies/mtimesx.zip` and add it to your Matlab path. <br />
    > add `/cisstICP/dependencies/MatlabDependencies/PLY_IO` to your Matlab path.  Note that the [PLY_IO library code](http://people.sc.fsu.edu/~jburkardt/m_src/ply_io/ply_io.html) included here containes a bug fix for `ply_read.m`, changing line 571 from

    if ( ( nargin > 1 & strcmpi(Str,'Tri') ) || nargout > 2 )

to

    if ( ( nargin > 1 && strcmpi(Str,'Tri') ) || nargout > 2 )

NOTE: If you do not have the Numerical Recipes code, then the Matlab interface for GIMLOP will fail to compile (since GIMLOP will not have been compiled in the C++ library), giving an error. This is not a problem, as the Matlab interfaces for the other algorithms won't be affected and can still be used. 

## Test Run

 - `cd` to 'cisstICP_App/cisstICP_App_build/Release/' and run `ICP_App` via command line
    > this should run the executable with default settings
 - Test whether everything is running correctly
    > run the following commands: <br />
    > - ICP_App --alg IMLP --out testingforrelease <br />
    > - ICP_App --alg DIMLP --out testingforrelease <br />
    > - ICP_App --alg IMLOP --out testingforrelease <br />
    > - ICP_App --alg DIMLOP --out testingforrelease <br />
    > - ICP_App --alg GIMLOP --out testingforrelease <br />
    > - ICP_App --alg GDIMLOP --out testingforrelease <br />
    > - ICP_App --alg PIMLOP --out testingforrelease <br />
    > `cd` into `cisstICP_App/tests` and run `compare_outputs_script.m` in Matlab. If all algorithms pass, then everything is running correctly. This is a good test to run when developing to make sure changes to the code have not broken anything
 - Test if the matlabICP libraries are working
    > `cd` to `matlabICP/TestApps/` and run `App_Test_IMLOP.m` in Matlab



