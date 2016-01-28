# See license at http://www.cisst.org/cisst/license.txt

cmake_minimum_required(VERSION 2.8)
project (cisstICP)

# NOTES:
#
# To get the last component of a path, use
# get_filename_component(... NAME)
# To get the parent directory, use
# get_filename_component(... PATH)

IF (CMAKE_CL_64 OR CMAKE_GENERATOR MATCHES Win64)
  set( X86_MODE 0 )  # 64-bit compile
ELSE (CMAKE_CL_64 OR CMAKE_GENERATOR MATCHES Win64)
  set( X86_MODE 1 )  # 32-bit compile
ENDIF (CMAKE_CL_64 OR CMAKE_GENERATOR MATCHES Win64)



# create a list of libraries needed for this project
set ( REQUIRED_CISST_LIBRARIES cisstCommon cisstVector cisstNumerical cisstOSAbstraction )

# find cisst and make sure the required libraries have been compiled
find_package (cisst REQUIRED ${REQUIRED_CISST_LIBRARIES})

if (cisst_FOUND_AS_REQUIRED)

  # load cisst configuration
  include (${CISST_USE_FILE})

  # IF (X86_MODE)
    # # GSL library
    # # (only have pre-compiled 32-bit versions of this library)
    # find_path( GSL_INCLUDE gsl
     # HINTS "C:/software/gsl-1.14-bin-msvc9/gsl-1.14/include"
    # )
    # find_library( GSL_LIB gsl.lib
     # HINTS "C:/software/gsl-1.14-bin-msvc9/gsl-1.14/lib"
    # )
    # find_library( GSL_CBLAS_LIB gslcblas.lib
     # HINTS "C:/software/gsl-1.14-bin-msvc9/gsl-1.14/lib"
    # )
    # include_directories(${GSL_INCLUDE})
  # ENDIF (X86_MODE)

  # dlib library
  #  (this library only consists of header files)
  find_path( DLIB_INCLUDE dlib
    NAMES "dlib/dlib_include_path_tutorial.txt"
    HINTS "C:/software/dlib-18.6"
    DOC  "Choose folder that CONTAINS the dlib folder"
    )
  include_directories(${DLIB_INCLUDE})

  # WildMagic5 Mathematics Library
  # Release
  find_path (WM5_BASE_DIR
    NAMES wm5path.txt
    PATHS "C:/software/WildMagic5"
    )
  find_path (WM5_LIB_PATH
    NAMES Wm5Core.lib
    #PATHS "C:/software/WildMagic5/SDK/Library/v120/x64/Release"
    HINTS ${WM5_BASE_DIR}"/SDK/Library/v120/x64/Release"
    )
  find_library( WM5_CORE_LIB Wm5Core.lib PATHS ${WM5_LIB_PATH})
  find_library( WM5_MATH_LIB Wm5Mathematics.lib PATHS ${WM5_LIB_PATH})
  # # Debug:
  # find_library( WM5MATH_LIB Wm5MathematicsD.lib Wm5CoreD.lib
    # HINTS "C:/software/WildMagic5/SDK/Library/v120/x64/Debug"
    # )
  find_path( WM5_LIB_INCLUDE
    NAMES "Wm5NoniterativeEigen3x3.h"
    PATHS
    HINTS ${WM5_BASE_DIR}"/SDK/Include"
    #HINTS "C:/software/WildMagic5/SDK/Include"
    DOC  ""
    )
  include_directories(${WM5_LIB_INCLUDE})  

  
  #=== ICP Library ===

  set (cisstICP_FILES
    cisstICP.cpp
    cisstICP.h
    # 3D Registration
    cisstCovTreeBase.cpp
    cisstCovTreeBase.h
    cisstCovTreeNode.cpp
    cisstCovTreeNode.h
    cisstCovTree_Mesh.cpp
    cisstCovTree_Mesh.h
    cisstCovTree_PointCloud.cpp
    cisstCovTree_PointCloud.h
    cisstAlgorithmCovTree.cpp
    cisstAlgorithmCovTree.h
    cisstAlgorithmCovTree_CP.cpp
    cisstAlgorithmCovTree_CP.h
    cisstAlgorithmCovTree_CP_Mesh.cpp
    cisstAlgorithmCovTree_CP_Mesh.h
    cisstAlgorithmCovTree_CP_PointCloud.cpp
    cisstAlgorithmCovTree_CP_PointCloud.h
    cisstAlgorithmCovTree_MLP.cpp
    cisstAlgorithmCovTree_MLP.h
    cisstAlgorithmCovTree_MLP_Mesh.cpp
    cisstAlgorithmCovTree_MLP_Mesh.h
    cisstAlgorithmCovTree_MLP_PointCloud.cpp
    cisstAlgorithmCovTree_MLP_PointCloud.h    
    cisstAlgorithmICP.cpp
    cisstAlgorithmICP.h
    cisstAlgorithmICP_StdICP.cpp
    cisstAlgorithmICP_StdICP.h
    cisstAlgorithmICP_StdICP_Mesh.h
    cisstAlgorithmICP_StdICP_PointCloud.h
    cisstAlgorithmICP_RobustICP.cpp
    cisstAlgorithmICP_RobustICP.h
    cisstAlgorithmICP_RobustICP_Mesh.h
    cisstAlgorithmICP_RobustICP_PointCloud.h
    cisstAlgorithmICP_IMLP.cpp
    cisstAlgorithmICP_IMLP.h
    cisstAlgorithmICP_IMLP_Mesh.cpp
    cisstAlgorithmICP_IMLP_Mesh.h
    cisstAlgorithmICP_IMLP_PointCloud.cpp
    cisstAlgorithmICP_IMLP_PointCloud.h
    cisstAlgorithmICP_IMLP_MahalDist.cpp
    cisstAlgorithmICP_IMLP_MahalDist.h 
    cisstAlgorithmICP_IMLP_MahalDist_PointCloud.cpp
    cisstAlgorithmICP_IMLP_MahalDist_PointCloud.h
    cisstAlgorithmICP_IMLP_MahalDist_Mesh.cpp
    cisstAlgorithmICP_IMLP_MahalDist_Mesh.h
    cisstAlgorithmICP_IMLP_ClosestPoint.cpp
    cisstAlgorithmICP_IMLP_ClosestPoint.h 
    cisstAlgorithmICP_IMLP_ClosestPoint_PointCloud.cpp
    cisstAlgorithmICP_IMLP_ClosestPoint_PointCloud.h
    cisstAlgorithmICP_IMLP_ClosestPoint_Mesh.cpp
    cisstAlgorithmICP_IMLP_ClosestPoint_Mesh.h
    cisstDirCovTreeBase.cpp
    cisstDirCovTreeBase.h
    cisstDirCovTreeNode.cpp
    cisstDirCovTreeNode.h
    cisstDirCovTree_Mesh.cpp
    cisstDirCovTree_Mesh.h
    cisstDirCovTree_PointCloud.cpp
    cisstDirCovTree_PointCloud.h
    cisstAlgorithmDirCovTree.cpp
    cisstAlgorithmDirCovTree.h
    cisstAlgorithmDirCovTree_vonMisesPrj.cpp
    cisstAlgorithmDirCovTree_vonMisesPrj.h
    cisstAlgorithmDirCovTree_vonMisesPrj_Mesh.cpp
    cisstAlgorithmDirCovTree_vonMisesPrj_Mesh.h    
    cisstAlgorithmDirICP.cpp
    cisstAlgorithmDirICP.h
    cisstAlgorithmDirICP_StdICP.cpp
    cisstAlgorithmDirICP_StdICP.h
    cisstAlgorithmDirICP_StdICP_Mesh.cpp
    cisstAlgorithmDirICP_StdICP_Mesh.h
    cisstAlgorithmDirICP_StdICP_PointCloud.cpp
    cisstAlgorithmDirICP_StdICP_PointCloud.h
    cisstAlgorithmDirICP_vMFG.cpp
    cisstAlgorithmDirICP_vMFG.h
    cisstAlgorithmDirICP_vMFG_Mesh.cpp
    cisstAlgorithmDirICP_vMFG_Mesh.h
    cisstAlgorithmDirICP_Kent.cpp
    cisstAlgorithmDirICP_Kent.h
    cisstAlgorithmDirICP_Kent_Mesh.cpp
    cisstAlgorithmDirICP_Kent_Mesh.h
    cisstAlgorithmDirICP_vonMisesPrj.cpp
    cisstAlgorithmDirICP_vonMisesPrj.h
    cisstAlgorithmDirICP_vonMisesPrj_Mesh.cpp
    cisstAlgorithmDirICP_vonMisesPrj_Mesh.h
    cisstAlgorithmDirICP_vonMisesPrj_dlibWrapper.cpp
    cisstAlgorithmDirICP_vonMisesPrj_dlibWrapper.h
    cisstMesh.cpp
    cisstMesh.h
    cisstTriangle.cpp
    cisstTriangle.h
    cisstTriangleClosestPointSolver.cpp
    cisstTriangleClosestPointSolver.h
    cisstBoundingBox.cpp
    cisstBoundingBox.h
    cisstPointCloud.h
    cisstPointCloud.cpp
    cisstPointCloudDir.h
    cisstPointCloudDir.cpp
    utilities.h
    utilities.cpp
    RegisterP2P.cpp
    RegisterP2P.h
    ./NumericalRecipes/mins.h
    ./NumericalRecipes/nr3.h
    wrapper_dlib.cpp
    wrapper_dlib.h
    Ellipsoid_OBB_Intersection_Solver.cpp
    Ellipsoid_OBB_Intersection_Solver.h
    cisstException.h
    # 2D Registration
    cisstDirCovTree2DBase.cpp
    cisstDirCovTree2DBase.h
    cisstDirCovTree2DNode.cpp
    cisstDirCovTree2DNode.h
    cisstDirCovTree2D_Edges.cpp
    cisstDirCovTree2D_Edges.h
    cisstAlgorithm2D_DirCovTree.cpp
    cisstAlgorithm2D_DirCovTree.h
    cisstAlgorithm2D_DirCovTree_CP.cpp
    cisstAlgorithm2D_DirCovTree_CP.h
    cisstAlgorithm2D_DirCovTree_CP_Edges.cpp
    cisstAlgorithm2D_DirCovTree_CP_Edges.h
    cisstAlgorithm2D_DirCovTree_vonMises.cpp
    cisstAlgorithm2D_DirCovTree_vonMises.h
    cisstAlgorithm2D_DirCovTree_vonMises_Edges.cpp
    cisstAlgorithm2D_DirCovTree_vonMises_Edges.h    
    # cisstAlgorithm2D_DirICP.cpp
    # cisstAlgorithm2D_DirICP.h
    # cisstAlgorithm2D_DirICP_StdICP.cpp
    # cisstAlgorithm2D_DirICP_StdICP.h
    # cisstAlgorithm2D_DirICP_StdICP_Edges.cpp
    # cisstAlgorithm2D_DirICP_StdICP_Edges.h
    cisstBoundingBox2D.cpp
    cisstBoundingBox2D.h
    cisstEdge2d.h
    cisstEdgeList2D.cpp
    cisstEdgeList2D.h
    cisstPointProjectionRoutines.cpp
    cisstPointProjectionRoutines.h
    utilities2D.h
    utilities2D.cpp
    cisstPointCloud2D.cpp
    cisstPointCloud2D.h
    )
  # IF (X86_MODE)
    # set (cisstICP_FILES
      # ${cisstICP_FILES}
      # wrapper_gsl.cpp
      # wrapper_gsl.h
      # )
  # ENDIF (X86_MODE)

  # cisstICP Library
  add_library (cisstICP #OBJECT
    ${cisstICP_FILES}    
    )
  include_directories(
    ${cisstICP_SOURCE_DIR}
    ${cisstICP_SOURCE_DIR}/NumericalRecipes
    )
  target_link_libraries( cisstICP
    ${WM5_MATH_LIB}
    ${WM5_CORE_LIB}
    )
  cisst_target_link_libraries ( cisstICP ${REQUIRED_CISST_LIBRARIES})
  # IF (X86_MODE)
    # target_link_libraries( cisstICP
      # ${GSL_LIB}
      # ${GSL_CBLAS_LIB}
      # )
  # ENDIF (X86_MODE)

  # TODO: how to adapt this for Release / Debug?
  add_custom_command( TARGET cisstICP
    POST_BUILD
    COMMAND lib.exe /OUT:./Release/cisstICPCombined.lib 
    #${CMAKE_CURRENT_BINARY_DIR}/cisstICP.lib   # only points to top-level of build directory
    ./Release/cisstICP.lib
    ${WM5_MATH_LIB}
    ${WM5_CORE_LIB}
    )
  # # this doesn't work for visual studio
  # add_custom_target(cisstICPCombined all
    # command ${cmake_cxx_archive_create} cissicpcombined.lib 
    # $<target_file:cissticp>
    # ${wm5_math_lib} ${wm5_core_lib}
    # )    
  
else (cisst_FOUND_AS_REQUIRED)
  message ("Information: code in ${CMAKE_CURRENT_SOURCE_DIR} will not be compiled, it requires      ${REQUIRED_CISST_LIBRARIES}")
endif (cisst_FOUND_AS_REQUIRED)
