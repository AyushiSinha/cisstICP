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
#include "cisstNumerical.h"
#include "utilities2D.h"

#define ENABLE_UTILITIES2D_DEBUG

// 2D method
void ComputeCovEigenDecomposition_SVD(const vct2x2 &C, vct2 &eigVal, vct2x2 &eigVct)
{
  // NOTE: for a covariance (positive-definite) matrix, the singular values are 
  //       synonymous with the eigenvalues since all eigenvalues must be positive
  
  // Compute SVD of C
  //   C = U*diag(S)*V'   where U = V
  static vctFixedSizeMatrix<double, 2, 2, VCT_COL_MAJOR> A;
  static vctFixedSizeMatrix<double, 2, 2, VCT_COL_MAJOR> U;
  static vctFixedSizeMatrix<double, 2, 2, VCT_COL_MAJOR> Vt;
  //static vct3 S;
  static nmrSVDFixedSizeData<2, 2, VCT_COL_MAJOR>::VectorTypeWorkspace workspace;
  try
  {
    A.Assign(C);  // must use "assign" rather than equals to properly transfer between different vector orderings
    nmrSVD(A, U, eigVal, Vt, workspace);
  }
  catch (...)
  {
    assert(0);
  }

#ifdef ENABLE_UTILITIES2D_DEBUG
  if (eigVal(1) < 0.0 || eigVal(1) > eigVal(0))
  {
    std::cout << std::endl << "========> ERROR: ComputeCovEigenDecomposition_SVD() eigenvalues misordered or less than zero!" << std::endl << std::endl;
    assert(0);
  }
#endif

  eigVct.Assign(U);
}