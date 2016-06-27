// ****************************************************************************
//
//    Copyright (c) 2014, Ayushi Sinha, Russell Taylor, Johns Hopkins University
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
#ifndef _DirPDTree2D_Points_h
#define _DirPDTree2D_Points_h


#include "DirPDTree2DBase.h"
#include "cisstPointCloud2D.h"


class DirPDTree2D_Points : public DirPDTree2DBase
{
    // NOTE: for function overrides, be sure "const" type is same as the base class otherwise,
    //       the base class function will be treated as a different function and not actually
    //       be overridden!

public:

    cisstPointCloud2D PointList(const vctDynamicVector<vct2> &points, 
		const vctDynamicVector<vct2> &orientations);


    //-- Methods --//

    // constructor
    //  nThresh   - min number of datums in subdivided node
    //  diagThresh - min physical size of subdivided node
    DirPDTree2D_Points(
        const vctDynamicVector<vct2> &points,
        const vctDynamicVector<vct2> &pointsNorm,
        int nThresh, double diagThresh, bool bUseOBB = true);

    // destructor
	~DirPDTree2D_Points();


    //-- Base Class Method Overrides --//

    //virtual vct2 DatumSortPoint(int datum) const;  // return sort point of this datum
    virtual vct2 DatumNorm(int datum) const;       // return normal orientation of this datum

    virtual void EnlargeBounds(const vctFrm2& F, int datum, BoundingBox2D& BB) const;
    virtual void EnlargeBounds(int datum, BoundingBox2D& BB) const;

};

#endif