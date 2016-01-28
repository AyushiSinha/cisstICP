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
#ifndef cisstEdgeList2D_h
#define cisstEdgeList2D_h


#include <cisstVector.h>
#include "cisstEdge2D.h"
#include "cisstException.h"


class cisstEdgeList2D 
{

public:

	vctDynamicVector<cisstEdge2D>  Edges;
  unsigned int numEdges;

  // constructors
  cisstEdgeList2D() {};
  //cisstEdgeList2D(const char *fn) { ReadEdgeFile(fn); }

  // destructor
  ~cisstEdgeList2D() {};

  void SetEdges( const vctDynamicVector<vct2> &V1,
                 const vctDynamicVector<vct2> &V2,
                 const vctDynamicVector<vct2> &Norm );

	// get mesh vertex indexes for the given triangle index
	inline void GetEdgeVertices( unsigned int edgeIndex, vct2 &v1, vct2 &v2 ) const
	{ v1=Edges[edgeIndex].V1;
		v2=Edges[edgeIndex].V2;
	}

  // get triangle normal for the given triangle index
  inline vct2 GetEdgeNorm(unsigned int edgeIndex) const
  { return Edges[edgeIndex].Norm;
  }

  //inline vct2 ProjectOnEdge( vct2 &x, unsigned int edgeIndex )
  //{ return Edges[edgeIndex].ProjectOnEdge( x );
  //}

private:

};

#endif