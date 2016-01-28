
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
#include "cisstPointProjectionRoutines.h"


vct2 ProjectOnLineSegment(const vct2 &x, const vct2 &a, const vct2 &b, double *pLambda)
{
    vct2 ax = x - a;
    vct2 ab = b - a;
    double lambda = (ax*ab) / (ab*ab);
    if (lambda <= 0.0)
    {
        if (pLambda) *pLambda = 0.0;
        return a;
    };
    if (lambda > 1.0)
    {
        if (pLambda) *pLambda = 1.0;
        return b;
    };
    if (pLambda) *pLambda = lambda;
    return a + ab*lambda;
};


vct3 ProjectOnLineSegment(const vct3 &x, const vct3 &a, const vct3 &b)
{
    vct3 ax = x - a;
    vct3 ab = b - a;
    double lambda = (ax*ab) / (ab*ab);
    if (lambda <= 0.0) { return a; };
    if (lambda > 1.0) { return b; };
    return a + ab*lambda;
};
