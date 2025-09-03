/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2016-2018 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::fft

Description
    Fast fourier transform using the fftw library.

    The complex transform field is returned in the field supplied.  The
    direction of transform is supplied as an argument (-1 = forward, 1 =
    reverse).  The dimensionality and organisation of the array of values
    in space is supplied in the nn indexing array.

Note
    The fftw library uses int only (no longs) for its dimensionality.

SourceFiles
    fft.C

\*---------------------------------------------------------------------------*/

#ifndef fft_H
#define fft_H

#include "complexFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fft
{

public:

    enum transformDirection
    {
        FORWARD_TRANSFORM = -1, //!< The sign -1 = FFTW_FORWARD
        REVERSE_TRANSFORM = 1   //!< The sign +1 = FFTW_BACKWARD
    };

    static void fftRenumberRecurse
    (
        List<complex>& data,
        List<complex>& renumData,
        const UList<int>& nn,
        label nnprod,
        label ii,
        label l1,
        label l2
    );

    //- fftRenumber: fold the n-d data array to get the fft components in
    //- the right places.
    static void fftRenumber(List<complex>& data, const UList<int>& nn);

    //- Transform real-value data
    //  - uses the fftw real to half-complex method
    //  - result size is field.size()/2 + 1
    static tmp<complexField> realTransform1D(const scalarField& field);


    //- Transform real-value data
    //  - uses the fftw real to half-complex method
    //  - result size is field.size()/2 + 1
    static tmp<complexField> realTransform1D(const tmp<scalarField>& field);


    //- Transform complex-value data
    static void transform
    (
        complexField& field,
        const UList<int>& nn,
        transformDirection fftDirection
    );


    static tmp<complexField> forwardTransform
    (
        const tmp<complexField>& field,
        const UList<int>& nn
    );


    static tmp<complexField> reverseTransform
    (
        const tmp<complexField>& field,
        const UList<int>& nn
    );


    static tmp<complexVectorField> forwardTransform
    (
        const tmp<complexVectorField>& field,
        const UList<int>& nn
    );


    static tmp<complexVectorField> reverseTransform
    (
        const tmp<complexVectorField>& field,
        const UList<int>& nn
    );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
