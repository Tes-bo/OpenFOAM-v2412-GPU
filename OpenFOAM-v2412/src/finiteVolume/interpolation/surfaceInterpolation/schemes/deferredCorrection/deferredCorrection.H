/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017 OpenCFD Ltd.
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
    Foam::deferredCorrection

Description
    Deferred correction interpolation scheme wrapper around a run-time
    selectable base scheme.

    The characteristics of the base scheme are recovered by applying an
    explicit correction to the upwind scheme weights.

Usage
    Example of the \c deferredCorrection scheme applied to the \c linear
    scheme:
    \verbatim
    divSchemes
    {
        .
        .
        div(phi,U)      Gauss deferredCorrection linear;
        .
        .
    }
    \endverbatim

SourceFiles
    deferredCorrection.C

SeeAlso
    Foam::upwind

\*---------------------------------------------------------------------------*/

#ifndef deferredCorrection_H
#define deferredCorrection_H

#include "upwind.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class deferredCorrection Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class deferredCorrection
:
    public upwind<Type>
{
    // Private data

        //- Base scheme
        tmp<surfaceInterpolationScheme<Type>> tbaseScheme_;


    // Private Member Functions

        //- No copy construct
        deferredCorrection(const deferredCorrection&) = delete;

        //- No copy assignment
        void operator=(const deferredCorrection&) = delete;


public:

    //- Runtime type information
    TypeName("deferredCorrection");


    // Constructors

        //- Construct from Istream.
        //  The name of the flux field is read from the Istream and looked-up
        //  from the mesh objectRegistry
        deferredCorrection
        (
            const fvMesh& mesh,
            Istream& is
        )
        :
            upwind<Type>(mesh, is),
            tbaseScheme_
            (
                surfaceInterpolationScheme<Type>::New(mesh, is)
            )
        {}

        //- Construct from faceFlux and Istream
        deferredCorrection
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            Istream& is
        )
        :
            upwind<Type>(mesh, faceFlux, is),
            tbaseScheme_
            (
                surfaceInterpolationScheme<Type>::New(mesh, faceFlux, is)
            )
        {}


    // Member Functions

        using upwind<Type>::weights;

        //- Return true if this scheme uses an explicit correction
        virtual bool corrected() const
        {
            return true;
        }

        //- Return the explicit correction to the face-interpolate
        virtual tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
        correction
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
