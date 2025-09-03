/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 OpenCFD Ltd.
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
    Foam::patchDistMethods::directionalMeshWave

Description
    Variant of \c meshWave distance-to-patch method.

    Ignores the component in the specified direction. Can be used e.g. to
    calculate the distance in the z-direction only.

Usage
    Example of specification in \c fvSchemes.wallDist:
    \verbatim
        wallDist
        {
            // Mandatory entries (unmodifiable)
            method      directionalMeshWave;
            n           (0 0 1);

            // Optional (inherited) entries (unmodifiable)
            ...
        }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                        | Type | Req'd | Dflt
      method       | Method name: directionalMeshWave   | word |  yes  | -
      n            | The direction component to ignore  | vector |  yes  | -
    \endtable

    The inherited entries are elaborated in:
     - \link meshWavePatchDistMethod.H \endlink

See also
    Foam::wallDist
    Foam::patchDistMethod::meshWave
    Foam::patchDistMethod::Poisson
    Foam::directionalWallPoint

SourceFiles
    directionalMeshWavePatchDistMethod.C

\*---------------------------------------------------------------------------*/

#ifndef directionalMeshWavePatchDistMethod_H
#define directionalMeshWavePatchDistMethod_H

#include "meshWavePatchDistMethod.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace patchDistMethods
{

/*---------------------------------------------------------------------------*\
                          Class directionalMeshWave Declaration
\*---------------------------------------------------------------------------*/

class directionalMeshWave
:
    public meshWave
{
    // Private Member Data

        //- The direction component to ignore
        const vector n_;

        //- No copy construct
        directionalMeshWave(const directionalMeshWave&) = delete;

        //- No copy assignment
        void operator=(const directionalMeshWave&) = delete;


public:

    //- Runtime type information
    TypeName("directionalMeshWave");


    // Constructors

        //- Construct from coefficients dictionary, mesh
        //- and fixed-value patch set
        directionalMeshWave
        (
            const dictionary& dict,
            const fvMesh& mesh,
            const labelHashSet& patchIDs
        );


    // Member Functions

        //- Correct the given distance-to-patch field
        virtual bool correct(volScalarField& y);

        //- Correct the given distance-to-patch and normal-to-patch fields
        virtual bool correct(volScalarField& y, volVectorField& n);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace patchDistMethods
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
