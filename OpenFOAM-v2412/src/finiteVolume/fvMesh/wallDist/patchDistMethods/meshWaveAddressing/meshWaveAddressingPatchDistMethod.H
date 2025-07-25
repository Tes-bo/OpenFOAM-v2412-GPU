/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023 OpenCFD Ltd.
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
    Foam::patchDistMethods::meshWaveAddressing

Description
    Variant of patchDistMethods::meshWave that stores nearest-wall addressing
    instead of directly transporting wall information.

    Example of the wallDist specification in fvSchemes:
    \verbatim
        wallDist
        {
            method meshWaveAddressing;

            // Optional entry enabling the calculation
            // of the normal-to-wall field
            nRequired true;
        }
    \endverbatim

See also
    Foam::patchDistMethod::meshWave
    Foam::wallDist

SourceFiles
    meshWaveAddressingPatchDistMethod.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_meshWaveAddressingPatchDistMethod_H
#define Foam_meshWaveAddressingPatchDistMethod_H

#include "patchDistMethod.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace patchDistMethods
{

/*---------------------------------------------------------------------------*\
                     Class meshWaveAddressing Declaration
\*---------------------------------------------------------------------------*/

class meshWaveAddressing
:
    public patchDistMethod
{
protected:

    // Protected Member Data

        //- Do accurate distance calculation for near-wall cells.
        const bool correctWalls_;


    // Protected Member Functions

        //- No copy construct
        meshWaveAddressing(const meshWaveAddressing&) = delete;

        //- No copy assignment
        void operator=(const meshWaveAddressing&) = delete;


public:

    //- Runtime type information
    TypeName("meshWaveAddressing");


    // Constructors

        //- Construct from coefficients dictionary, mesh
        //- and fixed-value patch set
        meshWaveAddressing
        (
            const dictionary& dict,
            const fvMesh& mesh,
            const labelHashSet& patchIDs
        );

        //- Construct from mesh, fixed-value patch set and flag specifying
        //- whether or not to correct wall.
        //  Calculate for all cells. correctWalls : correct wall (face&point)
        //  cells for correct distance, searching neighbours.
        meshWaveAddressing
        (
            const fvMesh& mesh,
            const labelHashSet& patchIDs,
            const bool correctWalls = true
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
