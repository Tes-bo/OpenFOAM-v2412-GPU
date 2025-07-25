/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2019 OpenCFD Ltd.
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
    Foam::patchDistMethods::exact

Description
    Calculation of exact distance to nearest patch for all cells and
    boundary by constructing a search tree for all patch faces.

See also
    Foam::patchDistMethod::meshWave
    Foam::wallDist

SourceFiles
    exactPatchDistMethod.C

\*---------------------------------------------------------------------------*/

#ifndef exactPatchDistMethod_H
#define exactPatchDistMethod_H

#include "patchDistMethod.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class distributedTriSurfaceMesh;

namespace patchDistMethods
{

/*---------------------------------------------------------------------------*\
                          Class exact Declaration
\*---------------------------------------------------------------------------*/

class exact
:
    public patchDistMethod
{
    // Private Member Data

        //- Cache surface+searching of patch
        mutable autoPtr<distributedTriSurfaceMesh> patchSurfPtr_;


    // Private Member Functions

        const distributedTriSurfaceMesh& patchSurface() const;

        //- No copy construct
        exact(const exact&) = delete;

        //- No copy assignment
        void operator=(const exact&) = delete;


public:

    //- Runtime type information
    TypeName("exactDistance");


    // Constructors

        //- Construct from coefficients dictionary, mesh
        //  and fixed-value patch set
        exact
        (
            const dictionary& dict,
            const fvMesh& mesh,
            const labelHashSet& patchIDs
        );

        //- Construct from mesh and fixed-value patch set
        exact
        (
            const fvMesh& mesh,
            const labelHashSet& patchIDs
        );


    // Member Functions

        //- Update cached geometry when the mesh moves
        virtual bool movePoints()
        {
            // ? Reconstruct patch surface?
            //patchSurfPtr_.clear();
            return true;
        }

        //- Update cached topology and geometry when the mesh changes
        virtual void updateMesh(const mapPolyMesh&)
        {
            // ? Reconstruct patch surface?
            //patchSurfPtr_.clear();
        }

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
