/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::fvBoundaryMeshMapper

Description
    Foam::fvBoundaryMeshMapper

\*---------------------------------------------------------------------------*/

#ifndef Foam_fvBoundaryMeshMapper_H
#define Foam_fvBoundaryMeshMapper_H

#include "PtrList.H"
#include "fvPatchMapper.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class fvBoundaryMeshMapper Declaration
\*---------------------------------------------------------------------------*/

class fvBoundaryMeshMapper
:
    public PtrList<fvPatchMapper>
{
public:

    // Public Types

        //- The patch mapper type
        typedef fvPatchMapper mapper_type;


    // Generated Methods

        //- No copy construct
        fvBoundaryMeshMapper(const fvBoundaryMeshMapper&) = delete;

        //- No copy assignment
        void operator=(const fvBoundaryMeshMapper&) = delete;


    // Constructors

        //- Construct from components
        fvBoundaryMeshMapper
        (
            const fvMesh& mesh,
            const faceMapper& faceMap
        )
        :
            PtrList<fvPatchMapper>(mesh.boundary().size())
        {
            const fvBoundaryMesh& patches = mesh.boundary();

            forAll(patches, patchi)
            {
                set
                (
                    patchi,
                    new fvPatchMapper
                    (
                        patches[patchi],
                        faceMap
                    )
                );
            }
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
