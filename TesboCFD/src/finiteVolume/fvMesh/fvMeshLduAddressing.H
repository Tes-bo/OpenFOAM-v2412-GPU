/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2021 OpenCFD Ltd.
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
    Foam::fvMeshLduAddressing

Description
    Foam::fvMeshLduAddressing

SourceFiles
    fvMeshLduAddressing.C

\*---------------------------------------------------------------------------*/

#ifndef fvMeshLduAddressing_H
#define fvMeshLduAddressing_H

#include "lduAddressing.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class fvMeshLduAddressing Declaration
\*---------------------------------------------------------------------------*/

class fvMeshLduAddressing
:
    public lduAddressing
{
    // Private Data

        //- Lower as a subList of allOwner
        const labelList::subList lowerAddr_;

        //- Upper as a reference to neighbour
        const labelList& upperAddr_;

        //- Patch addressing as a list of sublists
        const UPtrList<const labelUList> patchAddr_;

        //- Patch field evaluation schedule
        const lduSchedule& patchSchedule_;


    // Private Member Functions

        //- No copy construct
        fvMeshLduAddressing(const fvMeshLduAddressing&) = delete;

        //- No copy assignment
        void operator=(const fvMeshLduAddressing&) = delete;


public:

    // Constructors

        //- Construct from components
        explicit fvMeshLduAddressing(const fvMesh& mesh)
        :
            lduAddressing(mesh.nCells()),
            lowerAddr_
            (
                labelList::subList
                (
                    mesh.faceOwner(),
                    mesh.nInternalFaces()
                )
            ),
            upperAddr_(mesh.faceNeighbour()),
            patchAddr_(mesh.boundary().faceCells()),
            patchSchedule_(mesh.globalData().patchSchedule())
        {}


    //- Destructor
    ~fvMeshLduAddressing() = default;


    // Member Functions

        //- Return lower addressing (i.e. lower label = upper triangle)
        const labelUList& lowerAddr() const noexcept
        {
            return lowerAddr_;
        }

        //- Return upper addressing (i.e. upper label)
        const labelUList& upperAddr() const noexcept
        {
            return upperAddr_;
        }

        //- Return patch addressing
        const labelUList& patchAddr(const label patchi) const
        {
            return patchAddr_[patchi];
        }

        //- Return patch field evaluation schedule
        const lduSchedule& patchSchedule() const noexcept
        {
            return patchSchedule_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
