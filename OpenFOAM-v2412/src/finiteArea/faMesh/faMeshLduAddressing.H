/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 Wikki Ltd
    Copyright (C) 2019-2021 OpenCFD Ltd.
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
    Foam::faMeshLduAddressing

Description
    lduAddressing wrapper for faMesh

SourceFiles
    faMeshLduAddressing.C

Author
    Zeljko Tukovic, FMENA
    Hrvoje Jasak, Wikki Ltd.

\*---------------------------------------------------------------------------*/

#ifndef faMeshLduAddressing_H
#define faMeshLduAddressing_H

#include "lduAddressing.H"
#include "faMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class faMeshLduAddressing Declaration
\*---------------------------------------------------------------------------*/

class faMeshLduAddressing
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


public:

    // Generated Methods

        //- No copy construct
        faMeshLduAddressing(const faMeshLduAddressing&) = delete;

        //- No copy assignment
        void operator=(const faMeshLduAddressing&) = delete;


    // Constructors

        //- Construct from components
        faMeshLduAddressing(const faMesh& mesh)
        :
            lduAddressing(mesh.nFaces()),
            lowerAddr_
            (
                labelList::subList
                (
                    mesh.edgeOwner(),
                    mesh.nInternalEdges()
                )
            ),
            upperAddr_(mesh.edgeNeighbour()),
            patchAddr_(mesh.boundary().edgeFaces()),
            patchSchedule_(mesh.globalData().patchSchedule())
        {}


    //- Destructor
    virtual ~faMeshLduAddressing() = default;


    // Member Functions

        //- Return number of interfaces
        virtual label nPatches() const noexcept
        {
            return patchAddr_.size();
        }

        //- Return lower addressing (i.e. lower label = upper triangle)
        virtual const labelUList& lowerAddr() const noexcept
        {
            return lowerAddr_;
        }

        //- Return upper addressing (i.e. upper label)
        virtual const labelUList& upperAddr() const noexcept
        {
            return upperAddr_;
        }

        //- Return patch addressing
        virtual const labelUList& patchAddr(const label patchi) const
        {
            return patchAddr_[patchi];
        }

        //- Return patch field evaluation schedule
        virtual const lduSchedule& patchSchedule() const noexcept
        {
            return patchSchedule_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
