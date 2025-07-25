/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 Wikki Ltd
    Copyright (C) 2021-2023 OpenCFD Ltd.
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
    Foam::faPatchData

Description
    Helper class for holding data during faPatch construction.
    Most data members are exposed at the moment.

\*---------------------------------------------------------------------------*/

#ifndef Foam_faPatchData_H
#define Foam_faPatchData_H

#include "labelList.H"
#include "labelPair.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class faPatch;
class dictionary;

/*---------------------------------------------------------------------------*\
                         Class faPatchData Declaration
\*---------------------------------------------------------------------------*/

class faPatchData
{
public:

    // Data Members

        word name_;
        word type_;

        label ownerPolyPatchId_;
        label neighPolyPatchId_;

        //- The owner/neighbour for processor patches
        int ownerProcId_;
        int neighProcId_;

        // Storge (temporary or otherwise) for edge labels
        labelList edgeLabels_;


    // Constructors

        //- Default construct
        faPatchData();


    // Member Functions

    // Opaque read-only access

        //- Return the name
        const word& name() const noexcept;

        //- Contents transcribed into a patch dictionary,
        //- usually including the edge labels.
        dictionary dict(const bool withEdgeLabels = true) const;


    // Other Functions

        //- Has name/type etc
        bool good() const;

        //- Reset data
        void clear();

        //- Clear and populate with values from finiteArea patch
        void assign(const faPatch& fap);

        //- Set values consistent with a processor coupling
        bool assign_coupled(int ownProci, int neiProci);

        //- True if owner/neighbour processor ids are non-equal
        bool coupled() const noexcept
        {
            return (ownerProcId_ != neighProcId_);
        }

        //- Does this side own the patch? Also true for non-coupled patches
        bool owner() const noexcept
        {
            return (ownerProcId_ <= neighProcId_);
        }

        //- Does the other side own the patch?
        bool neighbour() const noexcept
        {
            return !owner();
        }

        //- Ordered match with owner/neighbour patchPair.
        //  \return
        //      - 0: matched none
        //      - 1: matched first only (owner)
        //      - 2: matched second only (neighbour)
        //      - 3: matched both (owner/neighbour)
        //      .
        int matchPatchPair(const labelPair& patchPair) const noexcept;

        //- Unordered comparison with owner/neighbour patchPair.
        int comparePatchPair(const labelPair& patchPair) const noexcept;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
