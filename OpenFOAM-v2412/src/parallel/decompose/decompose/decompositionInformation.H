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
    Foam::decompositionInformation

Description
    Abstract base class for decomposition

SourceFiles
    decompositionInformation.C

\*---------------------------------------------------------------------------*/

#ifndef decompositionInformation_H
#define decompositionInformation_H

#include "polyMesh.H"
#include "CompactListList.H"
#include "IOstreams.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class decompositionInformation Declaration
\*---------------------------------------------------------------------------*/

class decompositionInformation
{
    //- Simple storage for organizing min/max/median quantities
    struct stats
    {
        label min;
        label max;
        label median;

        stats() : min(0), max(0), median(0) {}

        void clear()
        {
            min = 0;
            max = 0;
            median = 0;
        }

        Ostream& print(Ostream& os) const;
    };

    // Private data

        labelListList distrib_;
        label nDomains_;

        stats cellsInfo_;
        stats neighInfo_;
        stats facesInfo_;


    // Private Member Functions

        //- Populate from cell decomposition list
        void populate
        (
            const labelUList& adjncy,
            const labelUList& xadj,
            const labelUList& decomp,
            const label nDomains
        );

        //- No copy construct
        decompositionInformation(const decompositionInformation&) = delete;

        //- No copy assignment
        void operator=(const decompositionInformation&) = delete;


public:

    // Constructors

        //- Construct
        decompositionInformation
        (
            const labelUList& adjncy,
            const labelUList& xadj,
            const labelUList& decomp,
            const label nDomains
        );

        //- Construct
        decompositionInformation
        (
            const CompactListList<label>& cellCells,
            const labelUList& decomp,
            const label nDomains
        );


    //- Destructor
    ~decompositionInformation() = default;


    // Member Functions

        void clear();

        label nDomains() const noexcept
        {
            return nDomains_;
        }

        void printSummary(Ostream& os) const;

        void printDetails(Ostream& os) const;

        void printAll(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
