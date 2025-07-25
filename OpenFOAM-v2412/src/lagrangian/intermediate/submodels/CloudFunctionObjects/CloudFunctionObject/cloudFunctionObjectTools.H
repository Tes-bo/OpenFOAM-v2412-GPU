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
    Foam::cloudFunctionObjectTools::collector

Description
    Implementation of template-invariant tools for
    various function objects such as Foam::ParticleHistogram.

SourceFiles
    cloudFunctionObjectTools.C

\*---------------------------------------------------------------------------*/

#ifndef cloudFunctionObjectTools_collector_H
#define cloudFunctionObjectTools_collector_H

#include "polyMesh.H"
#include "boundBox.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Namespace cloudFunctionObjectTools Declaration
\*---------------------------------------------------------------------------*/

namespace cloudFunctionObjectTools
{

/*---------------------------------------------------------------------------*\
                        Class collector Declaration
\*---------------------------------------------------------------------------*/

class collector
{
    // Private Data

        //- Flag to decide if the collector is patch based
        bool isPatch_;

        //- List of indices of collectors (zones or patches)
        labelList IDs_;

        //- List of names of collectors
        wordList names_;

        //- List of bounding-boxes of collectors
        List<boundBox> BBs_;


public:

    // Constructors

        //- Construct from dictionary
        collector
        (
            const dictionary& dict,
            const polyMesh& mesh
        );

        //- Copy construct
        collector(const collector& phc);

        //- No copy assignment
        void operator=(const collector&) = delete;


    //- Destructor
    virtual ~collector() = default;


    // Member Functions

    // Access

        //- Return the flag if the collector is patch based
        bool isPatch() const noexcept { return isPatch_; }

        //- Return number of collectors (zones or patches)
        label size() const noexcept { return IDs_.size(); }

        //- Return const reference to the indices of collectors
        const labelList& IDs() const noexcept { return IDs_; }

        //- Return const reference to the names of collectors
        const wordList& names() const noexcept { return names_; }

        //- Return const reference to the bounding-boxes of collectors
        const List<boundBox>& BBs() const noexcept { return BBs_; }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace cloudFunctionObjectTools
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
