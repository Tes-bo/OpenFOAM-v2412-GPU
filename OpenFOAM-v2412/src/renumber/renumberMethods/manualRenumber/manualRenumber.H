/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2022-2024 OpenCFD Ltd.
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
    Foam::manualRenumber

Description
    Renumber given a ordered-to-original cell association in a file

SourceFiles
    manualRenumber.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_manualRenumber_H
#define Foam_manualRenumber_H

#include "renumberMethod.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class manualRenumber Declaration
\*---------------------------------------------------------------------------*/

class manualRenumber
:
    public renumberMethod
{
    // Private Data

        fileName dataFile_;

public:

    //- Runtime type information
    TypeNameNoDebug("manual");


    // Constructors

        //- Construct with given data file
        explicit manualRenumber(const fileName& file);

        //- Construct given the renumber dictionary
        explicit manualRenumber(const dictionary& dict);


    //- Destructor
    virtual ~manualRenumber() = default;


    // Member Functions

        //- The renumbering method needs a polyMesh (for its IOobject)
        virtual bool needs_mesh() const { return true; }


    // With mesh topology

        //- Return the cell visit order (from ordered back to original cell id)
        //- using the mesh for its IOobject and instance.
        virtual labelList renumber(const polyMesh& mesh) const;


    // With explicit topology - Not implemented!

        //- Return the cell visit order (from ordered back to original cell id)
        //- Not implemented!
        virtual labelList renumber
        (
            const CompactListList<label>& cellCells
        ) const
        {
            NotImplemented;
            return labelList();
        }

        //- Return the cell visit order (from ordered back to original cell id)
        //- Not implemented!
        virtual labelList renumber
        (
            const labelListList& cellCells
        ) const
        {
            NotImplemented;
            return labelList();
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
