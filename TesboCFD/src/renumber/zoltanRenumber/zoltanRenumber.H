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
    Foam::zoltanRenumber

Description
    Use Zoltan

SourceFiles
    zoltanRenumber.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_zoltanRenumber_H
#define Foam_zoltanRenumber_H

#include "renumberMethod.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class zoltanRenumber Declaration
\*---------------------------------------------------------------------------*/

class zoltanRenumber
:
    public renumberMethod
{
    // Private Data

        const dictionary& coeffsDict_;


    // Private Member Functions

        //- No copy construct
        zoltanRenumber(const zoltanRenumber&) = delete;

        //- No copy assignment
        void operator=(const zoltanRenumber&) = delete;


public:

    //- Runtime type information
    TypeName("zoltan");


    // Constructors

        //- Construct given the renumber dictionary
        explicit zoltanRenumber(const dictionary& dict);


    //- Destructor
    virtual ~zoltanRenumber() = default;


    // Member Functions

        //- Renumbering method requires a polyMesh for its topology
        virtual bool needs_mesh() const { return true; }


    // With mesh topology

        //- Return the cell visit order (from ordered back to original cell id)
        //- uses the mesh for connectivity and global exchanges
        virtual labelList renumber
        (
            const polyMesh& mesh
        ) const;


    // With explicit topology - Not implemented!

        //- Return the cell visit order (from ordered back to original cell id).
        //- Not implemented!
        virtual labelList renumber
        (
            const CompactListList<label>& cellCells
        ) const
        {
            NotImplemented;
            return labelList();
        }

        //- Return the cell visit order (from ordered back to original cell id).
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
