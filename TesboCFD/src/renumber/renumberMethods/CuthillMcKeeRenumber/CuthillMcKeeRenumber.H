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
    Foam::CuthillMcKeeRenumber

Description
    Cuthill-McKee renumbering (CM or RCM)

SeeAlso
    Foam::meshTools::bandCompression

SourceFiles
    CuthillMcKeeRenumber.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_CuthillMcKeeRenumber_H
#define Foam_CuthillMcKeeRenumber_H

#include "renumberMethod.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class CuthillMcKeeRenumber Declaration
\*---------------------------------------------------------------------------*/

class CuthillMcKeeRenumber
:
    public renumberMethod
{
    // Private Data

        //- Use reverse indexing
        bool reverse_;


public:

    //- Runtime type information
    TypeNameNoDebug("CuthillMcKee");


    // Constructors

        //- Default construct, optionally with reverse
        explicit CuthillMcKeeRenumber(const bool reverse = false);

        //- Construct given the renumber dictionary
        explicit CuthillMcKeeRenumber(const dictionary& dict);

        //- Construct given the renumber dictionary (ignored)
        //- and specified reverse handling
        CuthillMcKeeRenumber(const dictionary& dict, const bool reverse);


    //- Destructor
    virtual ~CuthillMcKeeRenumber() = default;


    // Member Functions

        //- Toggle reverse on/off
        void reverse(bool on) noexcept { reverse_ = on; }


    // With mesh topology

        //- Return the cell visit order (from ordered back to original cell id)
        //- using the mesh to determine the connectivity.
        virtual labelList renumber
        (
            //! Mesh connectivity (see globalMeshData::calcCellCells)
            const polyMesh& mesh
        ) const;


    // With explicit topology

        //- Return the cell visit order (from ordered back to original cell id).
        virtual labelList renumber
        (
            //! Mesh connectivity
            const CompactListList<label>& cellCells
        ) const;

        //- Return the cell visit order (from ordered back to original cell id).
        virtual labelList renumber
        (
            //! Mesh connectivity
            const labelListList& cellCells
        ) const;
};


/*---------------------------------------------------------------------------*\
                   Class reverseCuthillMcKeeRenumber Declaration
\*---------------------------------------------------------------------------*/

//- Reverse Cuthill-McKee renumbering
class reverseCuthillMcKeeRenumber
:
    public CuthillMcKeeRenumber
{
public:

    //- Runtime type information
    TypeNameNoDebug("reverseCuthillMcKee");


    // Constructors

        //- Default construct
        reverseCuthillMcKeeRenumber();

        //- Construct given the renumber dictionary (ignored)
        explicit reverseCuthillMcKeeRenumber(const dictionary& dict);


    //- Destructor
    virtual ~reverseCuthillMcKeeRenumber() = default;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
