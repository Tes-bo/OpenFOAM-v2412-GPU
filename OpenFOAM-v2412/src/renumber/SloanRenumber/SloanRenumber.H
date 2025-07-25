/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2015 OpenFOAM Foundation
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
    Foam::SloanRenumber

Description
    Sloan renumbering algorithm

    E.g. http://www.apav.it/sito_ratio/file_pdf/ratio_4/capitolo_1.pdf

SourceFiles
    SloanRenumber.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_SloanRenumber_H
#define Foam_SloanRenumber_H

#include "renumberMethod.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class SloanRenumber Declaration
\*---------------------------------------------------------------------------*/

class SloanRenumber
:
    public renumberMethod
{
    // Private Data

        //- Use reverse indexing
        bool reverse_;


public:

    //- Runtime type information
    TypeName("Sloan");


    // Constructors

        //- Default construct, optionally with reverse
        explicit SloanRenumber(const bool reverse = false);

        //- Construct given the renumber dictionary
        explicit SloanRenumber(const dictionary& dict);


    //- Destructor
    virtual ~SloanRenumber() = default;


    // Member Functions

        //- Toggle reverse on/off
        void reverse(bool on) noexcept { reverse_ = on; }


    // With mesh topology

        //- Return the cell visit order (from ordered back to original cell id)
        //- using the mesh to determine the connectivity.
        virtual labelList renumber(const polyMesh& mesh) const;


    // With explicit topology

        //- Return the cell visit order (from ordered back to original cell id)
        virtual labelList renumber
        (
            const CompactListList<label>& cellCells
        ) const;

        //- Return the cell visit order (from ordered back to original cell id)
        virtual labelList renumber
        (
            const labelListList& cellCells
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
