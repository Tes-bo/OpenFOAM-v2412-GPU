/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2022 OpenCFD Ltd.
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
    Foam::graphWriters::xmgraceWriter

Description
    Output and \b agr file for \em xmgrace
    (http://plasma-gate.weizmann.ac.il/Grace/)

SourceFiles
    xmgraceGraphWriter.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_xmgraceGraphWriter_H
#define Foam_xmgraceGraphWriter_H

#include "graph.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace graphWriters
{

/*---------------------------------------------------------------------------*\
                        Class xmgraceWriter Declaration
\*---------------------------------------------------------------------------*/

class xmgraceWriter
:
    public graph::writer
{
public:

    //- Runtime type information
    TypeNameNoDebug("xmgr");

    //- Default construct
    xmgraceWriter() = default;

    //- Destructor
    ~xmgraceWriter() = default;


    // Member Functions

        //- The fileName extension for this graph format
        word ext() const { return "agr"; }

        //- Write
        void write(const graph&, Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace graphWriters
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
