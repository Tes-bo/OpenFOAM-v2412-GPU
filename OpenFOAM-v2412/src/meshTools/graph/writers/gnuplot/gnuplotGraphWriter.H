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
    Foam::graphWriters::gnuplotWriter

Description
    Output in gnuplot (http://www.gnuplot.info) format

SourceFiles
    gnuplotGraphWriter.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_gnuplotGraphWriter_H
#define Foam_gnuplotGraphWriter_H

#include "graph.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace graphWriters
{

/*---------------------------------------------------------------------------*\
                        Class gnuplotWriter Declaration
\*---------------------------------------------------------------------------*/

class gnuplotWriter
:
    public graph::writer
{
public:

    //- Runtime type information
    TypeNameNoDebug("gnuplot");

    //- Default construct
    gnuplotWriter() = default;

    //- Destructor
    ~gnuplotWriter() = default;


    // Member Functions

        //- The fileName extension for this graph format
        word ext() const { return "gplt"; }

        //- Write
        void write(const graph&, Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace graphWriters
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
