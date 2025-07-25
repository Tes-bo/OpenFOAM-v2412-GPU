/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2012 OpenFOAM Foundation
    Copyright (C) 2019-2022 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "gnuplotGraphWriter.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
typedef graph::writer graphWriter;

namespace graphWriters
{
    defineTypeName(gnuplotWriter);
    addToRunTimeSelectionTable(graphWriter, gnuplotWriter, word);
}
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::graphWriters::gnuplotWriter::write
(
    const graph& g,
    Ostream& os
) const
{
    os  << "set term pngcairo" << nl
        << "set output \"" << word(g.title()) << ".png\"" << nl
        << "set title " << g.title() << " 0,0" << nl << "show title" << nl
        << "set xlabel " << g.xName() << " 0,0" << nl << "show xlabel" << nl
        << "set ylabel " << g.yName() << " 0,0" << nl << "show ylabel" << nl;

    label nplots = 0;
    forAllConstIters(g, iter)
    {
        os  << (nplots++ ? ", \\" : "plot \\") << nl;
        os  << "'-' title " << iter()->name() << " with lines";
    }
    os << "; pause -1" << nl;


    forAllConstIters(g, iter)
    {
        os  << nl;
        writeXY(g.x(), *iter(), os);
    }
}


// ************************************************************************* //
