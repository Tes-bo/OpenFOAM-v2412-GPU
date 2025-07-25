/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "xmgraceGraphWriter.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

typedef graph::writer graphWriter;

namespace graphWriters
{
    defineTypeName(xmgraceWriter);
    addToRunTimeSelectionTable(graphWriter, xmgraceWriter, word);
}
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::graphWriters::xmgraceWriter::write
(
    const graph& g,
    Ostream& os
) const
{
    os  << "@title " << g.title() << nl
        << "@xaxis label " << g.xName() << nl
        << "@yaxis label " << g.yName() << nl;

    label nWritten_ = 0;

    forAllConstIters(g, iter)
    {
        os  << "@s" << nWritten_ << " legend "
            << iter()->name() << nl
            << "@target G0.S" << nWritten_ << nl
            << "@type xy" << nl;

        writeXY(g.x(), *iter(), os);
        os << endl;

        ++nWritten_;
    }
}


// ************************************************************************* //
