/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2012 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "indexedOctree.H"
#include "treeBoundBox.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(indexedOctreeBase, 0);

    scalar indexedOctreeBase::perturbTol_ = 10*SMALL;
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

void Foam::indexedOctreeBase::writeOBJ
(
    Ostream& os,
    const treeBoundBox& bb,
    label& vertIndex,
    const bool writeLinesOnly
)
{
    // Annotate with '#box' which can be grep'd for later
    os  << "#box" << nl;

    pointField pts(bb.points());

    for (const point& p : pts)
    {
        os  << "v " << p.x() << ' ' << p.y() << ' ' << p.z() << nl;
    }

    if (writeLinesOnly)
    {
        for (const edge& e : treeBoundBox::edges)
        {
            os  << "l " << e[0] + vertIndex + 1
                << ' ' << e[1] + vertIndex + 1 << nl;
        }
    }
    else
    {
        for (const face& f : treeBoundBox::faces)
        {
            os  << 'f';
            for (const label fpi : f)
            {
                os  << ' ' << fpi + vertIndex + 1;
            }
            os  << nl;
        }
    }

    vertIndex += pts.size();
}


// ************************************************************************* //
