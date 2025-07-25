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

\*---------------------------------------------------------------------------*/

#include "cloudFunctionObjectTools.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cloudFunctionObjectTools::collector::collector
(
    const dictionary& dict,
    const polyMesh& mesh
)
:
    isPatch_(false),
    IDs_(),
    names_(),
    BBs_()
{
    wordRes matcher;

    if (dict.readIfPresent("patches", matcher) && !matcher.empty())
    {
        isPatch_ = true;

        IDs_ = mesh.boundaryMesh().indices(matcher);

        names_.resize(IDs_.size());

        label count = 0;
        for (const label patchi : IDs_)
        {
            names_[count] = mesh.boundaryMesh()[patchi].name();
            ++count;
        }
    }
    else if (dict.readIfPresent("faceZones", matcher))
    {
        const faceZoneMesh& fzm = mesh.faceZones();

        IDs_ = fzm.indices(matcher);

        BBs_.resize(IDs_.size());
        names_.resize(IDs_.size());

        label count = 0;
        for (const label zonei : IDs_)
        {
            const faceZone& fz = fzm[zonei];
            names_[count] = fz.name();

            auto& bb = BBs_[count];
            ++count;

            bb.reset();

            const auto& faces = mesh.faces();
            const auto& points = mesh.points();
            for (const label facei : fz)
            {
                bb.add(points, faces[facei]);
            }
            bb.reduce();
            bb.inflate(0.05);
        }
    }

    if (matcher.empty() || IDs_.size() < 1)
    {
        FatalIOErrorInFunction(dict)
            << "No matching patches or face zones found: "
            << flatOutput(matcher) << nl
            << exit(FatalIOError);
    }
}


Foam::cloudFunctionObjectTools::collector::collector
(
    const collector& phc
)
:
    isPatch_(phc.isPatch_),
    IDs_(phc.IDs_),
    names_(phc.names_),
    BBs_(phc.BBs_)
{}


// ************************************************************************* //
