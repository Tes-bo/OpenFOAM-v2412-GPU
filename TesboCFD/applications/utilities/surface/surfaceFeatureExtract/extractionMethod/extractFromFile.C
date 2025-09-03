/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017 OpenCFD Ltd.
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

#include "extractFromFile.H"
#include "edgeMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace surfaceFeaturesExtraction
{
    addNamedToRunTimeSelectionTable
    (
        method,
        extractFromFile,
        dictionary,
        extractFromFile
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfaceFeaturesExtraction::extractFromFile::extractFromFile
(
    const dictionary& dict
)
:
    method()
{
    const dictionary& coeffDict =
        dict.optionalSubDict("extractFromFileCoeffs");

    coeffDict.readEntry("featureEdgeFile", featureEdgeFile_);
    coeffDict.readIfPresent("geometricTestOnly", geometricTestOnly_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::autoPtr<Foam::surfaceFeatures>
Foam::surfaceFeaturesExtraction::extractFromFile::features
(
    const triSurface& surf
) const
{
    edgeMesh eMesh(featureEdgeFile_);

    // Sometimes duplicate edges are present. Remove them.
    eMesh.mergeEdges();

    Info<< nl << "Reading existing feature edges from file "
        << featureEdgeFile_ << nl
        << "Selecting edges based purely on geometric tests: "
        << geometricTestOnly().c_str() << endl;

    return autoPtr<surfaceFeatures>::New
    (
        surf,
        eMesh.points(),
        eMesh.edges(),
        1e-6,  // mergeTol
        geometricTestOnly()
    );
}


// ************************************************************************* //
