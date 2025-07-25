/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2017-2023 OpenCFD Ltd.
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

#include "manualDecomp.H"
#include "labelIOList.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(manualDecomp, 0);
    addToRunTimeSelectionTable
    (
        decompositionMethod,
        manualDecomp,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::manualDecomp::manualDecomp
(
    const dictionary& decompDict,
    const word& regionName
)
:
    decompositionMethod(decompDict, regionName),
    dataFile_
    (
        findCoeffsDict(typeName + "Coeffs").get<fileName>("dataFile")
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::labelList Foam::manualDecomp::decompose
(
    const polyMesh& mesh,
    const pointField& points,
    const scalarField& pointWeights
) const
{
    labelList finalDecomp
    (
        labelIOList::readContents
        (
            IOobject
            (
                dataFile_,
                mesh.facesInstance(),
                mesh.thisDb(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            )
        )
    );

    // Check if the final decomposition is OK
    if (finalDecomp.size() != points.size())
    {
        FatalErrorInFunction
            << "Size of decomposition list does not correspond "
            << "to the number of points.  Size: "
            << finalDecomp.size() << " Number of points: "
            << points.size()
            << ".\n" << "Manual decomposition data read from file "
            << dataFile_ << "." << endl
            << exit(FatalError);
    }


    const MinMax<label> range = minMax(finalDecomp);

    if (range.min() < 0 || range.max() >= nDomains_)
    {
        FatalErrorInFunction
            << "According to the decomposition, cells assigned to "
            << "impossible processor numbers.  Min processor = "
            << range.min() << " Max processor = " << range.max()
            << ".\n" << "Manual decomposition data read from file "
            << dataFile_ << "." << endl
            << exit(FatalError);
    }

    return finalDecomp;
}


// ************************************************************************* //
