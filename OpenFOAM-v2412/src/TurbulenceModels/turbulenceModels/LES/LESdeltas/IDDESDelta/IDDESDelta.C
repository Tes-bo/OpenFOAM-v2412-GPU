/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2016-2020 OpenCFD Ltd.
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

#include "IDDESDelta.H"
#include "addToRunTimeSelectionTable.H"
#include "wallDist.H"
#include "maxDeltaxyz.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{
    defineTypeNameAndDebug(IDDESDelta, 0);
    addToRunTimeSelectionTable(LESdelta, IDDESDelta, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::LESModels::IDDESDelta::calcDelta()
{
    const volScalarField& hmax = hmaxPtr_();
    const fvMesh& mesh = turbulenceModel_.mesh();

    // Wall-normal vectors
    const volVectorField& n = wallDist::New(mesh).n();

    auto tfaceToFacenMax = volScalarField::New
    (
        "faceToFaceMax",
        IOobject::NO_REGISTER,
        mesh,
        dimensionedScalar(dimLength, Zero)
    );

    scalarField& faceToFacenMax = tfaceToFacenMax.ref().primitiveFieldRef();

    const cellList& cells = mesh.cells();
    const vectorField& faceCentres = mesh.faceCentres();

    forAll(cells, celli)
    {
        scalar maxDelta = 0.0;
        const labelList& cFaces = cells[celli];
        const vector nci = n[celli];

        forAll(cFaces, cFacei)
        {
            label facei = cFaces[cFacei];
            const point& fci = faceCentres[facei];

            forAll(cFaces, cFacej)
            {
                label facej = cFaces[cFacej];
                const point& fcj = faceCentres[facej];
                scalar ndfc = nci & (fcj - fci);

                if (ndfc > maxDelta)
                {
                    maxDelta = ndfc;
                }
            }
        }

        faceToFacenMax[celli] = maxDelta;
    }


    label nD = mesh.nGeometricD();

    if (nD == 2)
    {
        WarningInFunction
            << "Case is 2D, LES is not strictly applicable" << nl
            << endl;
    }
    else if (nD != 3)
    {
        FatalErrorInFunction
            << "Case must be either 2D or 3D" << exit(FatalError);
    }

    delta_.primitiveFieldRef() =
        min
        (
            max
            (
                max
                (
                    Cw_*wallDist::New(mesh).y(),
                    Cw_*hmax
                ),
                tfaceToFacenMax
            ),
            hmax
        );

    // Handle coupled boundaries
    delta_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::LESModels::IDDESDelta::IDDESDelta
(
    const word& name,
    const turbulenceModel& turbulence,
    const dictionary& dict
)
:
    LESdelta(name, turbulence),
    hmaxPtr_(nullptr),
    Cw_
    (
        dict.optionalSubDict(type() + "Coeffs").getOrDefault<scalar>
        (
            "Cw",
            0.15
        )
    )
{
    if (dict.optionalSubDict(type() + "Coeffs").found("hmax"))
    {
        // User-defined hmax
        hmaxPtr_ =
            LESdelta::New
            (
                IOobject::groupName("hmax", turbulence.U().group()),
                turbulence,
                dict.optionalSubDict(type() + "Coeffs"),
                "hmax"
            );
    }
    else
    {
        Info<< "Employing " << maxDeltaxyz::typeName << " for hmax" << endl;

        hmaxPtr_.reset
        (
            new maxDeltaxyz
            (
                IOobject::groupName("hmax", turbulence.U().group()),
                turbulence,
                dict.optionalSubDict(type() + "Coeffs")
            )
        );
    }

    calcDelta();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::LESModels::IDDESDelta::read(const dictionary& dict)
{
    const dictionary& coeffsDict(dict.optionalSubDict(type() + "Coeffs"));

    coeffsDict.readIfPresent<scalar>("Cw", Cw_);

    calcDelta();
}


void Foam::LESModels::IDDESDelta::correct()
{
    if (turbulenceModel_.mesh().changing())
    {
        hmaxPtr_->correct();
        calcDelta();
    }
}


// ************************************************************************* //
