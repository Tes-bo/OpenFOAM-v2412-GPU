/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 Wikki Ltd
    Copyright (C) 2021-2022 OpenCFD Ltd.
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

#include "wedgeFaPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "faMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(wedgeFaPatch, 0);
    addToRunTimeSelectionTable(faPatch, wedgeFaPatch, dictionary);
}

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::wedgeFaPatch::findAxisPoint() const
{
    // Find axis point

    const labelList& ptLabels = pointLabels();

    const labelListList& ptEdges = pointEdges();

    const vectorField& points = boundaryMesh().mesh().points();

    const scalarField& magL = magEdgeLengths();

    forAll(ptEdges, pointI)
    {
        if (ptEdges[pointI].size() == 1)
        {
            scalar r = mag((I-axis()*axis())&points[ptLabels[pointI]]);

            if (r < magL[ptEdges[pointI][0]])
            {
                axisPoint_ = ptLabels[pointI];
                break;
            }
        }
    }

    axisPointChecked_ = true;
}


// * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * * * //

Foam::wedgeFaPatch::wedgeFaPatch
(
    const word& name,
    const dictionary& dict,
    const label index,
    const faBoundaryMesh& bm,
    const word& patchType
)
:
    faPatch(name, dict, index, bm, patchType),
    wedgePolyPatchPtr_(nullptr),
    axisPoint_(-1),
    axisPointChecked_(false)
{
    if (ngbPolyPatchIndex() < 0)
    {
        FatalErrorInFunction
            << "Neighbour polyPatch index is not specified for faPatch "
            << this->name() << exit(FatalError);
    }

    const auto* wedgePtr = isA<wedgePolyPatch>
    (
        bm.mesh().mesh().boundaryMesh()[ngbPolyPatchIndex()]
    );

    if (wedgePtr)
    {
        wedgePolyPatchPtr_ = wedgePtr;
    }
    else
    {
        FatalErrorInFunction
            << "Neighbour polyPatch is not of type "
            << wedgePolyPatch::typeName
            << exit(FatalError);
    }
}


// ************************************************************************* //
