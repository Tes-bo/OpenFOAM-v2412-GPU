/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2017 OpenFOAM Foundation
    Copyright (C) 2015-2021 OpenCFD Ltd.
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

#include "decompositionConstraint.H"
#include "syncTools.H"
#include "cyclicAMIPolyPatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(decompositionConstraint, 1);
    defineRunTimeSelectionTable(decompositionConstraint, dictionary);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::decompositionConstraint::getMinBoundaryValue
(
    const polyMesh& mesh,
    const labelList& decomposition,
    labelList& destProc
) const
{
    destProc.setSize(mesh.nBoundaryFaces());
    destProc = labelMax;

    const polyBoundaryMesh& pbm = mesh.boundaryMesh();

    for (const polyPatch& pp : pbm)
    {
        const labelUList& faceCells = pp.faceCells();

        forAll(faceCells, i)
        {
            label bFacei = pp.offset()+i;
            destProc[bFacei] = decomposition[faceCells[i]];
        }
    }

    // Take minimum of coupled faces (over all patches!)
    syncTools::syncBoundaryFaceList(mesh, destProc, minEqOp<label>());

    // Do cyclicAMI ourselves
    for (const auto& pp : pbm)
    {
        const auto* ppp = isA<cyclicAMIPolyPatch>(pp);
        if (ppp)
        {
            const auto& cycPp = *ppp;
            const auto& nbrPp = cycPp.neighbPatch();
            const labelList nbrDecomp(decomposition, nbrPp.faceCells());
            labelList thisDecomp(decomposition, cycPp.faceCells());

            if (cycPp.owner())
            {
                cycPp.AMI().interpolateToSource
                (
                    nbrDecomp,
                    []
                    (
                        label& res,
                        const label facei,
                        const label& fld,
                        const scalar& w
                    )
                    {
                        res = min(res, fld);
                    },
                    thisDecomp,
                    thisDecomp      // used in case of low-weight-corr
                );
            }
            else
            {
                nbrPp.AMI().interpolateToTarget
                (
                    nbrDecomp,
                    []
                    (
                        label& res,
                        const label facei,
                        const label& fld,
                        const scalar& w
                    )
                    {
                        res = min(res, fld);
                    },
                    thisDecomp,
                    thisDecomp      // used in case of low-weight-corr
                );
            }

            forAll(thisDecomp, i)
            {
                label& proc = destProc[cycPp.offset()+i];
                proc = min(proc, thisDecomp[i]);
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::decompositionConstraint::decompositionConstraint
(
    const dictionary& constraintDict
)
:
    coeffDict_(constraintDict)
{}


Foam::decompositionConstraint::decompositionConstraint
(
    const dictionary& constraintDict,
    const word&
)
:
    coeffDict_(constraintDict)
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::decompositionConstraint>
Foam::decompositionConstraint::New
(
    const dictionary& dict
)
{
    return decompositionConstraint::New
    (
        dict,
        dict.get<word>("type")
    );
}


Foam::autoPtr<Foam::decompositionConstraint>
Foam::decompositionConstraint::New
(
    const dictionary& dict,
    const word& modelType
)
{
    Info<< "Selecting decompositionConstraint " << modelType << endl;

    auto* ctorPtr = dictionaryConstructorTable(modelType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "decompositionConstraint",
            modelType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<decompositionConstraint>(ctorPtr(dict));
}


// ************************************************************************* //
