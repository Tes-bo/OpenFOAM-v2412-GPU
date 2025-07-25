/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2024 OpenCFD Ltd.
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

#include "viewFactorHeatFlux.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(viewFactorHeatFlux, 0);
    addToRunTimeSelectionTable(functionObject, viewFactorHeatFlux, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::viewFactorHeatFlux::initialise()
{
    const auto& pbm = mesh_.boundaryMesh();
    const labelList selectedPatches = pbm.indices("viewFactorWall");


    // Initialise output file

    auto& os = file();

    writeHeader(os, "Radiation heat flux");
    writeCommented(os, "Time");

    for (const label patchi : selectedPatches)
    {
        writeTabbed(os, pbm[patchi].name());
    }

    os  << endl;


    // Create compact patch addressing

    label nFace = 0;
    forAll(selectedPatches, i)
    {
        const label patchi = selectedPatches[i];
        auto slice = compactPatchID_.slice(nFace, pbm[patchi].size());
        slice = i;
        nFace += slice.size();
    }


    if (Pstream::parRun())
    {
        mapDistPtr_.reset
        (
            new IOmapDistribute
            (
                IOobject
                (
                    "mapDist",
                    mesh_.facesInstance(),
                    mesh_,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE,
                    IOobject::NO_REGISTER
                )
            )
        );

        auto& mapDist = mapDistPtr_();

        mapDist.distribute(compactPatchID_);


        // Convert global addressing into compact form

        labelList compactGlobalIds(mapDist.constructSize(), Zero);

        globalIndex gi(nFace);

        SubList<label>(compactGlobalIds, nFace) = identity(gi.range());

        mapDist.distribute(compactGlobalIds);

        const label nTotalFaces = returnReduce(nFace, sumOp<label>());

        const labelList globalToCompact(invert(nTotalFaces, compactGlobalIds));

        for (labelList& visibleFaces : faceFaces_)
        {
            for (label& globalFacei : visibleFaces)
            {
                globalFacei = globalToCompact[globalFacei];
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::viewFactorHeatFlux::viewFactorHeatFlux
(
    const word& name,
    const Time& runTime,
    const dictionary& dict,
    const bool readFields
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeFile(mesh_, name, typeName, dict),
    qrName_(dict.getOrDefault<word>("qr", "qr")),
    mapDistPtr_(nullptr),
    faceFaces_
    (
        IOobject
        (
            "globalFaceFaces",
            mesh_.facesInstance(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        )
    ),
    Fij_
    (
        IOobject
        (
            "F",
            mesh_.facesInstance(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        )
    ),
    compactPatchID_(faceFaces_.size(), Zero)
{
    initialise();
}


Foam::functionObjects::viewFactorHeatFlux::viewFactorHeatFlux
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool readFields
)
:
    fvMeshFunctionObject(name, obr, dict),
    writeFile(mesh_, name, typeName, dict),
    qrName_(dict.getOrDefault<word>("qr", "qr")),
    mapDistPtr_(nullptr),
    faceFaces_
    (
        IOobject
        (
            "globalFaceFaces",
            mesh_.facesInstance(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        )
    ),
    Fij_
    (
        IOobject
        (
            "F",
            mesh_.facesInstance(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        )
    ),
    compactPatchID_(faceFaces_.size(), Zero)
{
   initialise();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::viewFactorHeatFlux::read(const dictionary& dict)
{
    if (fvMeshFunctionObject::read(dict))
    {
        dict.readIfPresent("qr", qrName_);

        return true;
    }

    return false;
}


bool Foam::functionObjects::viewFactorHeatFlux::execute()
{
    return true;
}


bool Foam::functionObjects::viewFactorHeatFlux::write()
{
    Log << type() << " " << name() <<  " write:" << nl;

    // Retrieve radiative heat flux field from the mesh database

    const auto* qrPtr = mesh_.cfindObject<volScalarField>(qrName_);

    if (qrPtr)
    {
        const auto& qr = *qrPtr;

        // const labelList selectedPatches =
        //     mesh_.boundaryMesh().indices
        //     (
        //         radiation::viewFactor::viewFactorWalls
        //     );
        const labelList selectedPatches =
            mesh_.boundaryMesh().indices("viewFactorWall");

        const auto& pbm = mesh_.boundaryMesh();
        const label nPatch = selectedPatches.size();

        // Accumulate qr per patch, and
        // assemble qr contributions in compact addressing
        scalarList qrPatch(nPatch, Zero);
        scalarList compactQr(faceFaces_.size());

        label compacti = 0;
        forAll(selectedPatches, i)
        {
            const label patchi = selectedPatches[i];
            const auto& qrp = qr.boundaryField()[patchi];
            forAll(qrp, facei)
            {
                compactQr[compacti] = qrp[facei];
                qrPatch[i] += qrp[facei];
                ++compacti;
            }
        }

        reduce(qrPatch, sumOp<scalarList>());


        if (Pstream::parRun())
        {
            mapDistPtr_->distribute(compactQr);
        }

        // Populate view factor radiation heat flux matrix
        scalarSquareMatrix qrVf(nPatch, Zero);

        forAll(Fij_, startFacei)
        {
            const scalar qr = compactQr[startFacei];
            const labelList& visibleSlots = faceFaces_[startFacei];
            const label i = compactPatchID_[startFacei];

            forAll(visibleSlots, visSloti)
            {
                const label sloti = visibleSlots[visSloti];
                const label j = compactPatchID_[sloti];

                qrVf[i][j] += Fij_[startFacei][visSloti]*qr;
            }
        }

        reduce(qrVf, sumOp<scalarSquareMatrix>());

        if (Pstream::master())
        {
            Log << "     Writing patch totals to " << file().name()
                << nl;

            // Write patch sums as time history
            writeCurrentTime(file());

            for (const auto& qrp : qrPatch)
            {
                file() << tab << qrp;
            }

            file() << endl;


            // Write view factor breakdown as matrix per write time
            auto osPtr = newFileAtTime("viewFactorQr", mesh_.time().value());
            auto& os = osPtr();

            Log << "     Writing view factor breakdown to " << os.name()
                << nl;


            forAll(selectedPatches, i)
            {
                if (i == 0)
                {
                    for (const label patchj : selectedPatches)
                    {
                        os  << tab << pbm[patchj].name();
                    }
                    os  << nl;
                }

                const label patchi = selectedPatches[i];

                os  << pbm[patchi].name();;

                forAll(selectedPatches, j)
                {
                    os  << tab << qrVf[i][j];
                }

                os  << nl;
            }

            os  << endl;
        }
    }

    Log << endl;

    return true;
}


void Foam::functionObjects::viewFactorHeatFlux::updateMesh(const mapPolyMesh&)
{
    NotImplemented;
}


void Foam::functionObjects::viewFactorHeatFlux::movePoints(const polyMesh&)
{
    NotImplemented;
}


// ************************************************************************* //
