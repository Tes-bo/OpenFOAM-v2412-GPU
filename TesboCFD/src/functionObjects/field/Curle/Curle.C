/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2022 OpenCFD Ltd.
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

#include "Curle.H"
#include "fvcDdt.H"
#include "mathematicalConstants.H"
#include "addToRunTimeSelectionTable.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(Curle, 0);
    addToRunTimeSelectionTable(functionObject, Curle, dictionary);
}
}

const Foam::Enum<Foam::functionObjects::Curle::modeType>
Foam::functionObjects::Curle::modeTypeNames_
({
    { modeType::point, "point" },
    { modeType::surface, "surface" },
});


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::Curle::Curle
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeFile(mesh_, name),
    pName_("p"),
    observerPositions_(),
    c0_(0),
    rawFilePtrs_(),
    inputSurface_(),
    surfaceWriterPtr_(nullptr)
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::Curle::read(const dictionary& dict)
{
    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

    if (!(fvMeshFunctionObject::read(dict) && writeFile::read(dict)))
    {
        return false;
    }

    dict.readIfPresent("p", pName_);

    patchIDs_ = pbm.patchSet(dict.get<wordRes>("patches")).sortedToc();

    if (patchIDs_.empty())
    {
        WarningInFunction
            << "No patches defined" << endl;
        return false;
    }

    const modeType inputMode = modeTypeNames_.get("input", dict);

    switch (inputMode)
    {
        case modeType::point:
        {
            observerPositions_ = dict.get<List<point>>("observerPositions");
            break;
        }
        case modeType::surface:
        {
            const fileName fName = (dict.get<fileName>("surface")).expand();
            inputSurface_ = MeshedSurface<face>(fName);

            observerPositions_ = inputSurface_.Cf();
        }
    }

    if (observerPositions_.empty())
    {
        WarningInFunction
            << "No observer positions defined" << endl;

        return false;
    }

    const auto outputMode = modeTypeNames_.get("output", dict);

    switch (outputMode)
    {
        case modeType::point:
        {
            rawFilePtrs_.setSize(observerPositions_.size());
            forAll(rawFilePtrs_, filei)
            {
                rawFilePtrs_.set
                (
                    filei,
                    newFileAtStartTime("observer" + Foam::name(filei))
                );

                if (rawFilePtrs_.set(filei))
                {
                    OFstream& os = rawFilePtrs_[filei];

                    writeHeaderValue(os, "Observer", filei);
                    writeHeaderValue(os, "Position", observerPositions_[filei]);
                    writeCommented(os, "Time");
                    writeTabbed(os, "p(Curle)");
                    os  << endl;
                }
            }
            break;
        }
        case modeType::surface:
        {
            if (inputMode != modeType::surface)
            {
                FatalIOErrorInFunction(dict)
                    << "Surface output is only available when input mode is "
                    << "type '" << modeTypeNames_[modeType::surface] << "'"
                    << abort(FatalIOError);
            }

            const word writerType = dict.get<word>("surfaceType");

            surfaceWriterPtr_ = surfaceWriter::New
            (
                writerType,
                surfaceWriter::formatOptions(dict, writerType)
            );

            // Use outputDir/TIME/surface-name
            surfaceWriterPtr_->useTimeDir(true);
        }
    }

    // Read the reference speed of sound
    dict.readEntry("c0", c0_);

    if (c0_ < VSMALL)
    {
        FatalErrorInFunction
            << "Reference speed of sound = " << c0_
            << " cannot be negative or zero."
            << abort(FatalError);
    }

    return true;
}


bool Foam::functionObjects::Curle::execute()
{
    if (!foundObject<volScalarField>(pName_))
    {
        return false;
    }

    const volScalarField& p = lookupObject<volScalarField>(pName_);
    const volScalarField::Boundary& pBf = p.boundaryField();
    const volScalarField dpdt(scopedName("dpdt"), fvc::ddt(p));
    const volScalarField::Boundary& dpdtBf = dpdt.boundaryField();
    const surfaceVectorField::Boundary& SfBf = mesh_.Sf().boundaryField();
    const surfaceVectorField::Boundary& CfBf = mesh_.Cf().boundaryField();

    scalarField pDash(observerPositions_.size(), 0);

    for (const label patchi : patchIDs_)
    {
        const scalarField& pp = pBf[patchi];
        const scalarField& dpdtp = dpdtBf[patchi];
        const vectorField& Sfp = SfBf[patchi];
        const vectorField& Cfp = CfBf[patchi];

        forAll(observerPositions_, pointi)
        {
            const vectorField r(observerPositions_[pointi] - Cfp);
            const scalarField invMagR(1/(mag(r) + ROOTVSMALL));

            pDash[pointi] +=
                sum((pp*sqr(invMagR) + invMagR/c0_*dpdtp)*(Sfp & (r*invMagR)));
        }
    }

    pDash /= 4*mathematical::pi;

    Pstream::listCombineReduce(pDash, plusEqOp<scalar>());

    if (surfaceWriterPtr_)
    {
        if (UPstream::master())
        {
            // Time-aware, with time spliced into the output path
            surfaceWriterPtr_->beginTime(time_);

            surfaceWriterPtr_->open
            (
                inputSurface_.points(),
                inputSurface_.surfFaces(),
                (baseFileDir()/name()/"surface"),
                false  // serial - already merged
            );

            surfaceWriterPtr_->write("p(Curle)", pDash);

            surfaceWriterPtr_->endTime();
            surfaceWriterPtr_->clear();
        }
    }
    else
    {
        forAll(observerPositions_, pointi)
        {
            if (rawFilePtrs_.set(pointi))
            {
                OFstream& os = rawFilePtrs_[pointi];
                writeCurrentTime(os);
                os  << pDash[pointi] << endl;
            }
        }
    }

    return true;
}


bool Foam::functionObjects::Curle::write()
{
    return true;
}


// ************************************************************************* //
