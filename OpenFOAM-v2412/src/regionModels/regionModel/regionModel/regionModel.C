/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2016-2022 OpenCFD Ltd.
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

#include "regionModel.H"
#include "fvMesh.H"
#include "Time.H"
#include "mappedWallPolyPatch.H"
#include "zeroGradientFvPatchFields.H"
#include "faceAreaWeightAMI.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
    defineTypeNameAndDebug(regionModel, 0);
}
}

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::regionModels::regionModel::constructMeshObjects()
{
    fvMesh* regionMeshPtr = time_.getObjectPtr<fvMesh>(regionName_);

    if (!regionMeshPtr)
    {
        regionMeshPtr = new fvMesh
        (
            IOobject
            (
                regionName_,
                time_.timeName(),
                time_,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                IOobject::REGISTER
            )
        );

        regionMeshPtr->objectRegistry::store();
    }
}


void Foam::regionModels::regionModel::initialise()
{
    if (debug)
    {
        Pout<< "regionModel::initialise()" << endl;
    }

    label nBoundaryFaces = 0;
    DynamicList<label> primaryPatchIDs;
    DynamicList<label> intCoupledPatchIDs;
    const polyBoundaryMesh& rbm = regionMesh().boundaryMesh();

    forAll(rbm, patchi)
    {
        const polyPatch& regionPatch = rbm[patchi];
        if (isA<mappedPatchBase>(regionPatch))
        {
            if (debug)
            {
                Pout<< "found " << mappedWallPolyPatch::typeName
                    <<  " " << regionPatch.name() << endl;
            }

            intCoupledPatchIDs.append(patchi);

            nBoundaryFaces += regionPatch.faceCells().size();

            const mappedPatchBase& mapPatch =
                refCast<const mappedPatchBase>(regionPatch);

            if
            (
                primaryMesh_.time().foundObject<polyMesh>
                (
                    mapPatch.sampleRegion()
                )
            )
            {

                const label primaryPatchi = mapPatch.samplePolyPatch().index();
                primaryPatchIDs.append(primaryPatchi);
            }
        }
    }

    primaryPatchIDs_.transfer(primaryPatchIDs);
    intCoupledPatchIDs_.transfer(intCoupledPatchIDs);

    if (!returnReduceOr(nBoundaryFaces))
    {
        WarningInFunction
            << "Region model has no mapped boundary conditions - transfer "
            << "between regions will not be possible" << endl;
    }

    if (!outputPropertiesPtr_)
    {
        const fileName uniformPath(word("uniform")/"regionModels");

        outputPropertiesPtr_.reset
        (
            new IOdictionary
            (
                IOobject
                (
                    regionName_ + "OutputProperties",
                    time_.timeName(),
                    uniformPath/regionName_,
                    primaryMesh_,
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE
                )
            )
        );
    }
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

bool Foam::regionModels::regionModel::read()
{
    if (regIOobject::read())
    {
        if (active_)
        {
            if (const dictionary* dictptr = findDict(modelName_ + "Coeffs"))
            {
                coeffs_ <<= *dictptr;
            }

            infoOutput_.readIfPresent("infoOutput", *this);
        }

        return true;
    }

    return false;
}


bool Foam::regionModels::regionModel::read(const dictionary& dict)
{
    if (active_)
    {
        if (const dictionary* dictptr = dict.findDict(modelName_ + "Coeffs"))
        {
            coeffs_ <<= *dictptr;
        }

        infoOutput_.readIfPresent("infoOutput", dict);
        return true;
    }

    return false;
}


const Foam::AMIPatchToPatchInterpolation&
Foam::regionModels::regionModel::interRegionAMI
(
    const regionModel& nbrRegion,
    const label regionPatchi,
    const label nbrPatchi,
    const bool flip
) const
{
    label nbrRegionID = interRegionAMINames_.find(nbrRegion.name());

    const fvMesh& nbrRegionMesh = nbrRegion.regionMesh();

    if (nbrRegionID != -1)
    {
        if (!interRegionAMI_[nbrRegionID].set(regionPatchi))
        {
            const polyPatch& p = regionMesh().boundaryMesh()[regionPatchi];
            const polyPatch& nbrP = nbrRegionMesh.boundaryMesh()[nbrPatchi];

            const int oldTag = UPstream::incrMsgType();

            interRegionAMI_[nbrRegionID].set
            (
                regionPatchi,
                AMIInterpolation::New
                (
                    faceAreaWeightAMI::typeName,
                    true, // requireMatch
                    flip
                )
            );

            interRegionAMI_[nbrRegionID][regionPatchi].calculate(p, nbrP);

            UPstream::msgType(oldTag);  // Restore tag
        }

        return interRegionAMI_[nbrRegionID][regionPatchi];
    }
    else
    {
        label nbrRegionID = interRegionAMINames_.size();

        interRegionAMINames_.append(nbrRegion.name());

        const polyPatch& p = regionMesh().boundaryMesh()[regionPatchi];
        const polyPatch& nbrP = nbrRegionMesh.boundaryMesh()[nbrPatchi];

        const label nPatch = regionMesh().boundaryMesh().size();


        interRegionAMI_.resize(nbrRegionID + 1);

        interRegionAMI_.set
        (
            nbrRegionID,
            new PtrList<AMIPatchToPatchInterpolation>(nPatch)
        );

        const int oldTag = UPstream::incrMsgType();

        interRegionAMI_[nbrRegionID].set
        (
            regionPatchi,
            AMIInterpolation::New
            (
                faceAreaWeightAMI::typeName,
                true, // requireMatch
                flip // reverse
            )
        );

        interRegionAMI_[nbrRegionID][regionPatchi].calculate(p, nbrP);

        UPstream::msgType(oldTag);  // Restore tag

        return interRegionAMI_[nbrRegionID][regionPatchi];
    }
}


Foam::label Foam::regionModels::regionModel::nbrCoupledPatchID
(
    const regionModel& nbrRegion,
    const label regionPatchi
) const
{
    label nbrPatchi = -1;

    // region
    const fvMesh& nbrRegionMesh = nbrRegion.regionMesh();

    // boundary mesh
    const polyBoundaryMesh& nbrPbm = nbrRegionMesh.boundaryMesh();

    const polyBoundaryMesh& pbm = regionMesh().boundaryMesh();

    if (regionPatchi > pbm.size() - 1)
    {
        FatalErrorInFunction
            << "region patch index out of bounds: "
            << "region patch index = " << regionPatchi
            << ", maximum index = " << pbm.size() - 1
            << abort(FatalError);
    }

    const polyPatch& pp = regionMesh().boundaryMesh()[regionPatchi];

    if (!isA<mappedPatchBase>(pp))
    {
        FatalErrorInFunction
            << "Expected a " << mappedPatchBase::typeName
            << " patch, but found a " << pp.type() << abort(FatalError);
    }

    const mappedPatchBase& mpb = refCast<const mappedPatchBase>(pp);

    // sample patch name on the primary region
    const word& primaryPatchName = mpb.samplePatch();

    // find patch on nbr region that has the same sample patch name
    forAll(nbrRegion.intCoupledPatchIDs(), j)
    {
        const label nbrRegionPatchi = nbrRegion.intCoupledPatchIDs()[j];

        const mappedPatchBase& mpb =
            refCast<const mappedPatchBase>(nbrPbm[nbrRegionPatchi]);

        if (mpb.samplePatch() == primaryPatchName)
        {
            nbrPatchi = nbrRegionPatchi;
            break;
        }
    }

    if (nbrPatchi == -1)
    {
        const polyPatch& p = regionMesh().boundaryMesh()[regionPatchi];

        FatalErrorInFunction
            << "Unable to find patch pair for local patch "
            << p.name() << " and region " << nbrRegion.name()
            << abort(FatalError);
    }

    return nbrPatchi;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionModels::regionModel::regionModel
(
    const fvMesh& mesh,
    const word& regionType
)
:
    IOdictionary
    (
        IOobject
        (
            regionType + "Properties",
            mesh.time().constant(),
            mesh.time(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    ),
    primaryMesh_(mesh),
    time_(mesh.time()),
    active_(false),
    infoOutput_(false),
    modelName_("none"),
    coeffs_(),
    outputPropertiesPtr_(nullptr),
    primaryPatchIDs_(),
    intCoupledPatchIDs_(),
    regionName_("none"),
    functions_(*this),
    interRegionAMINames_(),
    interRegionAMI_()
{}


Foam::regionModels::regionModel::regionModel
(
    const fvMesh& mesh,
    const word& regionType,
    const word& modelName,
    bool readFields
)
:
    IOdictionary
    (
        IOobject
        (
            regionType + "Properties",
            mesh.time().constant(),
            mesh.time(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    primaryMesh_(mesh),
    time_(mesh.time()),
    active_(get<Switch>("active")),
    infoOutput_(true),
    modelName_(modelName),
    coeffs_(subOrEmptyDict(modelName + "Coeffs")),
    outputPropertiesPtr_(nullptr),
    primaryPatchIDs_(),
    intCoupledPatchIDs_(),
    regionName_(lookup("region")),
    functions_(*this, subOrEmptyDict("functions"))
{
    if (active_)
    {
        constructMeshObjects();
        initialise();

        if (readFields)
        {
            read();
        }
    }
}


Foam::regionModels::regionModel::regionModel
(
    const fvMesh& mesh,
    const word& regionType,
    const word& modelName,
    const dictionary& dict,
    bool readFields
)
:
    IOdictionary
    (
        IOobject
        (
            regionType + "Properties",
            mesh.time().constant(),
            mesh.time(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::REGISTER
        ),
        dict
    ),
    primaryMesh_(mesh),
    time_(mesh.time()),
    active_(dict.get<Switch>("active")),
    infoOutput_(false),
    modelName_(modelName),
    coeffs_(dict.subOrEmptyDict(modelName + "Coeffs")),
    outputPropertiesPtr_(nullptr),
    primaryPatchIDs_(),
    intCoupledPatchIDs_(),
    regionName_(dict.lookup("region")),
    functions_(*this, subOrEmptyDict("functions"))
{
    if (active_)
    {
        constructMeshObjects();
        initialise();

        if (readFields)
        {
            read(dict);
        }
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::regionModels::regionModel::evolve()
{
    if (active_)
    {
        Info<< "\nEvolving " << modelName_ << " for region "
            << regionMesh().name() << endl;

        preEvolveRegion();

        evolveRegion();

        postEvolveRegion();

        // Provide some feedback
        if (infoOutput_)
        {
            Info<< incrIndent;
            info();
            Info<< endl << decrIndent;
        }

        if (time_.writeTime())
        {
            outputProperties().writeObject
            (
                IOstreamOption(IOstreamOption::ASCII, time_.writeCompression()),
                true
            );
        }
    }
}


void Foam::regionModels::regionModel::preEvolveRegion()
{
    functions_.preEvolveRegion();
}


void Foam::regionModels::regionModel::evolveRegion()
{}


void Foam::regionModels::regionModel::postEvolveRegion()
{
    functions_.postEvolveRegion();
}


void Foam::regionModels::regionModel::info()
{}


// ************************************************************************* //
