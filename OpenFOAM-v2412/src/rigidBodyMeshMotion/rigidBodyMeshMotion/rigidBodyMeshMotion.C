/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 OpenFOAM Foundation
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

#include "rigidBodyMeshMotion.H"
#include "addToRunTimeSelectionTable.H"
#include "polyMesh.H"
#include "pointPatchDist.H"
#include "pointConstraints.H"
#include "uniformDimensionedFields.H"
#include "forces.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(rigidBodyMeshMotion, 0);

    addToRunTimeSelectionTable
    (
        motionSolver,
        rigidBodyMeshMotion,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::rigidBodyMeshMotion::bodyMesh::bodyMesh
(
    const polyMesh& mesh,
    const word& name,
    const label bodyID,
    const dictionary& dict
)
:
    name_(name),
    bodyID_(bodyID),
    patches_(dict.get<wordRes>("patches")),
    patchSet_(mesh.boundaryMesh().patchSet(patches_)),
    di_(dict.get<scalar>("innerDistance")),
    do_(dict.get<scalar>("outerDistance")),
    weight_
    (
        IOobject
        (
            name_ + ".motionScale",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        pointMesh::New(mesh),
        dimensionedScalar(dimless, Zero)
    )
{}


Foam::rigidBodyMeshMotion::rigidBodyMeshMotion
(
    const polyMesh& mesh,
    const IOdictionary& dict
)
:
    displacementMotionSolver(mesh, dict, typeName),
    model_
    (
        mesh.time(),
        coeffDict(),
        IOobject
        (
            "rigidBodyMotionState",
            mesh.time().timeName(),
            "uniform",
            mesh
        ).typeHeaderOk<IOdictionary>(true)
      ? IOdictionary
        (
            IOobject
            (
                "rigidBodyMotionState",
                mesh.time().timeName(),
                "uniform",
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            )
        )
      : coeffDict()
    ),
    test_(coeffDict().getOrDefault("test", false)),
    rhoInf_(1.0),
    rhoName_(coeffDict().getOrDefault<word>("rho", "rho")),
    ramp_
    (
        Function1<scalar>::NewIfPresent("ramp", coeffDict(), word::null, &mesh)
    ),
    curTimeIndex_(-1),
    cOfGdisplacement_
    (
        coeffDict().getOrDefault<word>("cOfGdisplacement", "none")
    ),
    bodyIdCofG_(coeffDict().getOrDefault<label>("bodyIdCofG", -1))
{
    if (rhoName_ == "rhoInf")
    {
        readEntry("rhoInf", rhoInf_);
    }

    const dictionary& bodiesDict = coeffDict().subDict("bodies");

    for (const entry& dEntry : bodiesDict)
    {
        const keyType& bodyName = dEntry.keyword();
        const dictionary& bodyDict = dEntry.dict();

        if (bodyDict.found("patches"))
        {
            const label bodyID = model_.bodyID(bodyName);

            if (bodyID == -1)
            {
                FatalErrorInFunction
                    << "Body " << bodyName
                    << " has been merged with another body"
                       " and cannot be assigned a set of patches"
                    << exit(FatalError);
            }

            bodyMeshes_.append
            (
                new bodyMesh
                (
                    mesh,
                    bodyName,
                    bodyID,
                    bodyDict
                )
            );
        }
    }

    // Calculate scaling factor everywhere for each meshed body
    forAll(bodyMeshes_, bi)
    {
        const pointMesh& pMesh = pointMesh::New(mesh);

        pointPatchDist pDist(pMesh, bodyMeshes_[bi].patchSet_, points0());

        pointScalarField& scale = bodyMeshes_[bi].weight_;

        // Scaling: 1 up to di then linear down to 0 at do away from patches
        scale.primitiveFieldRef() =
            min
            (
                max
                (
                    (bodyMeshes_[bi].do_ - pDist.primitiveField())
                   /(bodyMeshes_[bi].do_ - bodyMeshes_[bi].di_),
                    scalar(0)
                ),
                scalar(1)
            );

        // Convert the scale function to a cosine
        scale.primitiveFieldRef() =
            min
            (
                max
                (
                    0.5
                  - 0.5
                   *cos(scale.primitiveField()
                   *Foam::constant::mathematical::pi),
                    scalar(0)
                ),
                scalar(1)
            );

        pointConstraints::New(pMesh).constrain(scale);
        //scale.write();
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::pointField>
Foam::rigidBodyMeshMotion::curPoints() const
{
    tmp<pointField> newPoints(points0() + pointDisplacement_.primitiveField());

    if (moveAllCells())
    {
        return newPoints;
    }
    else
    {
        auto ttransformedPts = tmp<pointField>::New(mesh().points());
        auto& transformedPts = ttransformedPts.ref();

        UIndirectList<point>(transformedPts, pointIDs()) =
            pointField(newPoints.ref(), pointIDs());

        return ttransformedPts;
    }
}


void Foam::rigidBodyMeshMotion::solve()
{
    const Time& t = mesh().time();

    if (mesh().nPoints() != points0().size())
    {
        FatalErrorInFunction
            << "The number of points in the mesh seems to have changed." << endl
            << "In constant/polyMesh there are " << points0().size()
            << " points; in the current mesh there are " << mesh().nPoints()
            << " points." << exit(FatalError);
    }

    // Store the motion state at the beginning of the time-step
    if (curTimeIndex_ != this->db().time().timeIndex())
    {
        model_.newTime();
        curTimeIndex_ = this->db().time().timeIndex();
    }

    const scalar ramp = (ramp_ ? ramp_->value(t.value()) : 1.0);

    if (t.foundObject<uniformDimensionedVectorField>("g"))
    {
        model_.g() =
            ramp*t.lookupObject<uniformDimensionedVectorField>("g").value();
    }

    vector oldPos(vector::uniform(GREAT));
    if (bodyIdCofG_ != -1)
    {
        oldPos = model_.cCofR(bodyIdCofG_);
    }

    if (test_)
    {
        const label nIter(coeffDict().get<label>("nIter"));

        for (label i=0; i<nIter; i++)
        {
            model_.solve
            (
                t.value(),
                t.deltaTValue(),
                scalarField(model_.nDoF(), Zero),
                Field<spatialVector>(model_.nBodies(), Zero)
            );
        }
    }
    else
    {
        const label nIter(coeffDict().getOrDefault<label>("nIter", 1));

        for (label i=0; i<nIter; i++)
        {
            Field<spatialVector> fx(model_.nBodies(), Zero);

            forAll(bodyMeshes_, bi)
            {
                const label bodyID = bodyMeshes_[bi].bodyID_;

                dictionary forcesDict;
                forcesDict.add("type", functionObjects::forces::typeName);
                forcesDict.add("patches", bodyMeshes_[bi].patches_);
                forcesDict.add("rhoInf", rhoInf_);
                forcesDict.add("rho", rhoName_);
                forcesDict.add("CofR", vector::zero);

                functionObjects::forces f("forces", db(), forcesDict);
                f.calcForcesMoments();

                fx[bodyID] = ramp*spatialVector(f.momentEff(), f.forceEff());
            }

            model_.solve
            (
                t.value(),
                t.deltaTValue(),
                scalarField(model_.nDoF(), Zero),
                fx
            );
        }

        if (cOfGdisplacement_ != "none")
        {
            if (bodyIdCofG_ != -1)
            {
                if
                (
                    db().time().foundObject<uniformDimensionedVectorField>
                    (
                        cOfGdisplacement_
                    )
                )
                {
                    auto& disp =
                        db().time().lookupObjectRef<uniformDimensionedVectorField>
                        (
                            cOfGdisplacement_
                        );

                    disp.value() += model_.cCofR(bodyIdCofG_) - oldPos;
                }
            }
            else
            {
                FatalErrorInFunction
                    << "CofGdisplacement is different to none." << endl
                    << "The model needs the entry body reference Id: bodyIdCofG."
                    << exit(FatalError);
            }
        }
    }

    if (Pstream::master() && model_.report())
    {
        forAll(bodyMeshes_, bi)
        {
            model_.status(bodyMeshes_[bi].bodyID_);
        }
    }

    // Update the displacements
    if (bodyMeshes_.size() == 1)
    {
        pointDisplacement_.primitiveFieldRef() = model_.transformPoints
        (
            bodyMeshes_[0].bodyID_,
            bodyMeshes_[0].weight_,
            points0()
        ) - points0();
    }
    else
    {
        labelList bodyIDs(bodyMeshes_.size());
        List<const scalarField*> weights(bodyMeshes_.size());
        forAll(bodyIDs, bi)
        {
            bodyIDs[bi] = bodyMeshes_[bi].bodyID_;
            weights[bi] = &bodyMeshes_[bi].weight_;
        }

         pointDisplacement_.primitiveFieldRef() =
             model_.transformPoints(bodyIDs, weights, points0()) - points0();

    }
    // Displacement has changed. Update boundary conditions
    pointConstraints::New
    (
        pointDisplacement_.mesh()
    ).constrainDisplacement(pointDisplacement_);
}


bool Foam::rigidBodyMeshMotion::writeObject
(
    IOstreamOption streamOpt,
    const bool writeOnProc
) const
{
    // Force ASCII writing
    streamOpt.format(IOstreamOption::ASCII);

    IOdictionary dict
    (
        IOobject
        (
            "rigidBodyMotionState",
            mesh().time().timeName(),
            "uniform",
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        )
    );

    model_.state().write(dict);
    return dict.regIOobject::writeObject(streamOpt, writeOnProc);
}


bool Foam::rigidBodyMeshMotion::read()
{
    if (displacementMotionSolver::read())
    {
        model_.read(coeffDict());

        return true;
    }

    return false;
}


// ************************************************************************* //
