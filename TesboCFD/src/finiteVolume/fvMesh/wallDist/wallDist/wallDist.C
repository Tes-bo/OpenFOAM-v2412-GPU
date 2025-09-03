/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2016 OpenFOAM Foundation
    Copyright (C) 2016-2021 OpenCFD Ltd.
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

#include "wallDist.H"
#include "wallPolyPatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(wallDist, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::wallDist::constructn() const
{
    n_.reset
    (
        new volVectorField
        (
            IOobject
            (
                "n" & patchTypeName_,
                mesh().time().timeName(),
                mesh().thisDb(),
                IOobjectOption::NO_REGISTER
            ),
            mesh(),
            dimensionedVector(dimless, Zero),
            patchDistMethod::patchTypes<vector>(mesh(), patchIDs_)
        )
    );

    const fvPatchList& patches = mesh().boundary();

    volVectorField::Boundary& nbf = n_.ref().boundaryFieldRef();

    for (const label patchi : patchIDs_)
    {
        nbf[patchi] == patches[patchi].nf();
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::wallDist::wallDist
(
    const fvMesh& mesh,
    const labelHashSet& patchIDs,
    const word& patchTypeName
)
:
    wallDist(mesh, word::null, patchIDs, patchTypeName)
{}


Foam::wallDist::wallDist
(
    const fvMesh& mesh,
    const word& defaultPatchDistMethod,
    const labelHashSet& patchIDs,
    const word& patchTypeName
)
:
    MeshObject_type(mesh),
    patchIDs_(patchIDs),
    patchTypeName_(patchTypeName),
    dict_
    (
        static_cast<const fvSchemes&>(mesh).subOrEmptyDict
        (
            patchTypeName_ & "Dist"
        )
    ),
    pdm_
    (
        patchDistMethod::New
        (
            dict_,
            mesh,
            patchIDs_,
            defaultPatchDistMethod
        )
    ),
    y_
    (
        IOobject
        (
            "y" & patchTypeName_,
            mesh.time().timeName(),
            mesh.thisDb(),
            IOobjectOption::NO_REGISTER
        ),
        mesh,
        dimensionedScalar("y" & patchTypeName_, dimLength, SMALL),
        patchDistMethod::patchTypes<scalar>(mesh, patchIDs_)
    ),
    n_(volVectorField::null()),
    updateInterval_(dict_.getOrDefault<label>("updateInterval", 1)),
    nRequired_(dict_.getOrDefault("nRequired", false)),
    requireUpdate_(true)
{
    if (nRequired_)
    {
        constructn();
    }

    movePoints();
}


Foam::wallDist::wallDist(const fvMesh& mesh, const word& patchTypeName)
:
    wallDist
    (
        mesh,
        mesh.boundaryMesh().findPatchIDs<wallPolyPatch>(),
        patchTypeName
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::wallDist::~wallDist()
{}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

bool Foam::wallDist::try_movePoints(const fvMesh& mesh)
{
    auto* ptr =
        mesh.getObjectPtr<UpdateableMeshObject<fvMesh>>("wallDist");

    if (ptr)
    {
        return ptr->movePoints();
    }

    return false;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::volVectorField& Foam::wallDist::n() const
{
    if (isNull(n_()))
    {
        WarningInFunction
            << "n requested but 'nRequired' not specified in the "
            << (patchTypeName_ & "Dist") << " dictionary" << nl
            << "    Recalculating y and n fields." << endl;

        nRequired_ = true;
        constructn();
        pdm_->correct(y_, n_.ref());
    }

    return n_();
}


bool Foam::wallDist::movePoints()
{
    if
    (
        (updateInterval_ > 0)
     && ((mesh_.time().timeIndex() % updateInterval_) == 0)
    )
    {
        requireUpdate_ = true;
    }

    if (requireUpdate_ && pdm_->movePoints())
    {
        DebugInfo<< "Updating wall distance" << endl;

        requireUpdate_ = false;

        if (nRequired_)
        {
            return pdm_->correct(y_, n_.ref());
        }
        else
        {
            return pdm_->correct(y_);
        }
    }

    return false;
}


void Foam::wallDist::updateMesh(const mapPolyMesh& mpm)
{
    pdm_->updateMesh(mpm);

    // Force update if performing topology change
    // Note: needed?
    // - field would have been mapped, so if using updateInterval option (!= 1)
    //   live with error associated of not updating and use mapped values?
    requireUpdate_ = true;
    movePoints();
}


// ************************************************************************* //
