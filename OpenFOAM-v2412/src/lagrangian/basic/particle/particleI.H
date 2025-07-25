/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2018, 2020 OpenFOAM Foundation
    Copyright (C) 2011-2022 OpenCFD Ltd.
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

#include "polyMesh.H"
#include "Time.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

inline Foam::barycentricTensor Foam::particle::stationaryTetTransform() const
{
    const tetPointRef tet = currentTetIndices().tet(mesh_);

    return barycentricTensor(tet.a(), tet.b(), tet.c(), tet.d());
}


inline void Foam::particle::movingTetGeometry
(
    const scalar fraction,
    Pair<vector>& centre,
    Pair<vector>& base,
    Pair<vector>& vertex1,
    Pair<vector>& vertex2
) const
{
    const triFace triIs(currentTetIndices().faceTriIs(mesh_));

    const pointField& ptsOld = mesh_.oldPoints();
    const pointField& ptsNew = mesh_.points();

    // !!! <-- We would be better off using mesh_.cellCentres() here. However,
    // we need to put a mesh_.oldCellCentres() method in for this to work. The
    // values obtained from the mesh and those obtained from the cell do not
    // necessarily match. See mantis #1993.
    //const vector ccOld = mesh_.cells()[celli_].centre(ptsOld, mesh_.faces());
    //const vector ccNew = mesh_.cells()[celli_].centre(ptsNew, mesh_.faces());

    const vector ccOld = mesh_.oldCellCentres()[celli_];
    const vector ccNew = mesh_.cellCentres()[celli_];

    // Old and new points and cell centres are not sub-cycled. If we are sub-
    // cycling, then we have to account for the timestep change here by
    // modifying the fractions that we take of the old and new geometry.
    const Pair<scalar> s = stepFractionSpan();
    const scalar f0 = s[0] + stepFraction_*s[1], f1 = fraction*s[1];

    centre[0] = ccOld + f0*(ccNew - ccOld);
    base[0] = ptsOld[triIs[0]] + f0*(ptsNew[triIs[0]] - ptsOld[triIs[0]]);
    vertex1[0] = ptsOld[triIs[1]] + f0*(ptsNew[triIs[1]] - ptsOld[triIs[1]]);
    vertex2[0] = ptsOld[triIs[2]] + f0*(ptsNew[triIs[2]] - ptsOld[triIs[2]]);

    centre[1] = f1*(ccNew - ccOld);
    base[1] = f1*(ptsNew[triIs[0]] - ptsOld[triIs[0]]);
    vertex1[1] = f1*(ptsNew[triIs[1]] - ptsOld[triIs[1]]);
    vertex2[1] = f1*(ptsNew[triIs[2]] - ptsOld[triIs[2]]);
}


inline Foam::Pair<Foam::barycentricTensor> Foam::particle::movingTetTransform
(
    const scalar fraction
) const
{
    Pair<vector> centre, base, vertex1, vertex2;
    movingTetGeometry(fraction, centre, base, vertex1, vertex2);

    return
        Pair<barycentricTensor>
        (
            barycentricTensor(centre[0], base[0], vertex1[0], vertex2[0]),
            barycentricTensor(centre[1], base[1], vertex1[1], vertex2[1])
        );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::label Foam::particle::getNewParticleID() const
{
    label id = particleCount_++;

    if (id == labelMax)
    {
        WarningInFunction
            << "Particle counter has overflowed. This might cause problems"
            << " when reconstructing particle tracks." << endl;
    }
    return id;
}


inline const Foam::polyMesh& Foam::particle::mesh() const noexcept
{
    return mesh_;
}


inline const Foam::barycentric& Foam::particle::coordinates() const noexcept
{
    return coordinates_;
}


inline Foam::label Foam::particle::cell() const noexcept
{
    return celli_;
}


inline Foam::label& Foam::particle::cell() noexcept
{
    return celli_;
}


inline Foam::label Foam::particle::tetFace() const noexcept
{
    return tetFacei_;
}


inline Foam::label& Foam::particle::tetFace() noexcept
{
    return tetFacei_;
}


inline Foam::label Foam::particle::tetPt() const noexcept
{
    return tetPti_;
}


inline Foam::label& Foam::particle::tetPt() noexcept
{
    return tetPti_;
}


inline Foam::label Foam::particle::face() const noexcept
{
    return facei_;
}


inline Foam::label& Foam::particle::face() noexcept
{
    return facei_;
}


inline Foam::scalar Foam::particle::stepFraction() const noexcept
{
    return stepFraction_;
}


inline Foam::scalar& Foam::particle::stepFraction() noexcept
{
    return stepFraction_;
}


inline Foam::label Foam::particle::origProc() const noexcept
{
    return origProc_;
}


inline Foam::label& Foam::particle::origProc() noexcept
{
    return origProc_;
}


inline Foam::label Foam::particle::origId() const noexcept
{
    return origId_;
}


inline Foam::label& Foam::particle::origId() noexcept
{
    return origId_;
}


inline Foam::Pair<Foam::scalar> Foam::particle::stepFractionSpan() const
{
    if (mesh_.time().subCycling())
    {
        const TimeState& tsNew = mesh_.time();
        const TimeState& tsOld = mesh_.time().prevTimeState();

        const scalar tFrac =
        (
            (tsNew.value() - tsNew.deltaTValue())
          - (tsOld.value() - tsOld.deltaTValue())
        )/tsOld.deltaTValue();

        const scalar dtFrac = tsNew.deltaTValue()/tsOld.deltaTValue();

        return Pair<scalar>(tFrac, dtFrac);
    }

    return Pair<scalar>(0, 1);
}


inline Foam::scalar Foam::particle::currentTimeFraction() const
{
    const Pair<scalar> s = stepFractionSpan();

    return s[0] + stepFraction_*s[1];
}


inline Foam::tetIndices Foam::particle::currentTetIndices() const noexcept
{
    return tetIndices(celli_, tetFacei_, tetPti_);
}


inline Foam::barycentricTensor Foam::particle::currentTetTransform() const
{
    if (mesh_.moving() && stepFraction_ != 1)
    {
        return movingTetTransform(0)[0];
    }

    return stationaryTetTransform();
}


inline Foam::vector Foam::particle::normal() const
{
    return currentTetIndices().faceTri(mesh_).unitNormal();
}


inline bool Foam::particle::onFace() const noexcept
{
    return facei_ >= 0;
}


inline bool Foam::particle::onInternalFace() const noexcept
{
    return onFace() && mesh_.isInternalFace(facei_);
}


inline bool Foam::particle::onBoundaryFace() const noexcept
{
    return onFace() && !mesh_.isInternalFace(facei_);
}


inline Foam::label Foam::particle::patch() const
{
    return onFace() ? mesh_.boundaryMesh().whichPatch(facei_) : -1;
}


inline Foam::vector Foam::particle::position() const
{
    return currentTetTransform() & coordinates_;
}


inline void Foam::particle::reset()
{
    stepFraction_ = 0;
    behind_ = 0;
    nBehind_ = 0;
}


void Foam::particle::patchData(vector& n, vector& U) const
{
    if (!onBoundaryFace())
    {
        FatalErrorInFunction
            << "Patch data was requested for a particle that isn't on a patch"
            << exit(FatalError);
    }

    if ((mesh_.moving() && stepFraction_ != 1))
    {
        Pair<vector> centre, base, vertex1, vertex2;
        movingTetGeometry(1, centre, base, vertex1, vertex2);

        n = triPointRef::unitNormal(base[0], vertex1[0], vertex2[0]);

        // Interpolate the motion of the three face vertices to the current
        // coordinates
        U =
            coordinates_.b()*base[1]
          + coordinates_.c()*vertex1[1]
          + coordinates_.d()*vertex2[1];

        // The movingTetGeometry method gives the motion as a displacement
        // across the time-step, so we divide by the time-step to get velocity
        U /= mesh_.time().deltaTValue();
    }
    else
    {
        n = currentTetIndices().faceTri(mesh_).unitNormal();

        U = Zero;
    }
}


// ************************************************************************* //
