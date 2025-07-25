/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2018-2020 OpenCFD Ltd.
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

#include "rotorDiskSource.H"
#include "volFields.H"
#include "unitConversion.H"

using namespace Foam::constant;

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class RhoFieldType>
void Foam::fv::rotorDiskSource::calculate
(
    const RhoFieldType& rho,
    const vectorField& U,
    const scalarField& thetag,
    vectorField& force,
    const bool divideVolume,
    const bool output
) const
{
    const scalarField& V = mesh_.V();

    // Logging info
    scalar dragEff = 0.0;
    scalar liftEff = 0.0;
    scalar AOAmin = GREAT;
    scalar AOAmax = -GREAT;

    // Cached position-dependent rotations available?
    const bool hasCache = bool(Rcyl_);

    forAll(cells_, i)
    {
        if (area_[i] > ROOTVSMALL)
        {
            const label celli = cells_[i];

            const scalar radius = x_[i].x();

            const tensor Rcyl =
            (
                hasCache
              ? (*Rcyl_)[i]
              : coordSys_.R(mesh_.C()[celli])
            );

            // Transform velocity into local cylindrical reference frame
            vector Uc = invTransform(Rcyl, U[celli]);

            // Transform velocity into local coning system
            Uc = transform(Rcone_[i], Uc);

            // Set radial component of velocity to zero
            Uc.x() = 0.0;

            // Set blade normal component of velocity
            Uc.y() = radius*omega_ - Uc.y();

            // Determine blade data for this radius
            // i2 = index of upper radius bound data point in blade list
            scalar twist = 0.0;
            scalar chord = 0.0;
            label i1 = -1;
            label i2 = -1;
            scalar invDr = 0.0;
            blade_.interpolate(radius, twist, chord, i1, i2, invDr);

            // Flip geometric angle if blade is spinning in reverse (clockwise)
            scalar alphaGeom = thetag[i] + twist;
            if (omega_ < 0)
            {
                alphaGeom = mathematical::pi - alphaGeom;
            }

            // Effective angle of attack
            scalar alphaEff = alphaGeom - atan2(-Uc.z(), Uc.y());
            if (alphaEff > mathematical::pi)
            {
                alphaEff -= mathematical::twoPi;
            }
            if (alphaEff < -mathematical::pi)
            {
                alphaEff += mathematical::twoPi;
            }

            AOAmin = min(AOAmin, alphaEff);
            AOAmax = max(AOAmax, alphaEff);

            // Determine profile data for this radius and angle of attack
            const label profile1 = blade_.profileID()[i1];
            const label profile2 = blade_.profileID()[i2];

            scalar Cd1 = 0.0;
            scalar Cl1 = 0.0;
            profiles_[profile1].Cdl(alphaEff, Cd1, Cl1);

            scalar Cd2 = 0.0;
            scalar Cl2 = 0.0;
            profiles_[profile2].Cdl(alphaEff, Cd2, Cl2);

            scalar Cd = invDr*(Cd2 - Cd1) + Cd1;
            scalar Cl = invDr*(Cl2 - Cl1) + Cl1;

            // Apply tip effect for blade lift
            scalar tipFactor = neg(radius/rMax_ - tipEffect_);

            // Calculate forces perpendicular to blade
            scalar pDyn = 0.5*rho[celli]*magSqr(Uc);

            scalar f = pDyn*chord*nBlades_*area_[i]/radius/mathematical::twoPi;
            vector localForce = vector(0.0, -f*Cd, tipFactor*f*Cl);

            // Accumulate forces
            dragEff += rhoRef_*localForce.y();
            liftEff += rhoRef_*localForce.z();

            // Transform force from local coning system into rotor cylindrical
            localForce = invTransform(Rcone_[i], localForce);

            // Transform force into global Cartesian coordinate system
            force[celli] = transform(Rcyl, localForce);

            if (divideVolume)
            {
                force[celli] /= V[celli];
            }
        }
    }

    if (output)
    {
        reduce(AOAmin, minOp<scalar>());
        reduce(AOAmax, maxOp<scalar>());
        reduce(dragEff, sumOp<scalar>());
        reduce(liftEff, sumOp<scalar>());

        Info<< type() << " output:" << nl
            << "    min/max(AOA)   = " << radToDeg(AOAmin) << ", "
            << radToDeg(AOAmax) << nl
            << "    Effective drag = " << dragEff << nl
            << "    Effective lift = " << liftEff << endl;
    }
}


template<class Type>
void Foam::fv::rotorDiskSource::writeField
(
    const word& name,
    const List<Type>& values,
    const bool writeNow
) const
{
    if (mesh_.time().writeTime() || writeNow)
    {
        if (cells_.size() != values.size())
        {
            FatalErrorInFunction
                << "Size mismatch. Number of cells "
                << cells_.size() << " != number of values "
                << values.size() << nl
                << abort(FatalError);
        }

        auto tfield = GeometricField<Type, fvPatchField, volMesh>::New
        (
            name,
            IOobject::NO_REGISTER,
            mesh_,
            dimensioned<Type>(dimless, Zero)
        );

        auto& field = tfield.ref().primitiveFieldRef();

        UIndirectList<Type>(field, cells_) = values;

        tfield().write();
    }
}


// ************************************************************************* //
