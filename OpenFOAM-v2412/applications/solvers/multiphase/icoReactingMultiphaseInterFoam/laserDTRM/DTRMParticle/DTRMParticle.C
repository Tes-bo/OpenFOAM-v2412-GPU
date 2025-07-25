/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2019 OpenCFD Ltd
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

#include "DTRMParticle.H"
#include "constants.H"
#include "physicoChemicalConstants.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::DTRMParticle::DTRMParticle
(
    const polyMesh& mesh,
    const vector& position,
    const vector& targetPosition,
    const scalar I,
    const label cellI,
    const scalar dA,
    const label transmissiveId
)
:
    particle(mesh, position, cellI),
    p0_(position),
    p1_(targetPosition),
    I0_(I),
    I_(I),
    dA_(dA),
    transmissiveId_(transmissiveId)
{}


Foam::DTRMParticle::DTRMParticle
(
    const polyMesh& mesh,
    const barycentric& coordinates,
    const label celli,
    const label tetFacei,
    const label tetPti,
    const vector& position,
    const vector& targetPosition,
    const scalar I,
    const scalar dA,
    const label transmissiveId
)
:
    particle(mesh, coordinates, celli, tetFacei, tetPti),
    p0_(position),
    p1_(targetPosition),
    I0_(I),
    I_(I),
    dA_(dA),
    transmissiveId_(transmissiveId)
{}


Foam::DTRMParticle::DTRMParticle(const DTRMParticle& p)
:
    particle(p),
    p0_(p.p0_),
    p1_(p.p1_),
    I0_(p.I0_),
    I_(p.I_),
    dA_(p.dA_),
    transmissiveId_(p.transmissiveId_)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::DTRMParticle::move
(
    Cloud<DTRMParticle>& spc,
    trackingData& td,
    const scalar trackTime
)
{
    td.switchProcessor = false;
    td.keepParticle = true;

    do
    {
        //Cache old data of particle to use for reflected particle
        const point pos0 = position();
        const label cell1 = cell();

        scalar f = 1 - stepFraction();
        const vector s = p1() - p0() - deviationFromMeshCentre();
        trackToAndHitFace(f*s, f, spc, td);

        const point p1 = position();
        vector dsv = p1 - pos0;
        scalar ds = mag(dsv);

        //const label cell1 = cell();

        //NOTE:
        // Under the new barocentric tracking alghorithm the newly
        // inserted particles are tracked to the nearest cell centre first,
        // then, given the direction, to a face. In both occasions the first call
        // to trackToAndHitFace returns ds = 0. In this case we do an extra
        // call to trackToAndHitFace to start the tracking.
        // This is a temporary fix until the tracking can handle it.
        if (ds == 0)
        {
            trackToAndHitFace(f*s, f, spc, td);
            dsv = p1 - position();
            ds = mag(dsv);
        }

        // Boltzman constant
        const scalar sigma = physicoChemical::sigma.value();

        label reflectedZoneId = td.relfectedCells()[cell1];

        if
        (
            (reflectedZoneId > -1)
         && (
                (transmissiveId_ == -1)
             || (transmissiveId_ != reflectedZoneId)
            )
        )
        {
            scalar rho(0);

            // Create a new reflected particle when the particles is not
            // transmissive and larger than an absolute I
            if (I_ > 0.01*I0_ && ds > 0)
            {
                vector pDir = dsv/ds;

                cellPointWeight cpw(mesh(), position(), cell1, face());
                vector nHat = td.nHatInterp().interpolate(cpw);

                nHat /= (mag(nHat) + ROOTSMALL);
                scalar cosTheta(-pDir & nHat);

                // Only new incoming rays
                if (cosTheta > SMALL)
                {
                    vector newDir =
                        td.reflection()
                        [
                            td.relfectedCells()[cell1]
                        ].R(pDir, nHat);

                    // reflectivity
                    rho =
                        min
                        (
                            max
                            (
                                td.reflection()
                                [
                                    td.relfectedCells()[cell1]
                                ].rho(cosTheta)
                                , 0.0
                            )
                            , 0.98
                        );

                    //scalar delaM = cbrt(mesh().cellVolumes()[cell0]);
                    scalar delaM = cbrt(mesh().cellVolumes()[cell1]);

                    const point insertP(position() - pDir*0.1*delaM);
                    label cellI = mesh().findCell(insertP);

                    if (cellI > -1)
                    {
                        DTRMParticle* pPtr = new DTRMParticle
                        (
                            mesh(),
                            insertP,
                            insertP + newDir*mesh().bounds().mag(),
                            I_*rho,
                            cellI,
                            dA_,
                            -1
                        );

                        // Add to cloud
                        spc.addParticle(pPtr);
                    }
                }
            }

            // Change transmissiveId of the particle
            transmissiveId_ = reflectedZoneId;

            scalar a = td.aInterp().interpolate(pos0, cell1);
            scalar e = td.eInterp().interpolate(pos0, cell1);
            scalar E = td.EInterp().interpolate(pos0, cell1);
            scalar T = td.TInterp().interpolate(pos0, cell1);


            // Left intensity after reflection
            const scalar Itran = I_*(1.0 - rho);
            const scalar I1 =
            (
                Itran
              + ds*(e*sigma*pow4(T)/mathematical::pi + E)
            ) / (1 + ds*a);

            td.Q(cell1) += (Itran - max(I1, 0.0))*dA_;

            I_ = I1;

            if (I_ <= 0.01*I0_)
            {
                stepFraction() = 1.0;
                break;
            }
        }
        else
        {
            scalar a = td.aInterp().interpolate(pos0, cell1);
            scalar e = td.eInterp().interpolate(pos0, cell1);
            scalar E = td.EInterp().interpolate(pos0, cell1);
            scalar T = td.TInterp().interpolate(pos0, cell1);

            const scalar I1 =
            (
                I_
                + ds*(e*sigma*pow4(T)/mathematical::pi + E)
            ) / (1 + ds*a);

            td.Q(cell1) += (I_ -  max(I1, 0.0))*dA_;

            I_ = I1;

            if ((I_ <= 0.01*I0_))
            {
                stepFraction() = 1.0;
                break;
            }
        }

    }while (td.keepParticle && !td.switchProcessor && stepFraction() < 1);

    return td.keepParticle;
}


void Foam::DTRMParticle::hitProcessorPatch
(
    Cloud<DTRMParticle>&,
    trackingData& td
)
{
    td.switchProcessor = true;
}


void Foam::DTRMParticle::hitWallPatch
(
    Cloud<DTRMParticle>&,
    trackingData& td
)
{
    td.keepParticle = false;
}


bool Foam::DTRMParticle::hitPatch
(
    Cloud<DTRMParticle>&,
    trackingData& td
)
{
    return false;
}


// ************************************************************************* //
