/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 Wikki Ltd
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

#include "boundedBackwardFaDdtScheme.H"
#include "faMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fa
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

scalar boundedBackwardFaDdtScheme::deltaT_() const
{
    return mesh().time().deltaTValue();
}


scalar boundedBackwardFaDdtScheme::deltaT0_() const
{
    return mesh().time().deltaT0Value();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

tmp<areaScalarField> boundedBackwardFaDdtScheme::facDdt
(
    const dimensionedScalar dt
)
{
    // No change compared to backward

    dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();

    const IOobject ddtIOobject
    (
        mesh().thisDb().newIOobject
        (
            "ddt("+dt.name()+')',
            { IOobject::REGISTER }
        )
    );

    scalar deltaT = deltaT_();
    scalar deltaT0 = deltaT0_();

    scalar coefft   = 1 + deltaT/(deltaT + deltaT0);
    scalar coefft00 = deltaT*deltaT/(deltaT0*(deltaT + deltaT0));
    scalar coefft0  = coefft + coefft00;

    if (mesh().moving())
    {
        tmp<areaScalarField> tdtdt
        (
            new areaScalarField
            (
                ddtIOobject,
                mesh(),
                dimensionedScalar(dt.dimensions()/dimTime, Zero)
            )
        );

        tdtdt.ref().primitiveFieldRef() = rDeltaT.value()*dt.value()*
        (
            coefft - (coefft0*mesh().S0() - coefft00*mesh().S00())/mesh().S()
        );

        return tdtdt;
    }

    return tmp<areaScalarField>::New
    (
        ddtIOobject,
        mesh(),
        dimensionedScalar(dt.dimensions()/dimTime, Zero),
        faPatchFieldBase::calculatedType()
    );
}


tmp<areaScalarField> boundedBackwardFaDdtScheme::facDdt0
(
    const dimensionedScalar dt
)
{
    // No change compared to backward

    dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();

    const IOobject ddtIOobject
    (
        mesh().thisDb().newIOobject
        (
            "ddt("+dt.name()+')',
            { IOobject::REGISTER }
        )
    );

    scalar deltaT = deltaT_();
    scalar deltaT0 = deltaT0_();

    scalar coefft   = 1 + deltaT/(deltaT + deltaT0);
    scalar coefft00 = deltaT*deltaT/(deltaT0*(deltaT + deltaT0));
    scalar coefft0  = coefft + coefft00;

    tmp<areaScalarField> tdtdt0
    (
        new areaScalarField
        (
            ddtIOobject,
            mesh(),
            -rDeltaT*(coefft0 - coefft00)*dt
        )
    );

    if (mesh().moving())
    {
        tdtdt0.ref().primitiveFieldRef() = (-rDeltaT.value()*dt.value())*
        (
            (coefft0*mesh().S0() - coefft00*mesh().S00())/mesh().S()
        );
    }

    return tdtdt0;
}


tmp<areaScalarField> boundedBackwardFaDdtScheme::facDdt
(
    const areaScalarField& vf
)
{
    dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();

    const IOobject ddtIOobject
    (
        mesh().thisDb().newIOobject
        (
            "ddt("+vf.name()+')',
            { IOobject::REGISTER }
        )
    );

    scalar deltaT = deltaT_();
    scalar deltaT0 = deltaT0_(vf);

    // Calculate unboundedness indicator
    // Note: all times moved by one because access to internal field
    // copies current field into the old-time level.
    areaScalarField phict
    (
        mag
        (
            vf.oldTime().oldTime()
          - vf.oldTime().oldTime().oldTime()
        )/
        (
            mag
            (
                vf.oldTime()
              - vf.oldTime().oldTime()
            )
          + dimensionedScalar("small", vf.dimensions(), SMALL)
        )
    );

    areaScalarField limiter(pos(phict) - pos(phict - scalar(1)));

    areaScalarField coefft(scalar(1) + limiter*deltaT/(deltaT + deltaT0));
    areaScalarField coefft00(limiter*sqr(deltaT)/(deltaT0*(deltaT + deltaT0)));
    areaScalarField coefft0(coefft + coefft00);

    if (mesh().moving())
    {
        return tmp<areaScalarField>
        (
            new areaScalarField
            (
                ddtIOobject,
                mesh(),
                rDeltaT.dimensions()*vf.dimensions(),
                rDeltaT.value()*
                (
                    coefft*vf.primitiveField() -
                    (
                        coefft0.primitiveField()
                        *vf.oldTime().primitiveField()*mesh().S0()
                      - coefft00.primitiveField()
                        *vf.oldTime().oldTime().primitiveField()
                       *mesh().S00()
                    )/mesh().S()
                ),
                rDeltaT.value()*
                (
                    coefft.boundaryField()*vf.boundaryField() -
                    (
                        coefft0.boundaryField()*
                        vf.oldTime().boundaryField()
                      - coefft00.boundaryField()*
                        vf.oldTime().oldTime().boundaryField()
                    )
                )
            )
        );
    }
    else
    {
        return tmp<areaScalarField>
        (
            new areaScalarField
            (
                ddtIOobject,
                rDeltaT*
                (
                    coefft*vf
                  - coefft0*vf.oldTime()
                  + coefft00*vf.oldTime().oldTime()
                )
            )
        );
    }
}


tmp<areaScalarField> boundedBackwardFaDdtScheme::facDdt0
(
    const areaScalarField& vf
)
{
    dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();

    const IOobject ddtIOobject
    (
        mesh().thisDb().newIOobject
        (
            "ddt0("+vf.name()+')',
            { IOobject::REGISTER }
        )
    );

    scalar deltaT = deltaT_();
    scalar deltaT0 = deltaT0_(vf);

    // Calculate unboundedness indicator
    // Note: all times moved by one because access to internal field
    // copies current field into the old-time level.
    areaScalarField phict
    (
        mag
        (
            vf.oldTime().oldTime()
          - vf.oldTime().oldTime().oldTime()
        )/
        (
            mag
            (
                vf.oldTime()
              - vf.oldTime().oldTime()
            )
          + dimensionedScalar("small", vf.dimensions(), SMALL)
        )
    );

    areaScalarField limiter(pos(phict) - pos(phict - scalar(1)));

    areaScalarField coefft(scalar(1) + limiter*deltaT/(deltaT + deltaT0));
    areaScalarField coefft00(limiter*sqr(deltaT)/(deltaT0*(deltaT + deltaT0)));
    areaScalarField coefft0(coefft + coefft00);

    if (mesh().moving())
    {
        return tmp<areaScalarField>
        (
            new areaScalarField
            (
                ddtIOobject,
                mesh(),
                rDeltaT.dimensions()*vf.dimensions(),
                rDeltaT.value()*
                (
                  - (
                        coefft0.primitiveField()*
                        vf.oldTime().primitiveField()*mesh().S0()
                      - coefft00.primitiveField()*
                        vf.oldTime().oldTime().primitiveField()
                       *mesh().S00()
                    )/mesh().S()
                ),
                rDeltaT.value()*
                (
                  - (
                        coefft0.boundaryField()*
                        vf.oldTime().boundaryField()
                      - coefft00.boundaryField()*
                        vf.oldTime().oldTime().boundaryField()
                    )
                )
            )
        );
    }
    else
    {
        return tmp<areaScalarField>
        (
            new areaScalarField
            (
                ddtIOobject,
                rDeltaT*
                (
                  - coefft0*vf.oldTime()
                  + coefft00*vf.oldTime().oldTime()
                )
            )
        );
    }
}


tmp<areaScalarField> boundedBackwardFaDdtScheme::facDdt
(
    const dimensionedScalar& rho,
    const areaScalarField& vf
)
{
    dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();

    const IOobject ddtIOobject
    (
        mesh().thisDb().newIOobject
        (
            "ddt("+rho.name()+','+vf.name()+')',
            { IOobject::REGISTER }
        )
    );

    scalar deltaT = deltaT_();
    scalar deltaT0 = deltaT0_(vf);

    // Calculate unboundedness indicator
    // Note: all times moved by one because access to internal field
    // copies current field into the old-time level.
    areaScalarField phict
    (
        mag
        (
            vf.oldTime().oldTime()
          - vf.oldTime().oldTime().oldTime()
        )/
        (
            mag
            (
                vf.oldTime()
              - vf.oldTime().oldTime()
            )
          + dimensionedScalar("small", vf.dimensions(), SMALL)
        )
    );

    areaScalarField limiter(pos(phict) - pos(phict - scalar(1)));

    areaScalarField coefft(scalar(1) + limiter*deltaT/(deltaT + deltaT0));
    areaScalarField coefft00(limiter*sqr(deltaT)/(deltaT0*(deltaT + deltaT0)));
    areaScalarField coefft0(coefft + coefft00);

    if (mesh().moving())
    {
        return tmp<areaScalarField>
        (
            new areaScalarField
            (
                ddtIOobject,
                mesh(),
                rDeltaT.dimensions()*rho.dimensions()*vf.dimensions(),
                rDeltaT.value()*rho.value()*
                (
                    coefft*vf.primitiveField() -
                    (
                        coefft0.primitiveField()*
                        vf.oldTime().primitiveField()*mesh().S0()
                      - coefft00.primitiveField()*
                        vf.oldTime().oldTime().primitiveField()
                       *mesh().S00()
                    )/mesh().S()
                ),
                rDeltaT.value()*rho.value()*
                (
                    coefft.boundaryField()*vf.boundaryField() -
                    (
                        coefft0.boundaryField()*
                        vf.oldTime().boundaryField()
                      - coefft00.boundaryField()*
                        vf.oldTime().oldTime().boundaryField()
                    )
                )
            )
        );
    }
    else
    {
        return tmp<areaScalarField>
        (
            new areaScalarField
            (
                ddtIOobject,
                rDeltaT*rho*
                (
                    coefft*vf
                  - coefft0*vf.oldTime()
                 + coefft00*vf.oldTime().oldTime()
                )
            )
        );
    }
}

tmp<areaScalarField> boundedBackwardFaDdtScheme::facDdt0
(
    const dimensionedScalar& rho,
    const areaScalarField& vf
)
{
    dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();

    const IOobject ddtIOobject
    (
        mesh().thisDb().newIOobject
        (
            "ddt0("+rho.name()+','+vf.name()+')',
            { IOobject::REGISTER }
        )
    );

    scalar deltaT = deltaT_();
    scalar deltaT0 = deltaT0_(vf);

    // Calculate unboundedness indicator
    // Note: all times moved by one because access to internal field
    // copies current field into the old-time level.
    areaScalarField phict
    (
        mag
        (
            vf.oldTime().oldTime()
          - vf.oldTime().oldTime().oldTime()
        )/
        (
            mag
            (
                vf.oldTime()
              - vf.oldTime().oldTime()
            )
          + dimensionedScalar("small", vf.dimensions(), SMALL)
        )
    );

    areaScalarField limiter(pos(phict) - pos(phict - scalar(1)));

    areaScalarField coefft(scalar(1) + limiter*deltaT/(deltaT + deltaT0));
    areaScalarField coefft00(limiter*sqr(deltaT)/(deltaT0*(deltaT + deltaT0)));
    areaScalarField coefft0(coefft + coefft00);

    if (mesh().moving())
    {
        return tmp<areaScalarField>
        (
            new areaScalarField
            (
                ddtIOobject,
                mesh(),
                rDeltaT.dimensions()*rho.dimensions()*vf.dimensions(),
                rDeltaT.value()*rho.value()*
                (
                   -(
                        coefft0.primitiveField()*
                        vf.oldTime().primitiveField()*mesh().S0()
                      - coefft00.primitiveField()*
                        vf.oldTime().oldTime().primitiveField()
                       *mesh().S00()
                    )/mesh().S()
                ),
                rDeltaT.value()*rho.value()*
                (
                   -(
                        coefft0.boundaryField()*
                        vf.oldTime().boundaryField()
                      - coefft00.boundaryField()*
                        vf.oldTime().oldTime().boundaryField()
                    )
                )
            )
        );
    }
    else
    {
        return tmp<areaScalarField>
        (
            new areaScalarField
            (
                ddtIOobject,
                rDeltaT*rho*
                (
                  - coefft0*vf.oldTime()
                 + coefft00*vf.oldTime().oldTime()
                )
            )
        );
    }
}


tmp<areaScalarField> boundedBackwardFaDdtScheme::facDdt
(
    const areaScalarField& rho,
    const areaScalarField& vf
)
{
    dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();

    const IOobject ddtIOobject
    (
        mesh().thisDb().newIOobject
        (
            "ddt("+rho.name()+','+vf.name()+')',
            { IOobject::REGISTER }
        )
    );

    scalar deltaT = deltaT_();
    scalar deltaT0 = deltaT0_(vf);

    // Calculate unboundedness indicator
    // Note: all times moved by one because access to internal field
    // copies current field into the old-time level.
    areaScalarField phict
    (
        mag
        (
            rho.oldTime().oldTime()*vf.oldTime().oldTime()
          - rho.oldTime().oldTime().oldTime()*vf.oldTime().oldTime().oldTime()
        )/
        (
            mag
            (
                rho.oldTime()*vf.oldTime()
              - rho.oldTime().oldTime()*vf.oldTime().oldTime()
            )
          + dimensionedScalar("small", rho.dimensions()*vf.dimensions(), SMALL)
        )
    );

    areaScalarField limiter(pos(phict) - pos(phict - scalar(1)));

    areaScalarField coefft(scalar(1) + limiter*deltaT/(deltaT + deltaT0));
    areaScalarField coefft00(limiter*sqr(deltaT)/(deltaT0*(deltaT + deltaT0)));
    areaScalarField coefft0(coefft + coefft00);

    if (mesh().moving())
    {
        return tmp<areaScalarField>
        (
            new areaScalarField
            (
                ddtIOobject,
                mesh(),
                rDeltaT.dimensions()*rho.dimensions()*vf.dimensions(),
                rDeltaT.value()*
                (
                    coefft*rho.primitiveField()*vf.primitiveField() -
                    (
                        coefft0.primitiveField()*
                        rho.oldTime().primitiveField()*
                        vf.oldTime().primitiveField()*mesh().S0()
                      - coefft00.primitiveField()*
                        rho.oldTime().oldTime().primitiveField()
                       *vf.oldTime().oldTime().primitiveField()*mesh().S00()
                    )/mesh().S()
                ),
                rDeltaT.value()*
                (
                    coefft.boundaryField()*vf.boundaryField() -
                    (
                        coefft0.boundaryField()*
                        rho.oldTime().boundaryField()*
                        vf.oldTime().boundaryField()
                      - coefft00.boundaryField()*
                        rho.oldTime().oldTime().boundaryField()*
                        vf.oldTime().oldTime().boundaryField()
                    )
                )
            )
        );
    }
    else
    {
        return tmp<areaScalarField>
        (
            new areaScalarField
            (
                ddtIOobject,
                rDeltaT*
                (
                    coefft*rho*vf
                  - coefft0*rho.oldTime()*vf.oldTime()
                  + coefft00*rho.oldTime().oldTime()*vf.oldTime().oldTime()
                )
            )
        );
    }
}


tmp<areaScalarField> boundedBackwardFaDdtScheme::facDdt0
(
    const areaScalarField& rho,
    const areaScalarField& vf
)
{
    dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();

    const IOobject ddtIOobject
    (
        mesh().thisDb().newIOobject
        (
            "ddt0("+rho.name()+','+vf.name()+')',
            { IOobject::REGISTER }
        )
    );

    scalar deltaT = deltaT_();
    scalar deltaT0 = deltaT0_(vf);

    // Calculate unboundedness indicator
    // Note: all times moved by one because access to internal field
    // copies current field into the old-time level.
    areaScalarField phict
    (
        mag
        (
            rho.oldTime().oldTime()*vf.oldTime().oldTime()
          - rho.oldTime().oldTime().oldTime()*vf.oldTime().oldTime().oldTime()
        )/
        (
            mag
            (
                rho.oldTime()*vf.oldTime()
              - rho.oldTime().oldTime()*vf.oldTime().oldTime()
            )
          + dimensionedScalar("small", rho.dimensions()*vf.dimensions(), SMALL)
        )
    );

    areaScalarField limiter(pos(phict) - pos(phict - scalar(1)));

    areaScalarField coefft(scalar(1) + limiter*deltaT/(deltaT + deltaT0));
    areaScalarField coefft00(limiter*sqr(deltaT)/(deltaT0*(deltaT + deltaT0)));
    areaScalarField coefft0(coefft + coefft00);

    if (mesh().moving())
    {
        return tmp<areaScalarField>
        (
            new areaScalarField
            (
                ddtIOobject,
                mesh(),
                rDeltaT.dimensions()*rho.dimensions()*vf.dimensions(),
                rDeltaT.value()*
                (
                  - (
                        coefft0.primitiveField()*
                        rho.oldTime().primitiveField()*
                        vf.oldTime().primitiveField()*mesh().S0()
                      - coefft00.primitiveField()*
                        rho.oldTime().oldTime().primitiveField()*
                        vf.oldTime().oldTime().primitiveField()*mesh().S00()
                    )/mesh().S()
                ),
                rDeltaT.value()*
                (
                  - (
                        coefft0.boundaryField()*
                        rho.oldTime().boundaryField()*
                        vf.oldTime().boundaryField()
                      - coefft00.boundaryField()*
                        rho.oldTime().oldTime().boundaryField()*
                        vf.oldTime().oldTime().boundaryField()
                    )
                )
            )
        );
    }
    else
    {
        return tmp<areaScalarField>
        (
            new areaScalarField
            (
                ddtIOobject,
                rDeltaT*
                (
                  - coefft0*rho.oldTime()*vf.oldTime()
                  + coefft00*rho.oldTime().oldTime()*vf.oldTime().oldTime()
                )
            )
        );
    }
}


tmp<faScalarMatrix> boundedBackwardFaDdtScheme::famDdt
(
    const areaScalarField& vf
)
{
    tmp<faScalarMatrix> tfam
    (
        new faScalarMatrix
        (
            vf,
            vf.dimensions()*dimArea/dimTime
        )
    );

    faScalarMatrix& fam = tfam.ref();

    scalar rDeltaT = 1.0/deltaT_();

    scalar deltaT = deltaT_();
    scalar deltaT0 = deltaT0_(vf);

    // Calculate unboundedness indicator
    // Note: all times moved by one because access to internal field
    // copies current field into the old-time level.
    scalarField phict
    (
        mag
        (
            vf.oldTime().oldTime().internalField()
          - vf.oldTime().oldTime().oldTime().internalField()
        )/
        (
            mag
            (
                vf.oldTime().internalField()
              - vf.oldTime().oldTime().internalField()
            )
            + SMALL
        )
    );

    scalarField limiter(pos(phict) - pos(phict - 1.0));

    scalarField coefft(1.0 + limiter*deltaT/(deltaT + deltaT0));
    scalarField coefft00(limiter*deltaT*deltaT/(deltaT0*(deltaT + deltaT0)));
    scalarField coefft0(coefft + coefft00);

    fam.diag() = (coefft*rDeltaT)*mesh().S();

    if (mesh().moving())
    {
        fam.source() = rDeltaT*
        (
            coefft0*vf.oldTime().primitiveField()*mesh().S0()
          - coefft00*vf.oldTime().oldTime().primitiveField()
           *mesh().S00()
        );
    }
    else
    {
        fam.source() = rDeltaT*mesh().S()*
        (
            coefft0*vf.oldTime().primitiveField()
          - coefft00*vf.oldTime().oldTime().primitiveField()
        );
    }

    return tfam;
}


tmp<faScalarMatrix> boundedBackwardFaDdtScheme::famDdt
(
    const dimensionedScalar& rho,
    const areaScalarField& vf
)
{
    tmp<faScalarMatrix> tfam
    (
        new faScalarMatrix
        (
            vf,
            rho.dimensions()*vf.dimensions()*dimArea/dimTime
        )
    );
    faScalarMatrix& fam = tfam.ref();

    scalar rDeltaT = 1.0/deltaT_();

    scalar deltaT = deltaT_();
    scalar deltaT0 = deltaT0_(vf);

    // Calculate unboundedness indicator
    // Note: all times moved by one because access to internal field
    // copies current field into the old-time level.
    scalarField phict
    (
        mag
        (
            vf.oldTime().oldTime().internalField()
          - vf.oldTime().oldTime().oldTime().internalField()
        )/
        (
            mag
            (
                vf.oldTime().internalField()
              - vf.oldTime().oldTime().internalField()
            )
            + SMALL
        )
    );

    scalarField limiter(pos(phict) - pos(phict - 1.0));

    scalarField coefft(1.0 + limiter*deltaT/(deltaT + deltaT0));
    scalarField coefft00(limiter*deltaT*deltaT/(deltaT0*(deltaT + deltaT0)));
    scalarField coefft0(coefft + coefft00);

    fam.diag() = (coefft*rDeltaT*rho.value())*mesh().S();

    if (mesh().moving())
    {
        fam.source() = rDeltaT*rho.value()*
        (
            coefft0*vf.oldTime().primitiveField()*mesh().S0()
          - coefft00*vf.oldTime().oldTime().primitiveField()
           *mesh().S00()
        );
    }
    else
    {
        fam.source() = rDeltaT*mesh().S()*rho.value()*
        (
            coefft0*vf.oldTime().primitiveField()
          - coefft00*vf.oldTime().oldTime().primitiveField()
        );
    }

    return tfam;
}


tmp<faScalarMatrix> boundedBackwardFaDdtScheme::famDdt
(
    const areaScalarField& rho,
    const areaScalarField& vf
)
{
    tmp<faScalarMatrix> tfam
    (
        new faScalarMatrix
        (
            vf,
            rho.dimensions()*vf.dimensions()*dimArea/dimTime
        )
    );
    faScalarMatrix& fam = tfam.ref();

    scalar rDeltaT = 1.0/deltaT_();

    scalar deltaT = deltaT_();
    scalar deltaT0 = deltaT0_(vf);

    // Calculate unboundedness indicator
    // Note: all times moved by one because access to internal field
    // copies current field into the old-time level.
    scalarField phict
    (
        mag
        (
            rho.oldTime().oldTime().internalField()*
            vf.oldTime().oldTime().internalField()
          - rho.oldTime().oldTime().oldTime().internalField()*
            vf.oldTime().oldTime().oldTime().internalField()
        )/
        (
            mag
            (
                rho.oldTime().internalField()*
                vf.oldTime().internalField()
              - rho.oldTime().oldTime().internalField()*
                vf.oldTime().oldTime().internalField()
            )
            + SMALL
        )
    );

    scalarField limiter(pos(phict) - pos(phict - 1.0));

    scalarField coefft(1.0 + limiter*deltaT/(deltaT + deltaT0));
    scalarField coefft00(limiter*deltaT*deltaT/(deltaT0*(deltaT + deltaT0)));
    scalarField coefft0(coefft + coefft00);

    fam.diag() = (coefft*rDeltaT)*rho.primitiveField()*mesh().S();

    if (mesh().moving())
    {
        fam.source() = rDeltaT*
        (
            coefft0*rho.oldTime().primitiveField()
           *vf.oldTime().primitiveField()*mesh().S0()
          - coefft00*rho.oldTime().oldTime().primitiveField()
           *vf.oldTime().oldTime().primitiveField()*mesh().S00()
        );
    }
    else
    {
        fam.source() = rDeltaT*mesh().S()*
        (
            coefft0*rho.oldTime().primitiveField()
           *vf.oldTime().primitiveField()
          - coefft00*rho.oldTime().oldTime().primitiveField()
           *vf.oldTime().oldTime().primitiveField()
        );
    }

    return tfam;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(boundedBackwardFaDdtScheme, 0);

faDdtScheme<scalar>::addIstreamConstructorToTable<boundedBackwardFaDdtScheme>
    addboundedBackwardFaDdtSchemeIstreamConstructorToTable_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fa

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
