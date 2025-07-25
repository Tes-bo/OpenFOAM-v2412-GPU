/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017 Volkswagen AG
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

#include "EulerFaD2dt2Scheme.H"
#include "faMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fa
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
scalar EulerFaD2dt2Scheme<Type>::deltaT_() const
{
    return mesh().time().deltaTValue();
}


template<class Type>
scalar EulerFaD2dt2Scheme<Type>::deltaT0_() const
{
    return mesh().time().deltaT0Value();
}


template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh>>
EulerFaD2dt2Scheme<Type>::facD2dt2
(
    const dimensioned<Type> dt
)
{

    scalar deltaT = deltaT_();
    scalar deltaT0 = deltaT0_();

    dimensionedScalar rDeltaT2 =
        4.0/sqr(mesh().time().deltaT() + mesh().time().deltaT0());

    scalar coefft   = (deltaT + deltaT0)/(2*deltaT);
    scalar coefft00 = (deltaT + deltaT0)/(2*deltaT0);

    const IOobject d2dt2IOobject
    (
        mesh().thisDb().newIOobject
        (
            "d2dt2("+dt.name()+')',
            { IOobject::REGISTER }
        )
    );

    if (mesh().moving())
    {
        scalar halfRdeltaT2 = 0.5*rDeltaT2.value();

        scalarField SS0(mesh().S() + mesh().S0());
        scalarField S0S00(mesh().S0() + mesh().S00());

        tmp<GeometricField<Type, faPatchField, areaMesh>> tdt2dt2
        (
            new GeometricField<Type, faPatchField, areaMesh>
            (
                d2dt2IOobject,
                mesh(),
                dimensioned<Type>(dt.dimensions()/dimTime/dimTime, Zero)
            )
        );

        tdt2dt2.ref().primitiveFieldRef() =
            halfRdeltaT2*dt.value()
           *(coefft*SS0 - (coefft*SS0 + coefft00*S0S00)  + coefft00*S0S00)
           /mesh().S();

        return tdt2dt2;
    }
    else
    {
        return tmp<GeometricField<Type, faPatchField, areaMesh>>
        (
            new GeometricField<Type, faPatchField, areaMesh>
            (
                d2dt2IOobject,
                mesh(),
                dimensioned<Type>(dt.dimensions()/dimTime/dimTime, Zero),
                faPatchFieldBase::calculatedType()
            )
        );
    }
}


template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh>>
EulerFaD2dt2Scheme<Type>::facD2dt2
(
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    dimensionedScalar rDeltaT2 =
        4.0/sqr(mesh().time().deltaT() + mesh().time().deltaT0());

    const IOobject d2dt2IOobject
    (
        mesh().thisDb().newIOobject
        (
            "d2dt2("+vf.name()+')',
            { IOobject::REGISTER }
        )
    );

    scalar deltaT = deltaT_();
    scalar deltaT0 = deltaT0_();

    scalar coefft   = (deltaT + deltaT0)/(2*deltaT);
    scalar coefft00 = (deltaT + deltaT0)/(2*deltaT0);
    scalar coefft0  = coefft + coefft00;

    if (mesh().moving())
    {
        scalar halfRdeltaT2 = 0.5*rDeltaT2.value();

        scalarField SS0(mesh().S() + mesh().S0());
        scalarField S0S00(mesh().S0() + mesh().S00());

        return tmp<GeometricField<Type, faPatchField, areaMesh>>
        (
            new GeometricField<Type, faPatchField, areaMesh>
            (
                d2dt2IOobject,
                mesh(),
                rDeltaT2.dimensions()*vf.dimensions(),
                halfRdeltaT2*
                (
                    coefft*SS0*vf.primitiveField()

                  - (coefft*SS0 + coefft00*S0S00)
                   *vf.oldTime().primitiveField()

                  + (coefft00*S0S00)*vf.oldTime().oldTime().primitiveField()
                )/mesh().S(),
                rDeltaT2.value()*
                (
                    coefft*vf.boundaryField()
                  - coefft0*vf.oldTime().boundaryField()
                  + coefft00*vf.oldTime().oldTime().boundaryField()
                )
            )
        );
    }
    else
    {
        return tmp<GeometricField<Type, faPatchField, areaMesh>>
        (
            new GeometricField<Type, faPatchField, areaMesh>
            (
                d2dt2IOobject,
                rDeltaT2*
                (
                    coefft*vf
                  - coefft0*vf.oldTime()
                  + coefft00*vf.oldTime().oldTime()
                )
            )
        );
    }
}


template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh>>
EulerFaD2dt2Scheme<Type>::facD2dt2
(
    const dimensionedScalar& rho,
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    dimensionedScalar rDeltaT2 =
        4.0/sqr(mesh().time().deltaT() + mesh().time().deltaT0());

    const IOobject d2dt2IOobject
    (
        mesh().thisDb().newIOobject
        (
            "d2dt2("+rho.name()+','+vf.name()+')',
            { IOobject::REGISTER }
        )
    );

    scalar deltaT = deltaT_();
    scalar deltaT0 = deltaT0_();

    scalar coefft   = (deltaT + deltaT0)/(2*deltaT);
    scalar coefft00 = (deltaT + deltaT0)/(2*deltaT0);
    scalar coefft0  = coefft + coefft00;

    if (mesh().moving())
    {
        scalar halfRdeltaT2 = 0.5*rDeltaT2.value();

        scalarField SS0rhoRho0
        (
            (mesh().S() + mesh().S0())*rho.value()
        );

        scalarField S0S00rho0Rho00
        (
            (mesh().S0() + mesh().S00())*rho.value()
        );

        return tmp<GeometricField<Type, faPatchField, areaMesh>>
        (
            new GeometricField<Type, faPatchField, areaMesh>
            (
                d2dt2IOobject,
                mesh(),
                rDeltaT2.dimensions()*rho.dimensions()*vf.dimensions(),
                halfRdeltaT2*
                (
                    coefft*SS0rhoRho0*vf.primitiveField()

                  - (coefft*SS0rhoRho0 + coefft00*S0S00rho0Rho00)
                   *vf.oldTime().primitiveField()

                  + (coefft00*S0S00rho0Rho00)
                   *vf.oldTime().oldTime().primitiveField()
                )/mesh().S(),
                rDeltaT2.value()*rho.value()*
                (
                    coefft*vf.boundaryField()
                  - coefft0*vf.oldTime().boundaryField()
                  + coefft00*vf.oldTime().oldTime().boundaryField()
                )
            )
        );
    }
    else
    {

        return tmp<GeometricField<Type, faPatchField, areaMesh>>
        (
            new GeometricField<Type, faPatchField, areaMesh>
            (
                d2dt2IOobject,
                rDeltaT2 * rho.value() *
                (
                    coefft*vf
                  - coefft0*vf.oldTime()
                  + coefft00*vf.oldTime().oldTime()
                )
            )
        );
    }
}


template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh>>
EulerFaD2dt2Scheme<Type>::facD2dt2
(
    const areaScalarField& rho,
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    dimensionedScalar rDeltaT2 =
        4.0/sqr(mesh().time().deltaT() + mesh().time().deltaT0());

    const IOobject d2dt2IOobject
    (
        mesh().thisDb().newIOobject
        (
            "d2dt2("+rho.name()+','+vf.name()+')',
            { IOobject::REGISTER }
        )
    );

    scalar deltaT = deltaT_();
    scalar deltaT0 = deltaT0_();

    scalar coefft   = (deltaT + deltaT0)/(2*deltaT);
    scalar coefft00 = (deltaT + deltaT0)/(2*deltaT0);

    if (mesh().moving())
    {
        scalar halfRdeltaT2 = 0.5*rDeltaT2.value();
        scalar quarterRdeltaT2 = 0.25*rDeltaT2.value();

        scalarField SS0rhoRho0
        (
            (mesh().S() + mesh().S0())
           *(rho.primitiveField() + rho.oldTime().primitiveField())
        );

        scalarField S0S00rho0Rho00
        (
            (mesh().S0() + mesh().S00())
           *(
               rho.oldTime().primitiveField()
             + rho.oldTime().oldTime().primitiveField()
            )
        );

        return tmp<GeometricField<Type, faPatchField, areaMesh>>
        (
            new GeometricField<Type, faPatchField, areaMesh>
            (
                d2dt2IOobject,
                mesh(),
                rDeltaT2.dimensions()*rho.dimensions()*vf.dimensions(),
                quarterRdeltaT2*
                (
                    coefft*SS0rhoRho0*vf.primitiveField()

                  - (coefft*SS0rhoRho0 + coefft00*S0S00rho0Rho00)
                   *vf.oldTime().primitiveField()

                  + (coefft00*S0S00rho0Rho00)
                   *vf.oldTime().oldTime().primitiveField()
                )/mesh().S(),
                halfRdeltaT2*
                (
                    coefft
                   *(rho.boundaryField() + rho.oldTime().boundaryField())
                   *vf.boundaryField()

                  - (
                        coefft
                       *(
                           rho.boundaryField()
                         + rho.oldTime().boundaryField()
                        )
                      + coefft00
                       *(
                           rho.oldTime().boundaryField()
                         + rho.oldTime().oldTime().boundaryField()
                        )
                    )*vf.oldTime().boundaryField()

                  + coefft00
                   *(
                       rho.oldTime().boundaryField()
                     + rho.oldTime().oldTime().boundaryField()
                    )*vf.oldTime().oldTime().boundaryField()
                )
            )
        );
    }
    else
    {
        dimensionedScalar halfRdeltaT2 = 0.5*rDeltaT2;

        areaScalarField rhoRho0(rho + rho.oldTime());
        areaScalarField rho0Rho00(rho.oldTime() + rho.oldTime().oldTime());

        return tmp<GeometricField<Type, faPatchField, areaMesh>>
        (
            new GeometricField<Type, faPatchField, areaMesh>
            (
                d2dt2IOobject,
                halfRdeltaT2*
                (
                    coefft*rhoRho0*vf
                  - (coefft*rhoRho0 + coefft00*rho0Rho00)*vf.oldTime()
                  + coefft00*rho0Rho00*vf.oldTime().oldTime()
                )
            )
        );
    }
}


template<class Type>
tmp<faMatrix<Type>>
EulerFaD2dt2Scheme<Type>::famD2dt2
(
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    tmp<faMatrix<Type>> tfam
    (
        new faMatrix<Type>
        (
            vf,
            vf.dimensions()*dimArea/dimTime/dimTime
        )
    );

    faMatrix<Type>& fam = tfam.ref();

    scalar deltaT = deltaT_();
    scalar deltaT0 = deltaT0_();

    scalar coefft = (deltaT + deltaT0)/(2*deltaT);
    scalar coefft00 = (deltaT + deltaT0)/(2*deltaT0);
    scalar coefft0 = coefft + coefft00;

    scalar rDeltaT2 = 4.0/sqr(deltaT + deltaT0);

    if (mesh().moving())
    {
        scalar halfRdeltaT2 = rDeltaT2/2.0;

        scalarField SS0(mesh().S() + mesh().S0());
        scalarField S0S00(mesh().S0() + mesh().S00());

        fam.diag() = (coefft*halfRdeltaT2)*SS0;

        fam.source() = halfRdeltaT2*
        (
            (coefft*SS0 + coefft00*S0S00)
           *vf.oldTime().primitiveField()

          - (coefft00*S0S00)*vf.oldTime().oldTime().primitiveField()
        );
    }
    else
    {
        fam.diag() = (coefft*rDeltaT2)*mesh().S();

        fam.source() = rDeltaT2*mesh().S()*
        (
            coefft0*vf.oldTime().primitiveField()
          - coefft00*vf.oldTime().oldTime().primitiveField()
        );
    }

    return tfam;
}


template<class Type>
tmp<faMatrix<Type>>
EulerFaD2dt2Scheme<Type>::famD2dt2
(
    const dimensionedScalar& rho,
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    tmp<faMatrix<Type>> tfam
    (
        new faMatrix<Type>
        (
            vf,
            rho.dimensions()*vf.dimensions()*dimArea
            /dimTime/dimTime
        )
    );

    faMatrix<Type>& fam = tfam.ref();

    scalar deltaT = deltaT_();
    scalar deltaT0 = deltaT0_();

    scalar coefft = (deltaT + deltaT0)/(2*deltaT);
    scalar coefft00 = (deltaT + deltaT0)/(2*deltaT0);

    scalar rDeltaT2 = 4.0/sqr(deltaT + deltaT0);

    if (mesh().moving())
    {
        scalar halfRdeltaT2 = 0.5*rDeltaT2;

        scalarField SS0(mesh().S() + mesh().S0());
        scalarField S0S00(mesh().S0() + mesh().S00());

        fam.diag() = rho.value()*(coefft*halfRdeltaT2)*SS0;

        fam.source() = halfRdeltaT2*rho.value()*
        (
            (coefft*SS0 + coefft00*S0S00)
           *vf.oldTime().primitiveField()

          - (coefft00*S0S00)*vf.oldTime().oldTime().primitiveField()
        );
    }
    else
    {
        fam.diag() = (coefft*rDeltaT2)*mesh().S()*rho.value();

        fam.source() = rDeltaT2*mesh().S()*rho.value()*
        (
            (coefft + coefft00)*vf.oldTime().primitiveField()
          - coefft00*vf.oldTime().oldTime().primitiveField()
        );
    }

    return tfam;
}


template<class Type>
tmp<faMatrix<Type>>
EulerFaD2dt2Scheme<Type>::famD2dt2
(
    const areaScalarField& rho,
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    tmp<faMatrix<Type>> tfam
    (
        new faMatrix<Type>
        (
            vf,
            rho.dimensions()*vf.dimensions()*dimArea/dimTime/dimTime
        )
    );
    faMatrix<Type>& fam = tfam.ref();

    scalar deltaT =deltaT_();
    scalar deltaT0 = deltaT0_();

    scalar coefft   = (deltaT + deltaT0)/(2*deltaT);
    scalar coefft00 = (deltaT + deltaT0)/(2*deltaT0);

    scalar rDeltaT2 = 4.0/sqr(deltaT + deltaT0);

    if (mesh().moving())
    {
        scalar quarterRdeltaT2 = 0.25*rDeltaT2;

        scalarField SS0rhoRho0
        (
            (mesh().S() + mesh().S0())
           *(rho.primitiveField() + rho.oldTime().primitiveField())
        );

        scalarField S0S00rho0Rho00
        (
            (mesh().S0() + mesh().S00())
           *(
               rho.oldTime().primitiveField()
             + rho.oldTime().oldTime().primitiveField()
            )
        );

        fam.diag() = (coefft*quarterRdeltaT2)*SS0rhoRho0;

        fam.source() = quarterRdeltaT2*
        (
            (coefft*SS0rhoRho0 + coefft00*S0S00rho0Rho00)
           *vf.oldTime().primitiveField()

          - (coefft00*S0S00rho0Rho00)
           *vf.oldTime().oldTime().primitiveField()
        );
    }
    else
    {
        scalar halfRdeltaT2 = 0.5*rDeltaT2;

        scalarField rhoRho0
        (
            rho.primitiveField() + rho.oldTime().primitiveField()
        );

        scalarField rho0Rho00
        (
            rho.oldTime().primitiveField()
          + rho.oldTime().oldTime().primitiveField()
        );

        fam.diag() = (coefft*halfRdeltaT2)*mesh().S()*rhoRho0;

        fam.source() = halfRdeltaT2*mesh().S()*
        (
            (coefft*rhoRho0 + coefft00*rho0Rho00)
           *vf.oldTime().primitiveField()

          - (coefft00*rho0Rho00)
           *vf.oldTime().oldTime().primitiveField()
        );
    }

    return tfam;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fa

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
