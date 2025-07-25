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

#include "EulerFaDdtScheme.H"
#include "faMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fa
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh>>
EulerFaDdtScheme<Type>::facDdt
(
    const dimensioned<Type> dt
)
{
    dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();

    const IOobject ddtIOobject
    (
        mesh().thisDb().newIOobject
        (
            "ddt("+dt.name()+')',
            { IOobject::REGISTER }
        )
    );

    if (mesh().moving())
    {
        tmp<GeometricField<Type, faPatchField, areaMesh>> tdtdt
        (
            new GeometricField<Type, faPatchField, areaMesh>
            (
                ddtIOobject,
                mesh(),
                dimensioned<Type>(dt.dimensions()/dimTime, Zero)
            )
        );

        tdtdt.ref().primitiveFieldRef() =
            rDeltaT.value()*dt.value()*(1.0 - mesh().S0()/mesh().S());

        return tdtdt;
    }


    return tmp<GeometricField<Type, faPatchField, areaMesh>>::New
    (
        ddtIOobject,
        mesh(),
        dimensioned<Type>(dt.dimensions()/dimTime, Zero),
        faPatchFieldBase::calculatedType()
    );
}


template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh>>
EulerFaDdtScheme<Type>::facDdt0
(
    const dimensioned<Type> dt
)
{
    dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();

    const IOobject ddtIOobject
    (
        mesh().thisDb().newIOobject
        (
            "ddt("+dt.name()+')',
            { IOobject::REGISTER }
        )
    );

    tmp<GeometricField<Type, faPatchField, areaMesh>> tdtdt0
    (
        new GeometricField<Type, faPatchField, areaMesh>
        (
            ddtIOobject,
            mesh(),
            -rDeltaT*dt
        )
    );

    if (mesh().moving())
    {
        tdtdt0.ref().primitiveFieldRef() =
            (-rDeltaT.value()*dt.value())*mesh().S0()/mesh().S();
    }

    return tdtdt0;
}


template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh>>
EulerFaDdtScheme<Type>::facDdt
(
    const GeometricField<Type, faPatchField, areaMesh>& vf
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

    if (mesh().moving())
    {
        return tmp<GeometricField<Type, faPatchField, areaMesh>>
        (
            new GeometricField<Type, faPatchField, areaMesh>
            (
                ddtIOobject,
                mesh(),
                rDeltaT.dimensions()*vf.dimensions(),
                rDeltaT.value()*
                (
                    vf()
                  - vf.oldTime()()*mesh().S0()/mesh().S()
                ),
                rDeltaT.value()*
                (
                    vf.boundaryField() - vf.oldTime().boundaryField()
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
                ddtIOobject,
                rDeltaT*(vf - vf.oldTime())
            )
        );
    }
}


template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh>>
EulerFaDdtScheme<Type>::facDdt0
(
    const GeometricField<Type, faPatchField, areaMesh>& vf
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

    if (mesh().moving())
    {
        return tmp<GeometricField<Type, faPatchField, areaMesh>>
        (
            new GeometricField<Type, faPatchField, areaMesh>
            (
                ddtIOobject,
                mesh(),
                rDeltaT.dimensions()*vf.dimensions(),
                (-rDeltaT.value())*vf.oldTime().internalField(),
                (-rDeltaT.value())*vf.oldTime().boundaryField()
            )
        );
    }
    else
    {
        return tmp<GeometricField<Type, faPatchField, areaMesh>>
        (
            new GeometricField<Type, faPatchField, areaMesh>
            (
                ddtIOobject,
                (-rDeltaT)*vf.oldTime()
            )
        );
    }
}


template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh>>
EulerFaDdtScheme<Type>::facDdt
(
    const dimensionedScalar& rho,
    const GeometricField<Type, faPatchField, areaMesh>& vf
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

    if (mesh().moving())
    {
        return tmp<GeometricField<Type, faPatchField, areaMesh>>
        (
            new GeometricField<Type, faPatchField, areaMesh>
            (
                ddtIOobject,
                mesh(),
                rDeltaT.dimensions()*rho.dimensions()*vf.dimensions(),
                rDeltaT.value()*rho.value()*
                (
                    vf()
                  - vf.oldTime()()*mesh().S0()/mesh().S()
                ),
                rDeltaT.value()*rho.value()*
                (
                    vf.boundaryField() - vf.oldTime().boundaryField()
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
                ddtIOobject,
                rDeltaT*rho*(vf - vf.oldTime())
            )
        );
    }
}

template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh>>
EulerFaDdtScheme<Type>::facDdt0
(
    const dimensionedScalar& rho,
    const GeometricField<Type, faPatchField, areaMesh>& vf
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

    if (mesh().moving())
    {
        return tmp<GeometricField<Type, faPatchField, areaMesh>>
        (
            new GeometricField<Type, faPatchField, areaMesh>
            (
                ddtIOobject,
                mesh(),
                rDeltaT.dimensions()*rho.dimensions()*vf.dimensions(),
                (-rDeltaT.value())*rho.value()*
                    vf.oldTime()()*mesh().S0()/mesh().S(),
                (-rDeltaT.value())*rho.value()*
                    vf.oldTime().boundaryField()
            )
        );
    }
    else
    {
        return tmp<GeometricField<Type, faPatchField, areaMesh>>
        (
            new GeometricField<Type, faPatchField, areaMesh>
            (
                ddtIOobject,
                (-rDeltaT)*rho*vf.oldTime()
            )
        );
    }
}


template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh>>
EulerFaDdtScheme<Type>::facDdt
(
    const areaScalarField& rho,
    const GeometricField<Type, faPatchField, areaMesh>& vf
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

    if (mesh().moving())
    {
        return tmp<GeometricField<Type, faPatchField, areaMesh>>
        (
            new GeometricField<Type, faPatchField, areaMesh>
            (
                ddtIOobject,
                mesh(),
                rDeltaT.dimensions()*rho.dimensions()*vf.dimensions(),
                rDeltaT.value()*
                (
                    rho()*vf()
                  - rho.oldTime()()
                   *vf.oldTime()()*mesh().S0()/mesh().S()
                ),
                rDeltaT.value()*
                (
                    rho.boundaryField()*vf.boundaryField()
                  - rho.oldTime().boundaryField()
                   *vf.oldTime().boundaryField()
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
                ddtIOobject,
                rDeltaT*(rho*vf - rho.oldTime()*vf.oldTime())
            )
        );
    }
}


template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh>>
EulerFaDdtScheme<Type>::facDdt0
(
    const areaScalarField& rho,
    const GeometricField<Type, faPatchField, areaMesh>& vf
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

    if (mesh().moving())
    {
        return tmp<GeometricField<Type, faPatchField, areaMesh>>
        (
            new GeometricField<Type, faPatchField, areaMesh>
            (
                ddtIOobject,
                mesh(),
                rDeltaT.dimensions()*rho.dimensions()*vf.dimensions(),
                rDeltaT.value()*
                (
                  - rho.oldTime()()
                   *vf.oldTime()()*mesh().S0()/mesh().S()
                ),
                rDeltaT.value()*
                (
                  - rho.oldTime().boundaryField()
                   *vf.oldTime().boundaryField()
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
                ddtIOobject,
                (-rDeltaT)*rho.oldTime()*vf.oldTime()
            )
        );
    }
}

template<class Type>
tmp<faMatrix<Type>>
EulerFaDdtScheme<Type>::famDdt
(
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    tmp<faMatrix<Type>> tfam
    (
        new faMatrix<Type>
        (
            vf,
            vf.dimensions()*dimArea/dimTime
        )
    );

    faMatrix<Type>& fam = tfam.ref();

    scalar rDeltaT = 1.0/mesh().time().deltaTValue();

    fam.diag() = rDeltaT*mesh().S();

    if (mesh().moving())
    {
        fam.source() = rDeltaT*vf.oldTime().primitiveField()*mesh().S0();
    }
    else
    {
        fam.source() = rDeltaT*vf.oldTime().primitiveField()*mesh().S();
    }

    return tfam;
}


template<class Type>
tmp<faMatrix<Type>>
EulerFaDdtScheme<Type>::famDdt
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
            rho.dimensions()*vf.dimensions()*dimArea/dimTime
        )
    );
    faMatrix<Type>& fam = tfam.ref();

    scalar rDeltaT = 1.0/mesh().time().deltaTValue();

    fam.diag() = rDeltaT*rho.value()*mesh().S();

    if (mesh().moving())
    {
        fam.source() = rDeltaT
            *rho.value()*vf.oldTime().primitiveField()*mesh().S0();
    }
    else
    {
        fam.source() = rDeltaT
            *rho.value()*vf.oldTime().primitiveField()*mesh().S();
    }

    return tfam;
}


template<class Type>
tmp<faMatrix<Type>>
EulerFaDdtScheme<Type>::famDdt
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
            rho.dimensions()*vf.dimensions()*dimArea/dimTime
        )
    );
    faMatrix<Type>& fam = tfam.ref();

    scalar rDeltaT = 1.0/mesh().time().deltaTValue();

    fam.diag() = rDeltaT*rho.primitiveField()*mesh().S();

    if (mesh().moving())
    {
        fam.source() = rDeltaT
            *rho.oldTime().primitiveField()
            *vf.oldTime().primitiveField()*mesh().S0();
    }
    else
    {
        fam.source() = rDeltaT
            *rho.oldTime().primitiveField()
            *vf.oldTime().primitiveField()*mesh().S();
    }

    return tfam;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fa

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
