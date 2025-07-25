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

#include "backwardFaDdtScheme.H"
#include "faMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fa
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
scalar backwardFaDdtScheme<Type>::deltaT_() const
{
    return mesh().time().deltaTValue();
}


template<class Type>
scalar backwardFaDdtScheme<Type>::deltaT0_() const
{
    return mesh().time().deltaT0Value();
}


template<class Type>
template<class GeoField>
scalar backwardFaDdtScheme<Type>::deltaT0_(const GeoField& vf) const
{
    if (vf.oldTime().timeIndex() == vf.oldTime().oldTime().timeIndex())
    {
        return GREAT;
    }
    else
    {
        return mesh().time().deltaT0Value();
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh>>
backwardFaDdtScheme<Type>::facDdt
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

    scalar deltaT = deltaT_();
    scalar deltaT0 = deltaT0_();

    scalar coefft   = 1 + deltaT/(deltaT + deltaT0);
    scalar coefft00 = deltaT*deltaT/(deltaT0*(deltaT + deltaT0));
    scalar coefft0  = coefft + coefft00;

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

        tdtdt.ref().primitiveFieldRef() = rDeltaT.value()*dt.value()*
        (
            coefft - (coefft0*mesh().S0() - coefft00*mesh().S00())/mesh().S()
        );

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
backwardFaDdtScheme<Type>::facDdt0
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

    scalar deltaT = deltaT_();
    scalar deltaT0 = deltaT0_();

    scalar coefft   = 1 + deltaT/(deltaT + deltaT0);
    scalar coefft00 = deltaT*deltaT/(deltaT0*(deltaT + deltaT0));
    scalar coefft0  = coefft + coefft00;

    tmp<GeometricField<Type, faPatchField, areaMesh>> tdtdt0
    (
        new GeometricField<Type, faPatchField, areaMesh>
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


template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh>>
backwardFaDdtScheme<Type>::facDdt
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

    scalar deltaT = deltaT_();
    scalar deltaT0 = deltaT0_(vf);

    scalar coefft   = 1 + deltaT/(deltaT + deltaT0);
    scalar coefft00 = deltaT*deltaT/(deltaT0*(deltaT + deltaT0));
    scalar coefft0  = coefft + coefft00;

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
                    coefft*vf() -
                    (
                        coefft0*vf.oldTime()()*mesh().S0()
                      - coefft00*vf.oldTime().oldTime()()
                       *mesh().S00()
                    )/mesh().S()
                ),
                rDeltaT.value()*
                (
                    coefft*vf.boundaryField() -
                    (
                        coefft0*vf.oldTime().boundaryField()
                      - coefft00*vf.oldTime().oldTime().boundaryField()
                    )
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


template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh>>
backwardFaDdtScheme<Type>::facDdt0
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

    scalar deltaT = deltaT_();
    scalar deltaT0 = deltaT0_(vf);

    scalar coefft   = 1 + deltaT/(deltaT + deltaT0);
    scalar coefft00 = deltaT*deltaT/(deltaT0*(deltaT + deltaT0));
    scalar coefft0  = coefft + coefft00;

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
                  - (
                        coefft0*vf.oldTime()()*mesh().S0()
                      - coefft00*vf.oldTime().oldTime()()
                       *mesh().S00()
                    )/mesh().S()
                ),
                rDeltaT.value()*
                (
                  - (
                        coefft0*vf.oldTime().boundaryField()
                      - coefft00*vf.oldTime().oldTime().boundaryField()
                    )
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
                rDeltaT*
                (
                  - coefft0*vf.oldTime()
                  + coefft00*vf.oldTime().oldTime()
                )
            )
        );
    }
}


template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh>>
backwardFaDdtScheme<Type>::facDdt
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

    scalar deltaT = deltaT_();
    scalar deltaT0 = deltaT0_(vf);

    scalar coefft   = 1 + deltaT/(deltaT + deltaT0);
    scalar coefft00 = deltaT*deltaT/(deltaT0*(deltaT + deltaT0));
    scalar coefft0  = coefft + coefft00;

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
                    coefft*vf.internalField() -
                    (
                        coefft0*vf.oldTime()()*mesh().S0()
                      - coefft00*vf.oldTime().oldTime()()
                       *mesh().S00()
                    )/mesh().S()
                ),
                rDeltaT.value()*rho.value()*
                (
                    coefft*vf.boundaryField() -
                    (
                        coefft0*vf.oldTime().boundaryField()
                      - coefft00*vf.oldTime().oldTime().boundaryField()
                    )
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

template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh>>
backwardFaDdtScheme<Type>::facDdt0
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

    scalar deltaT = deltaT_();
    scalar deltaT0 = deltaT0_(vf);

    scalar coefft   = 1 + deltaT/(deltaT + deltaT0);
    scalar coefft00 = deltaT*deltaT/(deltaT0*(deltaT + deltaT0));
    scalar coefft0  = coefft + coefft00;

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
                   -(
                        coefft0*vf.oldTime()()*mesh().S0()
                      - coefft00*vf.oldTime().oldTime()()
                       *mesh().S00()
                    )/mesh().S()
                ),
                rDeltaT.value()*rho.value()*
                (
                   -(
                        coefft0*vf.oldTime().boundaryField()
                      - coefft00*vf.oldTime().oldTime().boundaryField()
                    )
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
                rDeltaT*rho*
                (
                  - coefft0*vf.oldTime()
                 + coefft00*vf.oldTime().oldTime()
                )
            )
        );
    }
}


template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh>>
backwardFaDdtScheme<Type>::facDdt
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

    scalar deltaT = deltaT_();
    scalar deltaT0 = deltaT0_(vf);

    scalar coefft   = 1 + deltaT/(deltaT + deltaT0);
    scalar coefft00 = deltaT*deltaT/(deltaT0*(deltaT + deltaT0));
    scalar coefft0  = coefft + coefft00;

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
                    coefft*rho.internalField()*vf.internalField() -
                    (
                        coefft0*rho.oldTime()()
                       *vf.oldTime()()*mesh().S0()
                      - coefft00*rho.oldTime().oldTime()()
                       *vf.oldTime().oldTime()()*mesh().S00()
                    )/mesh().S()
                ),
                rDeltaT.value()*
                (
                    coefft*vf.boundaryField() -
                    (
                        coefft0*rho.oldTime().boundaryField()
                       *vf.oldTime().boundaryField()
                      - coefft00*rho.oldTime().oldTime().boundaryField()
                       *vf.oldTime().oldTime().boundaryField()
                    )
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


template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh>>
backwardFaDdtScheme<Type>::facDdt0
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

    scalar deltaT = deltaT_();
    scalar deltaT0 = deltaT0_(vf);

    scalar coefft   = 1 + deltaT/(deltaT + deltaT0);
    scalar coefft00 = deltaT*deltaT/(deltaT0*(deltaT + deltaT0));
    scalar coefft0  = coefft + coefft00;

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
                  - (
                        coefft0*rho.oldTime()()
                       *vf.oldTime()()*mesh().S0()
                      - coefft00*rho.oldTime().oldTime()()
                       *vf.oldTime().oldTime()()*mesh().S00()
                    )/mesh().S()
                ),
                rDeltaT.value()*
                (
                  - (
                        coefft0*rho.oldTime().boundaryField()
                       *vf.oldTime().boundaryField()
                      - coefft00*rho.oldTime().oldTime().boundaryField()
                       *vf.oldTime().oldTime().boundaryField()
                    )
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
                rDeltaT*
                (
                  - coefft0*rho.oldTime()*vf.oldTime()
                  + coefft00*rho.oldTime().oldTime()*vf.oldTime().oldTime()
                )
            )
        );
    }
}


template<class Type>
tmp<faMatrix<Type>>
backwardFaDdtScheme<Type>::famDdt
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

    scalar rDeltaT = 1.0/deltaT_();

    scalar deltaT = deltaT_();
    scalar deltaT0 = deltaT0_(vf);

    scalar coefft   = 1 + deltaT/(deltaT + deltaT0);
    scalar coefft00 = deltaT*deltaT/(deltaT0*(deltaT + deltaT0));
    scalar coefft0  = coefft + coefft00;

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


template<class Type>
tmp<faMatrix<Type>>
backwardFaDdtScheme<Type>::famDdt
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

    scalar rDeltaT = 1.0/deltaT_();

    scalar deltaT = deltaT_();
    scalar deltaT0 = deltaT0_(vf);

    scalar coefft   = 1 + deltaT/(deltaT + deltaT0);
    scalar coefft00 = deltaT*deltaT/(deltaT0*(deltaT + deltaT0));
    scalar coefft0  = coefft + coefft00;

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


template<class Type>
tmp<faMatrix<Type>>
backwardFaDdtScheme<Type>::famDdt
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

    scalar rDeltaT = 1.0/deltaT_();

    scalar deltaT = deltaT_();
    scalar deltaT0 = deltaT0_(vf);

    scalar coefft   = 1 + deltaT/(deltaT + deltaT0);
    scalar coefft00 = deltaT*deltaT/(deltaT0*(deltaT + deltaT0));
    scalar coefft0  = coefft + coefft00;

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

} // End namespace fa

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
