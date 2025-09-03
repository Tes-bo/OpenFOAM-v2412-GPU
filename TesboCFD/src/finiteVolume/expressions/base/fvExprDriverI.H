/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2010-2018 Bernhard Gschaider
    Copyright (C) 2019-2020 OpenCFD Ltd.
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

#include "Time.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline bool Foam::expressions::fvExprDriver::hasVariable
(
    const word& name
) const
{
    return delayedVariables_.found(name) || variables_.found(name);
}


inline const Foam::expressions::exprResult&
Foam::expressions::fvExprDriver::variable
(
    const word& name
) const
{
    if (delayedVariables_.found(name))
    {
        return delayedVariables_[name];
    }

    return variables_[name];
}


inline Foam::expressions::exprResult&
Foam::expressions::fvExprDriver::variable
(
    const word& name
)
{
    if (delayedVariables_.found(name))
    {
        return delayedVariables_[name];
    }

    return variables_[name];
}


template<class Type>
inline bool Foam::expressions::fvExprDriver::isVariable
(
    const word& name,
    const bool wantPointData,
    const label expectedSize
) const
{
    return
    (
        this->isLocalVariable<Type>(name, wantPointData, expectedSize)
     || this->isGlobalVariable<Type>(name, wantPointData, expectedSize)
    );
}


template<class Type>
inline bool Foam::expressions::fvExprDriver::isVariableOrField
(
    const word& name,
    const bool wantPointData,
    const label expectedSize
)
const
{
    return
    (
        this->isVariable<Type>(name, wantPointData, expectedSize)
     || this->isField<Type>(name, wantPointData)
    );
}


template<class GeomField>
inline Foam::tmp<GeomField>
Foam::expressions::fvExprDriver::getOrReadField
(
    const word& name,
    const bool mandatory,
    const bool getOldTime
)
{
    return this->getOrReadFieldImpl<GeomField>
    (
        name,
        this->mesh(),
        mandatory,
        getOldTime
    );
}


template<class GeomField>
inline Foam::tmp<GeomField>
Foam::expressions::fvExprDriver::getOrReadPointField
(
    const word& name,
    const bool mandatory,
    const bool getOldTime
)
{
    return this->getOrReadFieldImpl<GeomField>
    (
        name,
        pointMesh::New(this->mesh()),
        mandatory,
        getOldTime
    );
}


template<class GeomField, class Mesh>
inline Foam::tmp<GeomField>
Foam::expressions::fvExprDriver::readAndRegister
(
    const word& name,
    const Mesh& meshRef
)
{
    GeomField* ptr = new GeomField
    (
        IOobject
        (
            name,
            meshRef.thisDb().time().timeName(),
            meshRef.thisDb(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        meshRef
    );

    if (cacheReadFields())
    {
        DebugInfo
            << "Registering a copy of " << name << " with mesh" << nl;

        // This is clunky
        ptr->checkIn();
        return tmp<GeomField>(regIOobject::store(ptr));
    }

    return tmp<GeomField>(ptr);
}


// ************************************************************************* //
