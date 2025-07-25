/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2015-2023 OpenCFD Ltd.
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

#include "DimensionedField.H"
#include "dimensionedType.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Check that both fields use the same mesh
#undef  checkField
#define checkField(fld1, fld2, op)                                  \
if (&(fld1).mesh() != &(fld2).mesh())                               \
{                                                                   \
    FatalErrorInFunction                                            \
        << "Different mesh for fields "                             \
        << (fld1).name() << " and " << (fld2).name()                \
        << " during operation " << op                               \
        << abort(FatalError);                                       \
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type, class GeoMesh>
void Foam::DimensionedField<Type, GeoMesh>::checkFieldSize() const
{
    const label fieldSize = this->size();
    if (fieldSize)
    {
        const label meshSize = GeoMesh::size(this->mesh_);
        if (fieldSize != meshSize)
        {
            FatalErrorInFunction
                << "size of field = " << fieldSize
                << " is not the same as the size of mesh = "
                << meshSize
                << abort(FatalError);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type, class GeoMesh>
Foam::DimensionedField<Type, GeoMesh>::DimensionedField
(
    const IOobject& io,
    const Mesh& mesh,
    const dimensionSet& dims,
    const Field<Type>& field
)
:
    regIOobject(io),
    Field<Type>(field),
    mesh_(mesh),
    dimensions_(dims),
    oriented_()
{
    checkFieldSize();
}


template<class Type, class GeoMesh>
Foam::DimensionedField<Type, GeoMesh>::DimensionedField
(
    const IOobject& io,
    const Mesh& mesh,
    const dimensionSet& dims,
    Field<Type>&& field
)
:
    regIOobject(io),
    Field<Type>(std::move(field)),
    mesh_(mesh),
    dimensions_(dims)
{
    checkFieldSize();
}


template<class Type, class GeoMesh>
Foam::DimensionedField<Type, GeoMesh>::DimensionedField
(
    const IOobject& io,
    const Mesh& mesh,
    const dimensionSet& dims,
    List<Type>&& field
)
:
    regIOobject(io),
    Field<Type>(std::move(field)),
    mesh_(mesh),
    dimensions_(dims)
{
    checkFieldSize();
}


template<class Type, class GeoMesh>
Foam::DimensionedField<Type, GeoMesh>::DimensionedField
(
    const IOobject& io,
    const Mesh& mesh,
    const dimensionSet& dims,
    const tmp<Field<Type>>& tfield
)
:
    regIOobject(io),
    Field<Type>(tfield.constCast(), tfield.movable()),
    mesh_(mesh),
    dimensions_(dims),
    oriented_()
{
    tfield.clear();
    checkFieldSize();
}


template<class Type, class GeoMesh>
Foam::DimensionedField<Type, GeoMesh>::DimensionedField
(
    const IOobject& io,
    const Mesh& mesh,
    const dimensionSet& dims,
    const bool checkIOFlags
)
:
    regIOobject(io),
    Field<Type>(GeoMesh::size(mesh)),
    mesh_(mesh),
    dimensions_(dims),
    oriented_()
{
    if (checkIOFlags)
    {
        readIfPresent();
    }
}


template<class Type, class GeoMesh>
Foam::DimensionedField<Type, GeoMesh>::DimensionedField
(
    const IOobject& io,
    const Mesh& mesh,
    const Type& value,
    const dimensionSet& dims,
    const bool checkIOFlags
)
:
    regIOobject(io),
    Field<Type>(GeoMesh::size(mesh), value),
    mesh_(mesh),
    dimensions_(dims),
    oriented_()
{
    if (checkIOFlags)
    {
        readIfPresent();
    }
}


template<class Type, class GeoMesh>
Foam::DimensionedField<Type, GeoMesh>::DimensionedField
(
    const IOobject& io,
    const Mesh& mesh,
    const dimensioned<Type>& dt,
    const bool checkIOFlags
)
:
    DimensionedField<Type, GeoMesh>
    (
        io,
        mesh,
        dt.value(),
        dt.dimensions(),
        checkIOFlags
    )
{}


template<class Type, class GeoMesh>
Foam::DimensionedField<Type, GeoMesh>::DimensionedField
(
    const DimensionedField<Type, GeoMesh>& df
)
:
    regIOobject(df),
    Field<Type>(df),
    mesh_(df.mesh_),
    dimensions_(df.dimensions_),
    oriented_(df.oriented_)
{}


template<class Type, class GeoMesh>
Foam::DimensionedField<Type, GeoMesh>::DimensionedField
(
    DimensionedField<Type, GeoMesh>&& df
)
:
    DimensionedField<Type, GeoMesh>(df, true)
{}


template<class Type, class GeoMesh>
Foam::DimensionedField<Type, GeoMesh>::DimensionedField
(
    DimensionedField<Type, GeoMesh>& df,
    bool reuse
)
:
    regIOobject(df, reuse),
    Field<Type>(df, reuse),
    mesh_(df.mesh_),
    dimensions_(df.dimensions_),
    oriented_(df.oriented_)
{}


template<class Type, class GeoMesh>
Foam::DimensionedField<Type, GeoMesh>::DimensionedField
(
    const tmp<DimensionedField<Type, GeoMesh>>& tdf
)
:
    DimensionedField<Type, GeoMesh>(tdf.constCast(), tdf.movable())
{
    tdf.clear();
}


template<class Type, class GeoMesh>
Foam::DimensionedField<Type, GeoMesh>::DimensionedField
(
    const IOobject& io,
    const DimensionedField<Type, GeoMesh>& df
)
:
    regIOobject(io),
    Field<Type>(df),
    mesh_(df.mesh_),
    dimensions_(df.dimensions_),
    oriented_(df.oriented_)
{}


template<class Type, class GeoMesh>
Foam::DimensionedField<Type, GeoMesh>::DimensionedField
(
    const IOobject& io,
    DimensionedField<Type, GeoMesh>&& df
)
:
    DimensionedField<Type, GeoMesh>(io, df, true)
{}


template<class Type, class GeoMesh>
Foam::DimensionedField<Type, GeoMesh>::DimensionedField
(
    const IOobject& io,
    DimensionedField<Type, GeoMesh>& df,
    bool reuse
)
:
    regIOobject(io, df),
    Field<Type>(df, reuse),
    mesh_(df.mesh_),
    dimensions_(df.dimensions_),
    oriented_(df.oriented_)
{}


template<class Type, class GeoMesh>
Foam::DimensionedField<Type, GeoMesh>::DimensionedField
(
    const IOobject& io,
    const tmp<DimensionedField<Type, GeoMesh>>& tdf
)
:
    DimensionedField<Type, GeoMesh>(io, tdf.constCast(), tdf.movable())
{
    tdf.clear();
}


template<class Type, class GeoMesh>
Foam::DimensionedField<Type, GeoMesh>::DimensionedField
(
    const word& newName,
    const DimensionedField<Type, GeoMesh>& df
)
:
    regIOobject(newName, df, newName != df.name()),
    Field<Type>(df),
    mesh_(df.mesh_),
    dimensions_(df.dimensions_),
    oriented_(df.oriented_)
{}


template<class Type, class GeoMesh>
Foam::DimensionedField<Type, GeoMesh>::DimensionedField
(
    const word& newName,
    DimensionedField<Type, GeoMesh>&& df
)
:
    DimensionedField<Type, GeoMesh>(newName, df, true)
{}


template<class Type, class GeoMesh>
Foam::DimensionedField<Type, GeoMesh>::DimensionedField
(
    const word& newName,
    DimensionedField<Type, GeoMesh>& df,
    bool reuse
)
:
    regIOobject(newName, df, true),
    Field<Type>(df, reuse),
    mesh_(df.mesh_),
    dimensions_(df.dimensions_),
    oriented_(df.oriented_)
{}


template<class Type, class GeoMesh>
Foam::DimensionedField<Type, GeoMesh>::DimensionedField
(
    const word& newName,
    const tmp<DimensionedField<Type, GeoMesh>>& tdf
)
:
    DimensionedField<Type, GeoMesh>(newName, tdf.constCast(), tdf.movable())
{
    tdf.clear();
}


template<class Type, class GeoMesh>
Foam::tmp<Foam::DimensionedField<Type, GeoMesh>>
Foam::DimensionedField<Type, GeoMesh>::clone() const
{
    return tmp<DimensionedField<Type, GeoMesh>>::New(*this);
}


// * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * * * //

template<class Type, class GeoMesh>
Foam::DimensionedField<Type, GeoMesh>::~DimensionedField()
{
    // FUTURE: register cache field info
    // // this->db().cacheTemporaryObject(*this);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type, class GeoMesh>
Foam::tmp
<
    Foam::DimensionedField
    <
        typename Foam::DimensionedField<Type, GeoMesh>::cmptType, GeoMesh
    >
>
Foam::DimensionedField<Type, GeoMesh>::component
(
    const direction d
) const
{
    auto tresult = DimensionedField<cmptType, GeoMesh>::New
    (
        name() + ".component(" + ::Foam::name(d) + ')',
        mesh_,
        dimensions_
    );

    Foam::component(tresult.ref(), *this, d);

    return tresult;
}


template<class Type, class GeoMesh>
void Foam::DimensionedField<Type, GeoMesh>::replace
(
    const direction d,
    const DimensionedField
    <
        typename DimensionedField<Type, GeoMesh>::cmptType, GeoMesh
    >& df
)
{
    Field<Type>::replace(d, df);
}


template<class Type, class GeoMesh>
void Foam::DimensionedField<Type, GeoMesh>::replace
(
    const direction d,
    const tmp
    <
        DimensionedField
        <
            typename DimensionedField<Type, GeoMesh>::cmptType, GeoMesh
        >
    >& tdf
)
{
    replace(d, tdf());
    tdf.clear();
}


template<class Type, class GeoMesh>
Foam::tmp<Foam::DimensionedField<Type, GeoMesh>>
Foam::DimensionedField<Type, GeoMesh>::T() const
{
    auto tresult = DimensionedField<Type, GeoMesh>::New
    (
        name() + ".T()",
        mesh_,
        dimensions_
    );

    Foam::T(tresult.ref(), *this);

    return tresult;
}


template<class Type, class GeoMesh>
Foam::dimensioned<Type> Foam::DimensionedField<Type, GeoMesh>::average() const
{
    return
        dimensioned<Type>
        (
            this->name() + ".average()",
            this->dimensions(),
            gAverage(field())
        );
}


template<class Type, class GeoMesh>
Foam::dimensioned<Type> Foam::DimensionedField<Type, GeoMesh>::weightedAverage
(
    const DimensionedField<scalar, GeoMesh>& weightField
) const
{
    return
        dimensioned<Type>
        (
            this->name() + ".weightedAverage(weights)",
            this->dimensions(),
            gSum(weightField*field())/gSum(weightField)
        );
}


template<class Type, class GeoMesh>
Foam::dimensioned<Type> Foam::DimensionedField<Type, GeoMesh>::weightedAverage
(
    const tmp<DimensionedField<scalar, GeoMesh>>& tweightField
) const
{
    dimensioned<Type> result = weightedAverage(tweightField());
    tweightField.clear();
    return result;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type, class GeoMesh>
void Foam::DimensionedField<Type, GeoMesh>::operator=
(
    const DimensionedField<Type, GeoMesh>& df
)
{
    if (this == &df)
    {
        return;  // Self-assignment is a no-op
    }

    checkField(*this, df, "=");

    dimensions_ = df.dimensions();
    oriented_ = df.oriented();
    Field<Type>::operator=(df);
}


template<class Type, class GeoMesh>
void Foam::DimensionedField<Type, GeoMesh>::operator=
(
    const tmp<DimensionedField<Type, GeoMesh>>& tdf
)
{
    auto& df = tdf.constCast();

    if (this == &df)
    {
        return;  // Self-assignment is a no-op
    }

    checkField(*this, df, "=");

    dimensions_ = df.dimensions();
    oriented_ = df.oriented();
    this->transfer(df);
    tdf.clear();
}


template<class Type, class GeoMesh>
void Foam::DimensionedField<Type, GeoMesh>::operator=
(
    const dimensioned<Type>& dt
)
{
    dimensions_ = dt.dimensions();
    Field<Type>::operator=(dt.value());
}


#define COMPUTED_ASSIGNMENT(TYPE, op)                                          \
                                                                               \
template<class Type, class GeoMesh>                                            \
void Foam::DimensionedField<Type, GeoMesh>::operator op                        \
(                                                                              \
    const DimensionedField<TYPE, GeoMesh>& df                                  \
)                                                                              \
{                                                                              \
    checkField(*this, df, #op);                                                \
                                                                               \
    dimensions_ op df.dimensions();                                            \
    oriented_ op df.oriented();                                                \
    Field<Type>::operator op(df);                                              \
}                                                                              \
                                                                               \
template<class Type, class GeoMesh>                                            \
void Foam::DimensionedField<Type, GeoMesh>::operator op                        \
(                                                                              \
    const tmp<DimensionedField<TYPE, GeoMesh>>& tdf                            \
)                                                                              \
{                                                                              \
    operator op(tdf());                                                        \
    tdf.clear();                                                               \
}                                                                              \
                                                                               \
template<class Type, class GeoMesh>                                            \
void Foam::DimensionedField<Type, GeoMesh>::operator op                        \
(                                                                              \
    const dimensioned<TYPE>& dt                                                \
)                                                                              \
{                                                                              \
    dimensions_ op dt.dimensions();                                            \
    Field<Type>::operator op(dt.value());                                      \
}

COMPUTED_ASSIGNMENT(Type, +=)
COMPUTED_ASSIGNMENT(Type, -=)
COMPUTED_ASSIGNMENT(scalar, *=)
COMPUTED_ASSIGNMENT(scalar, /=)

#undef COMPUTED_ASSIGNMENT

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#undef checkField

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "DimensionedFieldIO.C"
#include "DimensionedFieldNew.C"
#include "DimensionedFieldFunctions.C"

// ************************************************************************* //
