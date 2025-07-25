/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019,2023 OpenCFD Ltd.
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

InNamespace
    Foam::FieldOps

Description
    Various utility functions to work on geometric fields

SourceFiles
    GeometricFieldOps.H

\*---------------------------------------------------------------------------*/

#ifndef Foam_GeometricFieldOps_H
#define Foam_GeometricFieldOps_H

#include "FieldOps.H"
#include "GeometricField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace FieldOps
{

/*---------------------------------------------------------------------------*\
                        Namespace FieldOps Declarations
\*---------------------------------------------------------------------------*/

//- Populate a geometric field as the result of a unary operation on an input.
//  It is permissible for inputs/outputs to refer to the same field(s),
//  but partially overlapping regions are ill-defined.
template
<
    class Tout, class T1, class UnaryOp,
    template<class> class PatchField, class GeoMesh
>
void assign
(
    GeometricField<Tout, PatchField, GeoMesh>& result,
    const GeometricField<T1, PatchField, GeoMesh>& a,
    const UnaryOp& op
)
{
    FieldOps::assign
    (
        result.primitiveFieldRef(),
        a.primitiveField(),
        op
    );

    auto& bfld = result.boundaryFieldRef();

    const label len = bfld.size();

    for (label i = 0; i < len; ++i)
    {
        FieldOps::assign
        (
            bfld[i],
            a.boundaryField()[i],
            op
        );
    }
    result.correctLocalBoundaryConditions();
    if (GeometricBoundaryField<Tout, PatchField, GeoMesh>::debug)
    {
        result.boundaryField().check();
    }
}


//- Populate a geometric field from the binary operation on two inputs.
//  It is permissible for inputs/outputs to refer to the same field(s),
//  but partially overlapping regions are ill-defined.
template
<
    class Tout, class T1, class T2, class BinaryOp,
    template<class> class PatchField, class GeoMesh
>
void assign
(
    GeometricField<Tout, PatchField, GeoMesh>& result,
    const GeometricField<T1, PatchField, GeoMesh>& a,
    const GeometricField<T2, PatchField, GeoMesh>& b,
    const BinaryOp& bop
)
{
    FieldOps::assign
    (
        result.primitiveFieldRef(),
        a.primitiveField(),
        b.primitiveField(),
        bop
    );

    auto& bfld = result.boundaryFieldRef();

    const label len = bfld.size();

    for (label i = 0; i < len; ++i)
    {
        FieldOps::assign
        (
            bfld[i],
            a.boundaryField()[i],
            b.boundaryField()[i],
            bop
        );
    }
    result.correctLocalBoundaryConditions();
    if (GeometricBoundaryField<Tout, PatchField, GeoMesh>::debug)
    {
        result.boundaryField().check();
    }
}


//- Emulate a ternary operation, selecting values from a or b
//- depending on the binary predicate.
template
<
    class T, class BinaryOp,
    template<class> class PatchField, class GeoMesh
>
void ternary
(
    GeometricField<T, PatchField, GeoMesh>& result,
    const GeometricField<T, PatchField, GeoMesh>& a,
    const GeometricField<T, PatchField, GeoMesh>& b,
    const BinaryOp& bop
)
{
    FieldOps::ternary
    (
        result.primitiveFieldRef(),
        a.primitiveField(),
        b.primitiveField(),
        bop
    );

    auto& bfld = result.boundaryFieldRef();

    const label len = bfld.size();

    for (label i = 0; i < len; ++i)
    {
        FieldOps::ternary
        (
            bfld[i],
            a.boundaryField()[i],
            b.boundaryField()[i],
            bop
        );
    }
    result.correctLocalBoundaryConditions();
    if (GeometricBoundaryField<T, PatchField, GeoMesh>::debug)
    {
        result.boundaryField().check();
    }
}


//- Emulate a ternary operation, selecting field values from a or b
//- depending on the conditional.
//
//  Since boolean fields are not normally used, a flip operation is
//  a general requirement.
template
<
    class T, class BoolType, class FlipOp,
    template<class> class PatchField, class GeoMesh
>
void ternarySelect
(
    GeometricField<T, PatchField, GeoMesh>& result,
    const GeometricField<BoolType, PatchField, GeoMesh>& cond,
    const GeometricField<T, PatchField, GeoMesh>& a,
    const GeometricField<T, PatchField, GeoMesh>& b,
    const FlipOp& flip
)
{
    FieldOps::ternarySelect
    (
        result.primitiveFieldRef(),
        cond.primitiveField(),
        a.primitiveField(),
        b.primitiveField(),
        flip
    );

    auto& bfld = result.boundaryFieldRef();

    const label len = bfld.size();

    for (label i = 0; i < len; ++i)
    {
        FieldOps::ternarySelect
        (
            bfld[i],
            cond.boundaryField()[i],
            a.boundaryField()[i],
            b.boundaryField()[i],
            flip
        );
    }
    result.correctLocalBoundaryConditions();
    if (GeometricBoundaryField<T, PatchField, GeoMesh>::debug)
    {
        result.boundaryField().check();
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FieldOps
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
