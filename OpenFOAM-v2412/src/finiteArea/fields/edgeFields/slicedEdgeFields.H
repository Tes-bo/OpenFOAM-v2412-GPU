/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2023 OpenCFD Ltd.
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

InClass
    Foam::slicedEdgeFields

\*---------------------------------------------------------------------------*/

#ifndef Foam_slicedEdgeFields_H
#define Foam_slicedEdgeFields_H

#include "SlicedGeometricField.H"
#include "slicedFaePatchField.H"
#include "edgeFaMesh.H"
#include "slicedEdgeFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<>
inline void
SlicedGeometricField<scalar, faePatchField, slicedFaePatchField, edgeMesh>::
correctBoundaryConditions()
{}

template<>
inline void
SlicedGeometricField<vector, faePatchField, slicedFaePatchField, edgeMesh>::
correctBoundaryConditions()
{}

template<>
inline void
SlicedGeometricField
<
    sphericalTensor,
    faePatchField,
    slicedFaePatchField,
    edgeMesh
>::
correctBoundaryConditions()
{}

template<>
inline void
SlicedGeometricField
<
    symmTensor,
    faePatchField,
    slicedFaePatchField,
    edgeMesh
>::
correctBoundaryConditions()
{}

template<>
inline void
SlicedGeometricField<tensor, faePatchField, slicedFaePatchField, edgeMesh>::
correctBoundaryConditions()
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
