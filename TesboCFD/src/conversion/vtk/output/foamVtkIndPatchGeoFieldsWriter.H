/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 OpenCFD Ltd.
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

Typedef
    Foam::vtk::indirectPatchGeoFieldsWriter

Description
    Write indirectPrimitivePatch faces/points
    as a vtp file or a legacy vtk file
    with support for geometric fields.

Typedef
    Foam::vtk::uindirectPatchGeoFieldsWriter

Description
    Write uindirectPrimitivePatch faces/points
    as a vtp file or a legacy vtk file
    with support for geometric fields.

See Also
    Foam::vtk::GenericPatchWriter

\*---------------------------------------------------------------------------*/

#ifndef Foam_vtk_IndPatchGeoFieldsWriter_H
#define Foam_vtk_IndPatchGeoFieldsWriter_H

#include "foamVtkGenericPatchGeoFieldsWriter.H"
#include "indirectPrimitivePatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace vtk
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

typedef
    GenericPatchGeoFieldsWriter<indirectPrimitivePatch>
    indirectPatchGeoFieldsWriter;

typedef
    GenericPatchGeoFieldsWriter<uindirectPrimitivePatch>
    uindirectPatchGeoFieldsWriter;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace vtk
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
