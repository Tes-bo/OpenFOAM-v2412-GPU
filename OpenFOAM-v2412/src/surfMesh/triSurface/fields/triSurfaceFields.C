/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2022 OpenCFD Ltd.
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

#include "triSurfaceFields.H"
#include "className.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

// Face fields
defineTemplateTypeName(triSurfaceLabelField);
defineTemplateTypeName(triSurfaceScalarField);
defineTemplateTypeName(triSurfaceVectorField);
defineTemplateTypeName(triSurfaceSphericalTensorField);
defineTemplateTypeName(triSurfaceSymmTensorField);
defineTemplateTypeName(triSurfaceTensorField);


// Point fields

defineTemplateTypeName(triSurfacePointLabelField);
defineTemplateTypeName(triSurfacePointScalarField);
defineTemplateTypeName(triSurfacePointVectorField);
defineTemplateTypeName(triSurfacePointSphericalTensorField);
defineTemplateTypeName(triSurfacePointSymmTensorField);
defineTemplateTypeName(triSurfacePointTensorField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// ************************************************************************* //
