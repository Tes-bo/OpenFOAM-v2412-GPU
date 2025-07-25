/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015 OpenFOAM Foundation
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
    Foam::globalIOFields

Description
    Typedefs for globalIOFields of primitive types. These are fully compatible
    with 'normal' IOFields except have global filePath() scope.

\*---------------------------------------------------------------------------*/

#ifndef Foam_globalIOFields_H
#define Foam_globalIOFields_H

#include "primitiveFields.H"
#include "GlobalIOField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    typedef GlobalIOField<label> labelGlobalIOField;
    typedef GlobalIOField<scalar> scalarGlobalIOField;
    typedef GlobalIOField<vector> vectorGlobalIOField;
    typedef GlobalIOField<sphericalTensor> sphericalTensorGlobalIOField;
    typedef GlobalIOField<symmTensor> symmTensorGlobalIOField;
    typedef GlobalIOField<tensor> tensorGlobalIOField;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
