/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2016 OpenFOAM Foundation
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

#include "flipOp.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

#undef  defineNegateOp
#define defineNegateOp(Type)                                                  \
                                                                              \
    template<> Type Foam::flipOp::operator()(const Type& val) const           \
    {                                                                         \
        return -val;                                                          \
    }


defineNegateOp(Foam::scalar);
defineNegateOp(Foam::vector);
defineNegateOp(Foam::sphericalTensor);
defineNegateOp(Foam::symmTensor);
defineNegateOp(Foam::tensor);
defineNegateOp(Foam::triad);

#undef defineNegateOp


// ************************************************************************* //
