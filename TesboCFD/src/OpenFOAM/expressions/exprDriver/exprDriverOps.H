/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019 OpenCFD Ltd.
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

Namespace
    Foam::expressions::exprDriverOps

Description
    Miscellaneous operations used in combination with the
    expressions::exprDriver

SourceFiles
    exprDriverOps.H

\*---------------------------------------------------------------------------*/

#ifndef expressions_exprDriverOps_H
#define expressions_exprDriverOps_H

#include "scalarOps.H"
#include "FieldOps.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace expressions
{
namespace exprDriverOps
{

/*---------------------------------------------------------------------------*\
                     Namespace exprDriverOps Declarations
\*---------------------------------------------------------------------------*/

//- A modulo operation similar to that used in swak4foam
//  This operation produces values in the -1/2 to +1/2 range.
struct swakModuloOp
{
    scalar operator()(const scalar& a, const scalar& b) const
    {
        if (Foam::mag(b) < ROOTVSMALL)
        {
            return 0;
        }

        const scalar val = std::fmod(a, b);

        if (Foam::mag(val) > (0.5*b))
        {
            if (val > 0)
            {
                return (val - b);
            }
            else
            {
                return (val + b);
            }
        }

        return val;
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace exprDriverOps
} // End namespace expressions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
