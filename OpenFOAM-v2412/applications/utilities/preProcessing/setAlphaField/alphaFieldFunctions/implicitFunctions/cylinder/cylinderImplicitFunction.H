/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019 OpenCFD Ltd.
    Copyright (C) 2019-2020 DLR
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

Class
    Foam::implicitFunctions::cylinderImplicitFunction

Description
    creates a infintite long cylinderImplicitFunction

    Original code supplied by Henning Scheufler, DLR (2019)

SourceFiles
    cylinderImplicitFunction.C

\*---------------------------------------------------------------------------*/

#ifndef implicitFunction_cylinderImplicitFunction_H
#define implicitFunction_cylinderImplicitFunction_H

#include "implicitFunction.H"
#include "point.H"
#include "tensor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace implicitFunctions
{

/*---------------------------------------------------------------------------*\
                   Class cylinderImplicitFunction Declaration
\*---------------------------------------------------------------------------*/

class cylinderImplicitFunction
:
    public implicitFunction
{
private:

    // Private Data

        //- Origin point
        const point origin_;

        //- Radius
        const scalar radius_;

        const scalar scale_;

        const vector direction_;

        const tensor project_;


public:

    //- Runtime type information
    TypeName("cylinder");


    // Constructors

        //- Construct from components
        cylinderImplicitFunction
        (
            const point& origin,
            const scalar radius,
            const scalar scale,
            const vector& direction
        );

        //- Construct from dictionary
        explicit cylinderImplicitFunction(const dictionary& dict);


    //- Destructor
    virtual ~cylinderImplicitFunction() = default;


    // Member Functions

    virtual scalar value(const vector& p) const
    {
        return (-mag(project_ & (p - origin_)) + radius_)*scale_;
    }

    virtual vector grad(const vector& p) const
    {
        return -(project_ & (p - origin_))*scale_;
    }

    virtual scalar distanceToSurfaces(const vector& p) const
    {
        return mag(mag(project_ & (p - origin_)) - radius_)*scale_;
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace implicitFunctions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
