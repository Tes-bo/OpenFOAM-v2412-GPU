/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2013 OpenFOAM Foundation
    Copyright (C) 2018-2024 OpenCFD Ltd.
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
    Foam::coordSystem::cartesian

Description
    A Cartesian coordinate system

    \heading Dictionary entries
    \table
        Property    | Description                           | Required | Default
        type        | Type name: cartesian                  | yes   |
    \endtable

SourceFiles
    cartesianCS.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_cartesianCS_H
#define Foam_cartesianCS_H

#include "coordinateSystem.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace coordSystem
{

/*---------------------------------------------------------------------------*\
                   Class coordSystem::cartesian Declaration
\*---------------------------------------------------------------------------*/

class cartesian
:
    public coordinateSystem
{
public:

    //- Runtime type information
    TypeNameNoDebug("cartesian");


    // Static Members

        //- Global (identity) cartesian coordinate system
        static const cartesian null;


    // Constructors

        //- Default construct. Identity coordinate system.
        cartesian() = default;

        //- Copy construct
        cartesian(const cartesian& csys) = default;

        //- Move construct
        cartesian(cartesian&& csys) = default;

        //- Copy construct from another coordinateSystem type
        explicit cartesian(const coordinateSystem& csys);

        //- Move construct from another coordinateSystem type
        explicit cartesian(coordinateSystem&& csys);

        //- Move construct from autoPtr of another coordinateSystem type
        explicit cartesian(autoPtr<coordinateSystem>&& csys);

        //- Copy construct from rotation with origin=0
        explicit cartesian(const coordinateRotation& crot);

        //- Move construct from rotation with origin=0
        explicit cartesian(coordinateRotation&& crot);

        //- Construct from origin and rotation
        cartesian(const point& origin, const coordinateRotation& crot);

        //- Construct from origin and 2 axes
        cartesian
        (
            const point& origin,
            const vector& axis,
            const vector& dirn
        );

        //- Construct named from origin and 2 axes
        cartesian
        (
            const word& name,
            const point& origin,
            const vector& axis,
            const vector& dirn
        );

        //- Construct from dictionary with optional
        //- read handling for the 'origin' entry (default: MUST_READ).
        //
        //  \note The readOrigin is downgraded to READ_IF_PRESENT
        //  if the dictionary itself is "coordinateSystem"
        explicit cartesian
        (
            const dictionary& dict,
            IOobjectOption::readOption readOrigin = IOobjectOption::MUST_READ
        );

        //- Construct from dictionary with optional subDict lookup and optional
        //- read handling for the 'origin' entry (default: MUST_READ).
        //
        //  \param dictName If non-empty, mandatory sub-dictionary to use.
        //
        //  \note The readOrigin is downgraded to READ_IF_PRESENT
        //  if the dictionary itself is "coordinateSystem"
        //  or if a sub-dictionary is being used
        cartesian
        (
            const dictionary& dict,
            const word& dictName,
            IOobjectOption::readOption readOrigin = IOobjectOption::MUST_READ
        );


        //- Return clone
        virtual autoPtr<coordinateSystem> clone() const
        {
            return coordinateSystem::Clone(*this);
        }


    //- Destructor
    virtual ~cartesian() = default;


    // Member Operators

        //- Copy assignment
        cartesian& operator=(const cartesian&) = default;

        //- Move assignment
        cartesian& operator=(cartesian&&) = default;

        //- Copy/move assignment from coordinateSystem, autoPtr
        using coordinateSystem::operator=;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace coordSystem

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Compatibility typedef 1806
typedef coordSystem::cartesian cartesianCS;

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
