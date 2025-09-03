/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2014 OpenFOAM Foundation
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
    Foam::coordSystem::cylindrical

Description
    A cylindrical coordinate system (r-theta-z).
    The coordinate system angle theta is always in radians.

    \heading Dictionary entries
    \table
        Property    | Description                           | Required | Default
        type        | type name: cylindrical                | yes   |
    \endtable

SourceFiles
    cylindricalCS.C
    cylindricalCSTransform.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_cylindricalCS_H
#define Foam_cylindricalCS_H

#include "coordinateSystem.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace coordSystem
{

/*---------------------------------------------------------------------------*\
                  Class coordSystem::cylindrical Declaration
\*---------------------------------------------------------------------------*/

class cylindrical
:
    public coordinateSystem
{
protected:

    // Protected Member Functions

        //- From local coordinate system to the global Cartesian system
        //- with optional translation for the origin
        virtual vector localToGlobal
        (
            const vector& local,
            bool translate
        ) const;

        //- From local coordinate system to the global Cartesian system
        //- with optional translation for the origin
        virtual tmp<vectorField> localToGlobal
        (
            const vectorField& local,
            bool translate
        ) const;

        //- From global Cartesian system to the local coordinate system
        //- with optional translation for the origin
        virtual vector globalToLocal
        (
            const vector& global,
            bool translate
        ) const;

        //- Convert from global Cartesian system to the local coordinate system
        //- with optional translation for the origin
        virtual tmp<vectorField> globalToLocal
        (
            const vectorField& global,
            bool translate
        ) const;


public:

    //- Runtime type information
    TypeNameNoDebug("cylindrical");


    // Static Members

        //- Global (identity) cylindrical coordinate system
        static const cylindrical null;


    // Constructors

        //- Default construct. Identity coordinate system.
        cylindrical() = default;

        //- Copy construct
        cylindrical(const cylindrical& csys) = default;

        //- Move construct
        cylindrical(cylindrical&& csys) = default;

        //- Copy construct from another coordinateSystem type
        explicit cylindrical(const coordinateSystem& csys);

        //- Move construct from another coordinateSystem type
        explicit cylindrical(coordinateSystem&& csys);

        //- Move construct from autoPtr of another coordinateSystem type
        explicit cylindrical(autoPtr<coordinateSystem>&& csys);

        //- Copy construct from rotation with origin=0
        explicit cylindrical(const coordinateRotation& crot);

        //- Construct from origin and rotation
        cylindrical(const point& origin, const coordinateRotation& crot);

        //- Construct from origin and single axis
        cylindrical(const point& origin, const vector& axis);

        //- Construct named from origin and single axis
        cylindrical(const word& name, const point& origin, const vector& axis);

        //- Construct from origin and two axes
        cylindrical
        (
            const point& origin,
            const vector& axis,
            const vector& dirn
        );

        //- Construct named from origin and two axes
        cylindrical
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
        explicit cylindrical
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
        cylindrical
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
    virtual ~cylindrical() = default;


    // Member Functions

        //- Treat the rotation tensor as non-uniform
        virtual bool uniform() const
        {
            return false;
        }


    // Rotations

        //- Position-dependent rotation tensors at multiple points
        using coordinateSystem::R;

        //- Position-dependent rotation tensor at a single point
        //- \return tensor
        virtual tensor R(const point& global) const;


    // Member Operators

        //- Copy assignment
        cylindrical& operator=(const cylindrical&) = default;

        //- Move assignment
        cylindrical& operator=(cylindrical&&) = default;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace coordSystem

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Compatibility typedef 1806
typedef coordSystem::cylindrical cylindricalCS;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
