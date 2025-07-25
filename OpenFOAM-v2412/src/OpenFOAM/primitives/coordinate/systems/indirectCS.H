/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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
    Foam::coordSystem::indirect

Description
    A coordinate system forward to a global coordinate system that is
    normally provided by the constant/coordinateSystems file.

    \heading Dictionary entries
    \table
        Property    | Description                           | Required | Default
        type        | Type name: indirect                   | yes   |
        name        | Name of the referenced system         | yes   |
    \endtable

SourceFiles
    indirectCS.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_indirectCS_H
#define Foam_indirectCS_H

#include "coordinateSystem.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace coordSystem
{

/*---------------------------------------------------------------------------*\
                    Class coordSystem::indirect Declaration
\*---------------------------------------------------------------------------*/

class indirect
:
    public coordinateSystem
{
    // Private Data

        //- The real coordinate system
        const coordinateSystem* backend_;


protected:

    // Protected Member Functions

        //- Convert from local coordinate system to the global Cartesian system
        //- with optional translation for the origin
        virtual vector localToGlobal
        (
            const vector& local,
            bool translate
        ) const
        {
            return backend_->localToGlobal(local, translate);
        }

        //- Convert from local coordinate system to the global Cartesian system
        //- with optional translation for the origin
        virtual tmp<vectorField> localToGlobal
        (
            const vectorField& local,
            bool translate
        ) const
        {
            return backend_->localToGlobal(local, translate);
        }

        //- Convert from global Cartesian system to the local coordinate system
        //- with optional translation for the origin
        virtual vector globalToLocal
        (
            const vector& global,
            bool translate
        ) const
        {
            return backend_->globalToLocal(global, translate);
        }

        //- Convert from global Cartesian system to the local coordinate system
        //- with optional translation for the origin
        virtual tmp<vectorField> globalToLocal
        (
            const vectorField& global,
            bool translate
        ) const
        {
            return backend_->globalToLocal(global, translate);
        }


public:

    //- Runtime type information
    TypeName("indirect");


    // Constructors

        //- Default construct is disallowed
        indirect() = delete;

        //- Copy construct
        indirect(const indirect& csys);

        //- Move construct
        indirect(indirect&& csys);

        //- Construct from global lookup
        indirect(const objectRegistry& obr, const word& name);

        //- Construct from global lookup.
        //  The readOrigin is unused.
        indirect
        (
            const objectRegistry& obr,
            const dictionary& dict,
            IOobjectOption::readOption readOrigin = IOobjectOption::NO_READ
        );

        //- Return clone
        virtual autoPtr<coordinateSystem> clone() const
        {
            return coordinateSystem::Clone(*this);
        }


    //- Destructor
    virtual ~indirect() = default;


    // Member Functions

    // Characteristics

        //- Is coordinate system good/valid?
        virtual bool good() const { return backend_ && backend_->good(); }

        //- True if the rotation tensor is uniform for all positions
        virtual bool uniform() const
        {
            return !backend_ || backend_->uniform();
        }

        //- Same as good() - 2023-07
        virtual bool valid() const { return this->good(); }


    // Access

        //- Reference to the underlying coordinate system
        virtual const coordinateSystem& cs() const
        {
            return *backend_;
        }

        //- Return origin
        virtual const point& origin() const
        {
            return backend_->origin();
        }

        //- The rotation specification
        virtual const coordinateRotation& rotation() const
        {
            return backend_->rotation();
        }

        //- Return the name
        virtual const word& name() const
        {
            return backend_->name_;
        }

        //- Return the optional note
        virtual const string& note() const
        {
            return backend_->note();
        }

        //- Return const reference to the rotation tensor
        virtual const tensor& R() const
        {
            return backend_->R();
        }

        //- The local Cartesian x-axis in global coordinates
        virtual const vector e1() const
        {
            return backend_->e1();
        }

        //- The local Cartesian y-axis in global coordinates
        virtual const vector e2() const
        {
            return backend_->e2();
        }

        //- The local Cartesian z-axis in global coordinates
        virtual const vector e3() const
        {
            return backend_->e3();
        }


    // Edit

        //- Rename (ignored)
        void rename(const word& newName) {}

        //- Provide non-constant access to the optional note
        string& note()
        {
            NotImplemented;
            return dummy_.note();
        }

        //- Edit access to origin (disallowed)
        virtual point& origin()
        {
            NotImplemented;
            return dummy_.origin();
        }

        //- Clear (ignored)
        virtual void clear() {}

        //- Change the rotation (disallowed)
        virtual void rotation(autoPtr<coordinateRotation>&& crot)
        {
            NotImplemented;
        }


    // Write

        //- Write
        virtual void write(Ostream& os) const;

        //- Write 'coordinateSystem' dictionary entry
        virtual void writeEntry(Ostream& os) const;

        //- Write dictionary entry
        virtual void writeEntry(const word& keyword, Ostream& os) const;


    // Member Operators

        //- No copy assignment
        void operator=(const coordinateSystem& csys) = delete;

        //- No move assignment
        void operator=(coordinateSystem&& csys) = delete;


    // Rotations

        //- Position-dependent rotation tensor at a single point
        virtual tensor R(const point& global) const
        {
            return backend_->R(global);
        }

        //- Position-dependent rotation tensors at multiple points
        virtual tmp<tensorField> R(const UList<point>& global) const
        {
            return backend_->R(global);
        }

        //- Position-dependent rotation tensors at multiple points
        virtual tmp<tensorField> R(const pointUIndList& global) const
        {
            return backend_->R(global);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace coordSystem
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
