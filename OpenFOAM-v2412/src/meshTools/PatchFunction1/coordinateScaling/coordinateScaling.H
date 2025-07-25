/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2022 OpenCFD Ltd.
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
    Foam::coordinateScaling

Description
    Helper class to wrap coordinate system and component-wise scaling

See also
    Foam::Function1Types

SourceFiles
    coordinateScaling.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_coordinateScaling_H
#define Foam_coordinateScaling_H

#include "coordinateSystem.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class coordinateScaling Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class coordinateScaling
{
    // Private Data

        //- Local coordinate system
        autoPtr<coordinateSystem> coordSys_;

        //- In local coordinate system component-wise scaling
        PtrList<Function1<Type>> scale_;

        //- Cache whether any scaling or coordinate system
        bool active_;


    // Private Member Functions

        //- No copy assignment
        void operator=(const coordinateScaling<Type>&) = delete;

public:

    // Constructors

        //- Default construct
        coordinateScaling();

        //- Construct from registry and dictionary
        coordinateScaling
        (
            const objectRegistry& obr,
            const dictionary& dict
        );

        //- Construct copy
        coordinateScaling(const coordinateScaling& rhs);


    //- Destructor
    virtual ~coordinateScaling() = default;


    // Member Functions

        //- Has any scaling or coordinate transformation
        bool active() const noexcept
        {
            return active_;
        }

        //- Access to optional coordinate system
        const autoPtr<coordinateSystem>& coordSys() const
        {
            return coordSys_;
        }

        //- Convert to local coordinates,
        //- pass-through if no coordinate system is active
        tmp<pointField> localPosition(const pointField& globalPos) const;

        //- Evaluate
        virtual tmp<Field<Type>> transform
        (
            const pointField& pos,
            const Field<Type>& local
        ) const;

        //- Write dictionary entry
        virtual void writeEntry(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Template specializations
template<>
tmp<Field<label>> coordinateScaling<label>::transform
(
    const pointField& pos,
    const Field<label>& local
) const;

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "coordinateScaling.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
