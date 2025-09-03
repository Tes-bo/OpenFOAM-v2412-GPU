/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2015-2022 OpenCFD Ltd.
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
    Foam::FieldMapper

Description
    Abstract base class to hold the Field mapping addressing and weights.

\*---------------------------------------------------------------------------*/

#ifndef Foam_FieldMapper_H
#define Foam_FieldMapper_H

#include "mapDistributeBase.H"
#include "nullObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class FieldMapper Declaration
\*---------------------------------------------------------------------------*/

class FieldMapper
{
public:

    // Constructors

        //- Default construct
        FieldMapper() = default;


    //- Destructor
    virtual ~FieldMapper() = default;


    // Member Functions

        //- The size of the mapper
        virtual label size() const = 0;

        //- Is it a direct (non-interpolating) mapper?
        virtual bool direct() const = 0;

        //- Does the mapper have remote contributions?
        virtual bool distributed() const
        {
            return false;
        }

        //- Any unmapped values?
        //  I.e. do all size() elements get value
        virtual bool hasUnmapped() const = 0;

        //- Return the direct addressing values
        virtual const labelUList& directAddressing() const
        {
            FatalErrorInFunction
                << "attempt to access null direct addressing"
                << abort(FatalError);

            return labelUList::null();
        }

        //- Return the distribution map
        virtual const mapDistributeBase& distributeMap() const
        {
            FatalErrorInFunction
                << "attempt to access null distributeMap"
                << abort(FatalError);

            return NullObjectRef<mapDistributeBase>();
        }

        //- Return the interpolation addressing
        virtual const labelListList& addressing() const
        {
            FatalErrorInFunction
                << "attempt to access null interpolation addressing"
                << abort(FatalError);

            return labelListList::null();
        }

        //- Return the interpolation weights
        virtual const scalarListList& weights() const
        {
            FatalErrorInFunction
                << "attempt to access null interpolation weights"
                << abort(FatalError);

            return scalarListList::null();
        }


    // Member Operators

        //- Perform mapping on the given field
        template<class Type>
        tmp<Field<Type>> operator()(const Field<Type>& fld) const
        {
            return tmp<Field<Type>>::New(fld, *this);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
