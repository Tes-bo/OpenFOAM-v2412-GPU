/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020-2021 OpenCFD Ltd.
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
    Foam::patchFunction1Base

Description
    Top level data entry class for use in dictionaries. Provides a mechanism
    to specify a variable as a certain type, e.g. constant or time varying, and
    provide functions to return the (interpolated) value, and integral between
    limits.

    Extends the Function1 class by adding autoMap and rMap functions

SourceFiles
    patchFunction1Base.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_patchFunction1Base_H
#define Foam_patchFunction1Base_H

#include "dictionary.H"
#include "polyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class objectRegistry;
class Time;

/*---------------------------------------------------------------------------*\
                     Class patchFunction1Base Declaration
\*---------------------------------------------------------------------------*/

class patchFunction1Base
:
    public refCount
{
protected:

    // Protected Data

        //- Name of entry
        const word name_;

        //- Reference to the patch
        const polyPatch& patch_;

        //- Generate face or point values on patch
        const bool faceValues_;


    // Protected Member Functions

        //- No copy assignment
        void operator=(const patchFunction1Base&) = delete;


public:

    // Constructors

        //- Construct from polyPatch and entry name
        patchFunction1Base
        (
            const polyPatch& pp,
            const word& entryName,
            const bool faceValues = true
        );

        //- Construct from polyPatch, dictionary and entry name
        patchFunction1Base
        (
            const polyPatch& pp,
            const word& entryName,
            const dictionary& dict,
            const bool faceValues = true
        );

        //- Copy construct
        explicit patchFunction1Base(const patchFunction1Base& rhs);

        //- Copy construct setting patch
        explicit patchFunction1Base
        (
            const patchFunction1Base& rhs,
            const polyPatch& pp
        );


    //- Destructor
    virtual ~patchFunction1Base() = default;


    // Member Functions

    // Access

        //- The name of the entry
        const word& name() const noexcept
        {
            return name_;
        }

        //- Reference to the patch
        const polyPatch& patch() const noexcept
        {
            return patch_;
        }

        //- Generate face or point values on patch?
        bool faceValues() const noexcept
        {
            return faceValues_;
        }

        //- Number of faces or points on the patch
        label size() const
        {
            return (faceValues_ ? patch_.size() : patch_.nPoints());
        }


        //- Return the associated registry (ie, the mesh)
        const objectRegistry* whichDb() const;

        //- Return the object registry (ie, the mesh)
        const objectRegistry& obr() const;

        //- Return false: function not created with time database
        /// bool isTime() const noexcept
        /// {
        ///     return false;
        /// }

        //- Return the time database
        const Time& time() const;


    // Manipulation

        //- Convert time
        virtual void userTimeToTime(const Time& t);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
