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

Class
    Foam::functionObjects::timeFunctionObject

Description
    Virtual base class for function objects with a reference to Time.

See also
    Foam::functionObject

SourceFiles
    timeFunctionObject.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_timeFunctionObject_H
#define functionObjects_timeFunctionObject_H

#include "functionObject.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                     Class timeFunctionObject Declaration
\*---------------------------------------------------------------------------*/

class timeFunctionObject
:
    public functionObject
{
protected:

    // Protected Member Data

        //- Reference to the time database
        const Time& time_;


    // Protected Member Functions

        //- Remove specified items from "functionObjectObjects"
        void clearOutputObjects(const wordList& objNames);


        //- No copy construct
        timeFunctionObject(const timeFunctionObject&) = delete;

        //- No copy assignment
        void operator=(const timeFunctionObject&) = delete;


public:

    // Constructors

        //- Construct from Time
        timeFunctionObject(const word& name, const Time& runTime);


    //- Destructor
    virtual ~timeFunctionObject() = default;


    // Member Functions

        //- Return time database
        const Time& time() const
        {
            return time_;
        }

        //- Write access to the output objects ("functionObjectObjects")
        //- registered on Time
        objectRegistry& storedObjects();

        //- Const access to the output objects ("functionObjectObjects")
        //- registered on Time
        const objectRegistry& storedObjects() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
