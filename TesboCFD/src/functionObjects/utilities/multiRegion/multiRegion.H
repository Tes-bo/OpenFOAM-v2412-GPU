/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 OpenCFD Ltd.
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
    Foam::functionObjects::multiRegion

Group
    grpUtilitiesFunctionObjects

Description
    Wrapper that clones the supplied function object for each region.

    Simplifies the setup of identical post-processing requirements for
    multi-region cases.
    Applies the supplied function to all regions by default.

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    multiRegion
    {
        // Mandatory entries
        type    multiRegion;
        libs    (utilityFunctionObjects);

        function
        {
            // Actual object specification
            type    fieldMinMax;
            libs    (fieldFunctionObjects);
            fields  (<field1> ... <fieldN>);
        }

        // Optional entries
        regions     (region1 region2);
    }
    \endverbatim

    where the entries mean:
    \table
      Property   | Description                          | Type | Reqd | Deflt
      type       | Type name: multiRegion               | word | yes  | -
      libs       | Library name: utilityFunctionObjects | word | yes  | -
      function   | Function object sub-dictionary       | dict | yes  | -
      regions    | List of region names                 | wordList | no  | all
    \endtable

See also
  - Foam::functionObject

SourceFiles
    multiRegion.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_multiRegion_H
#define functionObjects_multiRegion_H

#include "timeFunctionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                        Class multiRegion Declaration
\*---------------------------------------------------------------------------*/

class multiRegion
:
    public functionObjects::timeFunctionObject
{
protected:

    // Protected Data

        //- List of cloned objects (1 per region)
        PtrList<functionObject> functions_;


    // Protected Member Functions

        //- No copy construct
        multiRegion(const multiRegion&) = delete;

        //- No copy assignment
        void operator=(const multiRegion&) = delete;


public:

    //- Runtime type information
    TypeName("multiRegion");


    // Constructors

        //- Construct from Time and dictionary
        multiRegion
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~multiRegion() = default;


    // Member Functions

        //- Read the controls
        virtual bool read(const dictionary& dict);

        //- Execute, currently does nothing
        virtual bool execute();

        //- Write the multiRegion
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
