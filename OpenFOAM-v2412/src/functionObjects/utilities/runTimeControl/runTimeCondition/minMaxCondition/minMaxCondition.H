/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015 OpenFOAM Foundation
    Copyright (C) 2022 OpenCFD Ltd.
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
    Foam::minMaxCondition

Description
    Minimum/maximum run time conditions.  If the value type is not scalar,
    the magnitude of the value is used in the evaluation.

SourceFiles
    minMaxCondition.H
    minMaxCondition.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_runTimeControls_minMaxCondition_H
#define functionObjects_runTimeControls_minMaxCondition_H

#include "runTimeCondition.H"
#include "Enum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
namespace runTimeControls
{

/*---------------------------------------------------------------------------*\
                       Class minMaxCondition Declaration
\*---------------------------------------------------------------------------*/

class minMaxCondition
:
    public runTimeCondition
{
public:

    // Public enumerations

        // Mode type
        enum modeType
        {
            mdMin,          //!< Minimum
            mdMax           //!< Maximum
        };

        static const Enum<modeType> modeTypeNames_;


protected:

    // Protected data

        //- Name of function object to retrueve data from
        word functionObjectName_;

        //- Mode
        modeType mode_;

        //- Field names
        const wordList fieldNames_;

        //- Value to compare
        const scalar value_;

        //- Helper function to retrieve the value from the state dictionary
        template<class Type>
        void setValue
        (
            const word& valueType,
            const word& fieldName,
            scalar& value
        ) const;


public:

    //- Runtime type information
    TypeName("minMax");

    //- Constructor
    minMaxCondition
    (
        const word& name,
        const objectRegistry& obr,
        const dictionary& dict,
        stateFunctionObject& state
    );

    //- Destructor
    virtual ~minMaxCondition() = default;


    // Public Member Functions

        //- Apply the condition
        virtual bool apply();

        //- Write
        virtual void write();

        //- Reset
        virtual void reset();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
template<>
void minMaxCondition::setValue<Foam::scalar>
(
    const word& valueType,
    const word& fieldName,
    scalar& value
) const;

} // End namespace runTimeControls
} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "minMaxConditionTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
