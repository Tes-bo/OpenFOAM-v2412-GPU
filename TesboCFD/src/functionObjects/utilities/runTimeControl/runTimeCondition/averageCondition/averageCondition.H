/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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
    Foam::functionObjects::runTimeControls::averageCondition

Description
    Average run time condition - satisfied when average does not change by
    more than a given value.

See also
    - Foam::functionObjects::valueAverageBase

SourceFiles
    averageCondition.H
    averageCondition.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_runTimeControls_averageCondition_H
#define functionObjects_runTimeControls_averageCondition_H

#include "runTimeCondition.H"
#include "valueAverageBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
namespace runTimeControls
{

/*---------------------------------------------------------------------------*\
                      Class averageCondition Declaration
\*---------------------------------------------------------------------------*/

class averageCondition
:
    public runTimeCondition,
    public valueAverageBase
{
protected:

    // Protected Data

        //- Number of start-up iterations before allowing satisfied checks
        label nIterStartUp_;

        //- Current iteration count
        label iter_;


public:

    //- Runtime type information
    TypeName("average");

    //- Constructor
    averageCondition
    (
        const word& name,
        const objectRegistry& obr,
        const dictionary& dict,
        stateFunctionObject& state
    );

    //- Destructor
    virtual ~averageCondition() = default;


    // Public Member Functions

        //- Apply the condition
        virtual bool apply();

        //- Write
        virtual void write();

        //- Reset
        virtual void reset();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace runTimeControls
} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
