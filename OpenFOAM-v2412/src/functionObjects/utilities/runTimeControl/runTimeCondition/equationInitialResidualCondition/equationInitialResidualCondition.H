/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015 OpenFOAM Foundation
    Copyright (C) 2016-2022 OpenCFD Ltd.
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
    Foam::functionObjects::runTimeControls::equationInitialResidualCondition

Description
    Minimum or maximum initial residual run time condition

SourceFiles
    equationInitialResidualCondition.H
    equationInitialResidualCondition.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_runTimeControls_equationInitialResidualCondition_H
#define functionObjects_runTimeControls_equationInitialResidualCondition_H

#include "runTimeCondition.H"
#include "solverFieldSelection.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;

namespace functionObjects
{
namespace runTimeControls
{

/*---------------------------------------------------------------------------*\
              Class equationInitialResidualCondition Declaration
\*---------------------------------------------------------------------------*/

class equationInitialResidualCondition
:
    public runTimeCondition
{
public:

    enum operatingMode
    {
        omMin,          //!< Minimum
        omMax           //!< Maximum
    };

    static const Enum<operatingMode> operatingModeNames;


protected:

    // Protected data

        //- Field names
        solverFieldSelection fieldSelection_;

        //- Value to compare
        const scalar value_;

        //- Start checking from time - always skips first iteration
        scalar timeStart_;

        //- Operating mode
        operatingMode mode_;


    // Protected Member Functions

        //- Set the residual (scalar) value
        template<class Type>
        void setResidual
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const word& fieldName,
            const label componenti,
            bool& canSet,
            scalar& residual
        ) const;


public:

    //- Runtime type information
    TypeName("equationInitialResidual");

    //- Constructor
    equationInitialResidualCondition
    (
        const word& name,
        const objectRegistry& obr,
        const dictionary& dict,
        stateFunctionObject& state
    );

    //- Destructor
    virtual ~equationInitialResidualCondition() = default;


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

#ifdef NoRepository
    #include "equationInitialResidualConditionTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
