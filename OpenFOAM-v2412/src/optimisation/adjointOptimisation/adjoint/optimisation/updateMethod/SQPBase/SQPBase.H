/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2021 PCOpt/NTUA
    Copyright (C) 2013-2021 FOSS GP
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
    Foam::SQPBase

Description
    Base class for Sequantial Quadratic Programming (SQP) methods

SourceFiles
    SQPBase.C

\*---------------------------------------------------------------------------*/

#ifndef SQPBase_H
#define SQPBase_H

#include "constrainedOptimisationMethod.H"
#include "quasiNewton.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class SQPBase Declaration
\*---------------------------------------------------------------------------*/

class SQPBase
:
    public constrainedOptimisationMethod
{
protected:

    // Protected data

        //- Derivatives of the Lagrangian function
        scalarField LagrangianDerivatives_;

        //- The previous constraint derivatives
        List<scalarField> constraintDerivativesOld_;

        //- Lagrange multipliers
        scalarField lamdas_;

        //- Name of the objective folder
        fileName objFunctionFolder_;

        //- File including the l1 merit function
        autoPtr<OFstream> meritFunctionFile_;

        //- Penalty value for the merit function
        scalar mu_;

        //- Safety factor
        scalar delta_;


    // Protected Member Functions

        //- Get the part the merit function that depends on the constraints
        virtual scalar meritFunctionConstraintPart() const = 0;


private:

    // Private Member Functions

        //- No copy construct
        SQPBase(const SQPBase&) = delete;

        //- No copy assignment
        void operator=(const SQPBase&) = delete;

        //- Make folder holding the Lagrangian file
        void makeFolder();


public:

    //- Runtime type information
    TypeName("SQPBase");


    // Constructors

        //- Construct from components
        SQPBase
        (
            const fvMesh& mesh,
            const dictionary& dict,
            autoPtr<designVariables>& designVars,
            const updateMethod& UpdateMethod,
            const word& type
        );


    //- Destructor
    virtual ~SQPBase() = default;


    // Member Functions

        //- Write continuation info
        virtual bool addToFile(Ostream& os) const;

        //- Write info about the merit function
        virtual bool writeMeritFunction(const updateMethod& UpdateMethod);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
