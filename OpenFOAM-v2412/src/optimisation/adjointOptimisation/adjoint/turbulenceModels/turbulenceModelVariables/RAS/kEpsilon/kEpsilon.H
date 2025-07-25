/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2023 PCOpt/NTUA
    Copyright (C) 2013-2023 FOSS GP
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
    Foam::incompressible::RASVariables::kEpsilon

Description

SourceFiles
    kEpsilon.C

\*---------------------------------------------------------------------------*/

#ifndef incompressible_RASVariables_kEpsilon_H
#define incompressible_RASVariables_kEpsilon_H

#include "RASModelVariables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASVariables
{

/*---------------------------------------------------------------------------*\
                       Class kEpsilon Declaration
\*---------------------------------------------------------------------------*/

class kEpsilon
:
    public RASModelVariables
{

protected:

    // Protected data

        //- Average of the production term
        //  Avaraged separetely due the bi-zonal treatment next to the wall
        autoPtr<volScalarField::Internal> GMean_;


    // Protected Member functions

        virtual void allocateMeanFields();
        tmp<volScalarField::Internal> computeG();


public:

    //- Runtime type information
    TypeName("kEpsilon");


    // Constructors

        //- Construct from components
        kEpsilon
        (
            const fvMesh& mesh,
            const solverControl& SolverControl
        );


    //- Destructor
    virtual ~kEpsilon() = default;


    // Member Functions

        //- Return the turbulence production term
        virtual tmp<volScalarField::Internal> G();

        //- Compute mean fields on the fly
        virtual void computeMeanFields();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASVariables
} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
