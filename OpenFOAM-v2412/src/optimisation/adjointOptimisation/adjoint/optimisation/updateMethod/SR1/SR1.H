/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2021 PCOpt/NTUA
    Copyright (C) 2013-2021 FOSS GP
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
    Foam::SR1

Description
    The quasi-Newton Symmetric Rank One formula.
    The quasiNewton::Hessian corresponds to Hessian inverse for SR1

SourceFiles
    SR1.C

\*---------------------------------------------------------------------------*/

#ifndef SR1_H
#define SR1_H

#include "quasiNewton.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                             Class SR1 Declaration
\*---------------------------------------------------------------------------*/

class SR1
:
    public quasiNewton
{
protected:

    // Protected data

        //- For stability check
        scalar ratioThreshold_;


    // Protected Member Functions

        //- Update approximation of the inverse Hessian
        virtual void updateHessian();

        //- Update design variables
        virtual void update();


private:

    // Private Member Functions

        //- No copy construct
        SR1(const SR1&) = delete;

        //- No copy assignment
        void operator=(const SR1&) = delete;


public:

    //- Runtime type information
    TypeName("SR1");


    // Constructors

        //- Construct from components
        SR1
        (
            const fvMesh& mesh,
            const dictionary& dict,
            autoPtr<designVariables>& designVars,
            const label nConstraints,
            const word& type
        );


    //- Destructor
    virtual ~SR1() = default;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
