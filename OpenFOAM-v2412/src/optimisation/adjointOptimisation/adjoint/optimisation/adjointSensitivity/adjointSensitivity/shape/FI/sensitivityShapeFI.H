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
    Foam::sensitivityShapeFI

Description
    Class for computing Field Integral (FI)-based sensitivity derivatives

SourceFiles
    sensitivityShapeFI.C

\*---------------------------------------------------------------------------*/

#ifndef sensitivityShapeFI_H
#define sensitivityShapeFI_H

#include "ShapeSensitivitiesBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


/*---------------------------------------------------------------------------*\
                            Class sensitivityShapeFI Declaration
\*---------------------------------------------------------------------------*/

class sensitivityShapeFI
:
    public ShapeSensitivitiesBase
{
private:

    // Private Member Functions

        //- No copy construct
        sensitivityShapeFI(const sensitivityShapeFI&) = delete;

        //- No copy assignment
        void operator=(const sensitivityShapeFI&) = delete;


public:

    //- Runtime type information
    TypeName("shapeFI");


    // Constructors

        //- Construct from components
        sensitivityShapeFI
        (
            const fvMesh& mesh,
            const dictionary& dict,
            adjointSolver& adjointSolver
        );


    //- Destructor
    virtual ~sensitivityShapeFI() = default;


    // Member Functions

        //- Should the parameterization compute the internalField of dxdb
        virtual bool computeDxDbInternalField() const;

        //- Assemble sensitivities
        //  Solve the adjoint eikonal PDE, if needed, and calls the assembles
        //  the sensitivities
        virtual void assembleSensitivities
        (
            autoPtr<designVariables>& designVars
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
