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
    ESITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::sensitivityShapeESI

Description
    Class for computing sensitivity derivatives using the Enhanced Surface
    Integral (E-SI) formulation, when a parameterization scheme is inluded
    through the design variables. Sensitivity maps are implemented in class
    sensitivitySurfacePoints

SourceFiles
    sensitivityShapeESI.C

\*---------------------------------------------------------------------------*/

#ifndef sensitivityShapeESI_H
#define sensitivityShapeESI_H

#include "ShapeSensitivitiesBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class sensitivityShapeESI Declaration
\*---------------------------------------------------------------------------*/

class sensitivityShapeESI
:
    public ShapeSensitivitiesBase
{
protected:

    // Protected Member Functions

        //- Compute dxdbMult from its various components
        void computeDxDbMult();


private:

    // Private Member Functions

        //- No copy construct
        sensitivityShapeESI(const sensitivityShapeESI&) = delete;

        //- No copy assignment
        void operator=(const sensitivityShapeESI&) = delete;


public:

    //- Runtime type information
    TypeName("shapeESI");


    // Constructors

        //- Construct from components
        sensitivityShapeESI
        (
            const fvMesh& mesh,
            const dictionary& dict,
            adjointSolver& adjointSolver
        );


    //- Destructor
    virtual ~sensitivityShapeESI() = default;


    // Member Functions

        //- Read dict if changed
        virtual bool readDict(const dictionary& dict);

        //- Assemble sensitivities
        //  Solve the adjoint eikonal PDE and the adjoint grid displacement PDE,
        //  if needed, and assemble the sensitivities
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
