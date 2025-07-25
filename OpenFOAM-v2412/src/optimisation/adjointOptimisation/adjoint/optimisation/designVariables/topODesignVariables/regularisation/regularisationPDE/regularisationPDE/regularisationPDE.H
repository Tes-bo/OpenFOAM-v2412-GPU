
/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 PCOpt/NTUA
    Copyright (C) 2021 FOSS GP
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
    Foam::regularisationPDE

Description
    Base class for selecting the regulatisation PDE


SourceFiles
    regularisationPDE.C

\*---------------------------------------------------------------------------*/

#ifndef regularisationPDE_H
#define regularisationPDE_H

#include "fvMesh.H"
#include "topOZones.H"
#include "fvMatrices.H"
#include "regularisationRadius.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class regularisationPDE Declaration
\*---------------------------------------------------------------------------*/

class regularisationPDE
{

private:

    // Private Member Functions

        //- No copy construct
        regularisationPDE(const regularisationPDE&) = delete;

        //- No copy assignment
        void operator=(const regularisationPDE&) = delete;


protected:

    // Protected data

        const fvMesh& mesh_;

        const dictionary dict_;

        //- Cell zones related to topology optimisation
        const topOZones& zones_;

        //- Whether to apply a fixedValue BC or zeroGradient one to alphaTilda,
        //- when regularisation is performed
        bool growFromWalls_;


    // Protected Member Functions

        //- Set fixed bTilda values
        void setValues
        (
            fvScalarMatrix& bTildaEqn,
            const bool isTopoField,
            const scalar minSetValue = Zero
        );

        //- Gather the cells with constant values from all constrained zones
        void setValues
        (
            const fvMesh& mesh,
            DynamicList<label>& cells,
            DynamicList<scalar>& values,
            bool isTopoField,
            const scalar minSetValue = Zero
        );

        //- Compute smoothing radius, if not directly given
        scalar computeRadius();


public:

    //- Runtime type information
    TypeName("regularisationPDE");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            regularisationPDE,
            dictionary,
            (
                const fvMesh& mesh,
                const dictionary& dict,
                const topOZones& zones
            ),
            (mesh, dict, zones)
        );


    // Constructors

        //- Construct from components
        regularisationPDE
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const topOZones& zones
        );


    // Selectors

        //- Construct and return the selected regularisationPDE
        static autoPtr<regularisationPDE> New
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const topOZones& zones
        );


    //- Destructor
    virtual ~regularisationPDE() = default;


    // Member Functions

        //- Regularise field (or sensitivities)
        virtual void regularise
        (
            const volScalarField& aTilda,
            const scalarField& source,
            scalarField& result,
            const bool isTopoField,
            const regularisationRadius& radius,
            const scalar minSetValue = Zero,
            const bool fixATildaValues = true
        ) = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
