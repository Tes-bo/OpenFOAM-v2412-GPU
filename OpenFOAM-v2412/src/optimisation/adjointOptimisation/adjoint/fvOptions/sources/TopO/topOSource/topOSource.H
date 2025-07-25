/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020-2023 PCOpt/NTUA
    Copyright (C) 2020-2023 FOSS GP
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
    Foam::fv::topOSource

Group
    grpFvOptionsSources

Description
    Implements Brinkman penalisation terms for topology optimisation.
    Looks up the indicator field (beta) from the registry, through
    topOVariablesBase

SourceFiles
    topOSource.C

\*---------------------------------------------------------------------------*/

#ifndef topOSource_H
#define topOSource_H

#include "cellSetOption.H"
#include "topOInterpolationFunction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fv
{


/*---------------------------------------------------------------------------*\
                          Class topOSource Declaration
\*---------------------------------------------------------------------------*/

class topOSource
:
    public option
{

protected:

    // Protected data

        //- Interpolation function
        autoPtr<topOInterpolationFunction> interpolation_;

        //- Interpolation field name
        word interpolationFieldName_;

        //- Optional betaMax
        //  If not found, the one known by topOVariablesBase will be used.
        scalar betaMax_;

        //- Does this option apply to a Darcy flow model
        bool darcyFlow_;

        //- Dimensionless Darcy number
        autoPtr<scalar> Da_;


    // Protected Member Functions

        //- Compute the source term based on the indicator field
        virtual tmp<DimensionedField<scalar, volMesh>> getSource();


private:

    // Private Member Functions

        //- No copy construct
        topOSource(const topOSource&) = delete;

        //- No copy assignment
        void operator=(const topOSource&) = delete;


public:

    //- Runtime type information
    TypeName("topOSource");


    // Constructors

        //- Construct from components
        topOSource
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~topOSource() = default;


    // Member Functions

        //- Add implicit contribution to momentum equation
        virtual void addSup
        (
            fvMatrix<vector>& eqn,
            const label fieldi
        );

        //- Add implicit contribution to scalar equations
        //- (e.g. turbulence model)
        virtual void addSup
        (
            fvMatrix<scalar>& eqn,
            const label fieldi
        );

        //- Add implicit contribution to compressible momentum equation
        virtual void addSup
        (
            const volScalarField& rho,
            fvMatrix<vector>& eqn,
            const label fieldi
        );

        //- Add implicit contribution to compressible scalar equation
        virtual void addSup
        (
            const volScalarField& rho,
            fvMatrix<scalar>& eqn,
            const label fieldi
        );

        //- Multiply sensitivities with the derivative of the interpolation
        //- function
        virtual void postProcessSens
        (
            scalarField& sensField,
            const word& fieldName = word::null,
            const word& designValue = word::null
        );

        //- Read source dictionary
        virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
