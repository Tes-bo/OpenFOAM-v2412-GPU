/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2022 OpenCFD Ltd.
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
    Foam::interfaceCompositionModel

Description
    Generic base class for interface models. Mass transfer models are
    interface models between two thermos.
    Abstract class for mass transfer functions

SourceFiles
    interfaceCompositionModel.C

\*---------------------------------------------------------------------------*/

#ifndef interfaceCompositionModel_H
#define interfaceCompositionModel_H

#include "volFields.H"
#include "dictionary.H"
#include "Enum.H"
#include "runTimeSelectionTables.H"
#include "multiphaseInterSystem.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class phasePair;

namespace multiphaseInter
{

// Forward Declarations
class phaseModel;

/*---------------------------------------------------------------------------*\
                  Class interfaceCompositionModel Declaration
\*---------------------------------------------------------------------------*/

class interfaceCompositionModel
{
public:

    // Public Types

        //- Enumeration for variable based mass transfer models
        enum modelVariable
        {
            T,  /* temperature based */
            P,  /* pressure based */
            Y,   /* mass fraction based */
            alpha /* alpha source */
        };

protected:

    // Protected Data

        //- Selection names for the modelVariable
        static const Enum<modelVariable> modelVariableNames_;

        //- Enumeration for the model variable
        modelVariable modelVariable_;

        //- Add volume change in pEq
        bool includeVolChange_;

        //- Phase pair
        const phasePair& pair_;

        //- Names of the transferring specie
        word speciesName_;

        //- Reference to mesh
        const fvMesh& mesh_;


public:

    //- Runtime type information
    TypeName("interfaceCompositionModel");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            interfaceCompositionModel,
            dictionary,
            (
                const dictionary& dict,
                const phasePair& pair
            ),
            (dict, pair)
        );


    // Constructors

        //- Construct from a dictionary and a phase pair
        interfaceCompositionModel
        (
            const dictionary& dict,
            const phasePair& pair
        );


    //- Destructor
    virtual ~interfaceCompositionModel() = default;


    // Selectors

        static autoPtr<interfaceCompositionModel> New
        (
            const dictionary& dict,
            const phasePair& pair
        );


    // Member Functions

        //- Return the transferring species name
        const word transferSpecie() const;

        //- The phase pair
        const phasePair& pair() const;

        //- Return the multiphaseInterSystem this interface belongs to
        const multiphaseInterSystem& fluid() const;

        //- Interface mass fraction
        virtual tmp<volScalarField> Yf
        (
            const word& speciesName,
            const volScalarField& Tf
        ) const = 0;

        //- Mass fraction difference between the interface and the field
        virtual tmp<volScalarField> dY
        (
            const word& speciesName,
            const volScalarField& Tf
        ) const = 0;

        //- Specie mass diffusivity for pure mixture
        virtual tmp<volScalarField> Dfrom
        (
            const word& speciesName
        ) const = 0;

        //- Specie mass diffusivity for specie in a multicomponent
        virtual tmp<volScalarField> Dto
        (
             const word& speciesName
        ) const = 0;

        //- Latent heat (delta Hc)
        virtual tmp<volScalarField> L
        (
            const word& speciesName,
            const volScalarField& Tf
        ) const = 0;

        //- Explicit full mass transfer
        virtual tmp<volScalarField> Kexp
        (
            const volScalarField& field
        ) = 0;

        //- Implicit mass transfer
        virtual tmp<volScalarField> KSp
        (
            label modelVariable,
            const volScalarField& field
        ) = 0;

        //- Explicit mass transfer
        virtual tmp<volScalarField> KSu
        (
            label modelVariable,
            const volScalarField& field
        ) = 0;

        //- Reference value
        virtual const dimensionedScalar& Tactivate() const noexcept = 0;

        //- Add/subtract alpha*div(U) as a source term
        //- for alpha, substituting div(U) = mDot(1/rho1 - 1/rho2)
        virtual bool includeDivU() const noexcept;

        //- Add volume change in pEq
        bool includeVolChange();

        //- Returns the variable on which the model is based
        const word& variable() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam
} // End multiphaseInter

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
