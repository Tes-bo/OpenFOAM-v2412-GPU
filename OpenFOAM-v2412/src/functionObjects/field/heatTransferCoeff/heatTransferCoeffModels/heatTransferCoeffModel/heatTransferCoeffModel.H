/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2023 OpenCFD Ltd.
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

Namespace
    Foam::heatTransferCoeffModels

Description
    A namespace for various heat transfer coefficient model implementations.

Class
    Foam::heatTransferCoeffModel

Description
    A base class for heat transfer coefficient models.

SourceFiles
    heatTransferCoeffModel.C
    heatTransferCoeffModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_heatTransferCoeffModel_H
#define Foam_heatTransferCoeffModel_H

#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class fvMesh;

/*---------------------------------------------------------------------------*\
                   Class heatTransferCoeffModel Declaration
\*---------------------------------------------------------------------------*/

class heatTransferCoeffModel
{
protected:

    // Protected Data

        //- Const reference to the mesh
        const fvMesh& mesh_;

        //- List of (wall) patches to process (selected by name)
        labelList patchIDs_;

        //- Name of temperature field
        const word TName_;

        //- Name of radiative heat flux field
        word qrName_;


    // Protected Member Functions

        //- Set the heat transfer coefficient
        virtual void htc
        (
            volScalarField& htc,
            const FieldField<Field, scalar>& q
        ) = 0;


public:

    //- Runtime type information
    TypeName("heatTransferCoeffModel");


    // Declare runtime constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            heatTransferCoeffModel,
            dictionary,
            (
                const dictionary& dict,
                const fvMesh& mesh,
                const word& TName
            ),
            (dict, mesh, TName)
        );


    // Selectors

        //- Return a reference to the selected heat transfer coefficient model
        static autoPtr<heatTransferCoeffModel> New
        (
            const dictionary& dict,
            const fvMesh& mesh,
            const word& TName
        );


    // Constructors

        //- Construct from components
        heatTransferCoeffModel
        (
            const dictionary& dict,
            const fvMesh& mesh,
            const word& TName
        );

        //- No copy construct
        heatTransferCoeffModel(const heatTransferCoeffModel&) = delete;

        //- No copy assignment
        void operator=(const heatTransferCoeffModel&) = delete;


    //- Destructor
    virtual ~heatTransferCoeffModel() = default;


    // Member Functions

    // Access

        //- Return const reference to the mesh
        const fvMesh& mesh() const noexcept
        {
            return mesh_;
        }

        //- Return const reference to wall patches to process
        const labelList& patchIDs() const noexcept
        {
            return patchIDs_;
        }

        //- Return const reference to name of temperature field
        const word& TName() const noexcept
        {
            return TName_;
        }

        //- Return const reference to name of radiative heat-flux field
        const word& qrName() const noexcept
        {
            return qrName_;
        }


    // Evaluation

        //- Return boundary fields of heat-flux field
        tmp<FieldField<Field, scalar>> q() const;

        //- Calculate the heat transfer coefficient field and return true
        //- if successful
        virtual bool calc
        (
            volScalarField& result,
            const FieldField<Field, scalar>& q
        );


    // I-O

        //- Read from dictionary
        virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
