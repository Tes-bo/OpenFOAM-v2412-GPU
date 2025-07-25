/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017 IH-Cantabria
    Copyright (C) 2017 OpenCFD Ltd.
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
    Foam::fv::multiphaseMangrovesTurbulenceModel

Group
    grpFvOptionsSources

Description

Usage

SourceFiles
    multiphaseMangrovesTurbulenceModel.C

\*---------------------------------------------------------------------------*/

#ifndef multiphaseMangrovesTurbulenceModel_H
#define multiphaseMangrovesTurbulenceModel_H

#include "fvOption.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fv
{

/*---------------------------------------------------------------------------*\
              Class multiphaseMangrovesTurbulenceModel Declaration
\*---------------------------------------------------------------------------*/

class multiphaseMangrovesTurbulenceModel
:
    public fv::option
{
    // Private Member Functions

        //- No copy construct
        multiphaseMangrovesTurbulenceModel
        (
            const multiphaseMangrovesTurbulenceModel&
        ) = delete;

        //- No copy assignment
        void operator=(const multiphaseMangrovesTurbulenceModel&) = delete;



protected:

    // Protected data

        // Coefficients per cell zone

            //- Width of the vegetation element
            scalarList aZone_;

            //- Number of plants per unit of area
            scalarList NZone_;

            //- Ckp
            scalarList CkpZone_;

            //- Cep
            scalarList CepZone_;

            //- Drag coefficient
            scalarList CdZone_;

            //- Zone indices
            labelListList zoneIDs_;


        // Field properties

            //- Name of U; default = U
            word UName_;

            //- Name of k; default = k
            word kName_;

            //- Name of epsilon; default = epsilon
            word epsilonName_;


    // Protected Member Functions

        //- Return the k coefficient
        tmp<volScalarField> kCoeff(const volVectorField& U) const;

        //- Return the epsilon coefficient
        tmp<volScalarField> epsilonCoeff(const volVectorField& U) const;


public:

    //- Runtime type information
    TypeName("multiphaseMangrovesTurbulenceModel");


    // Constructors

        //- Construct from components
        multiphaseMangrovesTurbulenceModel
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~multiphaseMangrovesTurbulenceModel() = default;


    // Member Functions

        // Add explicit and implicit contributions

            //- Add implicit contribution to momentum equation
            virtual void addSup
            (
                fvMatrix<scalar>& eqn,
                const label fieldi
            );

            //- Add implicit contribution to compressible momentum equation
            virtual void addSup
            (
                const volScalarField& rho,
                fvMatrix<scalar>& eqn,
                const label fieldi
            );


        // IO

            //- Read dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
