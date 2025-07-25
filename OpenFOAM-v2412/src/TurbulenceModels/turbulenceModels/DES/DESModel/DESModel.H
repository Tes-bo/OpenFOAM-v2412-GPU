/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2015 OpenFOAM Foundation
    Copyright (C) 2022 Upstream CFD GmbH
    Copyright (C) 2016, 2022 OpenCFD Ltd.
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
    Foam::LESModels::DESModel

Description
    Templated abstract base class for DES models

SourceFiles
    DESModel.C

\*---------------------------------------------------------------------------*/

#ifndef LESModels_DESModel_H
#define LESModels_DESModel_H

#include "DESModelBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
                          Class DESModel Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class DESModel
:
    public DESModelBase,
    public LESeddyViscosity<BasicTurbulenceModel>
{
    // Private Member Functions

        //- No copy construct
        DESModel(const DESModel&) = delete;

        //- No copy assignment
        void operator=(const DESModel&) = delete;


protected:

    // Protected Data

        //- Model-specific transition constant
        dimensionedScalar Ctrans_;


public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;

    // Constructors

        //- Construct from components
        DESModel
        (
            const word& type,
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& propertiesName = turbulenceModel::propertiesName
        );


    //- Destructor
    virtual ~DESModel() = default;


    // Public Member Functions

        //- Re-read model coefficients if they have changed
        virtual bool read();

        //- Return the LES field indicator
        virtual tmp<volScalarField> LESRegion() const = 0;

        //- Return modified strain rate
        static tmp<volScalarField> Ssigma
        (
            const volTensorField& gradU,
            const dimensionedScalar& coeff
        );

        //- Return modified strain rate
        //  Note: uses Ctrans_ coefficient
        virtual tmp<volScalarField> Ssigma(const volTensorField& gradU) const;

        //- Return the shielding function
        virtual tmp<volScalarField> fd() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "DESModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
