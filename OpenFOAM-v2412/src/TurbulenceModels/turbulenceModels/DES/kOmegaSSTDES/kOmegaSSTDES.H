/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015 OpenFOAM Foundation
    Copyright (C) 2016-2019 OpenCFD Ltd.
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
    Foam::LESModels::kOmegaSSTDES

Group
    grpDESTurbulence

Description
    k-omega-SST DES turbulence model for incompressible and compressible flows.

    Reference:
    \verbatim
        Strelets, M. (2001).
        Detached Eddy Simulation of Massively Separated Flows.
        39th AIAA Aerospace Sciences Meeting and Exhibit, Reno, NV.
    \endverbatim

Note
    The default values of the DES constants implemented are code-specific
    values calibrated for OpenFOAM using decaying isotropic turbulence, and
    hence deviate slightly from the values suggested in the reference
    publication.

SourceFiles
    kOmegaSSTDES.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_kOmegaSSTDES_H
#define Foam_kOmegaSSTDES_H

#include "DESModel.H"
#include "kOmegaSSTBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
                        Class kOmegaSSTDES Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class kOmegaSSTDES
:
    public kOmegaSSTBase<DESModel<BasicTurbulenceModel>>
{
    // Private Member Functions

        //- No copy construct
        kOmegaSSTDES(const kOmegaSSTDES&) = delete;

        //- No copy assignment
        void operator=(const kOmegaSSTDES&) = delete;


protected:

    // Protected Data

        //- Switch to activate grey-area enhanced sigma-(D)DES
        Switch useSigma_;

        // Model coefficients

            dimensionedScalar kappa_;

            //- DES coefficients
            dimensionedScalar CDESkom_;
            dimensionedScalar CDESkeps_;


    // Protected Member Functions

        //- Blending for CDES parameter
        virtual tmp<volScalarField> CDES(const volScalarField& F1) const
        {
            return this->blend(F1, CDESkom_, CDESkeps_);
        }

        virtual void correctNut(const volScalarField& S2);
        virtual void correctNut();

        tmp<volScalarField> r
        (
            const volScalarField& nur,
            const volScalarField& magGradU
        ) const;

        //- Return square of strain rate
        virtual tmp<volScalarField> S2
        (
            const volTensorField& gradU
        ) const;

        //- Return length scale
        virtual tmp<volScalarField> dTilda
        (
            const volScalarField& magGradU,
            const volScalarField& CDES
        ) const;

        //- Return epsilon/k
        virtual tmp<volScalarField::Internal> epsilonByk
        (
            const volScalarField& F1,
            const volTensorField& gradU
        ) const;

        //- Return (G/nu)_0
        virtual tmp<volScalarField::Internal> GbyNu0
        (
            const volTensorField& gradU,
            const volScalarField& S2
        ) const;

        //- Return G/nu
        virtual tmp<volScalarField::Internal> GbyNu
        (
            const volScalarField::Internal& GbyNu0,
            const volScalarField::Internal& F2,
            const volScalarField::Internal& S2
        ) const;


public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;


    //- Runtime type information
    TypeName("kOmegaSSTDES");


    // Constructors

        //- Construct from components
        kOmegaSSTDES
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& propertiesName = turbulenceModel::propertiesName,
            const word& type = typeName
        );


    //- Destructor
    virtual ~kOmegaSSTDES() = default;


    // Member Functions

        //- Re-read model coefficients if they have changed
        virtual bool read();

        //- RAS length scale
        virtual tmp<volScalarField> lengthScaleRAS() const;

        //- LES length scale
        virtual tmp<volScalarField> lengthScaleLES
        (
            const volScalarField& CDES
        ) const;

        //- Return the LES field indicator
        virtual tmp<volScalarField> LESRegion() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#ifdef NoRepository
    #include "kOmegaSSTDES.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#endif

// ************************************************************************* //
