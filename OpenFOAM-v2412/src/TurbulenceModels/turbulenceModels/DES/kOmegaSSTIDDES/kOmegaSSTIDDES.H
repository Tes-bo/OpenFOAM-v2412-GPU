/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015 OpenFOAM Foundation
    Copyright (C) 2019-2022 OpenCFD Ltd.
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
    Foam::LESModels::kOmegaSSTIDDES

Group
    grpDESTurbulence

Description
    k-omega-SST IDDES turbulence model for
    incompressible and compressible flows.

    Reference:
    \verbatim
        Gritskevich, M. S., Garbaruk, A. V.,
        Schütze, J., & Menter, F. R. (2012).
        Development of DDES and IDDES formulations for
        the k-ω shear stress transport model.
        Flow, turbulence and combustion, 88(3), 431-449.
        DOI:10.1007/s10494-011-9378-4
    \endverbatim

SourceFiles
    kOmegaSSTIDDES.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_kOmegaSSTIDDES_H
#define Foam_kOmegaSSTIDDES_H

#include "kOmegaSSTDES.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
                       class kOmegaSSTIDDES Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class kOmegaSSTIDDES
:
    public kOmegaSSTDES<BasicTurbulenceModel>
{
    // Private Member Functions

        //- Check that the supplied delta is an IDDESDelta
        const IDDESDelta& setDelta() const;

        tmp<volScalarField> alpha() const;
        tmp<volScalarField> ft(const volScalarField& magGradU) const;
        tmp<volScalarField> fl(const volScalarField& magGradU) const;

        //- Delay function
        tmp<volScalarField> fdt(const volScalarField& magGradU) const;

        //- No copy construct
        kOmegaSSTIDDES(const kOmegaSSTIDDES&) = delete;

        //- No copy assignment
        void operator=(const kOmegaSSTIDDES&) = delete;


protected:

    // Protected Data

        // Model coefficients

            dimensionedScalar Cdt1_;
            dimensionedScalar Cdt2_;
            dimensionedScalar Cl_;
            dimensionedScalar Ct_;
            Switch fe_;


        //- IDDES delta
        const IDDESDelta& IDDESDelta_;


    // Protected Member Functions

        //- Return the length scale
        virtual tmp<volScalarField> dTilda
        (
            const volScalarField& magGradU,
            const volScalarField& CDES
        ) const;


public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;


    //- Runtime type information
    TypeName("kOmegaSSTIDDES");


    // Constructors

        //- Construct from components
        kOmegaSSTIDDES
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
    virtual ~kOmegaSSTIDDES() = default;


    // Member Functions

        //- Re-read model coefficients if they have changed
        virtual bool read();

        //- Return the shielding function
        virtual tmp<volScalarField> fd() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "kOmegaSSTIDDES.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
