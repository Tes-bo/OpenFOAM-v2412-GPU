/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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
    Foam::fv::PhaseLimitStabilization

Description
    Stabilization source for phase transport equations of \c Type,
    where \c \<Type\>=Scalar/Vector/SphericalTensor/SymmTensor/Tensor.

    Applies an implicit source to the phase transport equation for the specified
    \c field when the phase volume fraction is below \c residualAlpha.  The
    stabilization rate is provided by the registered
    uniformDimensionedScalarField \c rate, which could be extended to also
    support volScalarField and volScalarField::Internal field types.  The \c
    field is currently stabilized towards zero in the limit of the phase volume
    fraction approaching zero but this could be extended to support a
    specified value or a value or field looked-up from the database.

Usage
    Minimal example by using \c constant/fvOptions:
    \verbatim
    \<Type\>PhaseLimitStabilization1
    {
        // Mandatory entries (unmodifiable)
        type            \<Type\>PhaseLimitStabilization;
        field           \<fieldName\>;
        rate            \<rateName\>;

        // Mandatory entries (runtime modifiable)
        residualAlpha   1.0;

        // Mandatory/Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property  | Description                          | Type | Reqd | Dflt
      type      | Type name: \<Type\>PhaseLimitStabilization <!--
                -->                                    | word | yes  | -
      field     | Name of operand field                | word | yes  | -
      rate      | Name of rate field                   | word | yes  | -
      residualAlpha | Residual alpha value below which <!--
                --> stabilization is applied           | scalar | yes | -
    \endtable

    The inherited entries are elaborated in:
     - \link fvOption.H \endlink

    Example usage:
    \verbatim
    stabilization
    {
        type            symmTensorPhaseLimitStabilization;

        field           sigma.liquid;
        rate            rLambda.liquid;
        residualAlpha   1e-3;
    }
    \endverbatim

SourceFiles
    PhaseLimitStabilization.C

\*---------------------------------------------------------------------------*/

#ifndef PhaseLimitStabilization_H
#define PhaseLimitStabilization_H

#include "fvOption.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                Class PhaseLimitStabilization Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class PhaseLimitStabilization
:
    public fv::option
{
    // Private Data

        //- Name of operand field
        word fieldName_;

        //- Name of rate field
        word rateName_;

        //- Residual alpha value below which stabilization is applied
        scalar residualAlpha_;


public:

    //- Runtime type information
    TypeName("PhaseLimitStabilization");


    // Constructors

        //- Construct from components
        PhaseLimitStabilization
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- No copy construct
        PhaseLimitStabilization(const PhaseLimitStabilization&) = delete;

        //- No copy assignment
        void operator=(const PhaseLimitStabilization&) = delete;


    //- Destructor
    virtual ~PhaseLimitStabilization() = default;


    // Member Functions

        //- Source term to compressible phase equation
        virtual void addSup
        (
            const volScalarField& alpha,
            const volScalarField& rho,
            fvMatrix<Type>& eqn,
            const label fieldi
        );

        //- Read dictionary
        virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "PhaseLimitStabilization.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
