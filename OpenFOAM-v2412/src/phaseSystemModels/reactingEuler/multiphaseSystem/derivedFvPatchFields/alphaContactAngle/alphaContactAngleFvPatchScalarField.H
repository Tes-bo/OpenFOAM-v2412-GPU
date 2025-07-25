/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2018 OpenFOAM Foundation
    Copyright (C) 2021-2022 OpenCFD Ltd.
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
    Foam::alphaContactAngleFvPatchScalarField

Description
    Contact-angle boundary condition for multi-phase interface-capturing
    simulations.  Used in conjunction with \c multiphaseSystem.

Usage
    Example of the boundary condition specification:
    \verbatim
    <patch>
    {
        // Mandatory entries
        type                      alphaContactAngle;
        thetaProperties
        (
            (<phase1> <phase2>) <scalar1> <scalar2> <scalar3> <scalar4>
            (<phase3> <phase2>) <scalar1> <scalar2> <scalar3> <scalar4>
            ...
        );

        // Inherited entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property  | Description                             | Type | Reqd | Deflt
      type      | Type name: alphaContactAngle            | word | yes  | -
      thetaProperties | Contact-angle properties          | dict | yes  | -
      \<scalar1\>  | Equilibrium contact angle            | scalar | yes  |-
      \<scalar2\>  | Dynamic contact angle velocity scale | scalar | yes  |-
      \<scalar3\>  | Limiting advancing contact angle     | scalar | yes  |-
      \<scalar4\>  | Limiting receding contact angle      | scalar | yes  |-
    \endtable

    The inherited entries are elaborated in:
     - \link zeroGradientFvPatchFields.H \endlink

SourceFiles
    alphaContactAngleFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_reactingMultiphaseEuler_alphaContactAngleFvPatchScalarField_H
#define Foam_reactingMultiphaseEuler_alphaContactAngleFvPatchScalarField_H

#include "zeroGradientFvPatchFields.H"
#include "multiphaseSystem.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace reactingMultiphaseEuler
{

/*---------------------------------------------------------------------------*\
                    Class alphaContactAngleFvPatch Declaration
\*---------------------------------------------------------------------------*/

class alphaContactAngleFvPatchScalarField
:
    public zeroGradientFvPatchScalarField
{
public:

    class interfaceThetaProps
    {
        //- Equilibrium contact angle
        scalar theta0_;

        //- Dynamic contact angle velocity scale
        scalar uTheta_;

        //- Limiting advancing contact angle
        scalar thetaA_;

        //- Limiting receding contact angle
        scalar thetaR_;


    public:

        // Constructors

            //- Default construct
            interfaceThetaProps()
            {}

            //- Construct from Istream
            interfaceThetaProps(Istream&);


        // Member Functions

            //- Return the equilibrium contact angle theta0
            scalar theta0(bool matched=true) const
            {
                if (matched) return theta0_;
                else return 180.0 - theta0_;
            }

            //- Return the dynamic contact angle velocity scale
            scalar uTheta() const noexcept
            {
                return uTheta_;
            }

            //- Return the limiting advancing contact angle
            scalar thetaA(bool matched=true) const
            {
                if (matched) return thetaA_;
                else return 180.0 - thetaA_;
            }

            //- Return the limiting receding contact angle
            scalar thetaR(bool matched=true) const
            {
                if (matched) return thetaR_;
                else return 180.0 - thetaR_;
            }


        // IOstream operators

            friend Istream& operator>>(Istream&, interfaceThetaProps&);
            friend Ostream& operator<<(Ostream&, const interfaceThetaProps&);
    };

    typedef HashTable
    <
        interfaceThetaProps,
        phasePairKey,
        phasePairKey::hash
    > thetaPropsTable;


private:

    // Private Data

        //- Interface properties
        thetaPropsTable thetaProps_;


public:

    //- Runtime type information
    TypeName("reactingMultiphaseEuler::alphaContactAngle");


    // Constructors

        //- Construct from patch and internal field
        alphaContactAngleFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        alphaContactAngleFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given alphaContactAngleFvPatchScalarField
        //- onto a new patch
        alphaContactAngleFvPatchScalarField
        (
            const alphaContactAngleFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy setting internal field reference
        alphaContactAngleFvPatchScalarField
        (
            const alphaContactAngleFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Return a clone
        virtual tmp<fvPatchField<scalar>> clone() const
        {
            return fvPatchField<scalar>::Clone(*this);
        }

        //- Clone with an internal field reference
        virtual tmp<fvPatchField<scalar>> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return fvPatchField<scalar>::Clone(*this, iF);
        }


    // Member Functions

        //- Return the contact angle properties
        const thetaPropsTable& thetaProps() const noexcept
        {
            return thetaProps_;
        }

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace reactingMultiphaseEuler
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
