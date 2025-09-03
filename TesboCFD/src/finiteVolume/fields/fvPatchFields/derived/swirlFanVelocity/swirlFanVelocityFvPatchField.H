/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2022 OpenCFD Ltd.
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
    Foam::swirlFanVelocityFvPatchField

Group
    grpCoupledBoundaryConditions

Description
    This boundary condition provides a jump condition for \c U across a
    cyclic pressure jump condition and applies a transformation to \c U.

    The U-jump is specified with a swirl component as follows:
    \verbatim
         Utan = deltaP/rEff/fanEff/(rpm*pi/30.0);

         where

            deltaP : pressure drop across the cyclic.
            rEff   : effective radius
            fanEff : fan efficiency coefficient
            rpm    : RPM of the fan
    \endverbatim

    Alternatively an inner and outer radii can be used instead of \c rEff.
    The \c Utan is as follow for \c r > \c rInner and \c r < \c rOuter
    \verbatim
            Utan = deltaP/r/fanEff/(rpm/pi/30.0);

            where

                r : p - origin, p is the face center
    \endverbatim

    Outside \c rInner and \c rOuter, \c Utan=0. The input for this mode is:
    \verbatim
            useRealRadius   true;
            rInner          0.005;
            rOuter          0.01;
    \endverbatim

    The radial velocity is zero in the present model.

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        // Mandatory entries
        type            swirlFanVelocity;
        patchType       cyclic;
        rpm             <Function1>;

        // Optional entries
        phi             <word>;
        p               <word>;
        rho             <word>;
        origin          <vector>;
        fanEff          <scalar>;
        rEff            <scalar>;
        rInner          <scalar>;
        rOuter          <scalar>;
        useRealRadius   <bool>;

        // Inherited entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                       | Type | Reqd | Deflt
      type         | Type name: swirlFanVelocity       | word | yes  | -
      patchType    | Underlying patch type: cyclic     | word | yes  | -
      rpm          | RPM of the fan      | Function1\<scalar\> | yes | -
      phi          | Name of flux field                | word | no   | phi
      rho          | Name of density field             | word | no   | rho
      p            | Name of pressure field            | word | no   | p
      origin       | Fan centre               | vector | no   | calculated
      fanEff       | Fan efficiency                  | scalar | no   | 1
      rEff         | Effective radius                | scalar | no   | 0
      rInner       | Inner radius                    | scalar | no   | 0
      rOuter       | Outer radius                    | scalar | no   | 0
      useRealRadius| Flag to use inner/outer radii     | bool | no   | false
    \endtable

    The inherited entries are elaborated in:
      - \link fixedJumpFvPatchField.H \endlink
      - \link Function1.H \endlink

Note
  - Negative \c rpm will reverse the input tangential direction.
  - This boundary condition needs to be used with a pressure-jump (e.g. fan)
    condition with a non-zero dp, otherwise no swirl will be applied (dp=0).
  - Please ensure physical and complementary set-ups for the pressure-jump
    and \c swirlFanVelocity boundary conditions.

SourceFiles
    swirlFanVelocityFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef swirlFanVelocityFvPatchField_H
#define swirlFanVelocityFvPatchField_H

#include "fixedJumpFvPatchField.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class swirlFanVelocityFvPatchField Declaration
\*---------------------------------------------------------------------------*/

class swirlFanVelocityFvPatchField
:
     public fixedJumpFvPatchField<vector>
{
    // Private data

        //- Name of the flux field
        const word phiName_;

        //- Name of the pressure field
        const word pName_;

        //- Name of the rho field
        const word rhoName_;

        //- Origin of the rotation
        const vector origin_;

        //- Fan rpm
        autoPtr<Function1<scalar>> rpm_;

        //- Fan efficiency
        scalar fanEff_;

        //- Effective fan radius
        scalar rEff_;

        //- Inner radius
        scalar rInner_;

        //- Outer radius
        scalar rOuter_;

        //- Switch to use effective radius or inner and outer radius
        bool useRealRadius_;


    // Private Member Functions

        //- Calculate the fan pressure jump
        void calcFanJump();


public:

    //- Runtime type information
    TypeName("swirlFanVelocity");


    // Constructors

        //- Construct from patch and internal field
        swirlFanVelocityFvPatchField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        swirlFanVelocityFvPatchField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given swirlFanVelocityFvPatchField
        //- onto a new patch
        swirlFanVelocityFvPatchField
        (
            const swirlFanVelocityFvPatchField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        swirlFanVelocityFvPatchField
        (
            const swirlFanVelocityFvPatchField&
        );

        //- Construct as copy setting internal field reference
        swirlFanVelocityFvPatchField
        (
            const swirlFanVelocityFvPatchField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Return a clone
        virtual tmp<fvPatchField<vector>> clone() const
        {
            return fvPatchField<vector>::Clone(*this);
        }

        //- Clone with an internal field reference
        virtual tmp<fvPatchField<vector>> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return fvPatchField<vector>::Clone(*this, iF);
        }


    // Member functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //
