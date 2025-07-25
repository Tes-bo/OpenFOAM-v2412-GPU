/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2024 OpenCFD Ltd.
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
    Foam::fanFvPatchField

Group
    grpCoupledBoundaryConditions

Description
    This boundary condition provides a jump condition, using the \c cyclic
    condition as a base.

    The jump is specified as a \c Function1 type, to enable the use of, e.g.
    constant, polynomial, table values.

    The basis of the table is specified according to the \c mode:

    - velocity: deltap = F(velocity per face) \[DEFAULT\]
    - uniformVelocity: deltap = F(patch area-averaged velocity)
    - volumeFlowRate:  deltap = F(patch volume flow rate)
    - nonDimensional:  non-dim deltap = F(non-dim volume flow rate)

    Non-dimensional operation:

    As inputs it needs the fan RPM (rpm) and the mean diameter (dm).

    The non-dimensional U for the table is calculated as follows:

    \verbatim
        phi = 120*Un/(PI^3*dm^3*rpm)
        where:
            dm is the mean diameter.
            rpm is the RPM of the fan.
    \endverbatim

    The non-dimensional pressure:

    \verbatim
        Psi = 2 deltaP/(rho*(sqr(PI*omega*dm)))
        where:
            deltaP is the pressure drop
    \endverbatim

    The non-dimensional table should be given as Psi = F(phi).

Usage
    \table
        Property    | Description                           | Required | Default
        patchType   | underlying patch type should be \c cyclic | yes |
        mode        | jump table operating mode (see above) | no | velocity
        jumpTable   | jump data, e.g. \c csvFile            | yes |
        phi         | flux field name                       | no  | phi
        rho         | density field name                    | no  | rho
        rpm         | fan rpm (non-dimensional table)       | no |
        dm          | mean diameter (non-dimensional table) | no |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            fan;
        patchType       cyclic;
        jumpTable       csvFile;
        mode            velocity;

        jumpTableCoeffs
        {
            nHeaderLine     1;
            refColumn       0;
            componentColumns 1(1);
            separator       ",";
            mergeSeparators no;
            file            "<constant>/UvsPressure";
        }
        value           uniform 0;
    }
    \endverbatim

    The above example shows the use of a comma separated (CSV) file to specify
    the jump condition.

Note
    The underlying \c patchType should be set to \c cyclic

See also
    Foam::Function1Types

SourceFiles
    fanFvPatchField.C
    fanFvPatchFields.H
    fanFvPatchFields.C

\*---------------------------------------------------------------------------*/

#ifndef foam_fanFvPatchField_H
#define foam_fanFvPatchField_H

#include "uniformJumpFvPatchField.H"
#include "Function1.H"
#include "Enum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class fanFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class fanFvPatchField
:
    public uniformJumpFvPatchField<Type>
{
public:

    // Public Data Types

        //- Enumeration defining the operating modes
        enum class operatingMode
        {
            VELOCITY,           //!< velocity-based lookup
            UNIFORM_VELOCITY,   //!< uniform velocity-based lookup
            VOL_FLOW_RATE,      //!< volume-flow-rate-based lookup
            NON_DIMENSIONAL     //!< non-dimensional-based lookup
        };

        //- Names for the operating modes
        static const Enum<operatingMode> operatingModeNames_;


private:

    // Private Data

        //- Operating mode
        operatingMode operatingMode_;

        //- Name of the flux transporting the field
        word phiName_;

        //- Name of the density field for normalising the mass flux if necessary
        word rhoName_;

        //- Fan rpm (for non-dimensional curve)
        autoPtr<Function1<scalar>> rpm_;

        //- Fan mean diameter (for non-dimensional curve)
        autoPtr<Function1<scalar>> dm_;


    // Private Member Functions

        //- Calculate the fan pressure jump
        void calcFanJump();


public:

    //- Runtime type information
    TypeName("fan");


    // Constructors

        //- Construct from patch and internal field
        fanFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        fanFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given fanFvPatchField onto a new patch
        fanFvPatchField
        (
            const fanFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        fanFvPatchField
        (
            const fanFvPatchField<Type>&
        );

        //- Construct as copy setting internal field reference
        fanFvPatchField
        (
            const fanFvPatchField<Type>&,
            const DimensionedField<Type, volMesh>&
        );

        //- Return a clone
        virtual tmp<fvPatchField<Type>> clone() const
        {
            return fvPatchField<Type>::Clone(*this);
        }

        //- Clone with an internal field reference
        virtual tmp<fvPatchField<Type>> clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return fvPatchField<Type>::Clone(*this, iF);
        }


    // Member functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream& os) const;
};


//- Specialisation of the jump-condition for the pressure
template<>
void fanFvPatchField<scalar>::calcFanJump();


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fanFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
