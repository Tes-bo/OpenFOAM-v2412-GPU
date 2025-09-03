/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2021 OpenCFD Ltd
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
    Foam::fanPressureFvPatchScalarField

Group
    grpInletBoundaryConditions grpOutletBoundaryConditions

Description
    This boundary condition can be applied to assign either a pressure inlet
    or outlet total pressure condition for a fan.

    The switch nonDimensional can be used for a non-dimensional fan curve.
    It needs inputs rpm and dm of the fan.

    The nonDimensional flux is calculate as :

        phi = 4.0*mDot/(rho*sqr(PI)*dm^3*omega)
        where:
            dm is the mean diameter.
            omega is rad/sec.

    The nonDimensinal pressure :

        Psi = 2 deltaP/(rho*(sqr(PI*omega*dm)))
        where:
            deltaP is the pressure drop

    The non-dimensional table should be given as Psi = F(phi).

Usage
    \table
        Property     | Description                       | Required | Default
        fanCurve     | Pressure vs flow-rate             | yes |
        direction    | direction of flow through fan [in/out] | yes |
        p0           | environmental total pressure      | yes |
        nonDimensional | uses non-dimensional table      | no  | false
        rpm          | fan rpm for non-dimensional table | no  |
        dm           | mean diameter for non-dimensional table | no |
        file         | fan curve file name               | legacy |
        outOfBounds  | out of bounds handling            | legacy |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    inlet
    {
        type            fanPressure;
        direction       in;
        fanCurve
        {
            type table;
            file "<constant>/fanCurve";
            outOfBounds clamp;  // Optional out-of-bounds handling
        }
        p0              uniform 0;
        value           uniform 0;
    }

    // Legacy specification
    outlet
    {
        type            fanPressure;
        direction       out;
        file            "<constant>/fanCurve";
        outOfBounds     clamp;
        p0              uniform 0;
        value           uniform 0;
    }
    \endverbatim

Note
    For compatibility with older versions (OpenFOAM-v2006 and earlier),
    a missing \c fanCurve keyword is treated as a tableFile and makes
    the \c file keyword mandatory.

See also
    Foam::fanFvPatchField
    Foam::totalPressureFvPatchScalarField
    Foam::Function1

SourceFiles
    fanPressureFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef fanPressureFvPatchScalarField_H
#define fanPressureFvPatchScalarField_H

#include "totalPressureFvPatchScalarField.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class fanPressureFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class fanPressureFvPatchScalarField
:
    public totalPressureFvPatchScalarField
{
public:

        //- Fan flow direction
        enum fanFlowDirection : uint8_t
        {
            ffdIn,
            ffdOut
        };

        //- Fan flow direction names
        static const Enum<fanFlowDirection> fanFlowDirectionNames_;


private:

    // Private Data

        //- Run-time selectable fan curve
        autoPtr<Function1<scalar>> fanCurve_;

        //- Direction of flow through the fan relative to patch
        fanFlowDirection direction_;

        //- Use non-dimensional curve
        bool nonDimensional_;

        //- Fan rpm (for non-dimensional curve)
        autoPtr<Function1<scalar>> rpm_;

        //- Fan mean diameter (for non-dimensional curve)
        autoPtr<Function1<scalar>> dm_;


public:

    //- Runtime type information
    TypeName("fanPressure");


    // Constructors

        //- Construct from patch and internal field
        fanPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        fanPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  fanPressureFvPatchScalarField
        //  onto a new patch
        fanPressureFvPatchScalarField
        (
            const fanPressureFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        fanPressureFvPatchScalarField
        (
            const fanPressureFvPatchScalarField&
        );

        //- Construct as copy setting internal field reference
        fanPressureFvPatchScalarField
        (
            const fanPressureFvPatchScalarField&,
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


    // Member functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
