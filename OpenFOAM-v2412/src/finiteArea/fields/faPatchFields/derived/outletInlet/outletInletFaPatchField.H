/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023 OpenCFD Ltd.
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
    Foam::outletInletFaPatchField

Group
    grpInletBoundaryConditions

Description
    This boundary condition provides a generic inflow condition, with
    specified outflow for the case of reverse flow.

Usage
    \table
        Property     | Description             | Required    | Default value
        phi          | Flux field name         | no          | phi
        outletValue  | Outlet value for reverse flow | yes   |
        value        | initial field value     | optional |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            outletInlet;
        phi             phi;
        outletValue     uniform 0;
        value           uniform 0;
    }
    \endverbatim

    The mode of operation is determined by the sign of the flux across the
    patch edges.

Note
    Sign conventions:
    - Positive flux (out of domain): apply the "outletValue" fixed-value
    - Negative flux (into domain): apply zero-gradient condition

See also
    Foam::mixedFaPatchField
    Foam::zeroGradientFaPatchField
    Foam::inletOutletFaPatchField

SourceFiles
    outletInletFaPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_outletInletFaPatchField_H
#define Foam_outletInletFaPatchField_H

#include "mixedFaPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class outletInletFaPatch Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class outletInletFaPatchField
:
    public mixedFaPatchField<Type>
{
protected:

    // Protected Data

        //- Name of flux field (default: phi)
        word phiName_;


public:

    //- Runtime type information
    TypeName("outletInlet");


    // Constructors

        //- Construct from patch and internal field
        outletInletFaPatchField
        (
            const faPatch&,
            const DimensionedField<Type, areaMesh>&
        );

        //- Construct from patch, internal field and dictionary
        outletInletFaPatchField
        (
            const faPatch&,
            const DimensionedField<Type, areaMesh>&,
            const dictionary&
        );

        //- Construct by mapping given outletInletFaPatchField onto a new patch
        outletInletFaPatchField
        (
            const outletInletFaPatchField<Type>&,
            const faPatch&,
            const DimensionedField<Type, areaMesh>&,
            const faPatchFieldMapper&
        );

        //- Construct as copy
        outletInletFaPatchField
        (
            const outletInletFaPatchField<Type>&
        );

        //- Construct as copy setting internal field reference
        outletInletFaPatchField
        (
            const outletInletFaPatchField<Type>&,
            const DimensionedField<Type, areaMesh>&
        );

        //- Return clone
        virtual tmp<faPatchField<Type>> clone() const
        {
            return faPatchField<Type>::Clone(*this);
        }

        //- Clone with an internal field reference
        virtual tmp<faPatchField<Type>> clone
        (
            const DimensionedField<Type, areaMesh>& iF
        ) const
        {
            return faPatchField<Type>::Clone(*this, iF);
        }


    // Member Functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write includes "value" entry
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "outletInletFaPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
