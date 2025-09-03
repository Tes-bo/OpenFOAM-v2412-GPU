/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::turbulentInletFvPatchField

Group
    grpInletBoundaryConditions

Description
    This boundary condition produces spatiotemporal-variant field by summing
    a set of pseudo-random numbers and a given spatiotemporal-invariant mean
    field. The field can be any type, e.g. scalarField. At a single point and
    time, all components are summed by the same random number, e.g. velocity
    components (u, v, w) are summed by the same random number, p; thus, output
    is (u+p, v+p, w+p).

    The pseudo-random number generator obeys the probability density function
    of the uniform distribution constrained by the range [0:1]. The seed for
    the random number generator is hard-coded; therefore, it will produce the
    same sequence of random numbers at every execution.

    \f[
        x_p = (1 - \alpha) x_p^{n - 1} + \alpha (x_{ref} + c s R |x_{ref}|)
    \f]
    where

    \vartable
        x_p     | patch field
        x_{ref} | spatiotemporal-invariant patch scalar
        n       | time level
        \alpha  | a scalar attempting to build two-temporal-point correlations
                  by heuristically adding a fraction of the new random component
                  to the previous time patch field
        c       | a heuristic automatically calculated correction term
                  to compensate energy level losses due to the alpha scalar
        R       | pseudo-random number [HARD-CODED seed]
        s       | fluctuation scale (proportional to the xRef)
    \endvartable

Usage
    \table
        Property         | Description             | Required    | Default value
        fluctuationScale | RMS fluctuation scale (fraction of mean) | yes |
        referenceField   | reference (mean) field  | yes         |
        alpha | fraction of new random component added to previous  | no | 0.1
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        // Mandatory entries
        type             turbulentInlet;
        fluctuationScale 0.1;               // the term `s` above
        referenceField   uniform 10;        // the term `xRef` above

        // Optional entries
        alpha            0.1;               // the term `alpha` above
    }
    \endverbatim

Note
    This boundary condition should not be used for DES or LES computations as a
    turbulent velocity inflow condition, because the BC will not produce
    turbulence-alike time-series, and will decay almost immediately downstream
    of the inlet boundary although its historical name suggests the opposite.

    Nevertheless, the BC may be still used for other applications, e.g. as a
    uniform-random noise source in aeroacoustics.

SeeAlso
    Foam::fixedValueFvPatchField

SourceFiles
    turbulentInletFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef turbulentInletFvPatchField_H
#define turbulentInletFvPatchField_H

#include "Random.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class turbulentInletFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class turbulentInletFvPatchField
:
    public fixedValueFvPatchField<Type>
{
    // Private Data

        //- Random number generator
        Random ranGen_;

        //- Fluctuation scale
        Type fluctuationScale_;

        //- Reference field
        Field<Type> referenceField_;

        //- Fraction of RMS component to apply to last time step values
        scalar alpha_;

        //- Current time index (used for updating)
        label curTimeIndex_;


public:

    //- Runtime type information
    TypeName("turbulentInlet");


    // Constructors

        //- Construct from patch and internal field
        turbulentInletFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        turbulentInletFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given turbulentInletFvPatchField
        //- onto a new patch
        turbulentInletFvPatchField
        (
            const turbulentInletFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        turbulentInletFvPatchField
        (
            const turbulentInletFvPatchField<Type>&
        );

        //- Construct as copy setting internal field reference
        turbulentInletFvPatchField
        (
            const turbulentInletFvPatchField<Type>&,
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


    // Member Functions

        // Access

            //- Return the fluctuation scale
            const Type& fluctuationScale() const
            {
                return fluctuationScale_;
            }

            //- Return reference to the fluctuation scale to allow adjustment
            Type& fluctuationScale()
            {
                return fluctuationScale_;
            }

            //- Return the reference field
            const Field<Type>& referenceField() const
            {
                return referenceField_;
            }

            //- Return reference to the reference field to allow adjustment
            Field<Type>& referenceField()
            {
                return referenceField_;
            }


        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchField<Type>&,
                const labelList&
            );


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "turbulentInletFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
