/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 OpenCFD Ltd.
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
    Foam::fv::relaxedSnGrad

Group
    grpFvSnGradSchemes

Description
    Surface gradient scheme with under-/over-relaxed
    full or limited explicit non-orthogonal correction.

Usage
    Minimal example by using \c system/fvSchemes:
    \verbatim
    snGradSchemes
    {
        snGrad(<term>)       relaxed;
    }
    \endverbatim

    and by using \c system/fvSolution:
    \verbatim
    relaxationFactors
    {
        fields
        {
            snGrad(<term>)   <relaxation factor>;
        }
    }
    \endverbatim


SourceFiles
    relaxedSnGrad.C

\*---------------------------------------------------------------------------*/

#ifndef relaxedSnGrad_H
#define relaxedSnGrad_H

#include "correctedSnGrad.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                        Class relaxedSnGrad Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class relaxedSnGrad
:
    public snGradScheme<Type>
{
    // Private Data

        //- Type of correction scheme
        tmp<snGradScheme<Type>> correctedScheme_;


    // Private Member Functions

        //- No copy assignment
        void operator=(const relaxedSnGrad&) = delete;


public:

    //- Runtime type information
    TypeName("relaxed");


    // Constructors

        //- Construct from mesh
        relaxedSnGrad(const fvMesh& mesh)
        :
            snGradScheme<Type>(mesh),
            correctedScheme_(new correctedSnGrad<Type>(this->mesh()))
        {}

        //- Construct from mesh and data stream
        relaxedSnGrad(const fvMesh& mesh, Istream& schemeData)
        :
            snGradScheme<Type>(mesh),
            correctedScheme_(new correctedSnGrad<Type>(this->mesh()))
        {}


    //- Destructor
    virtual ~relaxedSnGrad() = default;


    // Member Functions

        //- Return the interpolation weighting factors for the given field
        virtual tmp<surfaceScalarField> deltaCoeffs
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const
        {
            return this->mesh().nonOrthDeltaCoeffs();
        }

        //- Return true if this scheme uses an explicit correction
        virtual bool corrected() const noexcept
        {
            return true;
        }

        //- Return the explicit correction to the relaxedSnGrad
        //- for the given field using the gradients of the field components
        virtual tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
        correction(const GeometricField<Type, fvPatchField, volMesh>&) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "relaxedSnGrad.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
