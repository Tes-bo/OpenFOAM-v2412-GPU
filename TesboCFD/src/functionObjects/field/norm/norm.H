/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 OpenCFD Ltd.
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
    Foam::functionObjects::norm

Group
    grpFieldFunctionObjects

Description
    Normalises an input field with a chosen norm,
    and outputs a new normalised field.

    Operands:
    \table
      Operand      | Type                       | Location
      input        | {vol,surface,polySurface}\<Type\>Field <!--
                   -->                          | \<case\>/\<time\>/\<inpField\>
      output file  | -                          | -
      output field | {vol,surface,polySurface}\<Type\>Field <!--
                   -->                          | \<case\>/\<time\>/\<outField\>
    \endtable

    where \c \<Type\>=Scalar/Vector/SphericalTensor/SymmTensor/Tensor.

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    norm1
    {
        // Mandatory entries
        type            norm;
        libs            (fieldFunctionObjects);
        field           <word>;
        norm            <word>;

        // Conditional entries

            // when norm == Lp
            p               <scalar>;

            // when norm == composite
            divisor         <Function1<scalar>>;

            // when norm == divisorField
            divisorField    <word>;

        // Inherited entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property | Description                        | Type | Reqd | Deflt
      type     | Type name: norm                    | word | yes  | -
      libs     | Library name: fieldFunctionObjects | word | yes  | -
      field    | Name of the operand field          | word | yes  | -
      norm     | Name of normalisation operation    | word | yes  | -
      p        | Norm exponent for the p-norm       | scalar | cndtnl | -
      divisor  | Norm divisor for the composite norm        <!--
               -->                     | Function1\<scalar\> | cndtnl | -
      divisorField | Divisor scalar field for the field norm <!--
               -->                                  | word | cndtnl | -
    \endtable

    Options for the \c norm entry:
    \verbatim
      L1            | L1/Taxicab norm
      L2            | L2/Euclidean norm
      Lp            | p norm
      max           | Maximum norm
      composite     | Composite norm comprising Function1 divisor
      divisorField  | Normalise by a given field
    \endverbatim

    The inherited entries are elaborated in:
      - \link functionObject.H \endlink
      - \link fieldExpression.H \endlink
      - \link Function1.H \endlink


    The normalisations are held according to the following expressions:

    - \c L1:

    \f[
        \mathbf{y} = \frac{\mathbf{x}}{\Sigma_{i=1}^n |x_i|}
    \f]

    - \c L2:

    \f[
        \mathbf{y} = \frac{\mathbf{x}}{\sqrt{x_1^2 + ... + x_n^2}}
    \f]

    - \c Lp:

    \f[
        \mathbf{y} = \frac{\mathbf{x}}{(\Sigma_{i=1}^n |x_i|^p)^{1/p}}
    \f]

    - \c max:

    \f[
        \mathbf{y} = \frac{\mathbf{x}}{max|x_i|}
    \f]

    - \c composite:

    \f[
        \mathbf{y} = \frac{\mathbf{x}}{f(t)}
    \f]

    - \c divisorField:

    \f[
        \mathbf{y} = \frac{\mathbf{x}}{\mathbf{z}}
    \f]

Note
  - Divisor quantity is perturbed by \c SMALL value to prevent any divisions
    by zero irrespective of whether the quantity is close to zero or not.

SourceFiles
    norm.C
    normImpl.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_norm_H
#define functionObjects_norm_H

#include "fieldExpression.H"
#include "Function1.H"
#include "polySurfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                            Class norm Declaration
\*---------------------------------------------------------------------------*/

class norm
:
    public fieldExpression
{
    // Private Enumerations

        //- Options for the norm algorithm
        enum normType : char
        {
            L1,             //!< "L1/Taxicab norm"
            L2,             //!< "L2/Euclidean norm"
            LP,             //!< "p norm"
            MAX,            //!< "Maximum norm"
            COMPOSITE,      //!< "Composite norm"
            FIELD           //!< "Normalise by a given field"
        };

        //- Names for normType
        static const Enum<normType> normTypeNames;


    // Private Data

        //- Norm algorithm
        enum normType norm_;

        //- Norm divisor for the composite norm
        autoPtr<Function1<scalar>> divisorPtr_;

        //- Divisor scalar field for the field norm
        word divisorFieldName_;

        //- Norm exponent for the p norm
        scalar p_;


    // Private Member Functions

        //- Calculate the chosen norm of the field and register the result
        template<class Type>
        bool calcNorm();

        //- Return the chosen norm of the field
        template<class GeoFieldType>
        tmp<GeoFieldType> calcNormType();

        //- Return the divisor volScalarField
        template<class Type>
        tmp<volScalarField> fieldNorm
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        //- Return the divisor surfaceScalarField
        template<class Type>
        tmp<surfaceScalarField> fieldNorm
        (
            const GeometricField<Type, fvsPatchField, surfaceMesh>&
        );

        //- Return the divisor polySurfaceScalarField
        template<class Type>
        tmp<polySurfaceScalarField> fieldNorm
        (
            const DimensionedField<Type, polySurfaceGeoMesh>&
        );

        //- Calculate the chosen norm of the field and return true if successful
        virtual bool calc();


public:

    //- Runtime type information
    TypeName("norm");


    // Constructors

        //- Construct from name, Time and dictionary
        norm
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- No copy construct
        norm(const norm&) = delete;

        //- No copy assignment
        void operator=(const norm&) = delete;


    //- Destructor
    virtual ~norm() = default;


    // Member Functions

        //- Read the dictionary data
        virtual bool read(const dictionary&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
