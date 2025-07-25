/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 OpenFOAM Foundation
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
    Foam::fv::FixedValueConstraint

Group
    grpFvOptionsConstraints

Description
    Constrain values of given fields of \c Type with
    a given maximum value within a specified region,
    where \c \<Type\>=Scalar/Vector/SphericalTensor/SymmTensor/Tensor.

Usage
    Minimal example by using \c constant/fvOptions:
    \verbatim
    \<Type\>FixedValueConstraint1
    {
        // Mandatory entries (unmodifiable)
        type            \<Type\>FixedValueConstraint;

        // Mandatory entries (runtime modifiable)
        fieldValues
        {
            \<fieldName1\>           \<value1\>;
            \<fieldName2\>           \<value2\>;
            ...
            \<fieldNameN\>           \<valueN\>;
        }

        // Mandatory/Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property  | Description                             | Type | Reqd  | Dflt
      type      | Type name: \<Type\>FixedValueConstraint | word | yes   | -
      fieldValues | Dictionary containing pairs of names and maximum <!--
           --> values of operand fields of \<Type\> | dictionary | yes   | -
    \endtable

    The inherited entries are elaborated in:
     - \link fvOption.H \endlink
     - \link cellSetOption.H \endlink


    For example to set the turbulence properties within a porous region:
    \verbatim
    porosityTurbulence
    {
        type            scalarFixedValueConstraint;
        active          yes;

        selectionMode   cellZone;
        cellZone        porosity;
        fieldValues
        {
            k           1;
            epsilon     150;
        }
    }
    \endverbatim

See also
  - Foam::fv::fixedTemperatureConstraint
  - Foam::fv::velocityDampingConstraint

SourceFiles
    FixedValueConstraint.C
    fixedValueConstraints.C

\*---------------------------------------------------------------------------*/

#ifndef FixedValueConstraint_H
#define FixedValueConstraint_H

#include "cellSetOption.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                       Class FixedValueConstraint Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class FixedValueConstraint
:
    public fv::cellSetOption
{
    // Private Data

        //- Field values
        List<Type> fieldValues_;


public:

    //- Runtime type information
    TypeName("FixedValueConstraint");


    // Constructors

        //- Construct from components
        FixedValueConstraint
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- No copy construct
        FixedValueConstraint(const FixedValueConstraint&) = delete;

        //- No copy assignment
        void operator=(const FixedValueConstraint&) = delete;


    //- Destructor
    virtual ~FixedValueConstraint() = default;


    // Member Functions

        //- Read source dictionary
        virtual bool read(const dictionary& dict);

        //- Set value on field
        virtual void constrain(fvMatrix<Type>& eqn, const label fieldi);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "FixedValueConstraint.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
