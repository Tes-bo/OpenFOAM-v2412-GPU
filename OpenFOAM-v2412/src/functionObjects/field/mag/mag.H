/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2016 OpenFOAM Foundation
    Copyright (C) 2019-2020 OpenCFD Ltd.
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
    Foam::functionObjects::mag

Group
    grpFieldFunctionObjects

Description
    Computes the magnitude of an input field.

    Operands:
    \table
      Operand      | Type                     | Location
      input      | {vol,surface}\<Type\>Field | $FOAM_CASE/\<time\>/\<inpField\>
      output file  | -                        | -
      output field | {vol,surface}ScalarField | $FOAM_CASE/\<time\>/\<outField\>
    \endtable

    where \c \<Type\>=Scalar/Vector/SphericalTensor/SymmTensor/Tensor.

    With the \c subRegion option, also supports fields on function object
    surface output (e.g. sampledSurfaces).

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    mag1
    {
        // Mandatory entries (unmodifiable)
        type            mag;
        libs            (fieldFunctionObjects);

        // Mandatory (inherited) entries (runtime modifiable)
        field           <field>;

        // Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                        | Type | Req'd | Dflt
      type         | Type name: mag                     | word |  yes  | -
      libs         | Library name: fieldFunctionObjects | word |  yes  | -
      field        | Name of the operand field          | word |  yes  | -
    \endtable

    The inherited entries are elaborated in:
     - \link functionObject.H \endlink
     - \link fieldExpression.H \endlink

    Minimal example by using the \c postProcess utility:
    \verbatim
        postProcess -func "mag(<field>)"
    \endverbatim

See also
    - Foam::functionObject
    - Foam::functionObjects::fvMeshFunctionObject
    - Foam::functionObjects::fieldExpression
    - ExtendedCodeGuide::functionObjects::field::mag

SourceFiles
    mag.C
    magTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_mag_H
#define functionObjects_mag_H

#include "fieldExpression.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                           Class mag Declaration
\*---------------------------------------------------------------------------*/

class mag
:
    public fieldExpression
{
    // Private Member Functions

        //- Calculate the magnitude of the field and register the result
        template<class Type>
        bool calcMag();

        //- Calculate the magnitude of the field and return true if successful
        virtual bool calc();


public:

    //- Runtime type information
    TypeName("mag");


    // Constructors

        //- Construct from Time and dictionary
        mag
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- No copy construct
        mag(const mag&) = delete;

        //- No copy assignment
        void operator=(const mag&) = delete;


    //- Destructor
    virtual ~mag() = default;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "magTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
