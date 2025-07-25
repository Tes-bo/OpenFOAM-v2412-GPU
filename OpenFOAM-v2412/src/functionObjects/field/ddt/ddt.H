/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017 OpenFOAM Foundation
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
    Foam::functionObjects::ddt

Group
    grpFieldFunctionObjects

Description
    Computes the Eulerian time derivative of an input volume field for
    time-variant simulations (not appropriate to steady-state simulations).

    Operands:
    \table
      Operand           | Type             | Location
      input             | vol\<Type\>Field | $FOAM_CASE/\<time\>/\<inpField\>
      output file       | -                | -
      output field      | vol\<Type\>Field | $FOAM_CASE/\<time\>/\<outField\>
    \endtable

    where \c \<Type\>=Scalar/Vector/SphericalTensor/SymmTensor/Tensor.

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    ddt1
    {
        // Mandatory entries (unmodifiable)
        type            ddt;
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
      type         | Type name: ddt                     | word |  yes  | -
      libs         | Library name: fieldFunctionObjects | word |  yes  | -
      field        | Name of operand field              | word |  yes  | -
    \endtable

    The inherited entries are elaborated in:
     - \link functionObject.H \endlink
     - \link fieldExpression.H \endlink

    Minimal example by using the \c postProcess utility:
    \verbatim
        postProcess -func "ddt(<field>)"
    \endverbatim

See also
    - Foam::functionObject
    - Foam::functionObjects::fvMeshFunctionObject
    - Foam::functionObjects::fieldExpression
    - ExtendedCodeGuide::functionObjects::field::ddt

SourceFiles
    ddt.C
    ddtTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_ddt_H
#define functionObjects_ddt_H

#include "fieldExpression.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                         Class ddt Declaration
\*---------------------------------------------------------------------------*/

class ddt
:
    public fieldExpression
{

    // Private Member Functions

        //- Calculate and store the result
        //  \return true if fieldName_ was known and had the correct type
        template<class Type>
        bool calcDdt();

        //- Calculate time derivative of the field
        //  \return true on success
        virtual bool calc();


public:

    //- Runtime type information
    TypeName("ddt");


    // Constructors

        //- Construct from Time and dictionary
        ddt
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- No copy construct
        ddt(const ddt&) = delete;

        //- No copy assignment
        void operator=(const ddt&) = delete;


    //- Destructor
    virtual ~ddt() = default;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ddtTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
