/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016 OpenFOAM Foundation
    Copyright (C) 2020-2021 OpenCFD Ltd.
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
    Foam::functionObjects::randomise

Group
    grpFieldFunctionObjects

Description
    Adds a random component to an input field,
    with a specified perturbation magnitude.

    Operands:
    \table
      Operand        | Type             | Location
      input          | vol\<Type\>Field | $FOAM_CASE/\<time\>/\<inpField\>
      output file    | -                | -
      output field   | vol\<Type\>Field | $FOAM_CASE/\<time\>/\<outField\>
    \endtable

    where \c \<Type\>=Scalar/Vector/SphericalTensor/SymmTensor/Tensor.

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    randomise1
    {
        // Mandatory entries (unmodifiable)
        type            randomise;
        libs            (fieldFunctionObjects);

        // Mandatory entries (runtime modifiable)
        magPerturbation 0.1;

        // Mandatory (inherited) entries (runtime modifiable)
        field           <field>;

        // Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                          | Type | Req'd | Dflt
      type         | Type name: randomise                 | word |  yes  | -
      libs         | Library name: fieldFunctionObjects   | word |  yes  | -
      regionType   | Selection type: all/cellSet/cellZone | word |  no  | all
      magPerturbation | The magnitude of the perturbation | scalar |  yes  | -
      field  | Name of the operand field            | word | yes  | -
      name   | Name of cellSet/cellZone if required | word | conditional  | -
    \endtable

    The inherited entries are elaborated in:
     - \link functionObject.H \endlink
     - \link fieldExpression.H \endlink

    Minimal example by using the \c postProcess utility:
    \verbatim
        postProcess -func "randomise(<field>)"
    \endverbatim

See also
    - Foam::functionObject
    - Foam::functionObjects::fieldExpression
    - Foam::functionObjects::fvMeshFunctionObject
    - ExtendedCodeGuide::functionObjects::field::randomise

SourceFiles
    randomise.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_randomise_H
#define functionObjects_randomise_H

#include "fieldExpression.H"
#include "volRegion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                           Class randomise Declaration
\*---------------------------------------------------------------------------*/

class randomise
:
    public fieldExpression,
    public volRegion
{
    // Private Data

        //- The magnitude of the perturbation
        scalar magPerturbation_;


    // Private Member Functions

        //- Calculate randomised field and register result
        //  Defined as file-local template implementation
        template<class Type>
        bool calcTemplate();


protected:

    // Protected Member Functions

        //- Calculate the randomised field and return true if successful
        virtual bool calc();

        //- No copy construct
        randomise(const randomise&) = delete;

        //- No copy assignment
        void operator=(const randomise&) = delete;


public:


    //- Runtime type information
    TypeName("randomise");


    // Constructors

        //- Construct from Time and dictionary
        randomise
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~randomise() = default;


    // Member Functions

        //- Read the randomise data
        virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
