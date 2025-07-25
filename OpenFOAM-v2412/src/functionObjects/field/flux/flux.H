/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2020 OpenCFD Ltd.
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
    Foam::functionObjects::flux

Group
    grpFieldFunctionObjects

Description
    Computes the flux of an input vector field.

    Operands:
    \table
      Operand      | Type                     | Location
      input        | {vol,surface}VectorField | $FOAM_CASE/\<time\>/\<inpField\>
      output file  | -                        | -
      output field | surfaceScalarField       | $FOAM_CASE/\<time\>/\<outField\>
    \endtable

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    flowType1
    {
        // Mandatory entries (unmodifiable)
        type            flux;
        libs            (fieldFunctionObjects);

        // Optional entries (runtime modifiable)
        rho             none;

        // Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
        Property     | Description                        | Type | Req'd | Dflt
        type         | Type name: flux                    | word |  yes  | -
        libs         | Library name: fieldFunctionObjects | word |  yes  | -
        rho          | Name of density field              | word |  no   | none
    \endtable

    The inherited entries are elaborated in:
     - \link functionObject.H \endlink
     - \link fieldExpression.H \endlink

    Usage by the \c postProcess utility is not available.

See also
    - Foam::functionObject
    - Foam::functionObjects::fieldExpression
    - Foam::functionObjects::fvMeshFunctionObject
    - ExtendedCodeGuide::functionObjects::field::flux

SourceFiles
    flux.C
    fluxTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_flux_H
#define functionObjects_flux_H

#include "fieldExpression.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                            Class flux Declaration
\*---------------------------------------------------------------------------*/

class flux
:
    public fieldExpression
{
    // Private Data

        //- Name of density field
        word rhoName_;


    // Private Member Functions

        //- Calculate the flux of a volVectorField and register the result
        template<class FieldType, class RhoFieldType>
        bool calcVolFlux(const RhoFieldType& rho);

        //- Calculate the flux of a surfaceVectorField and register the result
        template<class FieldType, class RhoFieldType>
        bool calcSurFlux(const RhoFieldType& rho);

        //- Calculate the flux field and return true if successful
        virtual bool calc();


public:

    //- Runtime type information
    TypeName("flux");


    // Constructors

        //- Construct from Time and dictionary
        flux
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- No copy construct
        flux(const flux&) = delete;

        //- No copy assignment
        void operator=(const flux&) = delete;


    //- Destructor
    virtual ~flux() = default;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fluxTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
