/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016 OpenFOAM Foundation
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
    Foam::functionObjects::XiReactionRate

Group
    grpFieldFunctionObjects

Description
    Writes the turbulent flame-speed and reaction-rate volScalarFields for the
    Xi-based combustion models.

    Operands:
    \table
      Operand        | Type           | Location
      input          | -              | -
      output file    | -              | -
      output field   | volScalarField | $FOAM_CASE/\<time\>/{St,wdot}
    \endtable

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    XiReactionRate1
    {
        // Mandatory entries (unmodifiable)
        type        XiReactionRate;
        libs        (fieldFunctionObjects);

        // Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property   | Description                        | Type | Req'd | Dflt
      type       | Type name: XiReactionRate          | word |  yes  | -
      libs       | Library name: fieldFunctionObjects | word |  yes  | -
    \endtable

    The inherited entries are elaborated in:
     - \link functionObject.H \endlink

    Usage by the \c postProcess utility is not available.

See also
    - Foam::functionObject
    - Foam::functionObjects::fvMeshFunctionObject
    - ExtendedCodeGuide::functionObjects::field::XiReactionRate

SourceFiles
    XiReactionRate.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_XiReactionRate_H
#define functionObjects_XiReactionRate_H

#include "fvMeshFunctionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class XiReactionRate Declaration
\*---------------------------------------------------------------------------*/

class XiReactionRate
:
    public fvMeshFunctionObject
{
public:

    //- Runtime type information
    TypeName("XiReactionRate");


    // Constructors

        //- Construct from Time and dictionary
        XiReactionRate
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- No copy construct
        XiReactionRate(const XiReactionRate&) = delete;

        //- No copy assignment
        void operator=(const XiReactionRate&) = delete;


    //- Destructor
    virtual ~XiReactionRate() = default;


    // Member Functions

        //- Read the reaction rate data
        virtual bool read(const dictionary&);

        //- Do nothing
        virtual bool execute();

        //- Write the reaction rate fields
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
