/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019 OpenCFD Ltd.
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
    Foam::functionObjects::filmFlux

Group
    grpFieldFunctionObjects

Description
    Calculates the mass flux for surface film models.

    Operands:
    \table
      Operand      | Type              | Location
      input        | dictionary        | constant/\<surfaceFilmProperties\>
      output file  | -                 | -
      output field | surfaceScalarField   <!--
                   -->                 | \<time\>/\<filmRegion\>/\<outField\>
    \endtable

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    filmFlux1
    {
        // Mandatory entries (unmodifiable)
        type        filmFlux;

        // Optional entries (runtime modifiable)
        film        <surfaceFilmPropertiesDictionary>;
        result      <resultName>;

        // Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property  | Description              | Type | Reqd | Dflt
      type      | Type name: filmFlux      | word | yes  | -
      film      | Name of surface film properties dictionary <!--
                -->                        | dict | no   | surfaceFilmProperties
      result    | Name of result field     | word | no   | filmFlux
    \endtable

    The inherited entries are elaborated in:
     - \link stateFunctionObject.H \endlink

    Usage by the \c postProcess utility is not available.

See also
  - Foam::functionObjects::stateFunctionObject

SourceFiles
    flux.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_filmFlux_H
#define functionObjects_filmFlux_H

#include "stateFunctionObject.H"
#include "kinematicSingleLayer.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                          Class filmFlux Declaration
\*---------------------------------------------------------------------------*/

class filmFlux
:
    public stateFunctionObject
{
    // Private Data

        typedef regionModels::surfaceFilmModels::kinematicSingleLayer filmType;

        //- Film model name
        word filmName_;

        //- Name of derived flux field
        word resultName_;


    // Private Member Functions

        //- Helper function to set the pointer to the film model
        const filmType& filmModel();


public:

    //- Runtime type information
    TypeName("filmFlux");


    // Constructors

        //- Construct from Time and dictionary
        filmFlux
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- No copy construct
        filmFlux(const filmFlux&) = delete;

        //- No copy assignment
        void operator=(const filmFlux&) = delete;


    //- Destructor
    virtual ~filmFlux() = default;


    // Member Functions

        //- Read the field data
        virtual bool read(const dictionary&);

        //- Execute
        virtual bool execute();

        //- Write the field
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
