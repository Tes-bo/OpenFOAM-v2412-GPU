/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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
    Foam::functionObjects::continuityError

Group
    grpFieldFunctionObjects

Description
    Computes local, global and cumulative continuity errors for a flux field.

    Local continuity error, \f$ \epsilon_{local} \f$:
    \f[
        \epsilon_{local} = \Delta_t \langle |x| \rangle
    \f]

    Global continuity error, \f$ \epsilon_{global} \f$:
    \f[
        \epsilon_{global} = \Delta_t \langle |x| \rangle
    \f]

    Cumulative continuity, \f$ \epsilon_{cum} \f$:
    \f[
        \epsilon_{cum} += \epsilon_{global}
    \f]

    where
    \vartable
      \Delta_t          | Time-step size
      \langle . \rangle | Cell-volume weighted average operator
      x                 | \f$ \div \phi \f$
      phi               | Flux field
    \endvartable

    Operands:
    \table
      Operand        | Type | Location
      input          | -    | -
      output file    | dat  | $FOAM_CASE/postProcessing/\<FO\>/\<time\>/\<file\>
      output field   | -    | -
    \endtable

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    continuityError1
    {
        // Mandatory entries (unmodifiable)
        type        continuityError;
        libs        (fieldFunctionObjects);

        // Optional entries (runtime modifiable)
        phi         phi;

        // Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                        | Type | Req'd | Dflt
      type         | Type name: continuityError         | word |  yes  | -
      libs         | Library name: fieldFunctionObjects | word |  yes  | -
      phi          | Name of flux field                 | word |  no   | phi
    \endtable

    The inherited entries are elaborated in:
     - \link functionObject.H \endlink
     - \link writeFile.H \endlink

    Usage by the \c postProcess utility is not available.

See also
    - Foam::functionObject
    - Foam::functionObjects::fvMeshFunctionObject
    - Foam::functionObjects::writeFile
    - ExtendedCodeGuide::functionObjects::field::continuityError

SourceFiles
    continuityError.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_continuityError_H
#define functionObjects_continuityError_H

#include "Switch.H"
#include "fvMeshFunctionObject.H"
#include "writeFile.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                         Class continuityError Declaration
\*---------------------------------------------------------------------------*/

class continuityError
:
    public fvMeshFunctionObject,
    public writeFile
{
protected:

    // Protected Data

        //- Name of flux field
        word phiName_;

        //- Cumulative error
        scalar cumulative_;


    // Protected Member Functions

        //- Output file header information
        virtual void writeFileHeader(Ostream& os);


public:

    //- Runtime type information
    TypeName("continuityError");


    // Constructors

        //- Construct from Time and dictionary
        continuityError
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- No copy construct
        continuityError(const continuityError&) = delete;

        //- No copy assignment
        void operator=(const continuityError&) = delete;


    //- Destructor
    virtual ~continuityError() = default;


    // Member Functions

        //- Read the field min/max data
        virtual bool read(const dictionary&);

        //- Execute, currently does nothing
        virtual bool execute();

        //- Write the continuityError
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
