/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2023 OpenCFD Ltd.
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
    Foam::functionObjects::limitFields

Group
    grpFieldFunctionObjects

Description
    Limits input fields to user-specified min and max bounds.

    Operands:
    \table
      Operand          | Type             | Location
      input            | vol\<Type\>Field | $FOAM_CASE/\<time\>/\<inpField\>
      output file      | -                | -
      output field     | vol\<Type\>Field | $FOAM_CASE/\<time\>/\<outField\>
    \endtable

    where \c \<Type\>=Scalar/Vector/SphericalTensor/SymmTensor/Tensor.

Note
    For non-scalar types of input field, the user limit is used to create a
    scaling factor using the field magnitudes.

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    limitFields1
    {
        // Mandatory entries (unmodifiable)
        type            limitFields;
        libs            (fieldFunctionObjects);

        // Mandatory entries (runtime modifiable)
        fields          (U);
        limit           max;
        max             100;

        // Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                        | Type | Req'd | Dflt
      type         | Type name: limitFields             | word |  yes  | -
      libs         | Library name: fieldFunctionObjects | word |  yes  | -
      fields       | List of fields to process          | wordList | yes | -
      limit        | Bound to limit - see below         | word |  yes  | -
      min          | Min limit value                  | scalar | conditional | -
      max          | Max limit value                  | scalar | conditional | -
    \endtable

    Options for the \c limit entry:
    \verbatim
        min  : specify a minimum value
        max  : specify a maximum value
        both : specify a minimum value and a maximum value
    \endverbatim

    The inherited entries are elaborated in:
     - \link functionObject.H \endlink

    Usage by the \c postProcess utility is not available.

See also
    - Foam::functionObject
    - Foam::functionObjects::fvMeshFunctionObject
    - ExtendedCodeGuide::functionObjects::field::limitFields

SourceFiles
    limitFields.C
    limitFieldsTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_limitFields_H
#define functionObjects_limitFields_H

#include "Enum.H"
#include "fvMeshFunctionObject.H"
#include "volFieldSelection.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                         Class limitFields Declaration
\*---------------------------------------------------------------------------*/

class limitFields
:
    public fvMeshFunctionObject
{
public:

    // Public Enumerations

        enum limitType : unsigned
        {
            CLAMP_NONE = 0,     //!< No limit
            CLAMP_MIN = 0x1,    //!< Clamp minimum value
            CLAMP_MAX = 0x2,    //!< Clamp maximum value
            CLAMP_RANGE = (CLAMP_MIN | CLAMP_MAX)  //!< Clamp min/max
        };


protected:

    // Protected Data

        //- Limit type names
        static const Enum<limitType> limitTypeNames_;

        //- Fields to limit
        volFieldSelection fieldSet_;

        //- Limiting type
        limitType withBounds_;

        //- Minimum limit
        scalar min_;

        //- Maximum limit
        scalar max_;


    // Protected Member Functions

        //- Limit a scalar field.
        //  \return true if field of this type was found.
        bool limitScalarField(const word& fieldName);

        //- Limit a field.
        //  \return true if field of this type was found.
        template<class Type>
        bool limitField(const word& fieldName);


public:

    //- Runtime type information
    TypeName("limitFields");


    // Constructors

        //- Construct from Time and dictionary
        limitFields
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- No copy construct
        limitFields(const limitFields&) = delete;

        //- No copy assignment
        void operator=(const limitFields&) = delete;


    //- Destructor
    virtual ~limitFields() = default;


    // Member Functions

        //- Read the field min/max data
        virtual bool read(const dictionary& dict);

        //- Execute, currently does nothing
        virtual bool execute();

        //- Write the limitFields
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "limitFieldsTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
