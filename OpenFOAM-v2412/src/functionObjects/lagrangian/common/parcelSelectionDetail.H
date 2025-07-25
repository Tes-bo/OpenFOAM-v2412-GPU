/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018 OpenCFD Ltd.
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
    Foam::Detail::parcelSelection

Description
    Selection of parcels based on their objectRegistry entries.
    Normally accessed via a dictionary entry.

    Example sub-dictionary entry
    \verbatim
    selection
    {
        stride
        {
            // every 10th parcelId
            action  use;
            source  stride;
            stride  10;
        }
        injector
        {
            // Only output from injectorID == 1
            action  subset;
            source  field;
            field   typeId;
            accept  (equal 1);
        }
        Umin
        {
            // Remove slow parcels
            action  subtract;
            source  field;
            field   U;
            accept  (less 1e-3);
        }
        diam
        {
            // Only particular diameter ranges
            action  subset;
            source  field;
            field   d;
            accept  (greater 1e-5) and (less 1e-3);
        }
    }
    \endverbatim

    \heading Entry type
    \table
        Property    | Description                           | Required | Default
        action      | all/clear/invert/ignore add/subtract/subset/use | yes |
        source      | field/stride                          | mostly  |
    \endtable

    \heading Stride source
    \table
        Property    | Description                           | Required | Default
        stride      | The stride for the parcel id          | yes |
    \endtable

    \heading Field source
    \table
        Property    | Description                           | Required | Default
        field       | The label/scalar/vector field name    | yes |
        accept      | Acceptance or test criterion          | yes |
    \endtable

    The \c accept criterion has two forms:
      -# single expression
         - (expr)
      -# composite expression
         - (expr) or (expr)
         - (expr) and (expr)

    The expressions are a (op scalar) pair that form a unary scalar
    predicate. The \a op is one of the following:
    - lt, less
    - le, lessEqual, lessEq
    - gt, greater
    - ge, greaterEqual, greaterEq
    - eq, equal
    - neq, notEqual

    For example,
    \verbatim
    accept  (less 10);
    accept  (less 10) or (greater 100);
    accept  (ge 10) and (le 20);
    \endverbatim

See also
    Foam::predicates::scalars
    Foam::functionObjects::vtkCloud

SourceFiles
    parcelSelectionDetail.C

\*---------------------------------------------------------------------------*/

#ifndef parcelSelectionDetail_H
#define parcelSelectionDetail_H

#include "bitSet.H"
#include "Enum.H"
#include "objectRegistry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Detail
{

/*---------------------------------------------------------------------------*\
                   Class Detail::parcelSelection Declaration
\*---------------------------------------------------------------------------*/

class parcelSelection
{
public:

    //- Enumeration defining the valid selection actions
    enum actionType
    {
        ALL,                //!< "all" - select all parcels
        CLEAR,              //!< "clear" - clear the selection
        INVERT,             //!< "invert" - invert the selection
        ADD,                //!< "add" - parcel selection
        SUBTRACT,           //!< "subtract" - remove parcel selection
        SUBSET,             //!< "subset" - subset parcel selection
        USE,                //!< "use" - same as clear + add
        IGNORE,             //!< "ignore" - dummy no-op
    };

    //- Names for the actionType
    static const Enum<actionType> actionNames;


    //- Enumeration defining the valid sources
    enum sourceType
    {
        FIELD,              //!< "field" - select based on field value
        STRIDE              //!< "stride" - select based on stride (parcel id)
    };

    //- Names for the sourceType
    static const Enum<sourceType> sourceNames;


    //- Enumeration defining and/or logic
    enum logicType { AND, OR };

    //- Names for the logicType
    static const Enum<logicType> logicNames;


protected:

    // Protected Data

        //- The filtered parcel addressing. Eg, for the current cloud.
        dictionary parcelSelect_;

        //- The filtered parcel addressing. Eg, for the current cloud.
        bitSet parcelAddr_;


    // Protected Member Functions

        //- Calculate parcel selection filter.
        //  \return True if the filter is applicable
        bool calculateFilter
        (
            const objectRegistry& obrTmp,
            const bool log = true
        );

public:

    // Constructors

        //- Default construct
        parcelSelection() = default;


    //- Destructor
    virtual ~parcelSelection() = default;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Detail
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
