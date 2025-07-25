/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2020 OpenCFD Ltd.
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
    Foam::functionObjects::columnAverage

Group
    grpFieldFunctionObjects

Description
    Computes the arithmetic average of given quantities along columns of cells
    in a given direction for structured-like layered meshes. It is, for example,
    useful for channel-like cases where spanwise average of a field is desired.
    However, the \c columnAverage function object does not operate on arbitrary
    unstructured meshes.

    For each patch face, calculates the average value of all cells attached in
    the patch face normal direction, and then pushes the average value back
    to all cells in the column.

    Operands:
    \table
      Operand           | Type               | Location
      input             | vol\<Type\>Field   | $FOAM_CASE/\<time\>/\<inpField\>
      output file       | -                  | -
      output field      | vol\<Type\>Field   | $FOAM_CASE/\<time\>/\<outField\>
    \endtable

    where \c \<Type\>=Scalar/Vector/SphericalTensor/SymmTensor/Tensor.

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    columnAverage1
    {
        // Mandatory entries (unmodifiable)
        type            columnAverage;
        libs            (fieldFunctionObjects);

        // Mandatory entries (runtime modifiable)
        patches         (<patch1> <patch2> ... <patchN>);
        fields          (<field1> <field2> ... <fieldN>);

        // Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                        | Type | Req'd | Dflt
      type         | Type name: columnAverage           | word |  yes  | -
      libs         | Library name: fieldFunctionObjects | word |  yes  | -
      patches      | Names of patches to collapse onto  | word |  yes  | -
      fields       | Names of the operand fields        | word |  yes  | -
    \endtable

    The inherited entries are elaborated in:
     - \link functionObject.H \endlink

    Usage by the \c postProcess utility is not available.

See also
    - Foam::functionObject
    - Foam::functionObjects::fvMeshFunctionObject
    - ExtendedCodeGuide::functionObjects::field::columnAverage

SourceFiles
    columnAverage.C
    columnAverageTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_columnAverage_H
#define functionObjects_columnAverage_H

#include "volFieldsFwd.H"
#include "fvMeshFunctionObject.H"
#include "volFieldSelection.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class globalIndex;
class meshStructure;

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                        Class columnAverage Declaration
\*---------------------------------------------------------------------------*/

class columnAverage
:
    public fvMeshFunctionObject
{
    // Private Data

        //- Patches on which to collapse the fields (in sorted order)
        labelList patchIDs_;

        //- Fields to collapse
        volFieldSelection fieldSet_;

        mutable autoPtr<globalIndex> globalFaces_;
        mutable autoPtr<globalIndex> globalEdges_;
        mutable autoPtr<globalIndex> globalPoints_;
        mutable autoPtr<meshStructure> meshStructurePtr_;


    // Private Member Functions

        //- Create the column average field name
        const word averageName(const word& fieldName) const;

        //- Return the column-based addressing
        const meshStructure& meshAddressing(const polyMesh&) const;

        //- Calculate the averaged field and return true if successful
        template<class Type>
        bool columnAverageField(const word& fieldName);


public:

    //- Runtime type information
    TypeName("columnAverage");


    // Constructors

        //- Construct from Time and dictionary
        columnAverage
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- No copy construct
        columnAverage(const columnAverage&) = delete;

        //- No copy assignment
        void operator=(const columnAverage&) = delete;


    //- Destructor
    virtual ~columnAverage() = default;


    // Member Functions

        //- Read the settings
        virtual bool read(const dictionary& dict);

        //- Execute, currently does nothing
        virtual bool execute();

        //- Write the results
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "columnAverageTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
