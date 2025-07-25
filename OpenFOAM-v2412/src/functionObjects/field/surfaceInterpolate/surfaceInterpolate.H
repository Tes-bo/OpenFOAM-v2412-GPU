/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::functionObjects::surfaceInterpolate

Group
    grpFieldFunctionObjects

Description
    Linearly interpolates volume fields to generate surface fields.

    Operands:
    \table
      Operand        | Type                  | Location
      input          | vol\<Type\>Field      | $FOAM_CASE/\<time\>/\<inpField\>
      output file    | -                     | -
      output field   | surface\<Type\>Field  | $FOAM_CASE/\<time\>/\<outField\>
    \endtable

    where \c \<Type\>=Scalar/Vector/SphericalTensor/SymmTensor/Tensor.

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    surfaceInterpolate1
    {
        // Mandatory entries (unmodifiable)
        type        surfaceInterpolate;
        libs        (fieldFunctionObjects);

        // Mandatory entries (runtime modifiable)
        fields      ((<inpField1> <outField1>) ... (<inpFieldN> <outFieldN>));

        // Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                           | Type | Req'd | Dflt
      type         | Type name: surfaceInterpolate         | word |  yes  | -
      libs         | Library name: fieldFunctionObjects    | word |  yes  | -
      fields       | Names of operand and output fields    | wordList | yes | -
    \endtable

    The inherited entries are elaborated in:
     - \link functionObject.H \endlink

    Usage by the \c postProcess utility is not available.

See also
    - Foam::functionObject
    - Foam::functionObjects::fvMeshFunctionObject
    - ExtendedCodeGuide::functionObjects::field::surfaceInterpolate

SourceFiles
    surfaceInterpolate.C
    surfaceInterpolateTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_surfaceInterpolate_H
#define functionObjects_surfaceInterpolate_H

#include "fvMeshFunctionObject.H"
#include "surfaceFieldsFwd.H"
#include "Tuple2.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class objectRegistry;
class dictionary;
class mapPolyMesh;

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                    Class surfaceInterpolate Declaration
\*---------------------------------------------------------------------------*/

class surfaceInterpolate
:
    public fvMeshFunctionObject
{
protected:

    // Protected Data

        //- Fields to process
        List<Tuple2<word, word>> fieldSet_;


    // Protected Member Functions

        //- Linearly interpolate volume fields to generate surface fields
        template<class Type>
        void interpolateFields();


public:

    //- Runtime type information
    TypeName("surfaceInterpolate");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        surfaceInterpolate
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- No copy construct
        surfaceInterpolate(const surfaceInterpolate&) = delete;

        //- No copy assignment
        void operator=(const surfaceInterpolate&) = delete;


    //- Destructor
    virtual ~surfaceInterpolate() = default;


    // Member Functions

        //- Read the controls
        virtual bool read(const dictionary&);

        //- Calculate the interpolated fields
        virtual bool execute();

        //- Write the interpolated fields
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "surfaceInterpolateTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
