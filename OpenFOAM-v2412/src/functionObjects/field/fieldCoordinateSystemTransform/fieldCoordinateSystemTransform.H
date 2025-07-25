/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2015-2020 OpenCFD Ltd.
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
    Foam::functionObjects::fieldCoordinateSystemTransform

Group
    grpFieldFunctionObjects

Description
    Transforms a user-specified selection of fields from global Cartesian
    coordinates to a local user-specified coordinate system.

    Operands:
    \table
      Operand      | Type                          | Location
      input        | {vol,surface}\<Type\>Field(s) <!--
               --> | $FOAM_CASE/\<time\>/\<inpField\>s
      output file  | -                             | -
      output field | {vol,surface}\<Type\>Field(s) <!--
               --> | $FOAM_CASE/\<time\>/\<outField\>s
    \endtable

    where \c \<Type\>=Scalar/Vector/SphericalTensor/SymmTensor/Tensor.

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    fieldCoordinateSystemTransform1
    {
        // Mandatory entries (unmodifiable)
        type    fieldCoordinateSystemTransform;
        libs    (fieldFunctionObjects);

        // Mandatory entries (runtime modifiable)
        fields  ( U UMean UPrime2Mean );
        coordinateSystem
        {
            origin      (0.001 0 0);
            rotation
            {
                type    axes;
                e1      (1 0.15 0);
                e3      (0 0 -1);
            }
        }

        // Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property  | Description                              | Type | Req'd | Dflt
      type      | Type name: fieldCoordinateSystemTransform | word |  yes  | -
      libs      | Library name: fieldFunctionObjects | word |  yes  | -
      fields    | Names of the operand fields               | wordList | yes | -
      coordinateSystem | Local coordinate system            | dict  | yes | -

    The inherited entries are elaborated in:
     - \link functionObject.H \endlink
     - \link coordinateSystem.H \endlink

    Usage by the \c postProcess utility is not available.

See also
    - Foam::functionObject
    - Foam::functionObjects::fvMeshFunctionObject
    - Foam::coordinateSystem
    - ExtendedCodeGuide::functionObjects::field::fieldCoordinateSystemTransform

SourceFiles
    fieldCoordinateSystemTransform.C
    fieldCoordinateSystemTransformTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_fieldCoordinateSystemTransform_H
#define functionObjects_fieldCoordinateSystemTransform_H

#include "fvMeshFunctionObject.H"
#include "coordinateSystem.H"
#include "volFieldSelection.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
               Class fieldCoordinateSystemTransform Declaration
\*---------------------------------------------------------------------------*/

class fieldCoordinateSystemTransform
:
    public fvMeshFunctionObject
{
protected:

    // Protected Data

        //- Fields to transform
        volFieldSelection fieldSet_;

        //- Coordinate system to transform to
        autoPtr<coordinateSystem> csysPtr_;

        //- Demand-driven non-uniform rotation field (surface fields)
        //  Eg, for cylindrical coordinates
        mutable autoPtr<surfaceTensorField> rotTensorSurface_;

        //- Demand-driven non-uniform rotation field (volume fields)
        //  Eg, for cylindrical coordinates
        mutable autoPtr<volTensorField> rotTensorVolume_;


    // Protected Member Functions

        //- Return the name of the transformed field
        word transformFieldName(const word& fieldName) const;

        //- Demand-driven non-uniform rotation field for surface fields
        const surfaceTensorField& srotTensor() const;

        //- Demand-driven non-uniform rotation field for volume fields
        const volTensorField& vrotTensor() const;


        //- Transform the given field
        template<class FieldType>
        void transformField(const FieldType& field);

        //- Transform the given field
        template<class FieldType, class RotationFieldType>
        void transformField
        (
            const RotationFieldType& rot,
            const FieldType& field
        );

        //- Transform the given field if has the specified element type
        template<class Type>
        void transform(const word& fieldName);


public:

    //- Runtime type information
    TypeName("fieldCoordinateSystemTransform");


    // Constructors

        //- Construct from Time and dictionary
        fieldCoordinateSystemTransform
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- No copy construct
        fieldCoordinateSystemTransform(const fieldCoordinateSystemTransform&)
            = delete;

        //- No copy assignment
        void operator=(const fieldCoordinateSystemTransform&) = delete;


    //- Destructor
    virtual ~fieldCoordinateSystemTransform() = default;


    // Member Functions

        //- Read the input data
        virtual bool read(const dictionary&);

        //- Calculate the transformed fields
        virtual bool execute();

        //- Write the transformed fields
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fieldCoordinateSystemTransformTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
