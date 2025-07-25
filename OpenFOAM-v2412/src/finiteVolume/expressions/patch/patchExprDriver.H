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
    Foam::expressions::patchExpr::parseDriver

Description
    Driver for patch expressions

    In addition to the standard mathematical functions, operations and
    logical and relational operations, the patch expressions support the
    following driver-specific functions:

    Functions
    \table
        Function    | Description                      | Number of arguments |
        pos         | The face centres                      | 0 |
        pts         | The face points                       | 0 |
        area        | The face area magnitudes              | 0 |
        weightAverage| Area weighted average                | 1 |
        weightSum   | Area weighted sum                     | 1 |
        face        | The face area normal vectors          | 0 |
        normal      | The face unit normal vectors          | 0 |
        point       | A point-field point value             | 1 |
        faceToPoint | Interpolate face values onto points   | 1 |
        pointToFace | Interpolate point values onto faces   | 1 |
        rand        | Random field                          | 0/1 |
        snGrad      | Surface normal field                  | 0 |
        internalField | Internal field next to patch        | 0 |
        neighbourField | patch field on opposite side of coupled patch | 0 |
    \endtable

Note
    Use namespace debug switch \c patchExpr for scanner (2), parser (4)
    or dictionary controls as per Foam::expressions::exprDriver.

SourceFiles
    patchExprDriver.C
    patchExprDriverFields.C
    patchExprDriverTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef expressions_patchExprDriver_H
#define expressions_patchExprDriver_H

#include "patchExprFwd.H"
#include "fvExprDriver.H"
#include "exprFieldAssociation.H"
#include "Enum.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "pointFields.H"
#include "genericRagelLemonDriver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace expressions
{
namespace patchExpr
{

/*---------------------------------------------------------------------------*\
                         Class parseDriver Declaration
\*---------------------------------------------------------------------------*/

class parseDriver
:
    public parsing::genericRagelLemonDriver,
    public expressions::fvExprDriver
{
    // Private Member Functions

        static const fvPatch& getFvPatch
        (
            const fvMesh& fvm,
            const dictionary& dict
        );


protected:

    // Protected Data

        //- The referenced patch
        const fvPatch& patch_;


    // Protected Member Functions

        //- Cell selections (as logical)
        tmp<boolField> field_cellSelection
        (
            const word& name,
            enum topoSetSource::sourceType setType
        ) const;

        //- Face selections (as logical)
        tmp<boolField> field_faceSelection
        (
            const word& name,
            enum topoSetSource::sourceType setType
        ) const;


public:

    ClassName("patchExpr::driver");

    // Generated Methods

        // No copy copy construct
        parseDriver(const parseDriver&) = delete;

        // No copy assignment
        void operator=(const parseDriver&) = delete;


    // Constructors

        //- Construct for specified patch, with dictionary information
        explicit parseDriver
        (
            const fvPatch& p,
            const dictionary& dict = dictionary::null
        );

        //- Construct for specified patch with copy of driver context
        parseDriver
        (
            const fvPatch& p,
            const parseDriver& driver,
            const dictionary& dict  // = dictionary::null
        );

        //- Construct with patchName for the given mesh
        parseDriver(const word& patchName, const fvMesh& mesh);

        //- Construct with "patch" (mandatory) and "region" (optional)
        //- specified in dictionary
        parseDriver(const dictionary& dict, const fvMesh& mesh);


    // Not generally clonable


    //- Destructor
    virtual ~parseDriver() = default;


    // Public Member Functions

        //- The mesh we are attached to
        virtual const fvMesh& mesh() const
        {
            return patch_.boundaryMesh().mesh();
        }

        //- The natural field size for the expression
        virtual label size() const
        {
            return patch_.patch().size();
        }

        //- The point field size for the expression
        virtual label pointSize() const
        {
            return patch_.patch().nPoints();
        }

        //- Field size associated with different geometric field types
        inline label size(const FieldAssociation geoType) const;


    // Evaluation

        //- Perform parsing on (sub) string
        using genericRagelLemonDriver::content;

        //- Execute the parser.
        //  The return value currently has no meaning.
        virtual unsigned parse
        (
            const std::string& expr,
            size_t pos = 0,
            size_t len = std::string::npos
        );


    // General

        //- Set result
        template<class Type>
        void setResult(Field<Type>* ptr, bool pointVal = false)
        {
            result().setResult<Type>(ptr, pointVal);
        }


        //- Retrieve variable as field if possible.
        //  Test tmp for validity to determine success of the operation.
        template<class Type>
        tmp<Field<Type>> getVariableIfAvailable(const word& fldName) const;


    // Field Retrieval

        //- Return named field
        template<class Type>
        tmp<Field<Type>> getField(const word& fldName);

        //- Retrieve field (vol field)
        template<class Type>
        tmp<Field<Type>> getVolField(const word& fldName);

        //- Retrieve field (surface field)
        template<class Type>
        tmp<Field<Type>> getSurfaceField(const word& fldName);

        //- Retrieve field (point field)
        template<class Type>
        tmp<Field<Type>> getPointField(const word& fldName);

        //- Return internal field next to patch
        template<class Type>
        tmp<Field<Type>> patchInternalField(const word& fldName);

        //- Return patchField on the opposite patch of a coupled patch
        template<class Type>
        tmp<Field<Type>> patchNeighbourField(const word& fldName);

        //- Return surface normal field (snGrad)
        template<class Type>
        tmp<Field<Type>> patchNormalField(const word& fldName);


    // Field "shape" conversions

        //- Interpolate face to point
        template<class Type>
        tmp<Field<Type>> faceToPoint(const Field<Type>& field) const;

        //- Interpolate point to face values
        template<class Type>
        tmp<Field<Type>> pointToFace(const Field<Type>& field) const;


    // Custom Field Functions

        //- The area-weighted average of a field
        template<class Type>
        Type areaAverage(const Field<Type>& fld) const
        {
            return weightedAverage(patch_.magSf(), fld);
        }

        //- The area-weighted sum of a field
        template<class Type>
        Type areaSum(const Field<Type>& fld) const
        {
            return weightedSum(patch_.magSf(), fld);
        }

        //- The face area magnitudes [magSf] - (swak = area)
        tmp<scalarField> field_faceArea() const;

        //- The face centres - (swak = pos)
        tmp<vectorField> field_faceCentre() const;

        //- The face areas with their vector direction [Sf] - (swak = face)
        tmp<vectorField> field_areaNormal() const;

        //- The face unit normal direction [nf] - (expression: normal)
        tmp<vectorField> field_unitNormal() const;

        //- The patch point locations - (swak = pts)
        tmp<vectorField> field_pointField() const;


        //- Cell selection (set)
        inline tmp<boolField> field_cellSet(const word& name) const;

        //- Cell selection (zone)
        inline tmp<boolField> field_cellZone(const word& name) const;

        //- Face selection (set)
        inline tmp<boolField> field_faceSet(const word& name) const;

        //- Face selection (zone)
        inline tmp<boolField> field_faceZone(const word& name) const;


        //- A uniform random field
        tmp<scalarField> field_rand(label seed=0, bool gaussian=false) const;

        //- A Gaussian random field
        tmp<scalarField> field_randGaussian(label seed=0) const
        {
            return field_rand(seed, true);
        }
};


// Template specializations

//- Retrieve field (surface field: bool)
template<>
tmp<Field<bool>> parseDriver::getSurfaceField<bool>(const word& fldName);

//- Retrieve field (point field: bool)
template<>
tmp<Field<bool>> parseDriver::getPointField<bool>(const word& fldName);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace patchExpr
} // End namespace expressions
} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "patchExprDriverI.H"

#ifdef NoRepository
    #include "patchExprDriverTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
