/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2022 OpenCFD Ltd.
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
    Foam::expressions::volumeExpr::parseDriver

Description
    Driver for volume, surface, point field expressions

    Additional Properties
    \table
        Property     | Description                          | Required | Default
        dimensions   | Dimensions for the expression result | no  |
    \endtable

    In addition to the standard mathematical functions, operations and
    logical and relational operations, the volume expressions support the
    following driver-specific functions:

    Functions
    \table
        Function    | Description                      | Number of arguments |
        vol         | The cell volumes                      | 0 |
        pos         | The cell centres                      | 0 |
        pts         | The cell points                       | 0 |
        area        | The face area magnitudes              | 0 |
        fpos        | The face centres                      | 0 |
        weightAverage| Volume or area weighted average      | 1 |
        weightSum   | Volume or area weighted sum           | 1 |
        face        | The face areaNormal vectors           | 0 |
        face        | A surface-field face value            | 1 |
        point       | A point-field point value             | 1 |
        cellToFace  | Interpolate cell values onto faces    | 1 |
        cellToPoint | Interpolate cell values onto points   | 1 |
        pointToCell | Interpolate point values onto cells   | 1 |
        reconstruct | Reconstruct cell vector from surface scalar | 1 |
        rand        | Random field                          | 0/1 |
    \endtable

    Selections
    \table
        Function| Description                           | Number of arguments |
        cset    | Logical vol field corresponding to cellSet    | 1 |
        fset    | Logical surf field corresponding to faceSet   | 1 |
        pset    | Logical point field corresponding to pointSet | 1 |
        czone   | Logical vol field corresponding to cellZone   | 1 |
        fzone   | Logical surf field corresponding to faceZone  | 1 |
        pzone   | Logical point field corresponding to pointZone| 1 |
    \endtable

Note
    Use namespace debug switch \c volumeExpr for scanner (2), parser (4)
    or dictionary controls as per Foam::expressions::exprDriver.

SourceFiles
    volumeExprDriver.C

\*---------------------------------------------------------------------------*/

#ifndef expressions_volumeExprDriver_H
#define expressions_volumeExprDriver_H

#include "volumeExprFwd.H"
#include "fvExprDriver.H"
#include "exprFieldAssociation.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "pointFields.H"
#include "genericRagelLemonDriver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace expressions
{
namespace volumeExpr
{

/*---------------------------------------------------------------------------*\
                         Class parseDriver Declaration
\*---------------------------------------------------------------------------*/

class parseDriver
:
    public parsing::genericRagelLemonDriver,
    public expressions::fvExprDriver
{
protected:

    // Protected Data

        //- The referenced mesh
        const fvMesh& mesh_;

        //- The results (volume, surface, point)
        autoPtr<regIOobject> resultField_;

        //- The result type-name.
        //  Normally volScalarField, surfaceVectorField etc,
        //  but Scalar is modified for logical as volScalarField etc
        word resultType_;

        //- A logical (bool-like) field (but actually a scalar)
        bool isLogical_;

        //- Requested use of dimensions
        bool hasDimensions_;

        //- A volume/surface/point field
        expressions::FieldAssociation fieldGeoType_;

        //- The result dimensions
        dimensionSet resultDimensions_;


    // Protected Member Functions

        //- Deep-copy the internalField as a result.
        //  Uses the isLogical() and isPointData() values to handle
        //  additional bookkeeping.
        //  For isLogical(), renames the resultType_ from '*Scalar*'
        //  to '*Logical*' (eg, volLogicalField)
        template<class Type>
        void setInternalFieldResult(const Field<Type>& fld);

        //- Cell selections (as logical)
        tmp<volScalarField> field_cellSelection
        (
            const word& name,
            enum topoSetSource::sourceType setType
        ) const;

        //- Face selections (as logical)
        tmp<surfaceScalarField> field_faceSelection
        (
            const word& name,
            enum topoSetSource::sourceType setType
        ) const;

        //- Point selections (as logical)
        tmp<pointScalarField> field_pointSelection
        (
            const word& name,
            enum topoSetSource::sourceType setType
        ) const;


public:

    ClassName("volumeExpr::driver");

    // Generated Methods

        // No copy copy construct
        parseDriver(const parseDriver&) = delete;

        // No copy assignment
        void operator=(const parseDriver&) = delete;


    // Constructors

        //- Construct for specified mesh, with dictionary information
        explicit parseDriver
        (
            const fvMesh& mesh,
            const dictionary& dict = dictionary::null
        );

        //- Construct for specified mesh with copy of driver context
        parseDriver
        (
            const fvMesh& mesh,
            const parseDriver& driver,
            const dictionary& dict
        );

        //- Construct with meshName for the given mesh
        parseDriver(const word& meshName, const fvMesh& mesh);

        //- Construct with patchName and region specified in dictionary
        parseDriver(const dictionary& dict, const fvMesh& mesh);


    // Not generally clonable


    //- Destructor
    virtual ~parseDriver() = default;


    // Public Member Functions

        //- The mesh we are attached to
        virtual const fvMesh& mesh() const
        {
            return mesh_;
        }

        //- The natural field size for the expression
        virtual label size() const
        {
            return mesh_.nCells();
        }

        //- The point field size for the expression
        virtual label pointSize() const
        {
            return mesh_.nPoints();
        }

        //- Field size associated with different geometric field types
        inline label size(const FieldAssociation geoType) const;

        //- Apply dimensions() to geometric fields
        inline bool hasDimensions() const noexcept;

        //- The preferred result dimensions (if any)
        inline const dimensionSet& dimensions() const noexcept;


        //- Clear out local copies of the field
        void clearField();


    // Reading

        //- Read variables, tables etc.
        //  Adds support for "dimensions"
        virtual bool readDict(const dictionary& dict);


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


    // Field Information

        //- The result type-name.
        //  Normally volScalarField, surfaceVectorField etc,
        //  but Scalar is modified for logical as volScalarField etc
        const word& resultType() const noexcept
        {
            return resultType_;
        }

        //- The geometric field association
        FieldAssociation fieldAssociation() const noexcept
        {
            return fieldGeoType_;
        }

        //- A logical (bool-like) field. Actually stored as a scalar.
        bool isLogical() const noexcept
        {
            return isLogical_;
        }

        //- A volume field
        bool isVolumeData() const noexcept
        {
            return fieldGeoType_ == FieldAssociation::VOLUME_DATA;
        }

        //- A surface field
        bool isFaceData() const noexcept
        {
            return fieldGeoType_ == FieldAssociation::FACE_DATA;
        }

        //- A point field
        bool isPointData() const noexcept
        {
            return fieldGeoType_ == FieldAssociation::POINT_DATA;
        }

        //- Test if stored result pointer is the specified type
        template<class GeoField>
        const GeoField* isResultType() const;

        //- Test if stored result pointer is the specified type
        //- and matches the specified logical type
        template<class GeoField>
        const GeoField* isResultType(bool logical, bool dieOnNull=false) const;

        //- A zero-initialized field with the same type as the result field.
        autoPtr<regIOobject> dupZeroField() const;


    // Set Fields

        //- Set result (vol field)
        template<class Type>
        void setResult(VolumeField<Type>* ptr, bool logical = false);

        //- Set result (surface field)
        template<class Type>
        void setResult(SurfaceField<Type>* ptr, bool logical = false);

        //- Set result (point field)
        template<class Type>
        void setResult(PointField<Type>* ptr, bool logical = false);


    // New Fields

        //- Return a new volume field with the mesh size
        template<class Type>
        tmp<VolumeField<Type>>
        newVolField(const Type& val = pTraits<Type>::zero) const;

        //- Return a new surface field with the mesh nInternalFaces size
        template<class Type>
        tmp<SurfaceField<Type>>
        newSurfaceField(const Type& val = pTraits<Type>::zero) const;

        //- Return a new point field with the mesh nPoints size
        template<class Type>
        tmp<PointField<Type>>
        newPointField(const Type& val = pTraits<Type>::zero) const;


        //- Retrieve field (vol field)
        template<class Type>
        tmp<VolumeField<Type>>
        getVolField(const word& fldName, bool getOldTime=false);

        //- Retrieve field (surface field)
        template<class Type>
        tmp<SurfaceField<Type>>
        getSurfaceField(const word& fldName, bool getOldTime=false);

        //- Retrieve field (surface field)
        template<class Type>
        tmp<PointField<Type>>
        getPointField(const word& fldName, bool getOldTime=false);


    // Field "shape" conversions

        //- Interpolate cell to face values
        template<class Type>
        tmp<SurfaceField<Type>>
        cellToFace(const VolumeField<Type>& field) const;

        //- Interpolate cell to point values
        template<class Type>
        tmp<PointField<Type>>
        cellToPoint(const VolumeField<Type>& field) const;

        //- Interpolate point to cell values
        template<class Type>
        tmp<VolumeField<Type>>
        pointToCell(const PointField<Type>& field) const;


    // Custom Field Functions

        //- The volume-weighted average of a field
        template<class Type>
        Type volAverage(VolumeField<Type>& fld) const
        {
            return weightedAverage(fld.mesh().V(), fld.primitiveField());
        }

        //- The volume-weighted sum of a field
        template<class Type>
        Type volSum(VolumeField<Type>& fld) const
        {
            return weightedSum(fld.mesh().V(), fld.primitiveField());
        }

        //- The area-weighted average of a field
        template<class Type>
        Type areaAverage(SurfaceField<Type>& fld) const
        {
            return weightedAverage
            (
                fld.mesh().magSf().primitiveField(),
                fld.primitiveField()
            );
        }

        //- The area-weighted sum of a field
        template<class Type>
        Type areaSum(SurfaceField<Type>& fld) const
        {
            return weightedSum
            (
                fld.mesh().magSf().primitiveField(),
                fld.primitiveField()
            );
        }


        //- The cell volumes - (swak = vol)
        tmp<volScalarField> field_cellVolume() const;

        //- The cell centres - (swak = pos)
        tmp<volVectorField> field_cellCentre() const;

        //- The face area magnitudes [magSf] - (swak = area)
        tmp<surfaceScalarField> field_faceArea() const;

        //- The face centres - (swak = fpos)
        tmp<surfaceVectorField> field_faceCentre() const;

        //- The face areas with their vector direction [Sf] - (swak = face)
        tmp<surfaceVectorField> field_areaNormal() const;

        //- The mesh point locations - (swak = pts)
        tmp<pointVectorField> field_pointField() const;


        //- Cell selection (set)
        inline tmp<volScalarField> field_cellSet(const word& name) const;

        //- Cell selection (zone)
        inline tmp<volScalarField> field_cellZone(const word& name) const;

        //- Face selection (set)
        inline tmp<surfaceScalarField> field_faceSet(const word& name) const;

        //- Face selection (zone)
        inline tmp<surfaceScalarField> field_faceZone(const word& name) const;

        //- Point selection (set)
        inline tmp<pointScalarField> field_pointSet(const word& name) const;

        //- Point selection (zone)
        inline tmp<pointScalarField> field_pointZone(const word& name) const;

        //- A uniform random field
        tmp<volScalarField> field_rand(label seed=0, bool gaussian=false) const;

        //- A Gaussian random field
        tmp<volScalarField> field_randGaussian(label seed=0) const
        {
            return field_rand(seed, true);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace volumeExpr
} // End namespace expressions
} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "volumeExprDriverI.H"

#ifdef NoRepository
    #include "volumeExprDriverTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
