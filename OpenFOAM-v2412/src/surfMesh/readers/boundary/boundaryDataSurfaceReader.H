/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 OpenCFD Ltd.
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
    Foam::boundaryDataSurfaceReader

Description
    A boundaryData format surface reader.
    However, the "surface" represented by boundaryData is actually only
    point data!

    // Points
    <case>/constant/region0/"boundaryData"/patchName/points

    // Values
    <case>/constant/region0/"boundaryData"/patchName/TIME/field

    \verbatim
    readOptions
    {
        boundaryData
        {
            points      points;
        }
    }
    \endverbatim

    Format options for boundaryData:
    \table
        Property | Description                              | Required | Default
        points   | Name of the "points" file                | no  | points
    \endtable

SourceFiles
    boundaryDataSurfaceReader.C
    boundaryDataSurfaceReaderTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_boundaryDataSurfaceReader_H
#define Foam_boundaryDataSurfaceReader_H

#include "surfaceReader.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class Time;

/*---------------------------------------------------------------------------*\
                  Class boundaryDataSurfaceReader Declaration
\*---------------------------------------------------------------------------*/

class boundaryDataSurfaceReader
:
    public surfaceReader
{
    // Private Data

        //- Base directory
        fileName baseDir_;

        //- Name of the 'points' file
        word pointsName_;

        //- Times
        instantList timeValues_;

        //- Field names
        mutable List<word> fieldNames_;

        //- Pointer to the surface
        autoPtr<meshedSurface> surfPtr_;


    // Private Member Functions

        //- Populate time instances
        void readCase();

        //- Read geometry
        void readGeometry(meshedSurface& surf, const label timeIndex);

        //- Read and return given field
        template<class Type>
        tmp<Field<Type>> readFieldTemplate
        (
            const label timeIndex,
            const label fieldIndex
        ) const;


public:

    //- Runtime type information
    TypeName("boundaryData");


    // Constructors

        //- Construct from fileName
        explicit boundaryDataSurfaceReader
        (
            const fileName& fName,
            const word& pointsName = "points"
        );

        //- Construct from fileName with reader options
        boundaryDataSurfaceReader
        (
            const fileName& fName,
            const dictionary& options,
            const word& pointsName = "points"
        );


    //- Destructor
    virtual ~boundaryDataSurfaceReader() = default;


    // Static Member Functions

        //- Read points file
        static pointField readPoints
        (
            const Time& runTime,
            const fileName& baseDir,
            const word& pointsName = "points"
        );

        //- Read points file
        static pointField readPoints
        (
            const fileName& baseDir,
            const word& pointsName = "points"
        );

        //- Read and return given field
        template<class Type>
        static tmp<Field<Type>> readField
        (
            const Time& runTime,
            const fileName& baseDir,
            const instant& timeDir,
            const word& fieldName,
            Type& avg
        );

        //- Read and return given field
        template<class Type>
        static tmp<Field<Type>> readField
        (
            const fileName& baseDir,
            const instant& timeDir,
            const word& fieldName,
            Type& avg
        );


    // Member Functions

        //- Return a reference to the surface geometry
        virtual const meshedSurface& geometry(const label timeIndex);

        //- Return a list of the available times
        virtual instantList times() const;

        //- Return a list of the available fields at a given time
        virtual wordList fieldNames(const label timeIndex) const;

        //- Return a scalar field at a given time
        virtual tmp<Field<scalar>> field
        (
            const label timeIndex,
            const label fieldIndex,
            const scalar& refValue = pTraits<scalar>::zero
        ) const;

        //- Return a vector field at a given time
        virtual tmp<Field<vector>> field
        (
            const label timeIndex,
            const label fieldIndex,
            const vector& refValue = pTraits<vector>::zero
        ) const;

        //- Return a sphericalTensor field at a given time
        virtual tmp<Field<sphericalTensor>> field
        (
            const label timeIndex,
            const label fieldIndex,
            const sphericalTensor& refValue = pTraits<sphericalTensor>::zero
        ) const;

        //- Return a symmTensor field at a given time
        virtual tmp<Field<symmTensor>> field
        (
            const label timeIndex,
            const label fieldIndex,
            const symmTensor& refValue = pTraits<symmTensor>::zero
        ) const;

        //- Return a tensor field at a given time
        virtual tmp<Field<tensor>> field
        (
            const label timeIndex,
            const label fieldIndex,
            const tensor& refValue = pTraits<tensor>::zero
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "boundaryDataSurfaceReaderTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
