/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2024 OpenCFD Ltd.
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
    Foam::ensightSurfaceReader

Description
    Ensight format surface reader

    \verbatim
    readOptions
    {
        ensight
        {
            debug       false;
            masterOnly  false;
        }
    }
    \endverbatim

    Format options for ensight:
    \table
        Property | Description                              | Required | Default
        debug    | Add debug flag                           | no  | false
        masterOnly | Read files on master and broadcast values | no  | true
    \endtable

SourceFiles
    ensightSurfaceReader.C
    ensightSurfaceReaderTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_ensightSurfaceReader_H
#define Foam_ensightSurfaceReader_H

#include "surfaceReader.H"
#include "ensightCase.H"
#include "ensightFaces.H"
#include "ensightReadFile.H"
#include "Pair.H"
#include "Tuple2.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class ensightSurfaceReader Declaration
\*---------------------------------------------------------------------------*/

class ensightSurfaceReader
:
    public surfaceReader
{
protected:

    // Data Types

        //- Handling of node/element id types (off, assign, ignore, given)
        enum idTypes : unsigned char
        {
            NONE = 0,       //!< "off", "assign"
            IGNORE = 1,     //!< Read but "ignore"
            GIVEN = 2       //!< Use "given" values (not supported)
        };


        //- Tuple of face type (tria3, quad4, nsided) and count
        typedef Tuple2<ensightFaces::elemType, label> faceInfoTuple;


    // Protected Data

        //- Read on master and broadcast (in parallel)
        bool masterOnly_;

        //- Format flag
        IOstreamOption::streamFormat readFormat_;

        //- Base directory
        fileName baseDir_;

        //- The timeset/fileset (if any) associated with the mesh
        labelPair meshTimeset_;

        //- Name of mesh file, including any subdirectory
        fileName meshFileName_;

        //- The timeset/fileset (if any) associated with fields
        List<labelPair> fieldTimesets_;

        //- Field names
        List<word> fieldNames_;

        //- Field file names
        List<string> fieldFileNames_;

        //- Start time index
        label timeStartIndex_;

        //- Time increment
        label timeIncrement_;

        //- Numbers for files
        labelList fileNumbers_;

        //- Times
        instantList timeValues_;

        //- Pointer to the surface
        autoPtr<meshedSurface> surfPtr_;

        //- List of face-type/count tuples
        List<faceInfoTuple> faceTypeInfo_;


    // Static Member Functions

        //- Helper function to read an ascii line from file,
        //- skipping blank lines and comments.
        //  \return True if reading was successful
        static bool readLine(ISstream& is, std::string& line);

        //- Check a section header
        static void checkSection
        (
            const word& expected,
            const string& buffer,
            const ISstream& is  // For errors
        );

        //- Read and check a section header
        static void debugSection(const word& expected, ISstream& is);

        //- Helper function to return Type from string
        template<class Type>
        static void readFrom(const std::string& buffer, Type& value);


    // Protected Member Functions

        //- Read (and discard) geometry file header.
        //  \return information about node/element id handling
        Pair<idTypes> readGeometryHeader(ensightReadFile& is) const;

        //- Read the case file
        void readCase(ISstream& is);

        //- Read and return surface geometry. Updates faceTypeInfo_
        meshedSurface readGeometry
        (
            const fileName& geometryFile,
            //! Optional index for transient single-file format
            const label timeIndex = 0
        );

        //- Helper function to return a field
        template<class Type>
        tmp<Field<Type>> readField
        (
            const fileName& dataFile,
            const word& fieldName,
            const label timeIndex = 0
        ) const;

        //- Helper function to return a field
        template<class Type>
        tmp<Field<Type>> readField
        (
            const label timeIndex,
            const label fieldIndex
        ) const;


public:

    // Generated Methods

        //- No copy construct
        ensightSurfaceReader(const ensightSurfaceReader&) = delete;

        //- No copy assignment
        void operator=(const ensightSurfaceReader&) = delete;


    //- Runtime type information
    TypeName("ensight");


    // Constructors

        //- Construct from fileName, with reader options
        explicit ensightSurfaceReader
        (
            const fileName& fName,
            const dictionary& options = dictionary()
        );


    //- Destructor
    virtual ~ensightSurfaceReader() = default;


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
    #include "ensightSurfaceReaderTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
