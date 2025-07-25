/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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
    Foam::adiosFoam::adiosReader

Description
    Base for reading ADIOS files.

SourceFiles
    adiosReader.C
    adiosReaderCloud.C
    adiosReaderField.C
    adiosReaderTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_adiosReader_H
#define Foam_adiosReader_H

#include "adiosCore.H"
#include "adiosTime.H"
#include "adiosFoamCloudInfo.H"
#include "adiosFoamFieldInfo.H"

#include "fileName.H"
#include "fvMesh.H"
#include "HashSet.H"
#include "Ostream.H"
#include <memory>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace adiosFoam
{

/*---------------------------------------------------------------------------*\
                         Class adiosReader Declaration
\*---------------------------------------------------------------------------*/

class adiosReader
{
private:

    // Private Data

        //- ADIOS instance
        adios2::ADIOS& adios_;

        //- ADIOS input IO definitions - mirror of InquireIO("read")
        mutable std::unique_ptr<adios2::IO> readIOPtr_;

        //- ADIOS file read - mirror of ...
        mutable std::unique_ptr<adios2::Engine> readFilePtr_;


    // Data populated by the scan method

        //- Available attributes
        std::map<std::string, adios2::Params> availableAttr_;

        //- Available variables
        std::map<std::string, adios2::Params> availableVars_;

        //- Label-size
        unsigned short labelByteSize_;

        //- Scalar-size
        unsigned short scalarByteSize_;

        //- The region names
        wordList regionNames_;

        //- The clouds information, per region
        HashTable<HashSet<cloudInfo>> clouds_;

        //- The fields information, per region
        HashTable<HashSet<fieldInfo>> fields_;


    // Private Member Functions

        //- Reset scanned values, does not open/close file
        void reset();

        //- scan opened file for variables/attributes
        void scan(bool verbose=false);


        //- Read varName of specified IOField Type and save into registry
        template<class Type>
        bool readIntoRegistry
        (
            objectRegistry& obr,
            const fileName& varName,
            word clsName = ""
        );


        //- Read into list (input type specified), with type conversion
        template<class Type, class SourceType>
        label readCopyList
        (
            List<Type>& list,
            const fileName& varName,
            const bool allowResize = false
        );

        //- Read into list (input type from Type), with type conversion
        template<class Type>
        label readList
        (
            List<Type>& list,
            const fileName& varName,
            const bool allowResize = false
        );

        //- Read into list (input type by name), with type conversion
        template<class Type>
        label readList
        (
            List<Type>& list,
            const fileName& varName,
            const std::string& varTypeName,
            const bool allowResize = false
        );

        template<class Type, class SourceType>
        bool readAttributeImpl
        (
            const string& name,
            Type& value,
            bool mandatory
        );

        template<class Type, class SourceType>
        bool readVariableImpl
        (
            const string& name,
            Type& value,
            bool mandatory
        );

        template<class Type, class SourceType>
        bool readListAttributeImpl
        (
            const string& name,
            List<Type>& value,
            bool mandatory
        );


        //- Disallow default bitwise copy construct
        adiosReader(const adiosReader&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const adiosReader&) = delete;


public:

    // Constructors

        //- Construct with reference to ADIOS backend
        explicit adiosReader(adios2::ADIOS& adios);

        //- Open specified filename and scan variables
        adiosReader(adios2::ADIOS& adios, const fileName& dataFile);


    //- Destructor
    ~adiosReader();


    // Member Functions

        //- ADIOS output IO definitions - mirror of InquireIO("read")
        adios2::IO& readerIO() { return *readIOPtr_; }


        //- Open specified filename, select local process and scan variables
        bool open(const fileName& dataFile);

        //- Close file and reset
        void close();

        //- True if the file has been opened
        bool isGood() const;


        //- The max buffer size for the process local variables
        inline size_t sizeOf() const;


        //- Cloud information (name, class) for region, or nullptr
        inline const HashSet<adiosFoam::cloudInfo>*
        getCloudInfo(const word& regName) const;

        //- Field information (name, class) for region, or nullptr
        inline const HashSet<adiosFoam::fieldInfo>*
        getFieldInfo(const word& regName) const;


        //- Extract cloud information as a registry
        std::unique_ptr<objectRegistry>
        getCloud(const word& regName, const word& cloudName);

        //- Extract cloud information as a registry
        std::unique_ptr<objectRegistry>
        getCloud(const adiosFoam::cloudInfo& cldInfo);


        //- Return string representation of attribute (if it exists)
        inline std::string inquireAttributeType(const std::string& name) const;

        //- Return string representation of variable (if it exists)
        inline std::string inquireVariableType(const std::string& name) const;


        //- True if the specified attribute exists, optionally with type-check
        template<class T = void>
        inline bool hasAttribute(const string& varName) const;

        //- True if the specified variable exists, optionally with type-check
        template<class T = void>
        inline bool hasVariable(const string& varName) const;


        //- Read integer attribute
        inline bool readIntAttribute
        (
            const string& attrName,
            label& value,
            bool mandatory = true
        );

        //- Read double attribute
        inline bool readScalarAttribute
        (
            const string& attrName,
            double& value,
            bool mandatory = true
        );

        //- Read string attribute
        inline bool readStringAttribute
        (
            const string& attrName,
            string& value,
            bool mandatory = true
        );

        //- Read integer attribute if it exists
        inline bool readIntAttributeIfPresent
        (
            const string& attrName,
            label& value
        );

        //- Read double attribute if it exists
        inline bool readScalarAttributeIfPresent
        (
            const string& attrName,
            double& value
        );

        //- Read string attribute if it exists
        inline bool readStringAttributeIfPresent
        (
            const string& attrName,
            string& value
        );

        //- Read int-list attribute if it exists
        inline bool readIntListAttributeIfPresent
        (
            const string& attrName,
            List<label>& lst
        );

        //- Read double-list attribute if it exists
        inline bool readScalarListAttributeIfPresent
        (
            const string& attrName,
            List<double>& lst
        );

        //- Read string-list attribute if it exists
        template<class StringType>
        bool readStringListAttributeIfPresent
        (
            const string& attrName,
            List<StringType>& lst
        );


        //- Get a mandatory integer attribute
        inline label getIntAttribute(const string& name);

        //- Get a mandatory scalar attribute
        inline double getDoubleAttribute(const string& name);

        //- Get a mandatory string attribute
        inline string getStringAttribute(const string& name);


        //- Get a mandatory integer list attribute
        inline List<label> getIntListAttribute(const string& name);

        //- Get a mandatory scalar list attribute
        inline List<double> getDoubleListAttribute(const string& name);

        //- Get a mandatory integer attribute
        template<class StringType>
        inline List<StringType> getStringListAttribute(const string& name);


        //- Read integer variable (single-value) if it exists
        inline bool readIntVariableIfPresent
        (
            const string& varName,
            label& value
        );

        //- Read scalar variable (single-value) if it exists
        inline bool readScalarVariableIfPresent
        (
            const string& varName,
            scalar& value
        );

        //- Get a mandatory integer variable (single-value)
        inline label getIntVariable(const string& name);

        //- Get a mandatory scalar variable (single-value)
        inline scalar getScalarVariable(const string& name);


#if 0
        //- Retrieve the variable directly into the memory specified
        bool getVariable(const string& varName, void* data) const;


        //- Retrieve variable into buffer, using the common adiosCore buffer.
        //  On success the buffer is resized according to the data.
        //  On failure the buffer is resized to zero.
        size_t getBuffered(const string& varName) const;


        //- Retrieve variable into buffer.
        //  On success the buffer is resized according to the data.
        //  On failure the buffer is resized to zero.
        template<class T>
        size_t getBuffered
        (
            const string& varName,
            DynamicList<T>& buffer
        ) const;
#endif

    // Functions for reading (fields, clouds)

        //- Read data for a given volume field type.
        //  Return the number of elements read.
        bool readFieldObject
        (
            regIOobject& obj,
            const adiosFoam::fieldInfo& src,
            const bool carp = false
        );

        //- Read GeometricField.
        //  Read internalField directly, boundaryField as a byte-stream
        //  Return True on success
        template<class Type, template<class> class PatchField, class GeoMesh>
        bool readGeometricField
        (
            GeometricField<Type, PatchField, GeoMesh>* pfield,
            const fileName& varName
        );

        //- Read a dimensioned field (eg, internalField) directly.
        //  Return True on success
        template<class Type, class GeoMesh>
        bool readField
        (
            DimensionedField<Type, GeoMesh>* pfield,
            const fileName& varName
        );

        //- Read GeometricField.
        //  Read internalField directly, boundaryField as a byte-stream
        //  Return the number of elements read.
        template<class Type, template<class> class PatchField, class GeoMesh>
        label readGeometricField
        (
            GeometricField<Type, PatchField, GeoMesh>& field,
            const fileName& varName
        );

        //- Read a dimensioned field (eg, internalField) directly.
        //  Return the number of elements read.
        template<class Type, class GeoMesh>
        label readField
        (
            DimensionedField<Type, GeoMesh>& field,
            const fileName& varName
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace adiosFoam
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "adiosReaderI.H"

#ifdef NoRepository
    #include "adiosReaderTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
