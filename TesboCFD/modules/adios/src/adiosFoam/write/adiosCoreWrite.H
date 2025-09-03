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
    Foam::adiosCoreWrite

Description
    Base for writing ADIOS files.

SourceFiles
    adiosCoreWrite.C
    adiosCoreWriteAttr.C
    adiosCoreWriteField.C
    adiosCoreWriteMesh.C
    adiosCoreWriteVar.C
    adiosCoreWriteTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_adiosCoreWrite_H
#define Foam_adiosCoreWrite_H

#include "adiosCore.H"
#include "adiosFoamCloudInfo.H"
#include "adiosFoamFieldInfo.H"
#include "dictionary.H"
#include "GeometricField.H"
#include <memory>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class cloud;
class polyMesh;
class TimeState;

namespace adiosFoam
{

/*---------------------------------------------------------------------------*\
                       Class adiosCoreWrite Declaration
\*---------------------------------------------------------------------------*/

class adiosCoreWrite
:
    public adiosCore
{
protected:

    // Private/Protected Data

        //- ADIOS instance
        std::unique_ptr<adios2::ADIOS> adiosPtr_;

        //- ADIOS output IO definitions - mirror of InquireIO("write")
        std::unique_ptr<adios2::IO> writeIOPtr_;

        //- ADIOS file write - mirror of ...
        std::unique_ptr<adios2::Engine> writeFilePtr_;


private:

        //- ADIOS input read method (default: BP3)
        word readMethod_;

        //- ADIOS output write method (default: MPI)
        word writeMethod_;

        //- ADIOS output write parameters (default: "")
        string writeParams_;

        // MPI parallel code variables
        //// MPI_Comm comm_;  // Pstream does not provide communicator


    // Private Member Functions

        //- Read the dictionary setup
        void read(const dictionary& dict);


protected:

    // Protected Member Functions

//         //- MPI communicator
//         MPI_Comm communicator() const
//         {
//             return comm_;
//         }


    // General functions

        //- Create and open output file
        bool open(const fileName& dataFile);

        bool beginWrite();
        bool endWrite();

        //- Close output file, flushing ADIOS buffer
        void close();

        //- Reset variable and attribute definitions
        void reset();


    // Member Functions

    // High-level methods

        //- Define and write base attributes
        void putBaseAttributes();

        //- Define and write attributes for OpenFOAM time information
        void putTimeAttributes(const TimeState& t);

        //- Define and write mesh patch attributes (name, type, etc)
        void putPatchAttributes(const polyMesh& mesh);

        //- Define/write mesh points
        void writeMeshPoints(const polyMesh& mesh);

        //- Define/write mesh face information
        void writeMeshFaces(const polyMesh& mesh);

        //- Define/write mesh proc-addressing
        void writeMeshAddressing(const polyMesh& mesh);

        //- Define/write list of mesh zones
        template<class ZoneMeshType>
        void writeMeshZones
        (
            const fileName& zonePath,
            const ZoneMeshType& zones
        );

        //- Define/write mesh points and face information
        void writeMesh(const polyMesh& mesh);


        //- Define/write adios variable for cloud fields
        //  \return number of fields handled
        label writeCloudRegistry(const cloud* cldPtr, objectRegistry& obr);


        //- Define/write adios variable for field object.
        //  \return true if object was handled, false otherwise
        bool writeFieldObject
        (
            const regIOobject& obj,
            const bool carp = false
        );

        //- Count non-processor boundaries
        template<class FvBoundaryMesh>
        static label nNonProcessor(const FvBoundaryMesh& bmesh);

        //- Define/write adios data for GeometricField.
        //  Writes internalField directly, boundaryField as a byte-stream
        //  \return True on success
        template<class Type, template<class> class PatchField, class GeoMesh>
        bool writeGeometricField
        (
            const GeometricField<Type, PatchField, GeoMesh>* field
        );

        //- Define/write adios data for GeometricField.
        //  Writes internalField directly, boundaryField as a byte-stream
        //  \return True on success
        template<class Type, template<class> class PatchField, class GeoMesh>
        bool writeGeometricField
        (
            const GeometricField<Type, PatchField, GeoMesh>& field
        );

        //- Define/write adios data for internal field
        //  \return True on success
        template<class Type, class GeoMesh>
        bool writeInternalField
        (
            const DimensionedField<Type, GeoMesh>* field
        );

        //- Define/write adios data for internal field.
        //  \return True on success
        template<class Type, class GeoMesh>
        bool writeInternalField
        (
            const DimensionedField<Type, GeoMesh>& field
        );

        //- Define/write adios data for field.
        //  \return True on success
        template<class Type, class GeoMesh>
        bool writeField
        (
            const fileName& varName,
            const word& local,
            const DimensionedField<Type, GeoMesh>& field
        );

        //- Define/write adios data for field.
        //  \return True on success
        template<class Type, class GeoMesh>
        bool writeField
        (
            const fileName& varName,
            const DimensionedField<Type, GeoMesh>& field
        );

        //- Write cloud fields (pre-collated in object registry)
        //  \return Names of fields written
        template<class Type>
        wordList writeCloudFields
        (
            const fileName& cloudVarName,
            const objectRegistry& obrTmp
        );


    // General ADIOS handling (lower-level methods)

    // Attributes

        //- Define and write a string attribute
        void putAttribute
        (
            const fileName& name,
            const std::string& value
        );

        //- Define and write an int attribute
        void putIntAttribute
        (
            const fileName& name,
            const int value
        );

        //- Define and write a double attribute
        void putDoubleAttribute
        (
            const fileName& name,
            const double value
        );

        //- Define and write an array attribute of integer values
        bool putListAttribute
        (
            const fileName& name,
            const UList<int>& list
        );

        //- Define and write an array attribute of double values
        bool putListAttribute
        (
            const fileName& name,
            const UList<double>& list
        );

        //- Define and write an array attribute of strings
        bool putListAttribute
        (
            const fileName& name,
            const UList<word>& list
        );

        //- Define and write an array attribute of strings
        bool putListAttribute
        (
            const fileName& name,
            const UList<string>& list
        );

        //- Define and write an "dimensions" attribute (an array of doubles).
        bool putDimensionsAttribute
        (
            const fileName& varName,
            const dimensionSet& dims
        );


    // Variables - lowest-level

        //- Generic define a variable, return adios variable
        template<class Type>
        adios2::Variable<typename pTraits<Type>::cmptType>
        definePutVariable
        (
            const fileName& name,
            const label local,
            const label global,
            const label offset
        );

        //- Define and write integer variable containing ONE element.
        //  To be stored as part of a global 1D array of this name.
        void putLabelVariable(const fileName& name, const label value);


        //- Define and write a variable from the list contents
        //  Uses local dimensions only.
        template<class Type>
        bool putListVariable(const fileName& name, const UList<Type>& list);


        //- Define and write an variable as a byte-stream.
        template<class Type>
        bool putStreamVariable(const fileName& name, const Type& obj);


public:

    //- No copy construct
    adiosCoreWrite(const adiosCoreWrite&) = delete;

    //- No copy assignment
    void operator=(const adiosCoreWrite&) = delete;


    // Constructors

        //- Construct for specified group name, using dictionary parameters
        explicit adiosCoreWrite(const dictionary& dict);


    //- Destructor
    virtual ~adiosCoreWrite();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


/*---------------------------------------------------------------------------*\
                       Class BoundaryWriter Declaration
\*---------------------------------------------------------------------------*/

template<class Type, template<class> class PatchField, class GeoMesh>
struct BoundaryWriter;


//- Helper class to stream a boundaryField as dictionary entries (for ADIOS),
//- calls internal writeEntries method
template<class Type, template<class> class PatchField, class GeoMesh>
struct BoundaryWriter
{
    using const_reference =
       const typename GeometricField<Type, PatchField, GeoMesh>::Boundary&;

    const_reference object;

    //- Constructor
    BoundaryWriter(const_reference& boundary)
    :
        object(boundary)
    {}


    // Ostream Operator

    inline friend Ostream& operator<<
    (
        Ostream& os,
        const BoundaryWriter& adapter
    )
    {
        adapter.object.writeEntries(os);
        return os;
    }
};

} // End namespace adiosFoam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "adiosCoreWriteTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
