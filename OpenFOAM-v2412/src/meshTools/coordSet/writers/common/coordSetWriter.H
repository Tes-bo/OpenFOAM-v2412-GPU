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
    Foam::coordSetWriter

Description
    Base class for writing coordSet(s) and tracks with fields.

    Example:
    \verbatim
        coordSet coords(...);

        // Construct writer of xmgr type
        autoPtr<coordSetWriter> writer(coordSetWriter::New("vtk"));

        writer.open(coords, path/name);

        writer.write("density", rho);
        writer.write("velocity", U);

    \endverbatim

SourceFiles
    coordSetWriterI.H
    coordSetWriter.C
    coordSetWriterBuffers.C
    coordSetWriterNew.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_coordSetWriter_H
#define Foam_coordSetWriter_H

#include "coordSet.H"
#include "typeInfo.H"
#include "vector.H"
#include "tensor.H"
#include "fileName.H"
#include "wordList.H"
#include "Field.H"
#include "DynamicList.H"
#include "PtrDynList.H"
#include "UPtrList.H"
#include "instant.H"
#include "InfoProxy.H"
#include "cartesianCS.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class coordSetWriter;
class Time;

Ostream& operator<<(Ostream& os, const InfoProxy<coordSetWriter>&);

/*---------------------------------------------------------------------------*\
                       Class coordSetWriter Declaration
\*---------------------------------------------------------------------------*/

class coordSetWriter
{
protected:

    // Protected Data

        //- Reference to coordinate set(s)
        UPtrList<const coordSet> coords_;

        //- Track times (eg, streamlines), one per coords_ entry
        List<scalarField> trackTimes_;

        //- The content is up-to-date?
        mutable bool upToDate_;

        //- Track if geometry has been written since the last open
        mutable bool wroteGeom_;

        /// //- Writing in parallel (via master)
        /// bool parallel_;

        //- Writer with buffering output
        mutable bool buffering_;

        //- Prefer tracks to points during single set writing
        bool useTracks_;

        //- Insert additional time sub-directory in the output path
        bool useTimeDir_;

        //- Additional output verbosity
        bool verbose_;

        //- The number of fields
        label nFields_;

        //- The current time value/name
        instant currTime_;

        //- The full output directory and file (coords) name
        fileName outputPath_;

        //- Output geometry scaling after rotate/translate
        scalar geometryScale_;

        //- The centre of rotation (untranslate, translate)
        point geometryCentre_;

        //- Local coordinate system transformation
        coordSystem::cartesian geometryTransform_;

        //- Field level to remove (on output)
        dictionary fieldLevel_;

        //- Field scaling (on output)
        dictionary fieldScale_;


    // Buffering

#undef  defineBufferMethod
#define defineBufferMethod(Type)                                              \
                                                                              \
        /* Names of Type fields */                                            \
        DynamicList<word> Type##Names_;                                       \
                                                                              \
        /* Values of Type fields */                                           \
        PtrDynList<Field<Type>> Type##Fields_;                                \
                                                                              \
        /* Add named Type field to buffering */                               \
        void appendField(const word& fieldName, const Field<Type>& vals)      \
        {                                                                     \
            Type##Names_.append(fieldName);                                   \
            Type##Fields_.append(vals.clone());                               \
        }

        defineBufferMethod(label);
        defineBufferMethod(scalar);
        defineBufferMethod(vector);
        defineBufferMethod(sphericalTensor);
        defineBufferMethod(symmTensor);
        defineBufferMethod(tensor);

#undef defineBufferMethod

    // Protected Member Functions

    // Buffering

        //- Write line contents (eg, buffered)
        static void writeLine(Ostream&, const UList<word>&, const char* sep);

        //- Write line contents (eg, buffered)
        static void writeLine(Ostream&, const UList<scalar>&, const char* sep);

        //- Clear out buffering
        void clearBuffers();

        //- The number of buffer data columns, after splitting into components
        label nDataColumns() const;

        //- Write buffered data
        virtual bool writeBuffered();

        //- Write buffered data
        void writeBufferContents
        (
            Ostream& os,
            const coordSet& coords,
            const char* sep
        ) const;

        //- Get buffered data line (components)
        void getBufferLine
        (
            DynamicList<scalar>& buf,
            const coordSet& coords,
            const label pointi
        ) const;


    // File Operations

        //- Get expected (characteristic) output file name - information only
        fileName getExpectedPath(const word& fileExt = word::null) const;

        //- Get field-prefixed output file name.
        //  Eg, dir/U_name.raw
        fileName getFieldPrefixedPath
        (
            const word& fieldName,
            const word& fileExt = word::null
        ) const;

        //- Verify that the outputPath_ has been set or FatalError
        void checkOpen() const;

        //- Perform any merging if not already upToDate (parallel)
        //- or simply mark as being up-to-date
        virtual bool merge() const;


    // Helpers

        template<class Type>
        tmp<Field<Type>> adjustFieldTemplate
        (
            const word& fieldName,
            const tmp<Field<Type>>& tfield
        ) const;

        //- Repackage field into a UPtrList
        template<class Type>
        static UPtrList<const Field<Type>>
        repackageFields
        (
            const Field<Type>& field
        );

        //- Repackage multiple fields into a UPtrList
        template<class Type>
        static UPtrList<const Field<Type>>
        repackageFields
        (
            const UList<Field<Type>>& fieldValues
        );

        //- Write coordinates and values
        template<class Type>
        static void writeTable
        (
            Ostream& os,
            const coordSet& coords,
            const UList<Type>& values,
            const char* sep
        );


    // Normal write templates

        //- Dummy templated write operation
        template<class Type>
        fileName writeTemplate
        (
            const word& fieldName,          //!< Name of field
            const Field<Type>& values       //!< Local field values to write
        )
        {
            /// if (!wroteGeom_)
            /// {
            ///     return this->write();
            /// }
            return fileName::null;
        }

        //- Dummy templated write operation. Multiple tracks
        template<class Type>
        fileName writeTemplate
        (
            const word& fieldName,
            const List<Field<Type>>& fieldValues
        )
        {
            return fileName::null;
        }


        //- No copy construct
        coordSetWriter(const coordSetWriter&) = delete;

        //- No copy assignment
        void operator=(const coordSetWriter&) = delete;


public:

    //- Runtime type information
    TypeName("coordSetWriter");

        // Declare run-time constructor selection table
        declareRunTimeSelectionTable
        (
            autoPtr,
            coordSetWriter,
            word,
            (),
            ()
        );

        declareRunTimeSelectionTable
        (
            autoPtr,
            coordSetWriter,
            wordDict,
            (
                const dictionary& writeOptions
            ),
            (writeOptions)
        );


    // Helpers

        //- Same as fileFormats::getFormatOptions
        static dictionary formatOptions
        (
            const dictionary& dict,
            const word& formatName,
            const word& entryName = "formatOptions"
        );

        //- Same as fileFormats::getFormatOptions
        static dictionary formatOptions
        (
            const dictionary& dict,
            const dictionary& setDict,
            const word& formatName,
            const word& entryName = "formatOptions"
        );


    // Selectors

        //- True if New is likely to succeed for this writeType
        static bool supportedType(const word& writeType);

        //- Return a reference to the selected writer
        static autoPtr<coordSetWriter> New(const word& writeFormat);

        //- Return a reference to the selected writer
        //  Select with extra write option
        static autoPtr<coordSetWriter> New
        (
            const word& writeFormat,
            const dictionary& writeOptions
        );


    // Constructors

        //- Default construct
        coordSetWriter();

        //- Default construct with specified options
        explicit coordSetWriter(const dictionary& options);


    //- Destructor. Calls close()
    virtual ~coordSetWriter();


    // Member Functions

    // Helpers

        //- Name suffix based on fieldName (underscore separator)
        static word suffix
        (
            const word& fldName,
            const word& fileExt = word::null
        );

        //- Name suffix based on fieldNames (underscore separator)
        static word suffix
        (
            const wordList& fieldNames,
            const word& fileExt = word::null
        );


    // Capability

        //- True if the format uses internal buffering (eg, column output)
        virtual bool buffering() const;

        //- Turn internal buffering on/off (only if supported by the writer)
        virtual bool buffering(const bool on);

        //- The writer is enabled. If the writer is not enabled, it may be
        //- possible for the caller to skip various preparatory operations.
        //  This method is primarily useful for the null writer
        virtual bool enabled() const
        {
            return true;
        }


    // Bookkeeping

        //- Does the writer need an update (eg, lagging behind other changes)
        virtual bool needsUpdate() const;

        //- Geometry or fields written since the last open?
        virtual bool wroteData() const;

        //- Mark that content changed and the writer will need an update,
        //- and set nFields = 0.
        //  May also free up unneeded data.
        //  Return false if it was previously already expired.
        virtual bool expire();

        //- Close any open output, remove coordSet associations and
        //- expire the writer.
        virtual void clear();


    // Content Association

        //- Set coordinates, can also be nullptr
        virtual void setCoordinates(const coordSet* coords);

        //- Set coordinates
        virtual void setCoordinates(const coordSet& coords);

        //- Set track coordinates
        virtual void setTracks(const UPtrList<coordSet>& tracks);

        //- Set track times
        virtual void setTrackTimes(const UList<scalarField>& times);


    // Queries, Access

        //- The number of associated points (local processor)
        label numPoints() const;

        //- The number of coordinate tracks
        label numTracks() const;

        //- Writer is associated with content
        bool hasCoords() const;

        //- Writer is not associated with content
        bool empty() const;

        //- Test if outputPath has been set
        inline bool is_open() const noexcept;

        //- The number of expected output fields.
        //  Currently only used by the legacy VTK format.
        inline label nFields() const noexcept;

        //- Set the number of expected output fields
        //  Currently only used by the legacy VTK format.
        //  \return old value
        inline label nFields(const label n) noexcept;

        //- Prefer tracks to points during single set writing
        inline bool useTracks() const noexcept;

        //- Enable/disable tracks preference
        //  \return old value
        inline bool useTracks(const bool on) noexcept;

        //- Should a time directory be spliced into the output path?
        inline bool useTimeDir() const noexcept;

        //- Enable/disable use of spliced output path
        //  \return old value
        inline bool useTimeDir(const bool on) noexcept;

        //- Get output verbosity
        inline bool verbose() const noexcept;

        //- Enable/disable verbose output
        //  \return old value
        inline bool verbose(bool on) noexcept;


    // Time

        //- True if there is a known time
        inline bool hasTime() const;

        //- The current time value/name
        inline const word& timeName() const;

        //- The current time value/name
        inline scalar timeValue() const;


        //- Set the current time
        void setTime(const instant& inst);

        //- Set current time from timeValue, auto generating the name
        void setTime(scalar timeValue);

        //- Set current time from timeValue and timeName
        void setTime(scalar timeValue, const word& timeName);

        //- Clear the current time
        void unsetTime();


        //- Begin a time-step
        virtual void beginTime(const Time& t);

        //- Begin a time-step
        virtual void beginTime(const instant& inst);

        //- End a time-step
        virtual void endTime();


    // Output

        //- Expected (characteristic) output file name - information only.
        //- Return empty when is_open() is false.
        virtual fileName path() const = 0;

        //- Write separate geometry to file.
        /// TBD:  virtual fileName write() = 0;

        //- Open for output on specified path, using existing content
        virtual void open(const fileName& outputPath);

        //- Open from components
        virtual void open
        (
            const coordSet& coords,
            const fileName& outputPath
        );

        //- Open from components
        virtual void open
        (
            const UPtrList<coordSet>& tracks,
            const fileName& outputPath
        );

        //- Finish output, performing any necessary cleanup
        //  Optional force disassociation with any coordSet(s)
        virtual void close(bool force = false);


    // Other

        //- Return info proxy,
        //- used to print information to a stream
        virtual InfoProxy<coordSetWriter> info() const
        {
            return *this;
        }

        //- Output info proxy
        friend Ostream& operator<<
        (
            Ostream& os,
            const InfoProxy<coordSetWriter>& iproxy
        );


    // Write methods

#undef  declareCoordSetWriterWriteMethod
#define declareCoordSetWriterWriteMethod(Type)                                 \
                                                                               \
        /*! \brief Write field of Type (point data) */                         \
        virtual fileName write                                                 \
        (                                                                      \
            const word& fieldName,          /*!< Name of field */              \
            const Field<Type>& field        /*!< Field values */               \
        ) = 0;                                                                 \
                                                                               \
        /*! \brief Write track fields of Type (point data) */                  \
        virtual fileName write                                                 \
        (                                                                      \
            const word& fieldName,                /*!< Name of field */        \
            const List<Field<Type>>& fieldValues  /*!< Per-track values */     \
        ) = 0;

        declareCoordSetWriterWriteMethod(label);
        declareCoordSetWriterWriteMethod(scalar);
        declareCoordSetWriterWriteMethod(vector);
        declareCoordSetWriterWriteMethod(sphericalTensor);
        declareCoordSetWriterWriteMethod(symmTensor);
        declareCoordSetWriterWriteMethod(tensor);


#undef  declareCoordSetWriterWriteMethod
#define declareCoordSetWriterWriteMethod(Type)                                 \
                                                                               \
        /*! \brief Write field of Type (point data) */                         \
        virtual fileName write                                                 \
        (                                                                      \
            const word& fieldName,          /*!< Name of field */              \
            const Field<Type>& values       /*!< Field values */               \
        );  /* override */                                                     \
                                                                               \
        /*! \brief Write track fields of Type (point data) */                  \
        virtual fileName write                                                 \
        (                                                                      \
            const word& fieldName,                /*!< Name of field */        \
            const List<Field<Type>>& fieldValues  /*!< Per-track values */     \
        ); /* override */
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "coordSetWriterI.H"

#ifdef NoRepository
    #include "coordSetWriterTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
