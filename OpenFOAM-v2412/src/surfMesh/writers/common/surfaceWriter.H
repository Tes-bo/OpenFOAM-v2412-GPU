/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2012 OpenFOAM Foundation
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

Namespace
    Foam::surfaceWriters

Description
    Namespace for surface writers

Class
    Foam::surfaceWriter

Description
    Base class for surface writers.

    The surfaceWriter interface is rather large since we need a writer that
    can either be initially defined without a surface association
    and have that added at a later stage, or be defined with a surface
    association.

    \verbatim
    formatOptions
    {
        default
        {
            verbose     true;
            commsType   scheduled;

            fieldLevel
            {
                "p.*"   1e5;    // Absolute -> gauge [Pa]
                T       273.15; // [K] -> [C]
                U       #eval{ 10/sqrt(3) };  // Uniform magU=10
            }
        }

        someFormat // Eg, ensight, vtk, etc
        {
            fieldScale
            {
                "p.*"   0.01;   // [Pa] -> [mbar]
            }

            scale   1000;  // [m] -> [mm]
            transform
            {
                origin  (0 0 0);
                rotationCentre  (0 0 0);
                rotation axisAngle;
                axis    (1 0 0);
                angle   45;
            }
        }
    }
    \endverbatim

    Format options:
    \table
        Property    | Description                           | Reqd | Default
        verbose     | Additional output verbosity           | no  | no
        commsType   | Communication type                    | no  | scheduled
        scale       | Output geometry scaling               | no  | 1
        transform   | Output coordinate transform           | no  |
        fieldLevel  | Subtract field level before scaling   | no  | empty dict
        fieldScale  | Output field scaling                  | no  | empty dict
    \endtable

Note
    The \c transform sub-dictionary also supports a \c rotationCentre
    keyword which applies \em untranslate by that amount prior to the rotation,
    and subsequently followed by a \em translate.

    For surface formats that require geometry in a separate file,
    it is the responsibility of the implementation (not the caller)
    to ensure that this occurs.

SourceFiles
    surfaceWriter.C
    surfaceWriterI.H

\*---------------------------------------------------------------------------*/

#ifndef Foam_surfaceWriter_H
#define Foam_surfaceWriter_H

#include "typeInfo.H"
#include "autoPtr.H"
#include "tmp.H"
#include "cartesianCS.H"
#include "Field.H"
#include "fileName.H"
#include "instant.H"
#include "mergedSurf.H"
#include "meshedSurfRef.H"
#include "UPstream.H"
#include "InfoProxy.H"
#include "runTimeSelectionTables.H"

// Deprecated methods
// #define Foam_surfaceWriter_directAccess

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class Time;
class surfaceWriter;

Ostream& operator<<(Ostream& os, const InfoProxy<surfaceWriter>& ip);


/*---------------------------------------------------------------------------*\
                        Class surfaceWriter Declaration
\*---------------------------------------------------------------------------*/

class surfaceWriter
{
protected:

    // Private Data

        //- Reference to surface or surface components
        meshedSurfRef surf_;

        //- Surface after merging (parallel)
        mutable mergedSurf mergedSurf_;

        //- The surface after point coordinate transforms and scaling
        mutable meshedSurfRef adjustedSurf_;

        //- Dimension for merging
        scalar mergeDim_;

        //- Output geometry scaling after rotate/translate
        scalar geometryScale_;

        //- The centre of rotation (untranslate, translate)
        point geometryCentre_;

        //- Local coordinate system transformation
        coordSystem::cartesian geometryTransform_;


protected:

    // Protected Data

        //- The topology/surface is up-to-date?
        mutable bool upToDate_;

        //- Track if geometry has been written since the last open
        mutable bool wroteGeom_;

        //- Writing in parallel (via master)
        bool parallel_;

        //- Insert additional time sub-directory in the output path
        bool useTimeDir_;

        //- Is point vs cell data
        bool isPointData_;

        //- Additional output verbosity
        bool verbose_;

        //- Communication type (for field merging)
        UPstream::commsTypes commType_;

        //- The number of fields
        label nFields_;

        //- The current time value/name
        instant currTime_;

        //- The full output directory and file (surface) name
        fileName outputPath_;

        //- Field level to remove (on output)
        dictionary fieldLevel_;

        //- Field scaling (on output)
        dictionary fieldScale_;


    // Protected Member Functions

        //- Verify that the outputPath_ has been set or FatalError
        void checkOpen() const;

        //- Merge surfaces if they are not already upToDate (parallel)
        //- or simply mark the surface as being up-to-date
        virtual bool merge() const;

        //- Merge surfaces (if not upToDate) and return merged (parallel)
        //- or regular surface (non-parallel)
        const meshedSurf& surface() const;

        //- Merge surfaces (if not upToDate) and return merged (parallel)
        //- or regular surface (non-parallel)
        //- and apply any coordinate system changes and/or output scaling.
        const meshedSurfRef& adjustSurface() const;

        //- Gather (merge) fields with renumbering and shrinking for point data
        template<class Type>
        tmp<Field<Type>> mergeFieldTemplate(const Field<Type>& fld) const;

        //- Apply refLevel and fieldScaling
        template<class Type>
        tmp<Field<Type>> adjustFieldTemplate
        (
            const word& fieldName,
            const tmp<Field<Type>>& tfield
        ) const;

#undef  declareSurfaceFieldMethods
#define declareSurfaceFieldMethods(Type)                                      \
                                                                              \
        tmp<Field<Type>> mergeField(const Field<Type>& fld) const;            \
                                                                              \
        tmp<Field<Type>> adjustField                                          \
        (                                                                     \
            const word& fieldName,                                            \
            const tmp<Field<Type>>& tfield                                    \
        ) const;

        declareSurfaceFieldMethods(label);
        declareSurfaceFieldMethods(scalar);
        declareSurfaceFieldMethods(vector);
        declareSurfaceFieldMethods(sphericalTensor);
        declareSurfaceFieldMethods(symmTensor);
        declareSurfaceFieldMethods(tensor);

        #undef declareSurfaceFieldMethods

        //- Dummy templated write operation
        template<class Type>
        fileName writeTemplate
        (
            const word& fieldName,          //!< Name of field
            const Field<Type>& localValues  //!< Local field values to write
        )
        {
            if (!wroteGeom_)
            {
                return this->write();
            }
            return fileName::null;
        }

public:

    // Public Data

        //- The default merge dimension (1e-8)
        static scalar defaultMergeDim;


    //- Runtime type information
    TypeName("surfaceWriter");


        // Declare run-time constructor selection table
        declareRunTimeSelectionTable
        (
            autoPtr,
            surfaceWriter,
            word,
            (),
            ()
        );

        declareRunTimeSelectionTable
        (
            autoPtr,
            surfaceWriter,
            wordDict,
            (
                const dictionary& writeOpts
            ),
            (writeOpts)
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
            const dictionary& surfDict,
            const word& formatName,
            const word& entryName = "formatOptions"
        );


    // Factory Methods / Selectors

        //- True if New is likely to succeed for this writeType
        static bool supportedType(const word& writeType);

        //- Optional select construct surfaceWriter.
        //  Return nullptr if the specified type is not supported.
        static autoPtr<surfaceWriter> TryNew(const word& writeType);

        //- Optional select construct surfaceWriter with extra write options.
        //  Return nullptr if the specified type is not supported.
        static autoPtr<surfaceWriter> TryNew
        (
            const word& writeType,
            const dictionary& writeOptions
        );

        //- Select construct a surfaceWriter
        static autoPtr<surfaceWriter> New(const word& writeType);

        //- Select construct a surfaceWriter with extra write options.
        static autoPtr<surfaceWriter> New
        (
            const word& writeType,
            const dictionary& writeOptions
        );


    // Constructors

        //- Default construct
        surfaceWriter();

        //- Default construct with specified options
        explicit surfaceWriter(const dictionary& options);


    //- Destructor. Calls close()
    virtual ~surfaceWriter();


    // Member Functions

    // Capability

        //- The writer is enabled. If the writer is not enabled, it may be
        //- possible for the caller to skip various preparatory operations.
        //  This method is primarily useful for the null writer
        virtual bool enabled() const
        {
            return true;
        }

        //- True if the surface format requires geometry in a separate file.
        virtual bool separateGeometry() const
        {
            return false;
        }

        //- True if the writer format uses faceIds as part of its output.
        //  Element ids are used by various CAE formats
        //  (abaqus, nastran, starcd, ...)
        virtual bool usesFaceIds() const
        {
            return false;
        }


    // Bookkeeping

        //- Does the writer need an update (eg, lagging behind surface changes)
        virtual bool needsUpdate() const;

        //- Geometry or fields written since the last open?
        virtual bool wroteData() const;

        //- Mark that surface changed and the writer will need an update,
        //- and set nFields = 0.
        //  May also free up unneeded data.
        //  Return false if it was previously already expired.
        virtual bool expire();

        //- Close any open output, remove association with a surface and
        //- expire the writer. The parallel flag remains untouched.
        virtual void clear();


    // Surface association

        //- Change association with a surface, expire the writer
        //- with defined parallel/serial treatment
        virtual void setSurface
        (
            const meshedSurf& surf,
            bool parallel
        );

        //- Change association with a surface, expire the writer
        //- with defined parallel/serial treatment
        virtual void setSurface
        (
            const pointField& points,
            const faceList& faces,
            bool parallel
        );

        //- Change association with a surface, expire the writer
        //- with the current parallel/serial treatment
        virtual void setSurface
        (
            const meshedSurf& surf
        );

        //- Change association with a surface, expire the writer
        //- with the current parallel/serial treatment
        virtual void setSurface
        (
            const pointField& points,
            const faceList& faces
        );


    // Queries, Access

        //- Test if outputPath has been set
        inline bool is_open() const noexcept;

        //- Writer is associated with a surface
        bool hasSurface() const;

        //- The surface to write is empty if the global number of faces is zero
        bool empty() const;

        //- The global number of faces for the associated surface
        label size() const;

        //- The number of expected output fields.
        //  Currently only used by the legacy VTK format.
        inline label nFields() const noexcept;

        //- Set the number of expected output fields
        //  Currently only used by the legacy VTK format.
        //  \return old value
        inline label nFields(const label n) noexcept;

        //- Are the field data to be treated as point data?
        inline bool isPointData() const noexcept;

        //- Set handling of field data to face/point data
        //  \return old value
        inline bool isPointData(const bool on) noexcept;

        //- Should a time directory be spliced into the output path?
        inline bool useTimeDir() const noexcept;

        //- Enable/disable use of spliced output path
        //  \return old value
        inline bool useTimeDir(const bool on) noexcept;

        //- Get output verbosity
        inline bool verbose() const noexcept;

        //- Enable/disable verbose output
        //  \return old value
        inline bool verbose(const bool on) noexcept;

        //- The current value of the point merge dimension (metre)
        inline scalar mergeDim() const noexcept;

        //- Change the point merge dimension (metre)
        //  \return old value
        inline scalar mergeDim(const scalar dist) noexcept;

        //- The current value of the geometry scaling
        inline scalar scale() const noexcept;

        //- Change the geometry scaling
        //  \return old value
        inline scalar scale(const scalar factor) noexcept;

        //- The current (cartesian) coordinate system transformation
        inline const coordSystem::cartesian& transform() const noexcept;


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

        //- Open for output on specified path, using existing surface
        virtual void open(const fileName& outputPath);

        //- Open from components
        virtual void open
        (
            const pointField& points,
            const faceList& faces,
            const fileName& outputPath,
            bool parallel
        );

        //- Open from components
        virtual void open
        (
            const meshedSurf& surf,
            const fileName& outputPath,
            bool parallel
        );

        //- Open from components, with the current parallel/serial treatment
        virtual void open
        (
            const pointField& points,
            const faceList& faces,
            const fileName& outputPath
        );

        //- Open from components, with the current parallel/serial treatment
        virtual void open
        (
            const meshedSurf& surf,
            const fileName& outputPath
        );

        //- Finish output, performing any necessary cleanup
        virtual void close();


    // Write

        //- Write separate surface geometry to file.
        virtual fileName write() = 0;


#undef  declareSurfaceWriterWriteMethod
#define declareSurfaceWriterWriteMethod(Type)                                  \
        /*! \brief Write field of Type (per face or vertex) */                 \
        virtual fileName write                                                 \
        (                                                                      \
            const word& fieldName,          /*!< Name of field */              \
            const Field<Type>& values       /*!< Field values to write */      \
        ) = 0

        declareSurfaceWriterWriteMethod(label);
        declareSurfaceWriterWriteMethod(scalar);
        declareSurfaceWriterWriteMethod(vector);
        declareSurfaceWriterWriteMethod(sphericalTensor);
        declareSurfaceWriterWriteMethod(symmTensor);
        declareSurfaceWriterWriteMethod(tensor);


#undef  declareSurfaceWriterWriteMethod
#define declareSurfaceWriterWriteMethod(Type)                                  \
                                                                               \
        /*! \brief Write field of Type (per face or vertex) */                 \
        virtual fileName write                                                 \
        (                                                                      \
            const word& fieldName,          /*!< Name of field */              \
            const Field<Type>& values       /*!< Field values to write */      \
        )   // override


    // Other IO

        //- Return info proxy,
        //- used to print information to a stream
        virtual InfoProxy<surfaceWriter> info() const noexcept
        {
            return *this;
        }

        //- Output info proxy
        friend Ostream& operator<<
        (
            Ostream& os,
            const InfoProxy<surfaceWriter>& iproxy
        );


    // Housekeeping

        #ifdef Foam_surfaceWriter_directAccess

        //- Access number of expected output fields.
        label& nFields() { return nFields_; }

        //- Access handling of face/point data
        bool& isPointData() { return isPointData_; }

        //- Access handling of spliced output path
        bool& useTimeDir() { return useTimeDir_; }

        //- Access output verbosity
        bool& verbose() { return verbose_; }

        //- Access value of the point merge dimension
        scalar& mergeDim() { return mergeDim_; }

        #endif
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "surfaceWriterI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
