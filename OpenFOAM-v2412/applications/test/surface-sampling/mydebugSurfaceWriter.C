/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022-2023 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "mydebugSurfaceWriter.H"
#include "globalIndex.H"
#include "argList.H"
#include "OFstream.H"
#include "OSspecific.H"
#include "IOmanip.H"
#include "Time.H"
#include "pointIOField.H"
#include "primitivePatch.H"
#include "profiling.H"
#include "surfaceWriterMethods.H"
#include "PrecisionAdaptor.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace surfaceWriters
{
    defineTypeName(mydebugWriter);
    addToRunTimeSelectionTable(surfaceWriter, mydebugWriter, word);
    addToRunTimeSelectionTable(surfaceWriter, mydebugWriter, wordDict);
}
}


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

// Type narrowing - base implementation is pass-through
template<class Type> struct narrowType
{
    typedef Type type;
};

template<> struct narrowType<double>
{
    typedef float type;
};

template<> struct narrowType<Vector<double>>
{
    typedef Vector<float> type;
};

template<> struct narrowType<SphericalTensor<double>>
{
    typedef SphericalTensor<float> type;
};

template<> struct narrowType<SymmTensor<double>>
{
    typedef SymmTensor<float> type;
};

// FIXME: Not sure why this one seems to be broken...
//
// template<> struct narrowType<Tensor<double>>
// {
//     typedef Tensor<float> type;
// };

} // End namespace Foam


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::surfaceWriters::mydebugWriter::mergeField
(
    const Field<Type>& fld
) const
{
    addProfiling(merge, "debugWriter::merge-field");

    // This is largely identical to surfaceWriter::mergeField()
    // but with narrowing for communication
    if (narrowTransfer_ && parallel_ && UPstream::parRun())
    {
        // The narrowed type
        typedef typename narrowType<Type>::type narrowedType;

        // Ensure geometry is also merged
        merge();

        // Gather all values
        auto tfield = tmp<Field<Type>>::New();
        auto& allFld = tfield.ref();

        // Update any expired global index (as required)

        const globalIndex& globIndex =
        (
            this->isPointData()
          ? mergedSurf_.pointGlobalIndex()
          : mergedSurf_.faceGlobalIndex()
        );

        ConstPrecisionAdaptor<narrowedType, Type> input(fld);
        PrecisionAdaptor<narrowedType, Type> output(allFld);

        globIndex.gather
        (
            input.cref(),   // fld,
            output.ref(),   // allFld,
            UPstream::msgType(),
            commType_,
            UPstream::worldComm
        );

        // Commit adapted content changes
        input.commit();
        output.commit();

        // Discard adaptors
        input.clear();
        output.clear();

        // Renumber (point data) to correspond to merged points
        if
        (
            UPstream::master()
         && this->isPointData()
         && mergedSurf_.pointsMap().size()
        )
        {
            inplaceReorder(mergedSurf_.pointsMap(), allFld);
            allFld.resize(mergedSurf_.points().size());
        }

        return tfield;
    }

    return surfaceWriter::mergeField(fld);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfaceWriters::mydebugWriter::mydebugWriter()
:
    surfaceWriter(),
    enableMerge_(true),
    enableWrite_(false),
    header_(true),
    narrowTransfer_(false),
    streamOpt_(IOstreamOption::BINARY)
{}


Foam::surfaceWriters::mydebugWriter::mydebugWriter
(
    const dictionary& options
)
:
    surfaceWriter(options),
    enableMerge_(options.getOrDefault("merge", true)),
    enableWrite_(options.getOrDefault("write", false)),
    header_(true),
    narrowTransfer_(options.getOrDefault("narrow", false)),
    streamOpt_(IOstreamOption::BINARY)
{
    Info<< "Using debug surface writer ("
        << (this->isPointData() ? "point" : "face") << " data):"
        << " commsType=" << UPstream::commsTypeNames[commType_]
        << " merge=" << Switch::name(enableMerge_)
        << " write=" << Switch::name(enableWrite_)
        << " narrow=" << Switch::name(narrowTransfer_)
        << endl;
}


Foam::surfaceWriters::mydebugWriter::mydebugWriter
(
    const meshedSurf& surf,
    const fileName& outputPath,
    bool parallel,
    const dictionary& options
)
:
    mydebugWriter(options)
{
    open(surf, outputPath, parallel);
}


Foam::surfaceWriters::mydebugWriter::mydebugWriter
(
    const pointField& points,
    const faceList& faces,
    const fileName& outputPath,
    bool parallel,
    const dictionary& options
)
:
    mydebugWriter(options)
{
    open(points, faces, outputPath, parallel);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::surfaceWriters::mydebugWriter::serialWriteGeometry
(
    const regIOobject& iopts,
    const meshedSurf& surf
)
{
    const pointField& points = surf.points();
    const faceList& faces = surf.faces();

    if (verbose_)
    {
        if (this->isPointData())
        {
            Info<< "Writing points: " << iopts.objectPath() << endl;
        }
        else
        {
            Info<< "Writing face centres: " << iopts.objectPath() << endl;
        }
    }

    // Like regIOobject::writeObject without instance() adaptation
    // since this would write to e.g. 0/ instead of postProcessing/

    autoPtr<primitivePatch> ppPtr;

    {
        OFstream os(iopts.objectPath(), streamOpt_);

        if (header_)
        {
            iopts.writeHeader(os);
        }

        if (this->isPointData())
        {
            // Just like writeData, but without copying beforehand
            os << points;
        }
        else
        {
            ppPtr.reset(new primitivePatch(SubList<face>(faces), points));

            // Just like writeData, but without copying beforehand
            os << ppPtr().faceCentres();
        }

        if (header_)
        {
            IOobject::writeEndDivider(os);
        }
    }
}


Foam::fileName Foam::surfaceWriters::mydebugWriter::write()
{
    checkOpen();

    // Geometry: rootdir/surfaceName/"points"
    // Field:    rootdir/surfaceName/<TIME>/field

    fileName surfaceDir = outputPath_;

    if (parallel_ && !enableMerge_)
    {
        if (verbose_)
        {
            Info<< "Not merging or writing" << nl;
        }

        // Pretend to have succeeded
        wroteGeom_ = true;
        return surfaceDir;
    }


    const meshedSurf& surf = surface();
    // const meshedSurfRef& surf = adjustSurface();

    // Dummy Time to use as objectRegistry
    autoPtr<Time> dummyTimePtr;

    if (enableWrite_)
    {
        #if (OPENFOAM <= 2306)
        dummyTimePtr = Time::New(argList::envGlobalPath());
        #else
        dummyTimePtr = Time::NewGlobalTime();
        #endif
    }
    else if (verbose_)
    {
        Info<< "Not writing: " << surf.faces().size() << " faces" << nl;
    }

    if (enableWrite_ && (UPstream::master() || !parallel_))
    {
        if (!isDir(surfaceDir))
        {
            mkDir(surfaceDir);
        }

        // Write sample locations
        pointIOField iopts
        (
            IOobject
            (
                surfaceDir/"points",
                *dummyTimePtr,
                IOobjectOption::NO_READ,
                IOobjectOption::NO_WRITE,
                IOobjectOption::NO_REGISTER
            )
        );
        iopts.note() = (this->isPointData() ? "point data" : "face data");

        serialWriteGeometry(iopts, surf);
    }

    wroteGeom_ = true;
    return surfaceDir;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Foam::fileName Foam::surfaceWriters::mydebugWriter::writeTemplate
(
    const word& fieldName,
    const Field<Type>& localValues
)
{
    checkOpen();

    // Geometry: rootdir/surfaceName/"points"
    // Field:    rootdir/surfaceName/<TIME>/field

    fileName surfaceDir = outputPath_;

    const fileName outputFile(surfaceDir/timeName()/fieldName);

    if (parallel_ && !enableMerge_)
    {
        if (verbose_)
        {
            Info<< "Not merging or writing" << nl;
        }

        // Pretend to have succeeded
        wroteGeom_ = true;
        return surfaceDir;
    }


    // Implicit geometry merge()
    tmp<Field<Type>> tfield = mergeField(localValues);

    // Dummy Time to use as objectRegistry
    autoPtr<Time> dummyTimePtr;

    if (enableWrite_)
    {
        #if (OPENFOAM <= 2306)
        dummyTimePtr = Time::New(argList::envGlobalPath());
        #else
        dummyTimePtr = Time::NewGlobalTime();
        #endif
    }
    else if (verbose_)
    {
        Info<< "Not writing: " << tfield().size()
            << ' ' << pTraits<Type>::typeName
            << " values" << nl;
    }

    const meshedSurf& surf = surface();
    // const meshedSurfRef& surf = adjustSurface();

    if (enableWrite_ && (UPstream::master() || !parallel_))
    {
        if (!isDir(outputFile.path()))
        {
            mkDir(outputFile.path());
        }

        // Write sample locations
        {
            pointIOField iopts
            (
                IOobject
                (
                    surfaceDir/"points",
                    *dummyTimePtr,
                    IOobjectOption::NO_READ,
                    IOobjectOption::NO_WRITE,
                    IOobjectOption::NO_REGISTER
                )
            );
            iopts.note() = (this->isPointData() ? "point data" : "face data");

            serialWriteGeometry(iopts, surf);
        }

        // Write field
        {
            IOField<Type> iofld
            (
                IOobject
                (
                    outputFile,
                    *dummyTimePtr,
                    IOobjectOption::NO_READ,
                    IOobjectOption::NO_WRITE,
                    IOobjectOption::NO_REGISTER
                )
            );
            iofld.note() = (this->isPointData() ? "point data" : "face data");

            OFstream os(iofld.objectPath(), streamOpt_);

            if (header_)
            {
                iofld.writeHeader(os);
            }

            // Just like writeData, but without copying beforehand
            os << tfield();

            if (header_)
            {
                IOobject::writeEndDivider(os);
            }
        }
    }

    wroteGeom_ = true;
    return surfaceDir;
}


// Field writing methods
defineSurfaceWriterWriteFields(Foam::surfaceWriters::mydebugWriter);


// ************************************************************************* //
