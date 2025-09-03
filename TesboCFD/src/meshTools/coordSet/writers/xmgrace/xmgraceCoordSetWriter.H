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
    Foam::coordSetWriters::xmgraceWriter

Description
    Write coordSet(s) in xmgrace format

    The formatOptions for xmgrace:
    \table
        Property    | Description                       | Required | Default
        buffer      | Use buffered output               | no  | true
        compression | Use file compression              | no  | false
        precision   | Write precision in ascii          | no  | same as IOstream
    \endtable

    When called with a single coordSet, will combine all output fields
    into a single file. Although it is nominally 'buffered', it actually
    just keeps track of the number of fields written into the file.

    When called as a track writer (eg, with multiple coordSets),
    will emit one file per field.

    \section Output file locations

    The \c rootdir normally corresponds to something like
    \c postProcessing/\<name\>

    \subsection Geometry and Fields
    \verbatim
    rootdir
    `-- timeName
        `-- setName.{agr}
    \endverbatim

SourceFiles
    xmgraceCoordSetWriter.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_coordSetWriters_xmgraceWriter_H
#define Foam_coordSetWriters_xmgraceWriter_H

#include "coordSetWriter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class OFstream;

namespace coordSetWriters
{

/*---------------------------------------------------------------------------*\
                        Class xmgraceWriter Declaration
\*---------------------------------------------------------------------------*/

class xmgraceWriter
:
    public coordSetWriter
{
    // Private Data

        //- Output stream option
        IOstreamOption streamOpt_;

        //- ASCII write precision
        unsigned precision_;

        //- Backend output file (when buffering)
        autoPtr<Foam::OFstream> ofile_;

        //- The number of fields written (when buffering)
        label nWritten_;


    // Private Member Functions

        //- Templated write operation
        template<class Type>
        fileName writeTemplate
        (
            const word& fieldName,
            const UPtrList<const Field<Type>>& fieldPtrs
        );

        //- Templated write operation
        template<class Type>
        fileName writeTemplate
        (
            const word& fieldName,    //!< Name of field
            const Field<Type>& vals   //!< Local field values to write
        );

        //- Templated write operation
        template<class Type>
        fileName writeTemplate
        (
            const word& fieldName,
            const List<Field<Type>>& fieldValues
        );


public:

    //- Runtime type information (no debug)
    TypeNameNoDebug("xmgr");


    // Constructors

        //- Default construct
        xmgraceWriter();

        //- Default construct with specified options
        explicit xmgraceWriter(const dictionary& options);

        //- Construct from components
        xmgraceWriter
        (
            const coordSet& coords,
            const fileName& outputPath,
            const dictionary& options = dictionary()
        );

        //- Construct from components
        xmgraceWriter
        (
            const UPtrList<coordSet>& tracks,
            const fileName& outputPath,
            const dictionary& options = dictionary()
        );


    //- Destructor. Calls close()
    virtual ~xmgraceWriter();


    // Member Functions

        //- Enable/disable internal buffering
        virtual bool buffering(const bool on);  // override

        //- Characteristic output file name - information only
        virtual fileName path() const;  // override

        //- Close and reset, clears backend.
        virtual void close(bool force = false);  // override

        //- Begin time step. Clears existing backend.
        virtual void beginTime(const Time& t);  // override

        //- Begin time step. Clears existing backend.
        virtual void beginTime(const instant& inst);  // override

        //- End time step. Clears existing backend.
        virtual void endTime();  // override


        declareCoordSetWriterWriteMethod(label);
        declareCoordSetWriterWriteMethod(scalar);
        declareCoordSetWriterWriteMethod(vector);
        declareCoordSetWriterWriteMethod(sphericalTensor);
        declareCoordSetWriterWriteMethod(symmTensor);
        declareCoordSetWriterWriteMethod(tensor);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace coordSetWriters
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
