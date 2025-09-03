/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2021-2023 OpenCFD Ltd.
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
    Foam::coordSetWriters::ensightCoordSetWriter

Description
    A coordSetWriter for Ensight format.

    \verbatim
    formatOptions
    {
        ensight
        {
            format          ascii;
            collateTimes    true;
        }
    }
    \endverbatim

    Format options:
    \table
        Property | Description                              | Required | Default
        format   | ascii/binary                             | no  | binary
        collateTimes | use common geometry for times        | no  | true
        timeFormat    | Time format (ensight case)          | no  | scientific
        timePrecision | Time precision (ensight case)       | no  | 5
     \endtable
    \endtable

SourceFiles
    ensightCoordSetWriter.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_coordSetWriters_ensightWriter_H
#define Foam_coordSetWriters_ensightWriter_H

#include "coordSetWriter.H"
#include "ensightCase.H"
#include "ensightWriterCaching.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class ensightGeoFile;

namespace coordSetWriters
{

/*---------------------------------------------------------------------------*\
                     Class ensightWriter Declaration
\*---------------------------------------------------------------------------*/

class ensightWriter
:
    public coordSetWriter
{
    // Data Types

        //- Types of element output
        enum elemOutputType
        {
            NO_ELEMENTS = 0,
            WRITE_POINTS = 1,
            WRITE_LINES = 2
        };


    // Private Data

        //- Ensight case options
        ensightCase::options caseOpts_;

        //- Collate times (default: true)
        bool collateTimes_;

        //- Cached information for times, geometry, fields (collated)
        ensightOutput::writerCaching caching_;


    // Private Member Functions

        //- Write geometry
        void writeGeometry
        (
            ensightGeoFile& os,
            elemOutputType elemOutput
        ) const;

        //- Write geometry
        fileName writeCollated(const bool writeTracks);

        //- Write geometry
        fileName writeUncollated(const bool writeTracks);

        //- Templated write operation - one file per timestep
        template<class Type>
        fileName writeCollated
        (
            const word& fieldName,
            const UPtrList<const Field<Type>>& fieldPtrs,
            elemOutputType elemOutput
        );

        //- Templated write operation - all time steps in single file
        template<class Type>
        fileName writeUncollated
        (
            const word& fieldName,
            const UPtrList<const Field<Type>>& fieldPtrs,
            elemOutputType elemOutput
        );

        //- Templated write operation
        template<class Type>
        fileName writeTemplate
        (
            const word& fieldName,
            const Field<Type>& values
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
    TypeNameNoDebug("ensight");


    // Constructors

        //- Default construct
        ensightWriter();

        //- Default construct with specified options
        explicit ensightWriter(const dictionary& options);

        //- Construct from components
        ensightWriter
        (
            const coordSet& coords,
            const fileName& outputPath,
            const dictionary& options = dictionary()
        );

        //- Construct from components
        ensightWriter
        (
            const UPtrList<coordSet>& tracks,
            const fileName& outputPath,
            const dictionary& options = dictionary()
        );


    //- Destructor. Calls close()
    virtual ~ensightWriter();


    // Member Functions

        //- Finish output, clears output times.
        //  Later reuse will rebuild times from fieldsDict file cache.
        virtual void close(bool force = false);  // override


    // Write

        //- Expected (characteristic) output file name - information only
        //  \warning incorrect for uncollated output!
        virtual fileName path() const;  // override

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
