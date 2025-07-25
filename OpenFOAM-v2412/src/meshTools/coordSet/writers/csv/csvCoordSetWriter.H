/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021-2022 OpenCFD Ltd.
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
    Foam::coordSetWriters::csvWriter

Description
    Write coordSet(s) in csv format

    The formatOptions for csv:
    \table
        Property    | Description                       | Required | Default
        buffer      | Use buffered output               | no  | true
        compression | Use file compression              | no  | false
        precision   | Write precision in ascii          | no  | same as IOstream
    \endtable

    When called with a single coordSet, will buffer the output of
    all fields and output together in the order of label/scalar/vector/...
    each of which is sorted alphabetically according to the field name.

    When called as a track writer (eg, with multiple coordSets),
    will emit one file per field.

SourceFiles
    csvCoorSetWriter.C
    csvCoorSetWriterImpl.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_coordSetWriters_csvWriter_H
#define Foam_coordSetWriters_csvWriter_H

#include "coordSetWriter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace coordSetWriters
{

/*---------------------------------------------------------------------------*\
                          Class csvWriter Declaration
\*---------------------------------------------------------------------------*/

class csvWriter
:
    public coordSetWriter
{
    // Private Data

        //- Output stream option
        IOstreamOption streamOpt_;

        //- ASCII write precision
        unsigned precision_;


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


protected:

    // Protected Member Functions

        //- Write buffered data
        virtual bool writeBuffered();


public:

    //- Runtime type information (no debug)
    TypeNameNoDebug("csv");


    // Constructors

        //- Default construct
        csvWriter();

        //- Default construct with specified options
        explicit csvWriter(const dictionary& options);

        //- Construct from components
        csvWriter
        (
            const coordSet& coords,
            const fileName& outputPath,
            const dictionary& options = dictionary()
        );

        //- Construct from components
        csvWriter
        (
            const UPtrList<coordSet>& tracks,
            const fileName& outputPath,
            const dictionary& options = dictionary()
        );


    //- Destructor. Calls close()
    virtual ~csvWriter();


    // Member Functions

        //- Enable/disable internal buffering
        virtual bool buffering(const bool on);  // override

        //- Characteristic output file name - information only
        //  \warning incorrect for unbuffered or track output!
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
