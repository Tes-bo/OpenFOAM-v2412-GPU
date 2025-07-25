/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::coordSetWriters::vtkWriter

Description

SourceFiles
    vtkCoordSetWriter.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_coordSetWriters_vtkWriter_H
#define Foam_coordSetWriters_vtkWriter_H

#include "coordSetWriter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class OFstream;

namespace vtk
{
// Forward Declarations
class outputOptions;
class coordSetWriter;
}

namespace coordSetWriters
{

/*---------------------------------------------------------------------------*\
                     Class vtkCoordSetWriter Declaration
\*---------------------------------------------------------------------------*/

class vtkWriter
:
    public coordSetWriter
{
    // Private Data

        //- The VTK output format type.
        //  Stored as a raw value to avoid a header dependency on fileFormats
        unsigned fmtType_;

        //- ASCII write precision
        unsigned precision_;

        //- Backend writer - master only
        autoPtr<Foam::vtk::coordSetWriter> writer_;


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
    TypeNameNoDebug("vtk");

    // Constructors

        //- Default construct
        vtkWriter();

        //- Construct with some output options
        explicit vtkWriter(const vtk::outputOptions& opts);

        //- Default construct with specified options
        explicit vtkWriter(const dictionary& options);

        //- Construct from components
        vtkWriter
        (
            const coordSet& coords,
            const fileName& outputPath,
            const dictionary& options = dictionary()
        );

        //- Construct from components
        vtkWriter
        (
            const UPtrList<coordSet>& tracks,
            const fileName& outputPath,
            const dictionary& options = dictionary()
        );


    //- Destructor. Calls close()
    virtual ~vtkWriter();


    // Member Functions

        //- Expected (characteristic) output file name - information only
        virtual fileName path() const;  // override

        //- Close and reset, clears backend.
        virtual void close(bool force = false);  // override

        //- Begin time step. Clears existing backend.
        virtual void beginTime(const Time& t);  // override

        //- Begin time step. Clears existing backend.
        virtual void beginTime(const instant& inst);  // override

        //- End time step. Clears existing backend.
        virtual void endTime();  // override


    // Write

        //- Write geometry to file.
        virtual fileName write(); // override

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
