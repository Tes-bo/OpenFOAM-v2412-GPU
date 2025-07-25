/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2022 OpenCFD Ltd.
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
    Foam::noiseModels::surfaceNoise

Description
    Perform noise analysis on surface-based pressure data.

    Input data is read from a dictionary, e.g.

    \verbatim
    // Pressure reference
    pRef            0;

    // Number of samples in sampling window, default = 2^16 (=65536)
    N               4096;

    // Lower frequency bounds
    fl              25;

    // Upper frequency bounds
    fu              25;

    // Start time
    startTime       0;

    windowModel     <modelType>
    <modelType>Coeffs
    {
        ...
    }

    // Input file
    file            "postProcessing/faceSource1/surface/patch/patch.case";
    //files           ("postProcessing/faceSource1/surface/patch/patch.case");

    // Write interval for FFT data, default = 1
    fftWriteInterval 100;

    // Area-weighted averaging switch, default = no (ensemble) for backwards
    // compatibility
    areaAverage     yes;

    // Surface reader
    reader          ensight;

    // Surface writer
    writer          ensight;  // none

    // Collate times for ensight output - ensures geometry is only written once
    writeOptions
    {
        ensight
        {
            collateTimes true;
        }

        // Write Prmsf; default = yes
        writePrmsf  no;

        // Write SPL; default = yes
        writeSPL    yes;

        // Write PSD; default = yes
        writePSD    yes;

        // Write PSDf; default = yes
        writePSDf   no;

        // Write writeOctaves; default = yes
        writeOctaves yes;
    }

    \endverbatim

    Communication Options
    \table
      Property     | Description                           | Type | Req'd | Dflt
      commsType    | Communication type for send/receive field | word | no  | scheduled
      broadcast    | (Experimental) broadcast all fields   | bool | no | false
    \endtable

SourceFiles
    surfaceNoise.C

SeeAlso
    noiseModel.H

\*---------------------------------------------------------------------------*/

#ifndef Foam_noiseModels_surfaceNoise_H
#define Foam_noiseModels_surfaceNoise_H

#include "noiseModel.H"
#include "scalarField.H"
#include "UPstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class globalIndex;
class surfaceReader;
class surfaceWriter;

namespace noiseModels
{

/*---------------------------------------------------------------------------*\
                        Class surfaceNoise Declaration
\*---------------------------------------------------------------------------*/

class surfaceNoise
:
    public noiseModel
{
protected:

    // Protected Data

        //- Input file names
        List<fileName> inputFileNames_;

        //- Name of pressure field
        word pName_;

        //- Index of pressure field in reader field list
        label pIndex_;

        //- Sample times
        scalarList times_;

        //- Time step (constant)
        scalar deltaT_;

        //- Start time index
        label startTimeIndex_;

        //- Global number of surface faces
        label nFaces_;

        //- Frequency data output interval, default = 1
        //  nSamples/2 data points are returned from the FFT, which can
        //  result in a very large number of output files (1 per frequency)
        label fftWriteInterval_;

        //- Apply area average; default = no (ensemble average) for backwards
        //- compatibility
        bool areaAverage_;

        //- Use broadcast to send entire field to sub-ranks
        bool useBroadcast_;

        //- Communication type (for sending/receiving fields)
        UPstream::commsTypes commType_;

        //- Reader type
        word readerType_;

        //- Pointer to the surface reader
        mutable autoPtr<surfaceReader> readerPtr_;

        //- Pointer to the surface writer
        mutable autoPtr<surfaceWriter> writerPtr_;


    // Protected Member Functions

        //- Initialise
        void initialise(const fileName& fName);

        //- Read surface data
        void readSurfaceData
        (
            const globalIndex& procFaceAddr,
            List<scalarField>& pData
        );

        //- Calculate the area average value
        scalar surfaceAverage
        (
            const scalarField& data,
            const globalIndex& procFaceAddr
        ) const;

        //- Write surface data to file
        //  Returns the area average value
        scalar writeSurfaceData
        (
            const fileName& outDirBase,
            const word& fName,
            const word& title,
            const scalar freq,
            const scalarField& data,
            const globalIndex& procFaceAddr,
            const bool writeSurface
        ) const;


public:

    //- Runtime type information
    TypeName("surfaceNoise");


    //- Constructor
    surfaceNoise
    (
        const dictionary& dict,
        const objectRegistry& obr,
        const word& name = typeName,
        const bool readFields = true
    );

    //- Destructor
    virtual ~surfaceNoise() = default;


    // Public Member Functions

        //- Read from dictionary
        virtual bool read(const dictionary& dict);

        //- Calculate
        virtual void calculate();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace noiseModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
