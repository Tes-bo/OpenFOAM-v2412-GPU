/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018 OpenCFD Ltd.
    Copyright (C) 2018 CINECA
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
    Foam::catalyst::coprocess

Description
    Low-level interface between OpenFOAM and ParaView Catalyst.

SourceFiles
    catalystCoprocess.C

\*---------------------------------------------------------------------------*/

#ifndef catalyst_coprocess_H
#define catalyst_coprocess_H

#include "className.H"
#include "wordList.H"
#include "stringList.H"
#include "catalystTools.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace catalyst
{

/*---------------------------------------------------------------------------*\
                     Class catalyst::coprocess Declaration
\*---------------------------------------------------------------------------*/

class coprocess
{
    // Private Data

        //- The currently active ParaView Catalyst co-process
        vtkCPProcessor* coproc_;


    // Private Member Functions

        //- Process multiple output channels
        template<class DataType>
        bool processImpl
        (
            const dataQuery& dataq,
            HashTable<vtkSmartPointer<DataType>>& outputs
        );

public:

    //- Declare type-name (with debug switch)
    ClassName("catalyst");


    // Constructors

        //- Default construct. Does not initialize catalyst.
        constexpr coprocess() noexcept
        :
            coproc_(nullptr)
        {}


    //- Destructor. Shutdown process.
    ~coprocess();


    // Member Functions

        //- \return True if the coprocess has been initialized.
        bool good() const;

        //- Reset/initialize pipeline with output directory, but without
        //- pipeline scripts.
        void reset(const fileName& outputDir);

        //- Reset/initialize pipeline with output directory and with
        //- pipeline scripts.
        //  The script names must have already been resolved prior calling this.
        void reset
        (
            const fileName& outputDir,
            const UList<string>& scripts
        );

        //- Query the coprocess pipelines if they should be executed at this
        //- iteration and possibly which fields they require.
        //
        // \param[in,out] dataq the data query for catalyst.
        //     On input it contains the published channel names with their
        //     published fields, as well as the current simulation time
        //     (index, value).
        //     On output the data description part of dataq will be filled
        //     with the requested field names added per channel.
        //
        // \return number of channels request at this time
        label query(dataQuery& dataq);

        //- Process the output channels
        bool process(const dataQuery& dataq, outputChannels& outputs);

        //- Finalize
        void stop();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace catalyst
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
