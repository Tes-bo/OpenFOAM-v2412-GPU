/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018 OpenCFD Ltd.
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
    Foam::catalyst::cloudInput

Description
    An input (source) for Paraview Catalyst from clouds (lagrangian).

    Produces a multi-block dataset with one block per cloud with pieces
    from each processor.

    Example specification:
    \verbatim
    someName
    {
        type        cloud;
        cloud       myCloud;
        fields      (U T rho);
    }
    \endverbatim

Usage
    \table
        Property    | Description                           | Required | Default
        type        | input type: \c cloud                  | yes   |
        region      | name for a single region              | no    | region0
        cloud       | name for a single cloud               | no  | defaultCloud
        clouds      | wordRe list of clouds                 | no    |
        fields      | wordRe list of fields                 | yes   |
    \endtable

    The output block structure:
    \verbatim
    |-- cloud0
    |   |-- piece0
    |   |-- ...
    |   \-- pieceN
    |-- ...
    \-- cloudN
        \-- ...
    \endverbatim

Note
    The channel name is that of the defining dictionary.

See also
    Foam::vtk::cloudAdaptor

SourceFiles
    catalystCloud.C
    catalystCloudTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef catalyst_cloudInput_H
#define catalyst_cloudInput_H

#include "wordRes.H"
#include "catalystInput.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace catalyst
{

/*---------------------------------------------------------------------------*\
                    Class catalyst::cloudInput Declaration
\*---------------------------------------------------------------------------*/

class cloudInput
:
    public catalystInput
{
protected:

    // Protected Data

        //- Reference to the time database
        const Time& time_;

        //- The region name for the clouds
        word regionName_;

        //- Requested names of clouds to process
        wordRes selectClouds_;

        //- Subset of cloud fields to process
        wordRes selectFields_;


    // Protected Member Functions

        //- No copy construct
        cloudInput(const cloudInput&) = delete;

        //- No copy assignment
        void operator=(const cloudInput&) = delete;

public:

    //- Declare type-name (with debug switch)
    ClassName("catalyst::cloud");


    // Constructors

        //- Construct from Time and dictionary
        cloudInput
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~cloudInput() = default;


    // Member Functions

        //- Read the specification
        virtual bool read(const dictionary& dict);

        //- Add available channels (with fields) to data query
        virtual label addChannels(dataQuery& dataq);

        //- Convert channels to vtkMultiBlockDataSet outputs
        virtual bool convert(dataQuery& dataq, outputChannels& outputs);

        //- Print information
        virtual Ostream& print(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace catalyst
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
