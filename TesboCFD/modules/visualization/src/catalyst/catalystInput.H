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
    Foam::catalyst::catalystInput

Description
    An abstract input (source) for Paraview Catalyst.

See also
    Foam::catalyst::coprocess
    Foam::catalyst::dataQuery
    Foam::catalyst::timeQuery

SourceFiles
    catalystInput.C

\*---------------------------------------------------------------------------*/

#ifndef catalyst_catalystInput_H
#define catalyst_catalystInput_H

#include "className.H"
#include "polyMesh.H"
#include "catalystTools.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class dictionary;
class mapPolyMesh;
class Time;

namespace catalyst
{

/*---------------------------------------------------------------------------*\
                   Class catalyst::catalystInput Declaration
\*---------------------------------------------------------------------------*/

class catalystInput
{
    // Private Data

        //- The Catalyst input channel name
        word name_;

public:

    //- Declare type-name (with debug switch)
    ClassName("catalyst::input");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            catalystInput,
            dictionary,
            (
                const word& name,
                const Time& runTime,
                const dictionary& dict
            ),
            (name, runTime, dict)
        );


    // Selectors

        //- Return a reference to the selected input
        static autoPtr<catalystInput> New
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    // Constructors

        //- Construct with given Catalyst input channel name
        explicit catalystInput(const word& channel);


    //- Destructor
    virtual ~catalystInput() = default;


    // Member Functions

        //- The name of the Catalyst input channel
        virtual const word& name() const
        {
            return name_;
        }

        //- Read the specification
        virtual bool read(const dictionary& dict);

        //- Update for changes of mesh or mesh point-motion
        virtual void update(polyMesh::readUpdateState state);

        //- Add channels (with fields) to data query
        virtual label addChannels(dataQuery& dataq) = 0;

        //- Convert channels to vtkMultiBlockDataSet outputs
        virtual bool convert(dataQuery& dataq, outputChannels& outputs) = 0;


        //- Print information
        virtual Ostream& print(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace catalyst
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
