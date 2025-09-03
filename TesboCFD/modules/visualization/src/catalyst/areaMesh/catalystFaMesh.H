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
    Foam::catalyst::faMeshInput

Description
    An input (source) for Paraview Catalyst from faMesh regions.

    Produces a multi-block dataset with one block per area mesh with pieces
    from each processor.

    Example specification:
    \verbatim
    someName
    {
        type            area;
        fields          (U p);
    }
    \endverbatim

Usage
    \table
        Property    | Description                           | Required | Default
        type        | input type: \c area                   | yes   |
        region      | name for a single region              | no    | region0
        area        | select a single area                  | no    |
        areas       | wordRe list of multiple areas         | no    |
        fields      | wordRe list of fields                 | yes   |
    \endtable

    The output block structure:
    \verbatim
    |-- area0
    |   |-- piece0
    |   |-- ...
    |   \-- pieceN
    |-- ...
    \-- areaN
        \-- ...
    \endverbatim

Note
    The channel name is that of the defining dictionary.

See also
    Foam::vtk::faMeshAdaptor

SourceFiles
    catalystFaMesh.C
    catalystFaMeshTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_catalyst_faMeshInput_H
#define Foam_catalyst_faMeshInput_H

#include "wordRes.H"
#include "HashPtrTable.H"
#include "catalystInput.H"
#include "foamVtkFaMeshAdaptor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace catalyst
{

/*---------------------------------------------------------------------------*\
                    Class catalyst::faMeshInput Declaration
\*---------------------------------------------------------------------------*/

class faMeshInput
:
    public catalystInput
{
protected:

    // Protected Data

        //- Reference to the time database
        const Time& time_;

        //- The polyMesh region for the area meshes
        word regionName_;

        //- Requested names of areas to process
        wordRes selectAreas_;

        //- Names of fields to process
        wordRes selectFields_;

        //- Pointers to the requested mesh regions
        HashTable<const faMesh*> meshes_;

        //- Backends for OpenFOAM to VTK translation (with internal caching)
        HashPtrTable<vtk::faMeshAdaptor> backends_;


    // Protected Member Functions

        //- Update/synchronize internals with catalyst backends
        void update();

        //- No copy construct
        faMeshInput(const faMeshInput&) = delete;

        //- No copy assignment
        void operator=(const faMeshInput&) = delete;

public:

    //- Declare type-name (with debug switch)
    ClassName("catalyst::faMesh");


    // Constructors

        //- Construct from Time and dictionary
        faMeshInput
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~faMeshInput() = default;


    // Member Functions

        //- Read the specification
        virtual bool read(const dictionary& dict);

        //- Update for changes of mesh or mesh point-motion
        virtual void update(polyMesh::readUpdateState state);

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
