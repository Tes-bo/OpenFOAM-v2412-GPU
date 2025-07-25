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
    Foam::catalyst::fvMeshInput

Description
    An input (source) for Paraview Catalyst from fvMesh regions.

    Produces a multi-block dataset with one block per region.
    Each region block contains up to two blocks corresponding
    to the internal (volume) mesh (block 0) and the boundary (block 1),
    which are further divided into sub-blocks for each patch.
    The lowest blocks contain pieces from each processor.

    Example specification:
    \verbatim
    someName
    {
        type        default;
        regions     ( ".*solid" )
        boundary    false;
        fields      (U p);
    }
    \endverbatim

Usage
    \table
        Property    | Description                           | Required | Default
        type        | input type: \c default                | no    | default
        region      | name for a single region              | no    | region0
        regions     | wordRe list of regions                | no    |
        patches     | explicit wordRe list of patches       | no    |
        fields      | wordRe list of fields                 | yes   |
        boundary    | convert boundary fields               | no    | true
        internal    | convert internal fields               | no    | true
        decompose   | decompose polyhedra (experimental)    | no    | false
    \endtable

    The output block structure:
    \verbatim
    |-- region0
    |   |-- internal
    |   |   |-- piece0
    |   |   |-- ...
    |   |   \-- pieceN
    |   \-- boundary
    |       |-- patch0
    |       |   |-- piece0
    |       |   |-- ...
    |       |   \-- pieceN
    |       |-- ...
    |       \-- patchN
    |           |-- piece0
    |           |-- ...
    |           \-- pieceN
    |-- ...
    \-- regionN
        |-- internal
        |   \-- ...
        \-- boundary
            \-- ...
    \endverbatim

Note
    The channel name is that of the defining dictionary.
    If the \c patches entry is missing or an empty list,
    all non-processor patches will be used for the boundary.
    When it is non-empty, only the explicitly specified (non-processor)
    patch names will be used.

See also
    Foam::vtk::fvMeshAdaptor

SourceFiles
    catalystFvMesh.C
    catalystFvMeshTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef catalyst_fvMeshInput_H
#define catalyst_fvMeshInput_H

#include "wordRes.H"
#include "HashPtrTable.H"
#include "catalystInput.H"
#include "foamVtkFvMeshAdaptor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace catalyst
{

/*---------------------------------------------------------------------------*\
                    Class catalyst::fvMeshInput Declaration
\*---------------------------------------------------------------------------*/

class fvMeshInput
:
    public catalystInput
{
protected:

    // Protected Data

        //- Channel enumeration from vtk::fvMeshAdaptor
        using channelType = vtk::fvMeshAdaptor::channelType;

        //- Reference to the time database
        const Time& time_;

        //- Sub-block selection (internal/boundary) mesh and fields
        channelType channelOpt_;

        //- Decompose polyhedra (experimental, perhaps of questionable use)
        bool decomposeOpt_;

        //- Requested names of regions to process
        wordRes selectRegions_;

        //- Requested names of patches to process
        wordRes selectPatches_;

        //- Names of fields to process
        wordRes selectFields_;

        //- Pointers to the requested mesh regions
        HashTable<const fvMesh*> meshes_;

        //- Backends for OpenFOAM to VTK translation (with internal caching)
        HashPtrTable<vtk::fvMeshAdaptor> backends_;


    // Protected Member Functions

        //- Update/synchronize internals with catalyst backends
        void update();

        //- No copy construct
        fvMeshInput(const fvMeshInput&) = delete;

        //- No copy assignment
        void operator=(const fvMeshInput&) = delete;

public:

    //- Declare type-name (with debug switch)
    ClassName("catalyst::fvMesh");


    // Constructors

        //- Construct from Time and dictionary
        fvMeshInput
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~fvMeshInput() = default;


    // Member Functions

        //- Read the specification
        virtual bool read(const dictionary& dict);

        //- Update for changes of mesh or mesh point-motion
        virtual void update(polyMesh::readUpdateState state);

        //- Add available channels (with fields) to data query
        //
        //  \return the number of channels added
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
