/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2023 OpenCFD Ltd.
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
    Foam::functionObjects::ensightWrite

Group
    grpUtilitiesFunctionObjects

Description
    Writes fields in ensight format.

    Example of function object specification:
    \verbatim
    ensight
    {
        type            ensightWrite;
        libs            (utilityFunctionObjects);
        writeControl    writeTime;
        writeInterval   1;
        format          binary;

        timeFormat      scientific;
        timePrecision   5;

        overwrite       true;
        width           12;

        fields          (U p);
        // excludeFields  ("force.*");

        selection
        {
            box
            {
                action  use;
                source  box;
                box     (-0.1 -0.01 -0.1) (0.1 0.30 0.1);
            }
            dome
            {
                action  add;
                shape   sphere;
                origin  (-0.1 -0.01 -0.1);
                radius  0.25;
            }
            centre
            {
                action  subtract;
                source  sphere;
                origin  (-0.1 -0.01 -0.1);
                radius  0.1;
            }
            blob
            {
                action  add;
                source  surface;
                surface triSurfaceMesh;
                name    blob.stl;
            }
        }
    }
    \endverbatim

    \heading Basic Usage
    \table
        Property    | Description                           | Required | Default
        type        | Type name: ensightWrite               | yes |
        fields      | Fields to output                      | yes |
        excludeFields | Exclude fields from output (wordRe list) | no  |
        boundary    | Convert boundary fields               | no  | true
        internal    | Convert internal fields               | no  | true
        nodeValues  | Write values at nodes                 | no  | false
    \endtable

    \heading Ensight Output Options
    \table
        Property    | Description                           | Required | Default
        format      | ascii or binary format                | no  | binary
        width       | Mask width for \c data/XXXX           | no  | 8
        directory   | The output directory name     | no | postProcessing/NAME
        overwrite   | Remove existing directory             | no  | false
        consecutive | Consecutive output numbering          | no  | false
        timeFormat  | Time format (ensight case)            | no  | scientific
        timePrecision | Time precision (ensight case)       | no  | 5
    \endtable

    \heading Output Selection
    \table
        Property    | Description                           | Required | Default
        region      | Name for a single region              | no  | region0
        faceZones   | Select faceZones to write             | no  |
        patches     | Limit to listed patches (wordRe list) | no  |
        excludePatches | Exclude specified patches          | no  |
        selection   | Cell selection (topoSet actions)      | no  | empty dict
    \endtable

Note
    The region of interest is defined by the selection dictionary
    as a set of actions (use,add,subtract,subset,invert).
    Omitting the selection dictionary is the same as specifying the
    conversion of all cells (in the selected regions).
    Omitting the patches entry is the same as specifying the conversion of all
    patches.

    Consecutive output numbering can be used in conjunction with \c overwrite.

See also
    Foam::functionObjects::vtkWrite
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::timeControl
    Foam::cellBitSet::select

SourceFiles
    ensightWrite.C
    ensightWriteImpl.C
    ensightWriteUpdate.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_functionObjects_ensightWrite_H
#define Foam_functionObjects_ensightWrite_H

#include "fvMeshFunctionObject.H"
#include "ensightCase.H"
#include "ensightMesh.H"
#include "ensightOutputFwd.H"

#include "interpolation.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvMeshSubsetProxy.H"
#include "searchableSurfaces.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class dictionary;

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                         Class ensightWrite Declaration
\*---------------------------------------------------------------------------*/

class ensightWrite
:
    public fvMeshFunctionObject
{
    // Private data

        //- Ensight writer options
        ensightMesh::options writeOpts_;

        //- Ensight case options
        ensightCase::options caseOpts_;

        //- The output directory
        fileName outputDir_;

        //- Consecutive output numbering
        bool consecutive_;

        //- Track changes in mesh geometry
        enum polyMesh::readUpdateState meshState_;

        //- Requested selection of fields to process
        wordRes selectFields_;

        //- Requested selection of fields to block
        wordRes blockFields_;

        //- Dictionary of volume selections
        dictionary selection_;

        //- Mesh subset handler
        fvMeshSubset meshSubset_;

        //- Ensight case handler
        autoPtr<ensightCase> ensCase_;

        //- Ensight mesh handler
        autoPtr<ensightMesh> ensMesh_;


    // Private Member Functions

        //- Ensight case handler
        ensightCase& ensCase()
        {
            return *ensCase_;
        }

        //- Ensight mesh handler
        ensightMesh& ensMesh()
        {
            return *ensMesh_;
        }


        //- Update mesh subset according to zones, geometry, bounds
        bool updateSubset(fvMeshSubset& subsetter) const;

        //- Read information for selections
        bool readSelection(const dictionary& dict);

        //- Update meshes, subMeshes etc.
        bool update();


    // Write

        //- Write selected volume fields.
        template<class Type>
        label writeVolFieldsImpl
        (
            ensightOutput::floatBufferType& scratch,
            const fvMeshSubset& proxy,
            const wordHashSet& candidateNames
        );

        //- Write all volume fields
        label writeAllVolFields
        (
            const fvMeshSubset& proxy,
            const wordHashSet& candidateNames
        );

        //- No copy construct
        ensightWrite(const ensightWrite&) = delete;

        //- No copy assignment
        void operator=(const ensightWrite&) = delete;


public:

    //- Runtime type information
    TypeName("ensightWrite");


    // Constructors

        //- Construct from runTime and dictionary.
        ensightWrite
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~ensightWrite() = default;


    // Member Functions

        //- Read the ensightWrite specification
        virtual bool read(const dictionary& dict);

        //- Do nothing
        virtual bool execute();

        //- Write fields, flush case file
        virtual bool write();

        //- Do nothing at the final time-loop
        virtual bool end();

        //- Update for changes of mesh
        virtual void updateMesh(const mapPolyMesh& mpm);

        //- Update for mesh point-motion
        virtual void movePoints(const polyMesh& mpm);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
