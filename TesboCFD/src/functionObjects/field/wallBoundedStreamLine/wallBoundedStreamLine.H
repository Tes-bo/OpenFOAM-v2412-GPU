/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2015-2020 OpenCFD Ltd.
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
    Foam::functionObjects::wallBoundedStreamLine

Group
    grpFieldFunctionObjects

Description
    Generates streamline data by sampling a set of user-specified fields along a
    particle track, transported by a user-specified velocity field, constrained
    to a patch.

    Operands:
    \table
      Operand       | Type  | Location
      input         | -     | -
      output file   | -     | $FOAM_CASE/postProcessing/\<FO\>/\<time\>/\<file\>
      output field  | -     | -
    \endtable

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    wallBoundedStreamLine1
    {
        // Mandatory entries (unmodifiable)
        type            wallBoundedStreamLine;
        libs            (fieldFunctionObjects);

        // Optional entries
        U               <velocity-name>;
        direction       forward;
        cloud           particleTracks;

        // Mandatory entries (runtime modifiable)
        fields          (<field1> ... <fieldN>);
        setFormat       vtk;
        lifeTime        10000;
        seedSampleSet
        {
            type        patchSeed;
            patches     (wall);
            axis        x;
            maxPoints   20000;
        }

        // Optional entries (runtime modifiable)
        bounds          (0.2 -10 -10)(0.22 10 10);
        trackLength     1e-3;
        nSubCycle       1;
        interpolationScheme cellPoint;

        // Deprecated
        // trackForward true;

        // Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                           | Type | Req'd | Dflt
      type         | Type name: wallBoundedStreamLine      | word |  yes  | -
      libs         | Library name: fieldFunctionObjects    | word |  yes  | -
      U            | Name of tracking velocity field       | word |  yes  | -
      fields       | Names of operand fields to sample     | wordList | yes | -
      setFormat    | Type of output data                   | word |  yes  | -
      direction    | Direction (enum) to track             | word |  no   | forward
      lifetime     | Maximum number of particle tracking steps | label | yes | -
      cloud        | Cloud name to use for streamlines     | word |  no  | typeName
      seedSampleSet| Seeding description (see below)       | dict |  yes  | -
      bounds       | Bounding box to trim tracks    | vector pair | no | -
      trackLength  | Tracking segment length              | scalar | no | VGREAT
      nSubCycle    | Number of tracking steps per cell    | label  | no | 1
      interpolationScheme | Interp. scheme for sample   | word | no | cellPoint
    \endtable

    Example types for the \c seedSampleSet sub-dict:
    \verbatim
      uniform   | uniform particle seeding
      cloud     | cloud of points
      patchSeed | seeding via patch faces
      triSurfaceMeshPointSet | points according to a tri-surface mesh
    \endverbatim

    Options for the \c setFormat entry:
    \verbatim
        csv
        ensight
        gnuplot
        nastran
        raw
        vtk
        xmgr
    \endverbatim

    Options for the \c direction entry:
    \verbatim
        forward
        backward
        bidirectional
    \endverbatim

    The inherited entries are elaborated in:
     - \link functionObject.H \endlink

    Usage by the \c postProcess utility is not available.

Note
    When specifying the track resolution, the \c trackLength OR \c nSubCycle
    option should be used.

See also
    - Foam::functionObject
    - Foam::functionObjects::streamLineBase
    - ExtendedCodeGuide::functionObjects::field::wallBoundedStreamLine

SourceFiles
    wallBoundedStreamLine.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_functionObjects_wallBoundedStreamLine_H
#define Foam_functionObjects_wallBoundedStreamLine_H

#include "streamLineBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                    Class wallBoundedStreamLine Declaration
\*---------------------------------------------------------------------------*/

class wallBoundedStreamLine
:
    public functionObjects::streamLineBase
{
protected:

    // Protected Member Functions

        //- Find wall tet on cell
        Tuple2<tetIndices, point> findNearestTet
        (
            const bitSet& isWallPatch,
            const point& seedPt,
            const label celli
        ) const;

        //- Push a point a tiny bit towards the centre of the triangle it is in
        //- to avoid tracking problems
        point pushIn
        (
            const triPointRef& tri,
            const point& pt
        ) const;


public:

    //- Runtime type information
    TypeName("wallBoundedStreamLine");


    // Constructors

        //- Construct from Time and dictionary
        wallBoundedStreamLine
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Construct from Time and dictionary and list of fields to sample
        wallBoundedStreamLine
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict,
            const wordList& fieldNames
        );

        //- No copy construct
        wallBoundedStreamLine(const wallBoundedStreamLine&) = delete;

        //- No copy assignment
        void operator=(const wallBoundedStreamLine&) = delete;


    //- Destructor
    virtual ~wallBoundedStreamLine() = default;


    // Member Functions

        //- Read settings
        virtual bool read(const dictionary&);

        //- Do the actual tracking to fill the track data
        virtual void track();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
