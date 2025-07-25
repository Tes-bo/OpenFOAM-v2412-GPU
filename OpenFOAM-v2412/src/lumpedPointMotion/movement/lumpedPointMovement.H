/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2020 OpenCFD Ltd.
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
    Foam::lumpedPointMovement

Description
    The movement \a driver that describes initial point locations,
    the current state of the points/rotations,
    and forwarding to the externalFileCoupler
    communication coordinator.

    The lumpedPointIOMovement class is simply a registered version of
    the same.

    See externalFileCoupler for more information about some of
    the communication parameters and setup.

    \heading Dictionary parameters
    \table
        Property    | Description                           | Required | Default
        points      | Initial locations of lumped points    | yes |
        pointLabels | The FEA ids for the points            | no  |
        origin      | Shift offset when reading points      | no  | (0 0 0)
        rotationOrder | The Euler rotation order            | no  | zxz
        degrees     | Input rotations in degrees            | no  | false
        relax       | Relaxation/scaling for updating positions | no  | 1
        controllers | Motion controllers (dictionary)       | yes |
        forces      | Force settings (dictionary)           | no  |
        communication | Communication settings (dictionary) | yes |
    \endtable

    \heading Parameters for communication dictionary
    \table
        Property    | Description                       | Required | Default
        inputName   | Name of positions input file      | yes |
        outputName  | Name of forces output file        | yes |
        logName     | Name of log file                  | no  | movement.log
        inputFormat | Input format: dictionary/plain    | yes |
        outputFormat | Output format: dictionary/plain  | yes |
        scaleInput  | Input scaling parameter dictionary | no |
        scaleOutput | Output scaling parameter dictionary | no  |
        calcFrequency | Calculation/coupling frequency    | no  | 1
    \endtable

    \heading Parameters for optional communication/scaleInput dictionary
    \table
        Property    | Description                       | Required | Default
        length      | Scaling for input positions       | no | 1
    \endtable

    \heading Parameters for optional communication/scaleOutput dictionary
    \table
        Property    | Description                       | Required | Default
        length      | Scaling for output positions      | no  | 1
        force       | Scaling for force                 | no  | 1
        moment      | Scaling for moment                | no  | 1
    \endtable

    \heading Parameters for optional forces dictionary
    \table
        Property    | Description                           | Required | Default
        p           | Name of the pressure field            | no  | p
        pRef        | Reference pressure in Pa              | no  | 0
        rhoRef      | Reference density for incompressible  | no  | 1
    \endtable

SourceFiles
    lumpedPointMovement.C

\*---------------------------------------------------------------------------*/

#ifndef lumpedPointMovement_H
#define lumpedPointMovement_H

#include "dictionary.H"
#include "scalarList.H"
#include "primitiveFields.H"
#include "IOobject.H"
#include "tmp.H"
#include "Tuple2.H"
#include "HashPtrTable.H"
#include "externalFileCoupler.H"
#include "lumpedPointController.H"
#include "lumpedPointInterpolator.H"
#include "lumpedPointState.H"
#include "Enum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class polyMesh;
class polyPatch;
class pointPatch;

/*---------------------------------------------------------------------------*\
                     Class lumpedPointMovement Declaration
\*---------------------------------------------------------------------------*/

class lumpedPointMovement
{
public:

    // Data Types

        //- Output format types
        enum class outputFormatType
        {
            PLAIN,      //!< "plain" is a simple ASCII format
            DICTIONARY  //!< "dictionary" is the OpenFOAM dictionary format
        };

        //- Output format types
        enum scalingType
        {
            LENGTH = 0, //!< The "length" scaling
            FORCE,      //!< The "force" scaling
            MOMENT      //!< The "moment" scaling
        };


    // Static Data

        //- Names for the output format types
        static const Enum<outputFormatType> formatNames;

        //- Names for the scaling types
        static const Enum<scalingType> scalingNames;


private:

    // Private Class

        //- The controller names and faceCentre mapping for a patch
        struct patchControl
        {
            wordList names_;
            labelList faceToPoint_;
            List<lumpedPointInterpolator> interp_;
        };


    // Private Data

        //- The offset for lumped points, used on input.
        point origin_;

        //- The initial state of positions with zero rotations.
        lumpedPointState state0_;

        //- The current state (positions, rotations)
        //  Eg, as response from external application
        lumpedPointState state_;

        //- The original point ids (optional information)
        labelList originalIds_;

        //- Connectivity for the controllers
        HashPtrTable<lumpedPointController> controllers_;

        //- The controller names and faceCentre mapping for patches
        Map<patchControl> patchControls_;

        //- The relaxation factor when moving points (default: 1)
        scalar relax_;

        //- Optional owner information (patch owner)
        label ownerId_;

        //- Dictionary of controls for force calculation
        dictionary forcesDict_;

        //- Communication control
        externalFileCoupler coupler_;

        //- File IO names
        word inputName_;
        word outputName_;
        word logName_;

        //- The input format for points, rotations
        lumpedPointState::inputFormatType inputFormat_;

        //- The output format for forces, moments
        outputFormatType outputFormat_;

        //- Scale factors for input/output files (optional)
        FixedList<scalar, 1> scaleInput_;
        FixedList<scalar, 3> scaleOutput_;

        //- Calculation frequency
        label calcFrequency_;

        //- The last timeIndex when coupling was triggered
        mutable label lastTrigger_;


    // Private Member Functions

        //- No copy construct
        lumpedPointMovement(const lumpedPointMovement&) = delete;

        //- No copy assignment
        void operator=(const lumpedPointMovement&) = delete;


public:

    // Static Data Members

        //- Debug switch
        static int debug;

        //- The canonical name ("lumpedPointMovement") for the dictionary
        static const word canonicalName;


    // Constructors

        //- Default construct
        lumpedPointMovement();

        //- Construct from dictionary, optionally with some owner information
        explicit lumpedPointMovement(const dictionary& dict, label ownerId=-1);


    //- Destructor
    virtual ~lumpedPointMovement() = default;


    // Member Functions

        //- Update settings from dictionary
        void readDict(const dictionary& dict);


        //- If no number of lumped points (locations) were specified
        inline bool empty() const;

        //- The number of lumped points (number of locations)
        inline label size() const;

        //- An owner Id, if needed for bookkeeping purposes
        inline label ownerId() const;

        //- Change the owner id, if needed for bookkeeping purposes
        inline void ownerId(label id);

        //- Communication control
        inline const externalFileCoupler& coupler() const;

        //- Communication control
        inline externalFileCoupler& coupler();

        //- Check if coupling is pending (according to the calcFrequency)
        bool couplingPending(const label timeIndex) const;

        //- Register that coupling is completed at this calcFrequency
        void couplingCompleted(const label timeIndex) const;

        //- The initial state (positions/rotations)
        inline const lumpedPointState& state0() const;

        //- The current state (positions/rotations)
        inline const lumpedPointState& state() const;

        //- The offset for lumped points, used on input.
        inline const point& origin() const;

        //- Scale the lumped points (on input).
        inline void scalePoints(lumpedPointState& state) const;

        //- The relaxation factor when changing states
        inline scalar relax() const;

        //- The relaxation factor when changing states
        inline scalar& relax();

        //- The input (state) file name
        inline const word& inputName() const;

        //- The output (forces) file name
        inline const word& outputName() const;

        //- The log file name
        inline const word& logName() const;

        //- The input (state) file format
        inline lumpedPointState::inputFormatType inputFormat() const;

        //- The output (forces) file format
        inline lumpedPointMovement::outputFormatType outputFormat() const;

        //- The Euler-angle rotation order
        inline quaternion::eulerOrder rotationOrder() const;

        //-  Rotation angles in degrees
        inline bool degrees() const;

        //- Check if patch control exists for specified patch
        inline bool hasPatchControl(const label patchIndex) const;

        //- Check if patch control exists for specified patch
        bool hasPatchControl(const polyPatch& pp) const;

        //- Check if patch control exists for specified patch
        bool hasInterpolator(const pointPatch& fpatch) const;

        //- Check if patch control exists for specified patch
        void checkPatchControl(const polyPatch& pp) const;

        //- Define pressure-zones mapping for faces in the specified patches.
        //  The face centres are compared to the controller points,
        //
        //  \param pp       The patch with a control
        //  \param ctrlNames The patch ids to be included in the mapping
        //  \param points0  The initial mesh points, prior to movement
        void setPatchControl
        (
            const polyPatch& pp,
            const wordList& ctrlNames,
            const pointField& points0
        );


        //- Define pressure-zones mapping for faces in the specified patches.
        //  The face centres are compared to the controller points,
        //
        //  \param mesh     The volume mesh reference
        //  \param patchIds The patch ids to be included in the mapping
        //  \param points0  The initial mesh points, prior to movement
        void setMapping
        (
            const polyMesh& mesh,
            const labelUList& patchIds,
            const pointField& points0
        );


        //- Check if patch control exists for specified patch
        void setInterpolator
        (
            const pointPatch& fpatch,
            const pointField& points0
        );

        //- True if the pressure-zones mapping has already been performed
        inline bool hasMapping() const;

        //- Check if patch interpolator exists for specified patch
        inline bool hasInterpolator(const label patchIndex) const;

        //- The areas for each pressure-zone.
        List<scalar> areas(const polyMesh& pmesh) const;

        //- The forces and moments acting on each pressure-zone.
        //  The zones must be previously defined via setMapping.
        bool forcesAndMoments
        (
            const polyMesh& pmesh,
            List<vector>& forces,
            List<vector>& moments
        ) const;


        //- Displace points according to the current state
        tmp<pointField> pointsDisplacement
        (
            const pointPatch& fpatch,
            const pointField& points0
        ) const;

        //- Displace points according to specified state
        tmp<pointField> pointsDisplacement
        (
            const lumpedPointState& state,
            const pointPatch& fpatch,
            const pointField& points0
        ) const;

        //- The points absolute position according to specified state
        tmp<pointField> pointsPosition
        (
            const lumpedPointState& state,
            const pointPatch& fpatch,
            const pointField& points0
        ) const;


        //- Write axis, locations, division as a dictionary
        void writeDict(Ostream& os) const;

        //- Write points, forces, moments. Only call from the master process
        bool writeData
        (
            Ostream& os,
            const UList<vector>& forces,
            const UList<vector>& moments,
            const outputFormatType fmt = outputFormatType::PLAIN,
            const Tuple2<scalar, scalar>* timesWritten = nullptr
        ) const;

        //- Write points, forces, moments
        bool writeData
        (
            const UList<vector>& forces,
            const UList<vector>& moments = List<vector>(),
            const Tuple2<scalar, scalar>* timesWritten = nullptr
        ) const;

        //- Read state from file, applying relaxation as requested
        bool readState();

        //- Write state as VTK PolyData format.
        void writeStateVTP
        (
            const lumpedPointState& state,
            const fileName& file
        ) const;

        //- Write state as VTK PolyData format.
        void writeStateVTP(const fileName& file) const;

        //- Write forces on points as VTK PolyData format.
        void writeForcesAndMomentsVTP
        (
            const fileName& file,
            const UList<vector>& forces,
            const UList<vector>& moments
        ) const;


        //- Write pressure-zones geometry, write as VTK PolyData format.
        void writeZonesVTP
        (
            const fileName& file,
            const polyMesh& mesh,
            const pointField& points0
        ) const;

        //- Write displaced geometry according to the current state,
        //  write as VTK PolyData format.
        void writeVTP
        (
            const fileName& file,
            const polyMesh& mesh,
            const pointField& points0
        ) const;

        //- Write displaced geometry according to the specified state,
        //  write as VTK PolyData format.
        void writeVTP
        (
            const fileName& file,
            const lumpedPointState& state,
            const polyMesh& mesh,
            const pointField& points0
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "lumpedPointMovementI.H"

#endif

// ************************************************************************* //
