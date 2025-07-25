/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2024 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify i
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
    Foam::functionObjects::externalCoupled

Group
    grpFieldFunctionObjects

Description
    Provides a simple file-based communication interface for explicit coupling
    with an external application, so that data is transferred to- and from
    OpenFOAM. The data exchange employs specialised boundary conditions to
    provide either one-way or two-way coupling models.

    The coupling is through plain text files where OpenFOAM boundary data
    is read/written as one line per face (data from all processors collated):
    \verbatim
        # Patch: <patch name>
        <fld1> <fld2> .. <fldn>             //face0
        <fld1> <fld2> .. <fldn>             //face1
        ..
        <fld1> <fld2> .. <fldn>             //faceN
    \endverbatim

    where the actual entries depend on the boundary condition type:
    - mixed: value, snGrad, refValue, refGrad, valueFraction
    - externalCoupledMixed: output of writeDataMaster
    - other: value, snGrad

    These text files are located in a user specified communications directory
    which gets read/written on the master processor only.

    In the communications directory the structure will be:
    \verbatim
        <regionsName>/<patchGroup>/<fieldName>.[in|out]
    \endverbatim

    (where \c regionsName is either the name of a single region or a composite
    of multiple region names)

    At start-up, the boundary creates a lock file, i.e.:
    \verbatim
        OpenFOAM.lock
    \endverbatim

    ... to signal the external source to wait. During the function object
    execution the boundary values are written to files (one per region,
    per patch(group), per field), e.g.
    \verbatim
        <regionsName>/<patchGroup>/<fieldName>.out
    \endverbatim

    The lock file is then removed, instructing the external source to take
    control of the program execution. When ready, the external program
    should create the return values, e.g. to files
    \verbatim
        <regionsName>/<patchGroup>/<fieldName>.in
    \endverbatim

    ... and then reinstate the lock file. The function object will then
    read these values, apply them to the boundary conditions and pass
    program execution back to OpenFOAM.

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    externalCoupled1
    {
        // Mandatory entries (unmodifiable)
        type            externalCoupled;
        libs            (fieldFunctionObjects);
        commsDir        "<case>/comms";
        regions
        {
            "(region1|region0)"         // Name of region(s)
            {
                TPatchGroup             // Name of patch(group)
                {
                    readFields  (p);    // List of fields to read
                    writeFields (T);    // List of fields to write
                }
            }
        }
        initByExternal  true;

        // Optional entries (runtime modifiable)
        waitInterval    1;
        timeOut         100;
        statusDone      done;  // Any arbitrary status=... value
        calcFrequency   1;

        // Optional (inherited) entries
        ...
    }
    \endverbatim

    This reads/writes (on the master processor) the directory:
    \verbatim
        comms/region0_region1/TPatchGroup/
    \endverbatim

    with contents:
    \verbatim
        patchPoints     (collected points)
        patchFaces      (collected faces)
        p.in            (input file of p, written by external application)
        T.out           (output file of T, written by OpenFOAM)
    \endverbatim

    The patchPoints/patchFaces files denote the (collated) geometry
    which will be written if it does not exist yet or can be written as
    a preprocessing step using the createExternalCoupledPatchGeometry
    application.

    The entries comprise:
    \table
      Property     | Description                        | Type | Req'd | Dflt
      type         | Type name: externalCoupled         | word |  yes  | -
      libs         | Library name: fieldFunctionObjects | word |  yes  | -
      commsDir     | Communication directory            | word |  yes  | -
      regions      | The regions to couple              | word |  yes  | -
      initByExternal | Initialization values supplied by external app  <!--
                                                    --> | bool | yes | -
      waitInterval | Wait interval in [s]               | label | no   | 1
      timeOut      | Timeout in [s]             | label | no | 100*waitInterval
      statusDone   | Lockfile status=... on termination | word  | no   | done
      calcFrequency | Calculation frequency             | label | no   | 1
    \endtable

    The inherited entries are elaborated in:
     - \link timeFunctionObject.H \endlink

    Usage by the \c postProcess utility is not available.

See also
    - Foam::functionObject
    - Foam::functionObjects::timeFunctionObject
    - Foam::externalFileCouple
    - ExtendedCodeGuide::functionObjects::field::externalCoupled

SourceFiles
    externalCoupled.C
    externalCoupledTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_functionObjects_externalCoupled_H
#define Foam_functionObjects_externalCoupled_H

#include "timeFunctionObject.H"
#include "externalFileCoupler.H"
#include "DynamicList.H"
#include "wordRes.H"
#include "scalarField.H"
#include "UPtrList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class fvMesh;
class IFstream;

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class externalCoupled Declaration
\*---------------------------------------------------------------------------*/

class externalCoupled
:
    public functionObjects::timeFunctionObject,
    public externalFileCoupler
{
    // Private Member Data

        //- Calculation frequency
        label calcFrequency_;

        //- The last timeIndex when coupling was triggered
        label lastTrigger_;

        //- Names of (composite) regions
        DynamicList<word> regionGroupNames_;

        // Per (composite) region the names of the regions
        DynamicList<wordList> regionGroupRegions_;

        // Per (composite) region the indices of the group information
        HashTable<labelList> regionToGroups_;

        // Per group the names of the patches/patchGroups
        DynamicList<wordRe> groupNames_;

        // Per group the names of the fields to read
        DynamicList<wordList> groupReadFields_;

        // Per group the names of the fields to write
        DynamicList<wordList> groupWriteFields_;

        //- Initialised coupling
        bool initialisedCoupling_;


    // Private Member Functions

        //- Return the file path to the communications directory for the region
        static fileName groupDir
        (
            const fileName& commsDir,
            const word& regionsName,
            const wordRe& groupName
        );

        //- Read data for a single region, single field
        template<class Type>
        bool readData
        (
            const UPtrList<const fvMesh>& meshes,
            const wordRe& groupName,
            const word& fieldName
        );

        //- Write data for a single region, single field
        template<class Type>
        bool writeData
        (
            const UPtrList<const fvMesh>& meshes,
            const wordRe& groupName,
            const word& fieldName
        ) const;

        void initCoupling();

        //- Read (and distribute) scalar columns from stream.
        //- Every processor gets nRows (= patch size) of these.
        //  Note: could make its argument
        //  ISstream& but then would need additional logic to construct valid
        //  stream on all processors.
        void readColumns
        (
            const label nRows,
            const label nColumns,
            autoPtr<IFstream>& masterFilePtr,
            List<scalarField>& data
        ) const;

        //- Read (and distribute) lines from stream.
        //- Every processor gets nRows (= patch size) of these
        //- (newline terminated).
        void readLines
        (
            const label nRows,
            autoPtr<IFstream>& masterFilePtr,
            //! [out] the nRows lines read (for this rank) - newline terminated
            std::string& lines
        ) const;

        static void checkOrder(const wordList& regionNames);

        //- Perform the coupling with necessary initialization etc.
        void performCoupling();


public:

    //- Runtime type information
    TypeName("externalCoupled");

    //- Name of patch key, e.g. '// Patch:' when looking for start of patch data
    static string patchKey;

    //- Inherited variable for logging
    using functionObject::log;


    // Constructors

        //- Construct given time and dictionary
        externalCoupled
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- No copy construct
        externalCoupled(const externalCoupled&) = delete;

        //- No copy assignment
        void operator=(const externalCoupled&) = delete;


    //- Destructor
    virtual ~externalCoupled() = default;


    // Member Functions

      // Function object control

        //- Called at each ++ or += of the time-loop
        virtual bool execute();

        //- Manual execute (sub-loop or when converged)
        virtual bool execute(const label subIndex);

        //- Called when Time::run() determines that the time-loop exits
        virtual bool end();

        //- Read and set the function object if its data have changed
        virtual bool read(const dictionary& dict);

        //- Write, currently a no-op
        virtual bool write();


      // File creation, removal

        //- Write data files (all regions, all fields) from master (OpenFOAM)
        virtual void writeDataMaster() const;

        //- Read data files (all regions, all fields) on master (OpenFOAM)
        virtual void readDataMaster();

        //- Remove data files written by master (OpenFOAM)
        virtual void removeDataMaster() const;

        //- Remove data files written by slave (external code)
        virtual void removeDataSlave() const;


      // Other

        //- Create single name by appending words (in sorted order),
        //- separated by '_'
        static word compositeName(const wordList&);

        //- Write geometry for the group as region/patch
        static void writeGeometry
        (
            const UPtrList<const fvMesh>& meshes,
            const fileName& commsDir,
            const wordRe& groupName
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "externalCoupledTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
