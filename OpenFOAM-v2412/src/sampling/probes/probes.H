/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2022 OpenCFD Ltd.
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
    Foam::probes

Group
    grpUtilitiesFunctionObjects

Description
    Set of locations to sample.

    Call write() to sample and write files.

    Example of function object specification:
    \verbatim
    probes
    {
        type            probes;
        libs            (sampling);

        // Name of the directory for probe data
        name            probes;

        // Write at same frequency as fields
        writeControl    writeTime;
        writeInterval   1;

        // Fields to be probed
        fields          (U "p.*");

        // Optional: do not recalculate cells if mesh moves
        fixedLocations  false;

        // Optional: interpolation scheme to use (default is cell)
        interpolationScheme cellPoint;

        probeLocations
        (
            ( 1e-06 0 0.01 )      // at inlet
            (0.21 -0.20999 0.01)  // at outlet1
            (0.21 0.20999 0.01)   // at outlet2
            (0.21 0 0.01)         // at central block
        );

        // Optional: filter out points that haven't been found. Default
        //           is to include them (with value -VGREAT)
        includeOutOfBounds  true;
    }
    \endverbatim

    Entries:
    \table
        Property | Description                              | Required | Default
        type     | Type-name: probes                        | yes |
        probeLocations | Probe locations                    | yes |
        fields   | word/regex list of fields to sample      | yes |
        interpolationScheme | scheme to obtain values       | no  | cell
        fixedLocations | Do not recalculate cells if mesh moves | no | true
        includeOutOfBounds | Include out-of-bounds locations | no | true
        sampleOnExecute | Sample on execution and store results | no | false
    \endtable

SourceFiles
    probes.C
    probesTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_probes_H
#define Foam_probes_H

#include "fvMeshFunctionObject.H"
#include "HashPtrTable.H"
#include "OFstream.H"
#include "polyMesh.H"
#include "pointField.H"
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "surfaceMesh.H"
#include "wordRes.H"
#include "IOobjectList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class Time;
class objectRegistry;
class dictionary;
class fvMesh;
class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                          Class probes Declaration
\*---------------------------------------------------------------------------*/

class probes
:
    public functionObjects::fvMeshFunctionObject,
    public pointField
{
protected:

    // Protected Classes

        //- Grouping of field names by GeometricField type
        template<class GeoField>
        struct fieldGroup : public DynamicList<word> {};


    // Data Types

        //- Local control for sampling actions
        enum sampleActionType : unsigned
        {
            ACTION_NONE  = 0,
            ACTION_WRITE = 0x1,
            ACTION_STORE = 0x2,
            ACTION_ALL = 0xF
        };


    // Protected Data

        //- Load fields from files (not from objectRegistry)
        bool loadFromFiles_;

        //- Fixed locations (default: true)
        //  Note: set to false for moving mesh calculations where locations
        //        should move with the mesh
        bool fixedLocations_;

        //- Include probes that were not found (default: true)
        bool includeOutOfBounds_;

        //- Output verbosity
        bool verbose_;

        //- Perform sample actions on execute as well
        bool onExecute_;

        //- Requested names of fields to probe
        wordRes fieldSelection_;

        //- Interpolation/sample scheme to obtain values at the points
        //  Note: only possible when fixedLocations_ is true
        word samplePointScheme_;


    // Calculated

        //- Current list of field names selected for sampling
        DynamicList<word> selectedFieldNames_;

        //- Categorized scalar/vector/tensor volume fields
        fieldGroup<volScalarField> scalarFields_;
        fieldGroup<volVectorField> vectorFields_;
        fieldGroup<volSphericalTensorField> sphericalTensorFields_;
        fieldGroup<volSymmTensorField> symmTensorFields_;
        fieldGroup<volTensorField> tensorFields_;

        //- Categorized scalar/vector/tensor surface fields
        fieldGroup<surfaceScalarField> surfaceScalarFields_;
        fieldGroup<surfaceVectorField> surfaceVectorFields_;
        fieldGroup<surfaceSphericalTensorField> surfaceSphericalTensorFields_;
        fieldGroup<surfaceSymmTensorField> surfaceSymmTensorFields_;
        fieldGroup<surfaceTensorField> surfaceTensorFields_;

        //- Cells to be probed (obtained from the locations)
        labelList elementList_;

        //- Faces to be probed
        labelList faceList_;

        //- Processor holding the cell or face (-1 if point not found
        //  on any processor)
        labelList processor_;

        //- Current open files (non-empty on master only)
        HashPtrTable<OFstream> probeFilePtrs_;

        //- Patch IDs on which the new probes are located (for patchProbes)
        labelList patchIDList_;

        //- Original probes location (only used for patchProbes)
        pointField oldPoints_;


    // Protected Member Functions

        //- Find cells and faces containing probes
        virtual void findElements(const fvMesh& mesh);

        //- Classify field types, close/open file streams
        //  \return number of fields to sample
        label prepare(unsigned request);

        //- Get from registry or load from disk
        template<class GeoField>
        tmp<GeoField> getOrLoadField(const word& fieldName) const;

        //- Store results: min/max/average/size
        template<class Type>
        void storeResults(const word& fieldName, const Field<Type>& values);


private:

    // Private Member Functions

        //- Create new streams as required
        void createProbeFiles(const wordList& fieldNames);

        //- Write field values
        template<class Type>
        void writeValues
        (
            const word& fieldName,
            const Field<Type>& values,
            const scalar timeValue
        );

        //- Sample and store/write all applicable sampled fields
        template<class GeoField>
        void performAction
        (
            const fieldGroup<GeoField>& fieldNames,  /* must be sorted */
            unsigned request
        );

        //- Perform sampling action with store/write
        bool performAction(unsigned request);

        //- No copy construct
        probes(const probes&) = delete;

        //- No copy assignment
        void operator=(const probes&) = delete;

public:

    //- Runtime type information
    TypeName("probes");


    // Constructors

        //- Construct from Time and dictionary
        probes
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict,
            const bool loadFromFiles = false,
            const bool readFields = true
        );


    //- Destructor
    virtual ~probes() = default;


    // Member Functions

        //- Enable/disable verbose output
        //  \return old value
        bool verbose(const bool on) noexcept;

        //- Return names of fields to probe
        virtual const wordRes& fieldNames() const noexcept
        {
            return fieldSelection_;
        }

        //- Return locations to probe
        virtual const pointField& probeLocations() const noexcept
        {
            return *this;
        }

        //- Return location for probe i
        virtual const point& probe(const label i) const
        {
            return operator[](i);
        }

        //- Cells to be probed (obtained from the locations)
        const labelList& elements() const noexcept
        {
            return elementList_;
        }

        //- Read the probes
        virtual bool read(const dictionary&);

        //- Sample and store result if the sampleOnExecute is enabled.
        virtual bool execute();

        //- Sample and write
        virtual bool write();

        //- Update for changes of mesh
        virtual void updateMesh(const mapPolyMesh&);

        //- Update for changes of mesh
        virtual void movePoints(const polyMesh&);

        //- Update for changes of mesh due to readUpdate
        virtual void readUpdate(const polyMesh::readUpdateState state)
        {}


    // Sampling

        //- Sample a volume field at all locations
        template<class Type>
        tmp<Field<Type>> sample(const VolumeField<Type>&) const;

        //- Sample a surface field at all locations
        template<class Type>
        tmp<Field<Type>> sample(const SurfaceField<Type>&) const;

        //- Sample a volume field at all locations
        template<class Type>
        tmp<Field<Type>> sample(const word& fieldName) const;

        //- Sample a surface field at all locations
        template<class Type>
        tmp<Field<Type>> sampleSurfaceField(const word& fieldName) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "probesTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
