/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2015-2024 OpenCFD Ltd.
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
    Foam::sampledSurfaces

Description
    Set of surfaces to sample.

    The write() method is used to sample and write files.

    Example of function object specification:

    \verbatim
    surfaces
    {
        type    surfaces;
        libs    (sampling);

        // Write at same frequency as fields
        writeControl    writeTime;
        writeInterval   1;

        // Fields to be sampled
        fields          (p U);

        // Scheme to obtain face centre value
        sampleScheme    cell;

        // Scheme to obtain node values
        // (only used if interpolate=true for the surfaces below)
        interpolationScheme cell;

        // Optional: registry storage
        store           true

        // Output surface format
        surfaceFormat   vtk;

        formatOptions
        {
            default
            {
                verbose     true;
            }
            vtk
            {
                precision   10;
            }
        }

        surfaces
        {
            f0surf
            {
                type        meshedSurface;
                surface     f0surf.obj;
                source      cells;

                // Optional: keep original regions
                keepIds     true;

                // Optional: generate values on points instead of faces
                interpolate true;

                // Optional: alternative output type
                surfaceFormat   ensight;

                // Optional: registry storage
                store       true
            }
        }
    }
    \endverbatim

    Entries:
    \table
        Property | Description                              | Required | Default
        type     | Type-name: surfaces                      | yes |
        surfaces | Dictionary or list of sample surfaces    | expected |
        fields   | word/regex list of fields to sample      | yes |
        sampleScheme | scheme to obtain face centre value   | no  | cell
        interpolationScheme | scheme to obtain node values  | no  | cellPoint
        surfaceFormat | output surface format               | yes |
        formatOptions | dictionary of format options        | no  |
        sampleOnExecute | Sample (store) on execution as well | no | false
        store    | Store surface/fields on registry         | no  | false
    \endtable

    Additional per-surface entries:
    \table
        Property | Description                              | Required | Default
        store    | Store surface/fields on registry         | no  |
        surfaceFormat | output surface format               | no  |
        formatOptions | dictionary of format options        | no  |
    \endtable

Note
    The interpolationScheme is only used if interpolate=true is used by any
    of the surfaces.

SourceFiles
    sampledSurfaces.C
    sampledSurfacesTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_sampledSurfaces_H
#define Foam_sampledSurfaces_H

#include "fvMeshFunctionObject.H"
#include "sampledSurface.H"
#include "surfaceWriter.H"
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "IOobjectList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class polySurface;

/*---------------------------------------------------------------------------*\
                       Class sampledSurfaces Declaration
\*---------------------------------------------------------------------------*/

class sampledSurfaces
:
    public functionObjects::fvMeshFunctionObject,
    public PtrList<sampledSurface>
{
    // Static Data Members

        //- Tolerance for merging points (fraction of mesh bounding box)
        static scalar mergeTol_;

        //- Local control for sampling actions
        enum sampleActionType : unsigned
        {
            ACTION_NONE  = 0,
            ACTION_WRITE = 0x1,
            ACTION_STORE = 0x2,
            ACTION_ALL = 0xF
        };


    // Private Data

        //- Load fields from files (not from objectRegistry)
        const bool loadFromFiles_;

        //- Output verbosity
        bool verbose_;

        //- Perform sample/store actions on execute as well
        bool onExecute_;

        //- Output path
        fileName outputPath_;


    // Read from dictionary

        //- Names of fields to sample
        wordRes fieldSelection_;

        //- Interpolation/sample scheme to obtain face values
        word sampleFaceScheme_;

        //- Interpolation/sample scheme to obtain node values
        word sampleNodeScheme_;


    // Output control

        //- Surface writers (one per surface)
        PtrList<surfaceWriter> writers_;

        //- Per-surface selection of store/write actions
        List<unsigned> actions_;

        //- Per-surface status of the surfaces
        List<bool> hasFaces_;


    // Private Member Functions

        //- Return the surfaces
        const PtrList<sampledSurface>& surfaces() const noexcept
        {
            return *this;
        }

        //- Return the surfaces
        PtrList<sampledSurface>& surfaces() noexcept
        {
            return *this;
        }

        //- A new surfaceWriter, with per-surface formatOptions
        static autoPtr<surfaceWriter> newWriter
        (
            word writerType,
            const dictionary& topDict,
            const dictionary& surfDict
        );


        //- Perform sampling action with store/write
        bool performAction(unsigned request);

        //- Count selected/sampled fields per surface
        //  Returns empty IOobjectList if loadFromFiles_ is not active
        IOobjectList preCheckFields();

        //- Write sampled fieldName on surface and on outputDir path
        template<class Type>
        void writeSurface
        (
            surfaceWriter& writer,
            const Field<Type>& values,
            const word& fieldName
        );

        //- Sample and store/write a specific volume field
        template<class Type>
        void performAction(const VolumeField<Type>& fld, unsigned request);

        //- Sample and store/write a specific surface field
        template<class Type>
        void performAction(const SurfaceField<Type>& fld, unsigned request);

        //- Sample and write all applicable sampled fields
        //  Only uses IOobjectList when loadFromFiles_ is active
        template<class GeoField>
        void performAction
        (
            const IOobjectList& objects,
            unsigned request
        );


        //- Put surface onto registry
        void storeRegistrySurface(const sampledSurface& s);

        //- Store sampled field onto surface registry (if surface exists)
        template<class Type, class GeoMeshType>
        void storeRegistryField
        (
            const sampledSurface& s,
            const word& fieldName,
            const dimensionSet& dims,
            Field<Type>&& values
        );

        //- Test surfaces for condition.
        //- Like std::any_of() but without any iterator requirements
        template<class Container, class Predicate>
        static bool testAny(const Container& items, const Predicate& pred);

        //- Do any surfaces need an update?
        virtual bool needsUpdate() const;

        //- Mark the surfaces as needing an update.
        //  Return false if all surfaces were already marked as expired.
        //  Optionally force expire, even if a surface has been marked as
        //  invariant.
        virtual bool expire(const bool force=false);

        //- Update the surfaces as required.
        //  Return false if no surfaces required an update.
        virtual bool update();

        //- No copy construct
        sampledSurfaces(const sampledSurfaces&) = delete;

        //- No copy assignment
        void operator=(const sampledSurfaces&) = delete;


public:

    //- Runtime type information
    TypeName("surfaces");


    // Constructors

        //- Construct from Time and dictionary
        sampledSurfaces
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Construct for given objectRegistry and dictionary
        //  allow the possibility to load fields from files
        sampledSurfaces
        (
            const word& name,
            const objectRegistry& obr,
            const dictionary& dict,
            const bool loadFromFiles = false
        );


    //- Destructor
    virtual ~sampledSurfaces() = default;


    // Member Functions

        //- Enable/disable verbose output
        //  \return old value
        bool verbose(bool on) noexcept;

        //- Return names of fields to sample
        const wordRes& fieldNames() const noexcept { return fieldSelection_; }

        //- Read the sampledSurfaces dictionary
        virtual bool read(const dictionary& dict);

        //- Sample and store if the sampleOnExecute is enabled.
        virtual bool execute();

        //- Sample and write
        virtual bool write();

        //- Update for changes of mesh - expires the surfaces
        virtual void updateMesh(const mapPolyMesh& mpm);

        //- Update for mesh point-motion - expires the surfaces
        virtual void movePoints(const polyMesh& mesh);

        //- Update for changes of mesh due to readUpdate - expires the surfaces
        virtual void readUpdate(const polyMesh::readUpdateState state);

        //- Get merge tolerance
        static scalar mergeTol() noexcept;

        //- Set merge tolerance and return old value
        static scalar mergeTol(scalar tol) noexcept;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "sampledSurfacesTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
