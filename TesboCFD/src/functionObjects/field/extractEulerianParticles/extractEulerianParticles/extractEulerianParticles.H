/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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
    Foam::functionObjects::extractEulerianParticles

Group
    grpFieldFunctionObjects

Description
    Generates particle size information from Eulerian calculations, e.g. \c VoF.

    Operands:
    \table
      Operand           | Type        | Location
      input             | -           | -
      output file       | -           | -
      output field 1    | scalarField | $OUTPUT/d
      output field 2    | scalarField | $OUTPUT/soi
      output field 3    | labelField  | $OUTPUT/tag
      output field 4    | vectorField | $OUTPUT/U
    \endtable

    where \c $OUTPUT=$FOAM_CASE/\<time\>/lagrangian/eulerianParticleCloud.

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    extractEulerianParticles1
    {
        // Mandatory entries (unmodifiable)
        type            extractEulerianParticles;
        libs            (fieldFunctionObjects);

        // Mandatory entries (runtime modifiable)
        faceZone        f0;
        alpha           alpha.water;

        // Optional entries (runtime modifiable)
        alphaThreshold  0.1;
        nLocations      0;
        U               U;
        rho             rho;
        phi             phi;
        minDiameter     1e-30;
        maxDiameter     1e30;

        // Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property   | Description                          | Type | Req'd | Dflt
      type       | Type name: extractEulerianParticles  | word |  yes  | -
      libs       | Library name: fieldFunctionObjects   | word |  yes  | -
      faceZone   | Name of faceZone used as collection surface | word  | yes | -
      alpha      | Name of phase indicator field        | word |  yes  | -
      alphaThreshold | Threshold for alpha field        | scalar | no  | 0.1
      nLocations | Number of injection bins to generate | label  | no  | 0
      U          | Name of velocity field               | word   | no  | U
      rho        | Name of density field                | word   | no  | rho
      phi        | Name of flux field                   | word   | no  | phi
      minDiameter | Minimum diameter                    | scalar | no  | SMALL
      maxDiameter | Maximum diameter                    | scalar | no  | GREAT
    \endtable

    The inherited entries are elaborated in:
     - \link functionObject.H \endlink
     - \link writeFile.H \endlink

    Usage by the \c postProcess utility is not available.

See also
    - Foam::functionObject
    - Foam::functionObjects::fvMeshFunctionObject
    - Foam::functionObjects::writeFile
    - Foam::eulerianParticle
    - ExtendedCodeGuide::functionObjects::field::extractEulerianParticles

SourceFiles
    extractEulerianParticles.C
    extractEulerianParticlesTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_extractEulerianParticles_H
#define functionObjects_extractEulerianParticles_H

#include "fvMeshFunctionObject.H"
#include "writeFile.H"
#include "runTimeSelectionTables.H"
#include "polyMesh.H"
#include "surfaceFieldsFwd.H"
#include "globalIndex.H"
#include "eulerianParticle.H"
#include "injectedParticleCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
/*---------------------------------------------------------------------------*\
           Class extractEulerianParticlesFunctionObject Declaration
\*---------------------------------------------------------------------------*/

class extractEulerianParticles
:
    public fvMeshFunctionObject,
    public writeFile
{
protected:

    // Protected Data

        //- Storage for collected particles
        injectedParticleCloud cloud_;


        // faceZone info

            //- Name of faceZone to sample
            word faceZoneName_;

            //- Index of the faceZone
            label zoneID_;

            //- Patch indices where faceZone face intersect patch
            labelList patchIDs_;

            //- Patch face indices where faceZone face intersect patch
            labelList patchFaceIDs_;


        // Field names

            //- Name of phase fraction field
            word alphaName_;

            //- Value of phase fraction used to identify particle boundaries
            scalar alphaThreshold_;

            //- Name of the velocity field, default = U
            word UName_;

            //- Name of the density field, default = rho
            word rhoName_;

            //- Name of the flux field, default ="rho"
            word phiName_;


        // Agglomeration

            //- Number of sample locations to generate
            label nInjectorLocations_;

            //- Agglomeration addressing from fine to coarse (local proc only)
            labelList fineToCoarseAddr_;

            //- Global coarse face addressing
            globalIndex globalCoarseFaces_;


        // Particle collection info

            //- Region indices in faceZone faces from last iteration
            labelList regions0_;

            //- Particle properties (partial, being accumulated)
            List<eulerianParticle> particles_;

            //- Map from region to index in particles_ list
            Map<label> regionToParticleMap_;

            //- Minimum diameter (optional)
            //  Can be used to filter out 'small' particles
            scalar minDiameter_;

            //- Maximum diameter (optional)
            //  Can be used to filter out 'large' particles
            scalar maxDiameter_;


        // Statistics

            //- Total number of collected particles
            label nCollectedParticles_;

            //- Total collected volume [m3]
            scalar collectedVolume_;

            //- Total number of discarded particles
            label nDiscardedParticles_;

            //- Total discarded volume [m3]
            scalar discardedVolume_;


    // Protected Member Functions

        //- Check that the faceZone is valid
        virtual void checkFaceZone();

        //- Initialise the particle collection bins
        virtual void initialiseBins();

        //- Return the volumetric flux
        virtual tmp<surfaceScalarField> phiU() const;

        //- Set the blocked faces, i.e. where alpha > alpha threshold value
        virtual void setBlockedFaces
        (
            const surfaceScalarField& alphaf,
            const faceZone& fz,
            boolList& blockedFaces
        );

        //- Calculate the addressing between regions between iterations
        //- Returns the number of active regions (particles)
        virtual void calculateAddressing
        (
            const label nRegionsNew,
            const scalar time,
            labelList& regionFaceIDs
        );

        //- Collect particles that have passed through the faceZone
        virtual void collectParticle
        (
            const scalar time,
            const label regioni
        );

        //- Process latest region information
        virtual void accumulateParticleInfo
        (
            const surfaceScalarField& alphaf,
            const surfaceScalarField& phi,
            const labelList& regionFaceIDs,
            const faceZone& fz
        );

        template<class Type>
        inline Type faceValue
        (
            const GeometricField<Type, fvsPatchField, surfaceMesh>& field,
            const label localFaceI,
            const label globalFaceI
        ) const;


public:

    // Static data members

        //- Runtime type information
        TypeName("extractEulerianParticles");


    // Constructors

        //- Construct from components
        extractEulerianParticles
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- No copy construct
        extractEulerianParticles(const extractEulerianParticles&) = delete;

        //- No copy assignment
        void operator=(const extractEulerianParticles&) = delete;


    //- Destructor
    virtual ~extractEulerianParticles() = default;


    // Member Functions

        //- Read the field min/max data
        virtual bool read(const dictionary&);

        //- Execute
        virtual bool execute();

        //- Write
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "extractEulerianParticlesTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
