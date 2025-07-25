/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 OpenCFD Ltd.
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
    Foam::zoneSubSet

Description
    Intermediate tool for handling \c cellZones
    for function objects (e.g. \c momentumError)
    wherein the sub-mesh option is available.

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    <functionObject>
    {
        // Mandatory and optional entries
        ...

        // Inherited entries
        cellZones         ( <cellZone1> <cellZone2> ... <cellZoneN> );
        nLayers           <label>;
    }
    \endverbatim

    where the entries mean:
    \table
      Property  | Description                             | Type  | Reqd | Deflt
      cellZones | Names of cellZones                    | wordRes | no   | -
      nLayers   | Number of cell layers around cellZones  | label | no   | 0
    \endtable

See also
  - Foam::functionObjects::momentumError
  - Foam::functionObjects::div

SourceFiles
    zoneSubSet.C
    zoneSubSetTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_zoneSubSet_H
#define Foam_zoneSubSet_H

#include "fvMeshSubset.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Detail
{

/*---------------------------------------------------------------------------*\
                           Class zoneSubSet Declaration
\*---------------------------------------------------------------------------*/

class zoneSubSet
{
    // Private Data

        //- Subsetting engine
        fvMeshSubset subsetter_;

        //- Matcher for zones
        wordRes zoneMatcher_;

        //- Number of layers around zones
        label nLayers_;

        //- The halo cells
        bitSet haloCells_;


    // Private Members Functions

        //- Initialise the sub-mesh
        bool correct();

        //- No copy construct
        zoneSubSet(const zoneSubSet&) = delete;

        //- No copy assignment
        void operator=(const zoneSubSet&) = delete;


public:

    //- Runtime type information
    TypeName("zoneSubSet");


    // Constructors

        //- Construct from components
        zoneSubSet
        (
            const fvMesh& mesh,
            const wordRes& zoneSelector,
            const label nZoneLayers = 0
        );

        //- Construct from mesh and dictionary
        zoneSubSet
        (
            const fvMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~zoneSubSet() = default;


    // Member Functions

    // Access

        //- The entire base mesh
        const fvMesh& baseMesh() const noexcept
        {
            return subsetter_.baseMesh();
        }

        //- The mesh subsetter
        const fvMeshSubset& subsetter() const
        {
            return subsetter_;
        }

        //- The mesh subsetter
        fvMeshSubset& subsetter()
        {
            return subsetter_;
        }

        //- True if sub-mesh is expected to be used
        bool useSubMesh() const noexcept
        {
            return !zoneMatcher_.empty();
        }

        //- Return the current zones selector
        const wordRes& zones() const noexcept
        {
            return zoneMatcher_;
        }

        //- Return number of layers around cell zones
        label nLayers() const noexcept
        {
            return nLayers_;
        }


    // Fields

        //- Map from the sub-mesh to original cell zones
        template<class Type>
        tmp<GeometricField<Type, fvPatchField, volMesh>> mapToZone
        (
            const GeometricField<Type, fvPatchField, volMesh>& subVolField
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Detail
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "zoneSubSetTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#endif

// ************************************************************************* //
