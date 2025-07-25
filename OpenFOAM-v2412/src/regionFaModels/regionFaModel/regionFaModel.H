/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2022 OpenCFD Ltd.
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
    Foam::regionModels::regionFaModel

Description
    Base class for area region models.

Usage
    Example of the model specification:
    \verbatim
    <patchName>
    {
        // Mandatory entries
        region          <word>;
        active          <bool>;

        <model>Coeffs
        {
            // subdictionary entries
        }

        // Optional entries
        infoOutput      <bool>;

        // Inherited entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property   | Description                         | Type  | Reqd | Deflt
      region     | Name of operand region              | word  | yes  | -
      active     | Flag to activate the model          | bool  | yes  | -
      infoOutput | Flag to activate information output | bool  | no   | false
    \endtable

SourceFiles
    regionFaModelI.H
    regionFaModel.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_regionFaModel_H
#define Foam_regionFaModel_H

#include "volMesh.H"
#include "volSurfaceMapping.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{

/*---------------------------------------------------------------------------*\
                         Class regionFaModel Declaration
\*---------------------------------------------------------------------------*/

class regionFaModel
:
    public IOdictionary
{
    // Private Member Functions

        //- Construct region mesh and fields
        void constructMeshObjects();

        //- Initialise the region
        void initialise();

        //- Read control parameters from dictionary
        bool init(const dictionary& dict);


protected:

    // Protected Data

        //- Reference to the primary mesh database
        const fvMesh& primaryMesh_;

        //- Reference to the time database
        const Time& time_;

        //- Active flag
        Switch active_;

        //- Active information output
        Switch infoOutput_;

        //- Model name
        const word modelName_;

        //- Pointer to the region mesh database
        autoPtr<faMesh> regionMeshPtr_;

        //- Model coefficients dictionary
        dictionary coeffs_;

        //- Dictionary of output properties
        autoPtr<IOdictionary> outputPropertiesPtr_;

        //- Volume/surface mapping
        autoPtr<volSurfaceMapping> vsmPtr_;

        //- Region name
        word regionName_;


public:

    //- Runtime type information
    TypeName("regionFaModel");


    //- Default name regionFaModel
    static const word regionFaModelName;

    // Constructors

        //- Construct from mesh and name and dict
        regionFaModel
        (
            const fvMesh& mesh,
            const word& regionType,
            const word& modelName,
            const dictionary& dict,
            bool readFields = true
        );

        //- No copy construct
        regionFaModel(const regionFaModel&) = delete;

        //- No copy assignment
        void operator=(const regionFaModel&) = delete;


    //- Destructor
    virtual ~regionFaModel() = default;


    // Member Functions

    // Access

        //- Return the reference to the primary mesh database
        const fvMesh& primaryMesh() const noexcept { return primaryMesh_; }

        //- Return the reference to the time database
        const Time& time() const noexcept { return time_; }

        //- Return the active flag
        Switch active() const noexcept { return active_; }

        //- Return the information flag
        Switch infoOutput() const noexcept { return infoOutput_; }

        //- Return the model name
        const word& modelName() const noexcept { return modelName_; }

        //- Return the region mesh database
        inline const faMesh& regionMesh() const;

        //- Return the region mesh database for manipulation
        inline faMesh& regionMesh();

        //- Return the model coefficients dictionary
        const dictionary& coeffs() const noexcept { return coeffs_; }

        //- Return const access to the output properties dictionary
        inline const IOdictionary& outputProperties() const;

        //- Return output properties dictionary
        inline IOdictionary& outputProperties();

        //- Return the solution dictionary
        inline const dictionary& solution() const;


    // Addressing

        //- List of patch IDs on the primary region coupled to this region
        inline const labelList& primaryPatchIDs() const;

        //- True if patchi on the primary region is coupled to this region
        inline bool isRegionPatch(const label patchi) const;


    // Helper Functions

        //- Return mapping between surface and volume fields
        const volSurfaceMapping& vsm() const;


    // Evolution

        //- Main driver routing to evolve the region - calls other evolves
        virtual void evolve();

        //- Pre-evolve region
        virtual void preEvolveRegion();

        //- Evolve the region
        virtual void evolveRegion();

        //- Post-evolve region
        virtual void postEvolveRegion();

        //- Courant number of the region
        virtual scalar CourantNumber() const;


    // IO

        //- Provide some feedback
        virtual void info() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace regionFaModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "regionFaModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
