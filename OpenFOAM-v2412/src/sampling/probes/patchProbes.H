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
    Foam::patchProbes

Description
    Set of locations to sample at patches

    Call write() to sample and write files.
    - find nearest location on nearest face
    - update *this with location (so header contains 'snapped' locations
    - use *this as the sampling location

    Example of function object specification:
    \verbatim
    patchProbes
    {
        type            patchProbes;
        libs            (sampling);

        // Name of the directory for probe data
        name            patchProbes;

        // Patches to sample (wildcards allowed)
        patches         (".*inl.*");

        // Write at same frequency as fields
        writeControl    writeTime;
        writeInterval   1;

        // Fields to be probed
        fields          (p U);

        // Locations to probe. These get snapped onto the nearest point
        // on the selected patches
        probeLocations
        (
            ( -100 0 0.01 )      // at inlet
        );
    }
    \endverbatim

SourceFiles
    patchProbes.C
    patchProbesTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_patchProbes_H
#define Foam_patchProbes_H

#include "probes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class patchProbes Declaration
\*---------------------------------------------------------------------------*/

class patchProbes
:
    public probes
{
protected:

    // Protected Data

        //- Patches to sample
        wordRes patchNames_;


    // Protected Member Functions

        //- Find elements containing patchProbes
        virtual void findElements(const fvMesh& mesh);  // override


private:

    // Private Member Functions

        //- Write field values
        template<class Type>
        void writeValues
        (
            const word& fieldName,
            const Field<Type>& values,
            const scalar timeValue
        );

        //- Sample and store/write applicable volume/surface fields
        template<class GeoField>
        void performAction
        (
            const fieldGroup<GeoField>& fieldNames,  /* must be sorted */
            unsigned request
        );

        //- Perform sampling action with store/write
        bool performAction(unsigned request);


        //- No copy construct
        patchProbes(const patchProbes&) = delete;

        //- No copy assignment
        void operator=(const patchProbes&) = delete;


public:

    //- Runtime type information
    TypeName("patchProbes");


    // Constructors

        //- Construct from Time and dictionary
        patchProbes
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict,
            const bool loadFromFiles = false,
            const bool readFields = true
        );

    //- Destructor
    virtual ~patchProbes() = default;


    // Member Functions

        //- Sample and store result if the sampleOnExecute is enabled.
        virtual bool execute();

        //- Sample and write
        virtual bool write();

        //- Read
        virtual bool read(const dictionary&);


    // Sampling

        //- Sample a volume field at all locations
        template<class Type>
        tmp<Field<Type>> sample(const VolumeField<Type>&) const;

        //- Sample a surface field at all locations
        template<class Type>
        tmp<Field<Type>> sample(const SurfaceField<Type>&) const;

        //- Sample a single field on all sample locations
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
    #include "patchProbesTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
