/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2024 OpenCFD Ltd.
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
    Foam::zoneBlended

Group
    grpFvSurfaceInterpolationSchemes

Description
    Multi-faceZone based blending differencing scheme.

    Schemes are set in dictonary format according to:

    \verbatim
    divSchemes
    {
        .
        .
        div(phi,U)      Gauss zoneBlended
        {
            default         defaultScheme;
            faceZone1       scheme1;
            faceZone2       scheme2;
            ...
            faceZoneN       schemeN;
        }
        .
        .
    }
    \endverbatim

    The \c default entry specifies the background scheme; additional schemes
    can be set per \c faceZone, e.g. \c scheme1 is applied to \c facZone1, 
    \c scheme2 is applied to \c facZone2 etc. 


Usage
    Example of the \c zoneBlended scheme to use \c linearUpwind as the
    background scheme and \c upwind in \c faceZone1:

    \verbatim
    divSchemes
    {
        .
        .
        div(phi,U)      Gauss zoneBlended 
        {
            default         linearUpwind grad(U); 
            faceZone1       upwind;
        };
        .
        .
    }
    \endverbatim

SourceFiles
    zoneBlended.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_zoneBlended_H
#define Foam_zoneBlended_H

#include "surfaceInterpolationScheme.H"
#include "blendedSchemeBase.H"
#include "surfaceInterpolate.H"
#include "UIndirectList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class zoneBlended Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class zoneBlended
:
    public surfaceInterpolationScheme<Type>,
    public blendedSchemeBase<Type>
{
    using SurfaceField = GeometricField<Type, fvsPatchField, surfaceMesh>;
    using VolumeField = GeometricField<Type, fvPatchField, volMesh>;

    // Private data

        //- Face zones
        wordList zoneNames_;

        //- Schemes
        //  Note: 0 index used to describe default/background scheme
        List<tmp<surfaceInterpolationScheme<Type>>> schemePtrs_;

        //- Corrected flag - true if any of the schemes has corrected() set
        bool corrected_;


    // Private Member Functions

        //- Set the lists of face zone names and schemes
        void setSchemes(const fvMesh& mesh, const dictionary& dict);

        //- Set the lists of face zone names and schemes
        void setSchemes
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            const dictionary& dict
        );

        //- Retrieve a scheme from the list
        const surfaceInterpolationScheme<Type>& scheme
        (
            const label schemei
        ) const
        {
            return schemePtrs_[schemei]();
        }

        //- Set destination values from source values for a face zone
        template<class FieldType>
        void setFaceZoneValues
        (
            FieldType& dest,
            const FieldType& src,
            const faceZone& fz
        ) const;

        //- Set destination values to zero for a face zone
        template<class FieldType>
        void zeroFaceZoneValues(FieldType& dest, const faceZone& fz) const;

        //- No copy construct
        zoneBlended(const zoneBlended&) = delete;

        //- No copy assignment
        void operator=(const zoneBlended&) = delete;


public:

    //- Runtime type information
    TypeName("zoneBlended");


    // Constructors

        //- Construct from mesh and Istream.
        //  The name of the flux field is read from the Istream and looked-up
        //  from the mesh objectRegistry
        zoneBlended(const fvMesh& mesh, Istream& is)
        :
            surfaceInterpolationScheme<Type>(mesh),
            zoneNames_(),
            schemePtrs_(),
            corrected_(false)
        {
            const dictionary dict(is);

            setSchemes(mesh, dict);
        }


        //- Construct from mesh, faceFlux and Istream
        zoneBlended
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            Istream& is
        )
        :
            surfaceInterpolationScheme<Type>(mesh),
            zoneNames_(),
            schemePtrs_(),
            corrected_(false)
        {
            const dictionary dict(is);

            setSchemes(mesh, faceFlux, dict);
        }


    // Member Functions

        //- Return the face-based blending factor
        virtual tmp<surfaceScalarField> blendingFactor
        (
            const VolumeField& vf
        ) const
        {
            const auto& mesh = vf.mesh();
            auto tbf = surfaceScalarField::New("blendingFactor", mesh, dimless);
            auto& bf = tbf.ref();
            auto& bbf = bf.boundaryFieldRef();
            bf = 0.0;

            const auto& pbm = mesh.boundaryMesh();
            const auto& zones = mesh.faceZones();

            // Use blending factor to show different zones
            for (label zonei=1; zonei<zoneNames_.size(); ++zonei)
            {
                const word& name = zoneNames_[zonei];
                const auto& fz = zones[name];

                for (const label facei : fz)
                {
                    if (mesh.isInternalFace(facei))
                    {
                        bf[facei] = zonei;
                    }
                    else
                    {
                        const labelPair pf = pbm.whichPatchFace(facei);
                        auto& pbf = bbf[pf.first()];
                        if (pbf.size())
                        {
                            pbf[pf.second()] = zonei;
                        }
                    }
                }
            }

            return tbf;
        }


        //- Return the interpolation weighting factors
        tmp<surfaceScalarField> weights(const VolumeField& vf) const
        {
            const auto& mesh = vf.mesh();
            auto tweights =
                surfaceScalarField::New("weights", vf.mesh(), dimless);
            auto& weights = tweights.ref();

            // Set default scheme weights
            weights = this->scheme(0).weights(vf);

            // Set face zone weights
            const auto& zones = mesh.faceZones();

            for (label schemei=1; schemei<schemePtrs_.size(); ++schemei)
            {
                const auto& scheme = this->scheme(schemei);

                auto tschemeWeights = scheme.weights(vf);
                const auto& schemeWeights = tschemeWeights();
                const auto& fz = zones[zoneNames_[schemei]];

                setFaceZoneValues(weights, schemeWeights, fz);
            }

            return tweights;
        }


        //- Return the face-interpolate of the given cell field
        //  with explicit correction
        tmp<SurfaceField> interpolate(const VolumeField& vf) const
        {
            return
                surfaceInterpolationScheme<Type>::interpolate
                (
                    vf,
                    weights(vf)
                );
        }


        //- Return true if this scheme uses an explicit correction
        virtual bool corrected() const
        {
            return corrected_;
        }


        //- Return the explicit correction to the face-interpolate
        //- for the given field
        virtual tmp<SurfaceField> correction
        (
            const VolumeField& vf
        ) const
        {
            const auto& mesh = vf.mesh();
            auto tcorr =
                SurfaceField::New("correction", vf.mesh(), vf.dimensions());
            auto& corr = tcorr.ref();
            corr = dimensioned<Type>(vf.dimensions(), Zero);

            // Set default scheme correction
            const auto& scheme0 = this->scheme(0);
            if (scheme0.corrected())
            {
                corr = scheme0.correction(vf);
            }

            // Only the default scheme exists - can exit early
            if (schemePtrs_.size() == 1) return tcorr;

            // Set correction field in faceZones
            const auto& zones = mesh.faceZones();

            for (label schemei=1; schemei<schemePtrs_.size(); ++schemei)
            {
                const auto& scheme = this->scheme(schemei);

                if (scheme.corrected())
                {
                    auto tschemeCorr = scheme.correction(vf);
                    const auto& schemeCorr = tschemeCorr();
                    const auto& fz = zones[zoneNames_[schemei]];

                    setFaceZoneValues(corr, schemeCorr, fz);
                }
                else
                {
                    if (scheme0.corrected())
                    {
                        // Remove correction from base scheme face zone faces
                        const auto& fz = zones[zoneNames_[schemei]];
                        zeroFaceZoneValues(corr, fz);
                    }
                }
            }

            return tcorr;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "zoneBlendedTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
