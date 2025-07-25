/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 Wikki Ltd
    Copyright (C) 2020 OpenCFD Ltd.
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
    Foam::skewCorrectedEdgeInterpolation

Description
    Linear/upwind blended differencing scheme

SourceFiles
    skewCorrectedEdgeInterpolationMake.C

\*---------------------------------------------------------------------------*/

#ifndef skewCorrectedEdgeInterpolation_H
#define skewCorrectedEdgeInterpolation_H

#include "edgeInterpolationScheme.H"
#include "gaussFaGrad.H"
#include "areaFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class skewCorrectedEdgeInterpolation Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class skewCorrectedEdgeInterpolation
:
    virtual public edgeInterpolationScheme<Type>
{
    // Private Data

        //- Edge-interpolation scheme
        tmp<edgeInterpolationScheme<Type>> tScheme_;


public:

    //- Runtime type information
    TypeName("skewCorrected");


    // Generated Methods

        //- No copy construct
        skewCorrectedEdgeInterpolation(const skewCorrectedEdgeInterpolation&)
            = delete;

        //- No copy assignment
        void operator=(const skewCorrectedEdgeInterpolation&) = delete;


    // Constructors

        //- Construct from Istream
        skewCorrectedEdgeInterpolation(const faMesh& mesh, Istream& is)
        :
            edgeInterpolationScheme<Type>(mesh),
            tScheme_
            (
                edgeInterpolationScheme<Type>::New(mesh, is)
            )
        {}

        //- Construct from mesh, faceFlux and blendingFactor
        skewCorrectedEdgeInterpolation
        (
            const faMesh& mesh,
            const edgeScalarField& faceFlux,
            Istream& is
        )
        :
            edgeInterpolationScheme<Type>(mesh),
            tScheme_
            (
                edgeInterpolationScheme<Type>::New(mesh, faceFlux, is)
            )
        {}


    // Member Functions

        //- Return the interpolation weighting factors
        virtual tmp<edgeScalarField> weights
        (
            const GeometricField<Type, faPatchField, areaMesh>& vf
        ) const
        {
            return tScheme_().weights(vf);
        }

        //- Return true if this scheme uses an explicit correction
        virtual bool corrected() const
        {
            return
                tScheme_().corrected() || (this->mesh()).skew();
        }

        tmp<GeometricField<Type, faePatchField, edgeMesh>>
        skewCorrection
        (
            const GeometricField<Type, faPatchField, areaMesh>& vf
        ) const
        {
            const faMesh& mesh = this->mesh();

            const edgeVectorField& scv = mesh.skewCorrectionVectors();

            tmp<GeometricField<Type, faePatchField, edgeMesh>> tsfCorr
            (
                new GeometricField<Type, faePatchField, edgeMesh>
                (
                    IOobject
                    (
                        "skewCorrected::skewCorrection(" + vf.name() + ')',
                        vf.instance(),
                        vf.db()
                    ),
                    mesh,
                    dimensioned<Type>(vf.dimensions(), Zero)
                )
            );

            GeometricField<Type, faePatchField, edgeMesh>& corr = tsfCorr.ref();

            for (direction cmpt=0; cmpt<pTraits<Type>::nComponents; ++cmpt)
            {
                corr.replace
                (
                    cmpt,
                    scv & linearEdgeInterpolation
                    <
                        typename outerProduct
                        <
                            vector,
                            typename pTraits<Type>::cmptType
                        >::type
                    >(mesh).interpolate
                    (
                        fa::gaussGrad<typename pTraits<Type>::cmptType>
                        (mesh).grad(vf.component(cmpt))
                    )
                );
            }

            return tsfCorr;
        }


        //- Return the explicit correction to the face-interpolate
        virtual tmp<GeometricField<Type, faePatchField, edgeMesh>>
        correction
        (
            const GeometricField<Type, faPatchField, areaMesh>& vf
        ) const
        {
            if
            (
                tScheme_().corrected()
             && (this->mesh()).skew()
            )
            {
                return tScheme_().correction(vf) + skewCorrection(vf);
            }
            else if (tScheme_().corrected())
            {
                return tScheme_().correction(vf);
            }
            else if ((this->mesh()).skew())
            {
                return skewCorrection(vf);
            }
            else
            {
                return
                    tmp<GeometricField<Type, faePatchField, edgeMesh>>
                    (
                        nullptr
                    );
            }
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
