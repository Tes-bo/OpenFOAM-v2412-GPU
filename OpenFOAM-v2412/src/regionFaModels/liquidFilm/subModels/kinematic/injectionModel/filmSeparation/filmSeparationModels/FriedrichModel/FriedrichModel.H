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
    Foam::filmSeparationModels::FriedrichModel

Description
    Computes film-separation properties from sharp edges for full separation
    (Friedrich et al., 2008) and partial separation (Zhang et al., 2018).

    The governing equations for full separation (Friedrich et al., 2008):

    \f[
        F_{ratio} =
            \frac{\rho \, |\vec{u}|^2 \, h \, sin(\theta)}
            {\sigma (1 + sin(\theta)) + \rho \, \mathbf{g}\, h\,L_b cos(\theta)}
    \f]

    with:

    \f[
        L_b = 0.0388 h^{0.5} \mathrm{Re}^{0.6} \mathrm{We}^{-0.5}
    \f]

    \f[
        \mathrm{Re} = \frac{h \, |\vec{u}| \, \rho}{\mu}
    \f]

    \f[
        \mathrm{We} = \frac{h \, \rho_p (\vec{u}_p - \vec{u})^2}{2 \sigma}
    \f]

    where:

    \vartable
      F_{ratio}     | Force ratio
      h             | Film thickness
      \rho          | Film density
      \rho_p        | Primary-phase (gas) density
      \sigma        | Film surface tension
      \mu           | Film dynamic viscosity
      \mathbf{u}    | Film velocity
      \mathbf{g}    | Gravitational acceleration
      \theta        | Sharp-corner angle
      L_b           | Characteristic breakup length
    \endvartable

    The onset of film separation is triggered and the film is assumed
    fully separated when \f$F_{ratio}>1\f$; otherwise, it remains attached.


    The governing equations for partial separation (Zhang et al., 2018):

    \f[
        m_{ratio} = C_0 + C_1 \, exp\left(-\frac{F_{ratio}}{C_2}\right)
    \f]

    where:

    \vartable
      m_{ratio}     | Mass fraction of separated film mass
      C_0           | Empirical constant (0.882)
      C_1           | Empirical constant (-1.908)
      C_2           | Empirical constant (1.264)
    \endvartable

    With the above model modification, the film separation begins when
    \f$F_{ratio}>1\f$; however, only the portion with \f$m_{ratio}\f$
    separates while the rest stays attached.


    Reference:
    \verbatim
    Governing equations for the full film-separation model (tag:FLW):
        Friedrich, M. A., Lan, H., Wegener, J. L.,
        Drallmeier, J. A., & Armaly, B. F. (2008).
        A separation criterion with experimental
        validation for shear-driven films in separated flows.
        J. Fluids Eng. May 2008, 130(5): 051301.
        DOI:10.1115/1.2907405

    Governing equations for the partial film-separation model (tag:ZJD):
        Zhang, Y., Jia, M., Duan, H., Wang, P.,
        Wang, J., Liu, H., & Xie, M. (2018).
        Numerical and experimental study of spray impingement and liquid
        film separation during the spray/wall interaction at expanding corners.
        International Journal of Multiphase Flow, 107, 67-81.
        DOI:10.1016/j.ijmultiphaseflow.2018.05.016
    \endverbatim

Usage
    Minimal example in boundary-condition files:
    \verbatim
    filmSeparationCoeffs
    {
        // Mandatory entries
        model               Friedrich;
        rhop                <scalar>;

        // Optional entries
        separationType      <word>;
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                        | Type | Reqd | Deflt
      model        | Model name: Friedrich              | word | yes  | -
      rhop         | Primary-phase density            | scalar | yes  | -
      separationType | Film separation type             | word | no   | full
    \endtable

    Options for the \c separationType entry:
    \verbatim
      full          | Full film separation (Friedrich et al., 2008)
      partial       | Partial film separation (Zhang et al., 2018)
    \endverbatim

SourceFiles
    FriedrichModel.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_filmSeparationModels_FriedrichModel_H
#define Foam_filmSeparationModels_FriedrichModel_H

#include "filmSeparationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace filmSeparationModels
{

/*---------------------------------------------------------------------------*\
                      Class FriedrichModel Declaration
\*---------------------------------------------------------------------------*/

class FriedrichModel
:
    public filmSeparationModel
{
    // Private Enumerations

        //- Options for the film separation type
        enum separationType : char
        {
            FULL = 0,  //!< "Full film separation"
            PARTIAL    //!< "Partial film separation"
        };

        //- Names for separationType
        static const Enum<separationType> separationTypeNames;


    // Private Data

        //- Film separation type
        enum separationType separation_;

        //- Approximate uniform mass density of primary phase
        scalar rhop_;

        //- Magnitude of the gravitational acceleration
        scalar magG_;

        //- Empirical constant for the partial separation model
        scalar C0_;

        //- Empirical constant for the partial separation model
        scalar C1_;

        //- Empirical constant for the partial separation model
        scalar C2_;

        //- List of flags identifying sharp-corner edges where separation
        //- may occur
        bitSet cornerEdges_;

        //- Corner angles of sharp-corner edges where separation may occur
        scalarList cornerAngles_;


    // Private Member Functions

        //- Return Boolean list of identified sharp-corner edges
        bitSet calcCornerEdges() const;

        //- Return true if the given edge is identified as sharp
        bool isCornerEdgeSharp
        (
            const vector& faceCentreO,
            const vector& faceCentreN,
            const vector& faceNormalO,
            const vector& faceNormalN
        ) const;

        //- Return the list of sharp-corner angles for each edge
        scalarList calcCornerAngles() const;

        //- Return the sharp-corner angle for a given edge
        scalar calcCornerAngle
        (
            const vector& faceNormalO,
            const vector& faceNormalN
        ) const;

        //- Return Boolean list of identified separation faces
        bitSet calcSeparationFaces() const;

        //- Return true if the given face is identified as a separation face
        void isSeparationFace
        (
            bitSet& separationFaces,
            const scalar phiEdge,
            const label faceO,
            const label faceN = -1
        ) const;

        //- Return the list of sharp-corner angles for each face
        scalarList calcSeparationAngles(const bitSet& separationFaces) const;

        //- Return the film-separation force ratio per face
        tmp<scalarField> Fratio() const;


public:

    //- Runtime type information
    TypeName("Friedrich");


    // Constructors

        //- Construct from the base film model and dictionary
        FriedrichModel
        (
            const regionModels::areaSurfaceFilmModels::liquidFilmBase& film,
            const dictionary& dict
        );


    // Destructor
    virtual ~FriedrichModel() = default;


    // Member Functions

    // Evaluation

        //- Calculate the mass ratio of film separation
        virtual tmp<scalarField> separatedMassRatio() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace filmSeparationModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
