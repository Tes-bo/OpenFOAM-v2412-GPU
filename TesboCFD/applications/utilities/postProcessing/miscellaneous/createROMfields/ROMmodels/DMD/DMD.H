/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023 OpenCFD Ltd.
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
    Foam::ROMmodels::DMD

Description
    Field creation model using the streaming total dynamic mode
    decomposition method (STDMD).

    The governing equations of the field creation are as follows:

    \f[
        \mathbf{x}_\tau \approx
        \tilde{\mathbf{x}_\tau} =
        \left( \sum_{i=0}^{N-1} \phi_i \alpha_i v_i^\tau \right)_{real}​
    \f]

    with

    \f[
        \tau = \frac{t - t_o}{\Delta t}​
    \f]

    where:

    \vartable
        \mathbf{x}    | Field snapshot at time t
        \tilde{\mathbf{x}} | Reconstructed field snapshot at time t (complex)
        N             | Number of modes
        i             | Mode index
        \tau          | Nondimensional time
        t             | Time [s]
        t_o           | Start time (of mode decomposition calculations) [s]
        \Delta t      | Time-step size of mode decomposition [s]
        \phi          | Mode (complex)
        \alpha        | Mode amplitude (complex)
        v             | Mode eigenvalue (complex)
    \endvartable

    References:
    \verbatim
        Governing equations (tag:K):
            Kiewat, M. (2019).
            Streaming modal decomposition approaches for vehicle aerodynamics.
            PhD thesis. Munich: Technical University of Munich.
            URL:mediatum.ub.tum.de/doc/1482652/1482652.pdf
    \endverbatim

    Operands:
    \table
      Operand      | Type                 | Location
      input        | {vol,surface}\<Type\>Field | \<time\>/\<inpField\>
      output file  | -                    | -
      output field | {vol,surface}\<Type\>Field | \<time\>/\<outField\>
    \endtable

    where \c \<Type\>=Scalar/Vector/SphericalTensor/SymmTensor/Tensor.

Usage
    Minimal example by using \c system/ROMfieldsDict:
    \verbatim
        // Mandatory entries
        field           <word>;
        object          <word>;
        deltaT          <scalar>;
        time            <scalar>;
        modes           <labelList>;
        amplitudes      <complexList>;
        eigenvalues     <complexList>;

        // Optional entries
        startTime       <scalar>;
        dimensions      <dimensionSet>;
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                        | Type | Reqd  | Deflt
      field        | Name of reconstructed field        | word | yes   | -
      object       | Name of operand function object    | word | yes   | -
      deltaT       | Time-step size of mode decomposition | scalar | yes | -
      time         | Time instant where modes are located | scalar | yes | -
      modes        | List of mode indices               | labelList | yes | -
      amplitudes   | Amplitude coefficients             | complexList | yes | -
      eigenvalues  | Eigenvalues                        | complexList | yes | -
      startTime | Start time for mode-information collection | scalar | no | 0
      dimensions | Dimensions of reconstructed fields | dimensionSet | no | -
    \endtable

SourceFiles
    DMD.C
    DMDImpl.C

\*---------------------------------------------------------------------------*/

#ifndef ROMmodels_DMD_H
#define ROMmodels_DMD_H

#include "ROMmodels/ROMmodel/ROMmodel.H"
#include "dimensionSet.H"
#include "complex.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace ROMmodels
{
/*---------------------------------------------------------------------------*\
                            Class DMD Declaration
\*---------------------------------------------------------------------------*/

class DMD
:
    public ROMmodel
{
    // Private Typedefs

        typedef List<complex> complexList;


    // Private Data

        //- Name of reconstructed field
        word fieldName_;

        //- Name of operand function object
        word objectName_;

        //- Time-step size of mode decomposition (not that of the simulation)
        scalar deltaT_;

        //- Time instant where modes are located
        scalar time_;

        //- Start time for mode-information collection
        scalar startTime_;

        //- List of mode indices
        labelList modes_;

        //- Dimensions of reconstructed fields
        dimensionSet dims_;

        //- Amplitude coefficients
        complexList amps_;

        //- Eigenvalues
        complexList evals_;


    // Private Member Functions

        //- Return names of mode fields
        wordList modeNames(const word& modeType) const;

        //- Implementation for creating and writing fields
        template<class GeoType>
        bool createAndWriteImpl() const;


public:

    //- Runtime type information
    TypeName("DMD");


    // Constructors

        //- Construct from components
        DMD
        (
            Time& runTime,
            fvMesh& mesh,
            const dictionary& dict,
            const instantList& times
        );

        //- No copy construct
        DMD(const DMD&) = delete;

        //- No copy assignment
        void operator=(const DMD&) = delete;


    //- Destructor
    virtual ~DMD() = default;


    // Member Functions

        //- Read model settings
        virtual bool read(const dictionary& dict);

        //- Create and write fields
        virtual bool createAndWrite();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace ROMmodels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
