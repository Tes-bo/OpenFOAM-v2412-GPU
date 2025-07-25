/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 OpenCFD Ltd.
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
    Foam::fa::limitHeight

Description
    Limits the film height to a specified \c max value.

    The \c limitHeight applies an explicit correction to limit the max film
    height to the user-specified value. This will ensure that the max value
    is adhered to, but may become unstable if the limiting instigates an
    oscillatory/destabilising effect into neighbouring faces that then spreads.

Usage
    Minimal example by using \c constant/faOptions:
    \verbatim
    <faOption>
    {
        // Mandatory entries
        type            limitHeight;
        active          yes;
        selectionMode   all;
        max             <scalar>;

        // Optional entries
        h               <word>;

        // Inherited entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property  | Description                   | Type   | Reqd  | Deflt
      type      | Type name: limitHeight        | word   | yes   | -
      max       | Maximum height limit [m]      | scalar | yes   | -
      h         | Name of operand height field  | word   | no    | h
    \endtable

    The inherited entries are elaborated in:
      - \link faceSetOption.H \endlink

SourceFiles
    limitHeight.C

\*---------------------------------------------------------------------------*/

#ifndef fa_limitHeight_H
#define fa_limitHeight_H

#include "faceSetOption.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fa
{

/*---------------------------------------------------------------------------*\
                        Class limitHeight Declaration
\*---------------------------------------------------------------------------*/

class limitHeight
:
    public faceSetOption
{
protected:

    // Protected Data

        //- Name of operand height field
        word hName_;

        //- Maximum height [m]
        scalar max_;


    // Protected Member Functions

        //- No copy construct
        limitHeight(const limitHeight&) = delete;

        //- No copy assignment
        void operator=(const limitHeight&) = delete;


public:

    //- Runtime type information
    TypeName("limitHeight");


    // Constructors

        //- Construct from components
        limitHeight
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~limitHeight() = default;


    // Member Functions

        //- Read dictionary
        virtual bool read(const dictionary& dict);

        //- Correct the height field
        virtual void correct(areaScalarField& h);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
