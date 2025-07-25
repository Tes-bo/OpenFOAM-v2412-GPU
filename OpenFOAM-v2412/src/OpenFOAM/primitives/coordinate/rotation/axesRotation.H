/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2017-2024 OpenCFD Ltd.
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
    Foam::coordinateRotations::axes

Description
    A coordinateRotation specified using global axes.

    The rotation is defined by a combination of vectors (e1/e2), (e2/e3)
    or (e3/e1). Any nonorthogonality is absorbed into the second vector.

    \verbatim
    rotation
    {
        type        axes;
        e1          (1 0 0);
        e2          (0 1 0);
    }
    \endverbatim

    \heading Dictionary entries
    \table
        Property    | Description                           | Reqd  | Default
        type        | Type name: axes                       | yes   |
        e1          | local x-axis                          | partly |
        e2          | local y-axis                          | partly |
        e3          | local z-axis                          | partly |
    \endtable

Note
    It is also possible to specify in terms of \c axis and \c direction.
    For cylindrical coordinates, the \c axis would correspond to the
    \a z-axis and the \c direction to the \a r-axis.
    Also accepts "axesRotation" (OpenFOAM-v1806) for the type.

SourceFiles
    axesRotation.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_coordinateRotations_axes_H
#define Foam_coordinateRotations_axes_H

#include "coordinateRotation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace coordinateRotations
{

/*---------------------------------------------------------------------------*\
                  Class coordinateRotations::axes Declaration
\*---------------------------------------------------------------------------*/

class axes
:
    public coordinateRotation
{
public:

    //- The order/combination of local axes for the axes-rotation definition
    //  Note that these follow the right-hand rule.
    enum axisOrder
    {
        E1_E2,  //!< The axis1 (dominant) is local X, axis2 is local Y.
        E2_E3,  //!< The axis1 (dominant) is local Y, axis2 is local Z.
        E3_E1,  //!< The axis1 (dominant) is local Z, axis2 is local X.
        E3_E1_COMPAT  //!< E3_E1 specified as axis/direction.
    };


protected:

    // Protected Data

        //- The primary axis
        vector axis1_;

        //- The secondary axis
        vector axis2_;

        //- The axis order
        axisOrder order_;


    // Protected Member Functions

        //- Read from dictionary
        void read(const dictionary& dict);


public:

    //- Runtime type information
    TypeNameNoDebug("axes");


    // Constructors

        //- Default construct - an identity transform
        axes();

        //- Copy construct
        axes(const axes& crot);

        //- Construct from two axes
        axes(const vector& axis1, const vector& axis2, axisOrder order=E3_E1);

        //- Construct from a single axis (as e3) using a best-guess for the
        //- second axis.
        //  The largest component and its sign are used when guessing
        //  an appropriate orientation (direction).
        explicit axes(const vector& axis);

        //- Construct from dictionary
        explicit axes(const dictionary& dict);

        //- Return clone
        autoPtr<coordinateRotation> clone() const
        {
            return coordinateRotation::Clone(*this);
        }


    //- Destructor
    virtual ~axes() = default;


    // Static Member Functions

        //- The rotation tensor calculated from two axes and their order.
        //  The input axes will be normalised.
        //  If axis2 is zero, an axis orthogonal to axis1 will be guessed.
        //  The largest component and its sign are used when guessing
        //  an appropriate orientation (direction).
        static tensor rotation
        (
            const vector& axis1,
            const vector& axis2,
            axisOrder order = E3_E1
        );


    // Member Functions

        //- Reset specification
        virtual void clear();

        //- The rotation tensor calculated from the specified axes and order
        virtual tensor R() const;

        //- Write information
        virtual void write(Ostream& os) const;

        //- Write dictionary entry
        virtual void writeEntry(const word& keyword, Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace coordinateRotations

// Compatibility typedef 1806
typedef coordinateRotations::axes axesRotation;

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
