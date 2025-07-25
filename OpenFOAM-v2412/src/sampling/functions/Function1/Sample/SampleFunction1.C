/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021-2022 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "SampleFunction1.H"
#include "volFields.H"
#include "interpolation.H"
#include "pointIOField.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::Function1Types::Sample<Type>::setSampleCell() const
{
    const polyMesh& mesh = this->template mesh<polyMesh>();

    const auto& points = static_cast<const pointIOField&>(mesh.points());

    if (pointEventNo_ < points.eventNo())
    {
        pointEventNo_ = points.eventNo();

        celli_ = mesh.findCell(position_);

        if (returnReduceAnd(celli_ < 0))
        {
            FatalErrorInFunction
                << "Sample cell could not be found at position "
                << position_ << nl
                << exit(FatalError);
        }

        if (debug)
        {
            Pout<< "Position: " << position_
                << " celli:" << celli_
                << " eventNo:" << pointEventNo_
                << " points eventNo:" << points.eventNo()
                << endl;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::Function1Types::Sample<Type>::Sample
(
    const word& entryName,
    const dictionary& dict,
    const objectRegistry* obrPtr
)
:
    Function1<Type>(entryName, dict, obrPtr),
    fieldName_(dict.get<word>("field")),
    position_(dict.get<point>("position")),
    interpolationScheme_
    (
        dict.getOrDefault<word>("interpolationScheme", "cell")
    ),
    celli_(-1),
    pointEventNo_(-1)
{}


template<class Type>
Foam::Function1Types::Sample<Type>::Sample(const Sample& s)
:
    Function1<Type>(s),
    fieldName_(s.fieldName_),
    position_(s.position_),
    interpolationScheme_(s.interpolationScheme_),
    celli_(s.celli_),
    pointEventNo_(s.pointEventNo_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Type Foam::Function1Types::Sample<Type>::value(const scalar x) const
{
    const auto& mesh = this->template mesh<fvMesh>();

    const auto* fieldPtr =
        mesh.template cfindObject<VolumeField<Type>>(fieldName_);

    if (!fieldPtr)
    {
        FatalErrorInFunction
            << "Unable to find field " << fieldName_ << " on the mesh database"
            << ". Valid " << VolumeField<Type>::typeName << " fields are:"
            << mesh.template sortedNames<VolumeField<Type>>()
            << exit(FatalError);
    }


    // Might trigger parallel comms (e.g. volPointInterpolation, if
    // result is not yet cached) so have all processors do it
    autoPtr<interpolation<Type>> interpolator
    (
        interpolation<Type>::New(interpolationScheme_, *fieldPtr)
    );

    Type result = pTraits<Type>::min;

    setSampleCell();

    if (celli_ != -1)
    {
        result = interpolator().interpolate(position_, celli_, -1);
    }

    reduce(result, maxOp<Type>());

    DebugInfo << "sampled value: " << result << endl;

    return result;
}


template<class Type>
Type Foam::Function1Types::Sample<Type>::integrate
(
    const scalar x1,
    const scalar x2
) const
{
    NotImplemented;

    return Zero;
}


template<class Type>
void Foam::Function1Types::Sample<Type>::writeEntries(Ostream& os) const
{
    os.writeEntry("field", fieldName_);
    os.writeEntry("position", position_);

    os.writeEntryIfDifferent<word>
    (
        "interpolationScheme", "cell", interpolationScheme_
    );
}


template<class Type>
void Foam::Function1Types::Sample<Type>::writeData(Ostream& os) const
{
    Function1<Type>::writeData(os);
    os.endEntry();

    os.beginBlock(word(this->name() + "Coeffs"));
    writeEntries(os);
    os.endBlock();
}


// ************************************************************************* //
