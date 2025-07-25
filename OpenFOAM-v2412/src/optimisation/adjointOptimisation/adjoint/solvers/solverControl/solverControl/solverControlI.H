/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2019 PCOpt/NTUA
    Copyright (C) 2013-2019 FOSS GP
    Copyright (C) 2019 OpenCFD Ltd.
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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::dictionary Foam::solverControl::solverDict() const
{
    return solver_.dict();
}


inline const Foam::dictionary Foam::solverControl::solutionDict() const
{
    return solverDict().subDict("solutionControls");
}


inline bool Foam::solverControl::printMaxMags() const
{
    return printMaxMags_;
}


inline bool Foam::solverControl::storeInitValues() const
{
    return storeInitValues_;
}


inline Foam::label Foam::solverControl::iter() const
{
    return iter_;
}


inline Foam::label& Foam::solverControl::averageIter()
{
    return averageIter_;
}


inline Foam::label Foam::solverControl::averageIter() const
{
    return averageIter_;
}


inline Foam::label Foam::solverControl::averageStartIter() const
{
    return averageStartIter_;
}


inline bool Foam::solverControl::doAverageIter() const
{
    if (average_ && iter_ >= averageStartIter_)
    {
        return true;
    }
    else
    {
        return false;
    }
}


inline bool Foam::solverControl::useAveragedFields() const
{
    if (average_ && averageIter_)
    {
        return true;
    }
    else
    {
        return false;
    }
}


inline bool Foam::solverControl::average() const
{
    return average_;
}


inline const Foam::solver& Foam::solverControl::getSolver() const
{
    return solver_;
}


// ************************************************************************* //
