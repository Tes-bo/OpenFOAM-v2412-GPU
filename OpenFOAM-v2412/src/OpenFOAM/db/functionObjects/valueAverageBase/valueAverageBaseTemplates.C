/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015 OpenFOAM Foundation
    Copyright (C) 2015-2022 OpenCFD Ltd.
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

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class Type, class Type2>
bool Foam::functionObjects::valueAverageBase::calc
(
    const label fieldi,
    bool& converged,
    dictionary& dict
)
{
    const word& fieldName = fieldNames_[fieldi];

    const word valueType =
        state_.objectResultType(functionObjectName_, fieldName);

    if (pTraits<Type>::typeName != valueType)
    {
        return false;
    }

    const scalar dt = state_.time().deltaTValue();

    const Type2 currentValue =
        state_.getObjectResult<Type>(functionObjectName_, fieldName);

    // Current mean value
    const word meanName(fieldName + "Mean");
    Type2 meanValue = state_.getResult<Type2>(meanName);

    switch (windowType_)
    {
        case windowType::NONE:
        {
            const scalar Dt = totalTime_[fieldi];
            const scalar beta = dt/Dt;
            meanValue = (1 - beta)*meanValue + beta*currentValue;

            break;
        }
        case windowType::APPROXIMATE:
        {
            const scalar Dt = totalTime_[fieldi];
            scalar beta = dt/Dt;
            if (Dt - dt >= window_)
            {
                beta = dt/window_;
            }
            else
            {
                converged = false;
            }

            meanValue = (1 - beta)*meanValue + beta*currentValue;

            break;
        }
        case windowType::EXACT:
        {
            FIFOStack<scalar> windowTimes;
            FIFOStack<Type2> windowValues;
            dictionary& fieldDict = dict.subDict(fieldName);
            fieldDict.readIfPresent("windowTimes", windowTimes);
            fieldDict.readIfPresent("windowValues", windowValues);

            // Increment time for all existing values
            for (scalar& dti : windowTimes)
            {
                dti += dt;
            }

            // Remove any values outside the window
            bool removeValue = true;
            while (removeValue && windowTimes.size())
            {
                removeValue = windowTimes.first() > window_;

                if (removeValue)
                {
                    windowTimes.pop();
                    windowValues.pop();
                }
            }

            // Add the current value
            windowTimes.push(dt);
            windowValues.push(currentValue);

            // Calculate the window average
            auto timeIter = windowTimes.cbegin();
            auto valueIter = windowValues.cbegin();

            meanValue = Zero;
            Type2 valueOld(Zero);

            for
            (
                label i = 0;
                timeIter.good();
                ++i, ++timeIter, ++valueIter
            )
            {
                const Type2& value = valueIter();
                const scalar dt = timeIter();

                meanValue += dt*value;

                if (i)
                {
                    meanValue -= dt*valueOld;
                }

                valueOld = value;
            }

            meanValue /= windowTimes.first();

            // Store the state information for the next step
            fieldDict.set("windowTimes", windowTimes);
            fieldDict.set("windowValues", windowValues);

            break;
        }
    }

    scalar delta = mag(meanValue - currentValue);

    Log << indent << "    " << meanName << ": " << meanValue
        << ", delta: " << delta << nl;

    file() << tab << meanValue;

    state_.setResult(meanName, meanValue);

    if ((tolerance_ > 0) && (delta > tolerance_))
    {
        converged = false;
    }

    return true;
}


// ************************************************************************* //
