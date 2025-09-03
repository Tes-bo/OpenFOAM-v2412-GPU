/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2019 OpenCFD Ltd.
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

#include "objectRegistry.H"
#include "Time.H"

template<class Type>
bool Foam::functionObjects::fieldAverageItem::calculateMeanField
(
    const objectRegistry& obr
) const
{
    if (!mean_)
    {
        return false;
    }

    const Type* baseFieldPtr = obr.findObject<Type>(fieldName_);

    if (!baseFieldPtr)
    {
        return false;
    }

    const Type& baseField = *baseFieldPtr;
    Type& meanField = obr.lookupObjectRef<Type>(meanFieldName_);

    switch (windowType_)
    {
        case windowType::NONE:
        {
            scalar dt = this->dt(obr.time().deltaTValue());
            scalar Dt = this->Dt();
            scalar beta = dt/Dt;

            meanField = (1 - beta)*meanField + beta*baseField;

            break;
        }
        case windowType::APPROXIMATE:
        {
            scalar dt = this->dt(obr.time().deltaTValue());
            scalar Dt = this->Dt();
            scalar beta = dt/Dt;

            if (Dt - dt >= window_)
            {
                beta = dt/window_;
            }

            meanField = (1 - beta)*meanField + beta*baseField;

            break;
        }
        case windowType::EXACT:
        {
            switch (base_)
            {
                case baseType::ITER:
                {
                    // Uniform time step - can use simplified algorithm
                    // Note: stores an additional old time field, but only
                    //       needs to do 1 field lookup

                    label n = windowTimes_.size();
                    const Type& lastField =
                        obr.lookupObject<Type>(windowFieldNames_.first());

                    if (n <= round(window_))
                    {
                        scalar beta = 1.0/scalar(n);
                        meanField = (1 - beta)*meanField + beta*baseField;
                    }
                    else
                    {
                        meanField += (baseField - lastField)/scalar(n - 1);
                    }

                    break;
                }
                case baseType::TIME:
                {
                    // Assuming non-uniform time step
                    // Note: looks up all window fields from the registry

                    meanField = 0*baseField;

                    auto timeIter = windowTimes_.cbegin();
                    auto nameIter = windowFieldNames_.cbegin();

                    const Type* wOld = nullptr;

                    for
                    (
                        ;
                        timeIter.good();
                        ++timeIter, ++nameIter
                    )
                    {
                        const word& fieldName = nameIter();
                        const scalar dt = timeIter();
                        const Type* w = obr.findObject<Type>(fieldName);

                        meanField += dt*(*w);

                        if (wOld)
                        {
                            meanField -= dt*(*wOld);
                        }

                        wOld = w;
                    }

                    meanField /= windowTimes_.first();

                    break;
                }
                default:
                {
                    FatalErrorInFunction
                        << "Unhandled baseType enumeration "
                        << baseTypeNames_[base_]
                        << abort(FatalError);
                }
            }

            break;
        }
        default:
        {
            FatalErrorInFunction
                << "Unhandled windowType enumeration "
                << windowTypeNames_[windowType_]
                << abort(FatalError);
        }
    }

    return true;
}


template<class Type1, class Type2>
bool Foam::functionObjects::fieldAverageItem::calculatePrime2MeanField
(
    const objectRegistry& obr
) const
{
    if (!prime2Mean_)
    {
        return false;
    }

    const Type1* baseFieldPtr = obr.findObject<Type1>(fieldName_);

    if (!baseFieldPtr)
    {
        return false;
    }

    const Type1& baseField = *baseFieldPtr;
    const Type1& meanField = obr.lookupObject<Type1>(meanFieldName_);

    Type2& prime2MeanField =
        obr.lookupObjectRef<Type2>(prime2MeanFieldName_);

    switch (windowType_)
    {
        case windowType::NONE:
        {
            scalar dt = this->dt(obr.time().deltaTValue());
            scalar Dt = this->Dt();
            scalar beta = dt/Dt;

            prime2MeanField =
                (1 - beta)*prime2MeanField
              + beta*sqr(baseField)
              - sqr(meanField);

            break;
        }
        case windowType::APPROXIMATE:
        {
            scalar dt = this->dt(obr.time().deltaTValue());
            scalar Dt = this->Dt();
            scalar beta = dt/Dt;

            if (Dt - dt >= window_)
            {
                beta = dt/window_;
            }

            prime2MeanField =
                (1 - beta)*prime2MeanField
              + beta*sqr(baseField)
              - sqr(meanField);

            break;
        }
        case windowType::EXACT:
        {
            // Not storing old time mean fields - treat all as TIME (integrated)
            prime2MeanField = 0*prime2MeanField;

            auto timeIter = windowTimes_.cbegin();
            auto nameIter = windowFieldNames_.cbegin();

            switch (base_)
            {
                case baseType::ITER:
                {
                    // ITER method stores an additional entry compared to TIME
                    ++timeIter;
                    ++nameIter;

                    if (!timeIter.good()) return false;

                    break;
                }
                default:
                {}
            }


            scalar windowLength = timeIter();

            const Type1* wOld = nullptr;

            for
            (
                ;
                timeIter.good();
                ++timeIter, ++nameIter
            )
            {
                const word& fieldName = nameIter();
                const scalar dt = timeIter();
                const Type1* w = obr.findObject<Type1>(fieldName);

                prime2MeanField += dt*(sqr((*w) - meanField));

                if (wOld)
                {
                    prime2MeanField -= dt*(sqr((*wOld) - meanField));
                }

                wOld = w;
            }

            prime2MeanField /= windowLength;

            break;
        }
        default:
        {
            FatalErrorInFunction
                << "Unhandled windowType enumeration "
                << windowTypeNames_[windowType_]
                << abort(FatalError);
        }
    }

    return true;
}


// ************************************************************************* //
