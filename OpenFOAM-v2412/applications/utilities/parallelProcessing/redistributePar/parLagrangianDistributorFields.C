/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022-2023 OpenCFD Ltd.
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

#include "parLagrangianDistributor.H"
#include "unmappedPassivePositionParticleCloud.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::parLagrangianDistributor::readAllFields
(
    const passivePositionParticleCloud& cloud,
    const bool haveCloud,
    const IOobjectList& cloudObjs,
    const wordRes& selectedFields
)
{
    label nTotal = 0;

    do
    {
        #undef  doLocalCode
        #define doLocalCode(Type)                                             \
        {                                                                     \
            nTotal += parLagrangianDistributor::readFields                    \
            <IOField<Type>>                                                   \
            (                                                                 \
                cloud,                                                        \
                haveCloud,                                                    \
                cloudObjs,                                                    \
                selectedFields                                                \
            );                                                                \
                                                                              \
            nTotal += parLagrangianDistributor::readFields                    \
            <IOField<Field<Type>>>                                            \
            (                                                                 \
                cloud,                                                        \
                haveCloud,                                                    \
                cloudObjs,                                                    \
                selectedFields                                                \
            );                                                                \
                                                                              \
            nTotal += parLagrangianDistributor::readFields                    \
            <CompactIOField<Field<Type>, Type>>                               \
            (                                                                 \
                cloud,                                                        \
                haveCloud,                                                    \
                cloudObjs,                                                    \
                selectedFields                                                \
            );                                                                \
        }

        doLocalCode(label);
        doLocalCode(scalar);
        doLocalCode(vector);
        doLocalCode(sphericalTensor);
        doLocalCode(symmTensor);
        doLocalCode(tensor);

        #undef doLocalCode
    }
    while (false);

    return nTotal;
}


//Foam::label Foam::parLagrangianDistributor::readAllFields
//(
//    const passivePositionParticleCloud& cloud,
//    const wordRes& selectedFields
//)
//{
//    IOobjectList cloudObjs(cloud, cloud.time().timeName());
//    return readAllFields(cloud, cloudObjs, selectedFields);
//}


Foam::label Foam::parLagrangianDistributor::distributeAllFields
(
    const mapDistributeBase& lagrangianMap,
    const word& cloudName,
    const bool haveCloud,
    const IOobjectList& cloudObjs,
    const wordRes& selectedFields
) const
{
    label nTotal = 0;

    do
    {
        #undef  doLocalCode
        #define doLocalCode(Type)                                             \
        {                                                                     \
            nTotal += this->distributeFields<Type>                            \
            (                                                                 \
                lagrangianMap,                                                \
                cloudName,                                                    \
                haveCloud,                                                    \
                cloudObjs,                                                    \
                selectedFields                                                \
            );                                                                \
                                                                              \
            nTotal += this->distributeFieldFields<Type>                       \
            (                                                                 \
                lagrangianMap,                                                \
                cloudName,                                                    \
                haveCloud,                                                    \
                cloudObjs,                                                    \
                selectedFields                                                \
            );                                                                \
        }

        doLocalCode(label);
        doLocalCode(scalar);
        doLocalCode(vector);
        doLocalCode(sphericalTensor);
        doLocalCode(symmTensor);
        doLocalCode(tensor);

        #undef doLocalCode
    }
    while (false);

    return nTotal;
}


Foam::label Foam::parLagrangianDistributor::distributeAllStoredFields
(
    const mapDistributeBase& lagrangianMap,
    passivePositionParticleCloud& cloud
) const
{
    label nTotal = 0;

    do
    {
        #undef  doLocalCode
        #define doLocalCode(Type)                                             \
        {                                                                     \
            nTotal += this->distributeStoredFields                            \
            <IOField<Type>>                                                   \
            (                                                                 \
                lagrangianMap,                                                \
                cloud                                                         \
            );                                                                \
                                                                              \
            nTotal += this->distributeStoredFields                            \
            <IOField<Field<Type>>>                                            \
            (                                                                 \
                lagrangianMap,                                                \
                cloud                                                         \
            );                                                                \
                                                                              \
            nTotal += this->distributeStoredFields                            \
            <CompactIOField<Field<Type>, Type>>                               \
            (                                                                 \
                lagrangianMap,                                                \
                cloud                                                         \
            );                                                                \
        }

        doLocalCode(label);
        doLocalCode(scalar);
        doLocalCode(vector);
        doLocalCode(sphericalTensor);
        doLocalCode(symmTensor);
        doLocalCode(tensor);

        #undef doLocalCode
    }
    while (false);

    return nTotal;
}


// ************************************************************************* //
