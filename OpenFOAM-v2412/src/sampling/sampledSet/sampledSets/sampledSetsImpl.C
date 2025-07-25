/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2015-2023 OpenCFD Ltd.
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

#include "sampledSets.H"
#include "globalIndex.H"
#include "interpolation.H"
#include "volFields.H"
#include "volPointInterpolation.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class GeoField>
Foam::tmp<GeoField>
Foam::sampledSets::getOrLoadField(const word& fieldName) const
{
    tmp<GeoField> tfield;

    if (loadFromFiles_)
    {
        tfield.emplace
        (
            IOobject
            (
                fieldName,
                mesh_.time().timeName(),
                mesh_.thisDb(),
                IOobjectOption::MUST_READ,
                IOobjectOption::NO_WRITE,
                IOobjectOption::NO_REGISTER
            ),
            mesh_
        );
    }
    else
    {
        tfield.cref(mesh_.cfindObject<GeoField>(fieldName));
    }

    return tfield;
}


template<class Type>
void Foam::sampledSets::writeCoordSet
(
    coordSetWriter& writer,
    const Field<Type>& values,
    const word& fieldName
)
{
    fileName outputName;
    if (Pstream::master())
    {
        outputName = writer.write(fieldName, values);
    }
    Pstream::broadcast(outputName);

    if (outputName.size())
    {
        // Case-local file name with "<case>" to make relocatable

        dictionary propsDict;
        propsDict.add
        (
            "file",
            time_.relativePath(outputName, true)
        );
        setProperty(fieldName, propsDict);
    }
}


template<class Type>
void Foam::sampledSets::performAction
(
    const VolumeField<Type>& fld,
    unsigned request
)
{
    const word& fieldName = fld.name();
    const scalar timeValue = fld.time().timeOutputValue();

    // The interpolator for this field
    autoPtr<interpolation<Type>> interpPtr;

    if (!samplePointScheme_.empty() && samplePointScheme_ != "cell")
    {
        interpPtr.reset(interpolation<Type>::New(samplePointScheme_, fld));
    }

    const unsigned int width(IOstream::defaultPrecision() + 7);
    OFstream* osptr = nullptr;

    if (writeAsProbes_ && (request & ACTION_WRITE))
    {
        osptr = createProbeFile(fieldName);

        if (Pstream::master() && osptr)
        {
            (*osptr) << setw(width) << timeValue;
        }
    }

    // Ensemble min/max/avg values
    Type avgEnsemble = Zero;
    label sizeEnsemble = 0;
    MinMax<Type> limitsEnsemble;

    forAll(*this, seti)
    {
        const sampledSet& s = (*this)[seti];
        const globalIndex& globIdx = globalIndices_[seti];
        const labelList& globOrder = gatheredSorting_[seti];

        const word& setName = s.name();
        Field<Type> values(s.size());

        if (interpPtr)
        {
            forAll(s, samplei)
            {
                const point& p = s[samplei];
                const label celli = s.cells()[samplei];
                const label facei = s.faces()[samplei];

                if (celli == -1 && facei == -1)
                {
                    // Special condition for illegal sampling points
                    values[samplei] = pTraits<Type>::max;
                }
                else
                {
                    values[samplei] = interpPtr().interpolate(p, celli, facei);
                }
            }
        }
        else
        {
            forAll(s, samplei)
            {
                const label celli = s.cells()[samplei];

                if (celli == -1)
                {
                    values[samplei] = pTraits<Type>::max;
                }
                else
                {
                    values[samplei] = fld[celli];
                }
            }
        }

        // Collect data from all processors
        globIdx.gatherInplace(values);

        // Local min/max/avg values - calculate on master
        Type avgValue = Zero;
        label sizeValue = 0;
        MinMax<Type> limits;

        if (Pstream::master())
        {
            avgValue = sum(values);
            sizeValue = values.size();
            limits = MinMax<Type>(values);

            // Ensemble values
            avgEnsemble += avgValue;
            sizeEnsemble += sizeValue;
            limitsEnsemble += limits;

            if (sizeValue)
            {
                avgValue /= sizeValue;
            }

            // Use sorted order
            values = UIndirectList<Type>(values, globOrder)();
        }
        Pstream::broadcasts(UPstream::worldComm, avgValue, sizeValue, limits);

        // Store results: min/max/average/size with the name of the set
        // for scoping.
        // Eg, average(lines,T) ...
        const word resultArg('(' + setName + ',' + fieldName + ')');

        this->setResult("average" + resultArg, avgValue);
        this->setResult("min" + resultArg, limits.min());
        this->setResult("max" + resultArg, limits.max());
        this->setResult("size" + resultArg, sizeValue);

        if (verbose_)
        {
            Info<< name() << ' ' << setName << " : " << fieldName << nl
                << "    avg: " << avgValue << nl
                << "    min: " << limits.min() << nl
                << "    max: " << limits.max() << nl << nl;
        }

        if ((request & ACTION_WRITE) != 0)
        {
            if (writeAsProbes_)
            {
                if (osptr)
                {
                    for (const Type& val : values)
                    {
                        (*osptr) << ' ' << setw(width) << val;
                    }
                }
            }
            else
            {
                writeCoordSet<Type>(writers_[seti], values, fieldName);
            }
        }
    }

    // Finish probes write
    if (Pstream::master() && osptr)
    {
        (*osptr) << endl;
    }

    if (sizeEnsemble)
    {
        avgEnsemble /= sizeEnsemble;
    }

    if (size())
    {
        Pstream::broadcasts
        (
            UPstream::worldComm,
            avgEnsemble,
            sizeEnsemble,
            limitsEnsemble
        );

        // Store results: min/max/average/size for the ensemble
        // Eg, average(T) ...
        const word resultArg('(' + fieldName + ')');

        this->setResult("average" + resultArg, avgEnsemble);
        this->setResult("min" + resultArg, limitsEnsemble.min());
        this->setResult("max" + resultArg, limitsEnsemble.max());
        this->setResult("size" + resultArg, sizeEnsemble);
    }
}


template<class GeoField>
void Foam::sampledSets::performAction
(
    const IOobjectList& objects,
    unsigned request
)
{
    wordList fieldNames;
    if (loadFromFiles_)
    {
        fieldNames = objects.sortedNames<GeoField>(fieldSelection_);
    }
    else
    {
        fieldNames = mesh_.thisDb().sortedNames<GeoField>(fieldSelection_);
    }

    for (const word& fieldName : fieldNames)
    {
        tmp<GeoField> tfield = getOrLoadField<GeoField>(fieldName);

        if (tfield)
        {
            performAction<typename GeoField::value_type>(tfield(), request);
        }
    }
}


// ************************************************************************* //
