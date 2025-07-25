/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2019 OpenCFD Ltd.
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

#include "fileName.H"

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

inline const Foam::HashSet<Foam::adiosFoam::cloudInfo>*
Foam::adiosFoam::adiosReader::getCloudInfo(const word& regName) const
{
    const auto iter = clouds_.cfind(regName);

    if (iter.found())
    {
        return &(iter.val());
    }

    return nullptr;
}


inline const Foam::HashSet<Foam::adiosFoam::fieldInfo>*
Foam::adiosFoam::adiosReader::getFieldInfo(const word& regName) const
{
    const auto iter = fields_.cfind(regName);

    if (iter.found())
    {
        return &(iter.val());
    }

    return nullptr;
}


inline std::string
Foam::adiosFoam::adiosReader::inquireAttributeType
(
    const std::string& name
) const
{
    const auto iter = availableAttr_.find(name);

    if (iter != availableAttr_.end())
    {
        return iter->second.find("Type")->second;
    }

    return "";
}


inline std::string
Foam::adiosFoam::adiosReader::inquireVariableType
(
    const std::string& name
) const
{
    const auto iter = availableVars_.find(name);

    if (iter != availableVars_.end())
    {
        return iter->second.find("Type")->second;
    }

    return "";
}


template<class T>
inline bool Foam::adiosFoam::adiosReader::hasAttribute
(
    const string& name
) const
{
    const auto iter = availableAttr_.find(name);

    return
    (
        iter != availableAttr_.end()
     &&
        (
            std::is_void<T>::value
         || (adios2::GetType<T>() == iter->second.find("Type")->second)
        )
    );
}


template<class T>
inline bool Foam::adiosFoam::adiosReader::hasVariable
(
    const string& name
) const
{
    const auto iter = availableVars_.find(name);

    return
    (
        iter != availableVars_.end()
     &&
        (
            std::is_void<T>::value
         || (adios2::GetType<T>() == iter->second.find("Type")->second)
        )
    );
}


inline bool Foam::adiosFoam::adiosReader::readIntAttribute
(
    const string& name,
    label& value,
    bool mandatory
)
{
    // Probably don't distinguish between signed/unsigned
    // mostly just use signed anyhow

    return readAttributeImpl<label, int>(name, value, mandatory);
}


inline bool Foam::adiosFoam::adiosReader::readScalarAttribute
(
    const string& name,
    double& value,
    bool mandatory
)
{
    return readAttributeImpl<double, double> (name, value, mandatory);
}


inline bool Foam::adiosFoam::adiosReader::readStringAttribute
(
    const string& name,
    string& value,
    bool mandatory
)
{
    return readAttributeImpl<string, std::string>(name, value, mandatory);
}


inline bool Foam::adiosFoam::adiosReader::readIntAttributeIfPresent
(
    const string& name,
    label& value
)
{
    return readIntAttribute(name, value, false);
}


inline bool Foam::adiosFoam::adiosReader::readScalarAttributeIfPresent
(
    const string& name,
    double& value
)
{
    return readScalarAttribute(name, value, false);
}


inline bool Foam::adiosFoam::adiosReader::readStringAttributeIfPresent
(
    const string& name,
    string& value
)
{
    return readStringAttribute(name, value, false);
}


inline bool Foam::adiosFoam::adiosReader::readIntListAttributeIfPresent
(
    const string& name,
    List<label>& list
)
{
    return readListAttributeImpl<label, int>(name, list, false);
}


inline bool Foam::adiosFoam::adiosReader::readScalarListAttributeIfPresent
(
    const string& name,
    List<double>& list
)
{
    return readListAttributeImpl<double, double>(name, list, false);
}


inline Foam::label
Foam::adiosFoam::adiosReader::getIntAttribute(const string& name)
{
    label value;

    if (!readIntAttributeIfPresent(name, value))
    {
        FatalErrorInFunction
            << "integer attribute missing: " << name
            << exit(FatalIOError);
    }

    return value;
}


inline double
Foam::adiosFoam::adiosReader::getDoubleAttribute(const string& name)
{
    double value;

    if (!readScalarAttributeIfPresent(name, value))
    {
        FatalErrorInFunction
            << "double attribute missing: " << name
            << exit(FatalIOError);
    }

    return value;
}


inline Foam::string
Foam::adiosFoam::adiosReader::getStringAttribute(const string& name)
{
    string value;

    if (!readStringAttributeIfPresent(name, value))
    {
        FatalErrorInFunction
            << "string attribute missing: " << name
            << exit(FatalIOError);
    }

    return value;
}


inline Foam::List<Foam::label>
Foam::adiosFoam::adiosReader::getIntListAttribute(const string& name)
{
    List<label> value;

    if (!readIntListAttributeIfPresent(name, value))
    {
        FatalErrorInFunction
            << "int-list attribute missing: " << name
            << exit(FatalIOError);
    }

    return value;
}


inline Foam::List<double>
Foam::adiosFoam::adiosReader::getDoubleListAttribute(const string& name)
{
    List<double> value;

    if (!readScalarListAttributeIfPresent(name, value))
    {
        FatalErrorInFunction
            << "double-list attribute missing: " << name
            << exit(FatalIOError);
    }

    return value;
}


template<class StringType>
inline Foam::List<StringType> Foam::adiosFoam::adiosReader::getStringListAttribute
(
    const string& name
)
{
    List<StringType> value;

    if (!readStringListAttributeIfPresent(name, value))
    {
        FatalErrorInFunction
            << "string-list attribute missing: " << name
            << exit(FatalIOError);
    }

    return value;
}



inline bool Foam::adiosFoam::adiosReader::readIntVariableIfPresent
(
    const string& name,
    label& value
)
{
    return readVariableImpl<label, int>(name, value, false);
}


inline bool Foam::adiosFoam::adiosReader::readScalarVariableIfPresent
(
    const string& name,
    scalar& value
)
{
    return readVariableImpl<scalar, double>(name, value, false);
}


inline Foam::label
Foam::adiosFoam::adiosReader::getIntVariable(const string& name)
{
    label value;

    if (!readIntVariableIfPresent(name, value))
    {
        FatalErrorInFunction
            << "integer variable missing: " << name
            << exit(FatalIOError);
    }

    return value;
}


inline Foam::scalar
Foam::adiosFoam::adiosReader::getScalarVariable(const string& name)
{
    scalar value;

    if (!readScalarVariableIfPresent(name, value))
    {
        FatalErrorInFunction
            << "scalar value missing: " << name
            << exit(FatalIOError);
    }

    return value;
}


// ************************************************************************* //
