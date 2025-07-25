/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018 OpenCFD Ltd.
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

#include "DirLister.H"
#include <dirent.h>

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

bool Foam::DirLister::const_iterator::open(const fileName& dir)
{
    if (!dir.empty())
    {
        dirptr_ = ::opendir(dir.c_str());
    }

    return dirptr_;
}


void Foam::DirLister::const_iterator::close()
{
    if (dirptr_)
    {
        ::closedir(dirptr_);
        dirptr_ = nullptr;
    }
}


Foam::word Foam::DirLister::next(DIR* dirPtr) const
{
    // Read next entry in the directory

    struct dirent *list;
    while (dirPtr && (list = ::readdir(dirPtr)) != nullptr)
    {
        const std::string item(list->d_name);

        // Ignore files/directories beginning with "."
        // These are the ".", ".." directories and any hidden files/dirs
        if (item.empty() || item[0] == '.')
        {
            continue;
        }

        // Validate filename without spaces, quotes, etc in the name.
        // No duplicate slashes to strip - dirent will not have them anyhow.

        word name(fileName::validate(item));
        if (name != item)
        {
            // ++nFailed_;
            continue;
        }

        bool ok = false;
        fileName::Type fType = fileName::UNDEFINED;

        if
        (
            (requestedType_ == fileName::DIRECTORY)
         || (requestedType_ == fileName::FILE && !fileName::isBackup(name))
        )
        {
            fType = (dirName_/name).type(followLink_);

            // A DIRECTORY or FILE was request, so only accept the same type
            ok = (requestedType_ == fType);
        }
        else if (requestedType_ == fileName::UNDEFINED)
        {
            fType = (dirName_/name).type(followLink_);

            // An UNDEFINED type was requested, so accept DIRECTORY or FILE
            ok =
            (
                (fType == fileName::DIRECTORY)
             || (fType == fileName::FILE && !fileName::isBackup(name))
            );
        }

        if (ok)
        {
            if (fType == fileName::FILE && stripgz_ && name.has_ext("gz"))
            {
                name.remove_ext();
            }

            if (!name.empty() && accept(name))
            {
                return name;
            }
        }
    }

    return word();
}


// * * * * * * * * * * * * * * * Const Iterator  * * * * * * * * * * * * * * //

bool Foam::DirLister::const_iterator::next()
{
    if (lister_ && dirptr_)
    {
        name_ = lister_->next(dirptr_);
        if (name_.empty())
        {
            close();
        }
    }

    return dirptr_;
}


// ************************************************************************* //
