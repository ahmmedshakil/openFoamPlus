/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2016 OpenFOAM Foundation
    Modified code Copyright (C) 2016 OpenCFD Ltd.
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

InClass
    Foam::particle

Description
    Macros for adding to particle property lists

\*---------------------------------------------------------------------------*/

#ifndef particleMacros_H
#define particleMacros_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Define a static 'propertyList' for particle properties
//  The property list is space-delimited with brackets for vector groupings
// \sa AddToPropertyList
#define DefinePropertyList(str)                                                \
                                                                               \
    static string propertyList_;                                               \
                                                                               \
    static string propertyList()                                               \
    {                                                                          \
        return str;                                                            \
    }


//- Add to existing static 'propertyList' for particle properties
//  The property list is space-delimited with brackets for vector groupings
// \sa DefinePropertyList
#define AddToPropertyList(ParcelType, str)                                     \
                                                                               \
    static string propertyList_;                                               \
                                                                               \
    static string propertyList()                                               \
    {                                                                          \
        return ParcelType::propertyList() + str;                               \
    }


//- Define a static 'propertyTypes' for the types of particle properties
//  Brace brackets are used to delimit binary write groups
// \sa AddToPropertyTypes
#define DefinePropertyTypes(str)                                               \
                                                                               \
    static string propertyTypes_;                                              \
                                                                               \
    static string propertyTypes()                                              \
    {                                                                          \
        return str;                                                            \
    }


//- Add to existing static 'propertyTypes' for the types of particle properties
//  Brace brackets are used to delimit binary write groups
// \sa AddToPropertyTypes
#define AddToPropertyTypes(ParcelType, str)                                    \
                                                                               \
    static string propertyTypes_;                                              \
                                                                               \
    static string propertyTypes()                                              \
    {                                                                          \
        return ParcelType::propertyTypes() + str;                              \
    }


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //
