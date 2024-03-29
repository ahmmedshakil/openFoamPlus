/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Released 2008-2011 OpenCFD Ltd.
    Copyright (C) 2011 OpenFOAM Foundation
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
    Foam::fileFormats::surfaceFormatsCore

Description
    A collection of helper functions for reading/writing surface formats.

SourceFiles
    surfaceFormatsCore.C

\*---------------------------------------------------------------------------*/

#ifndef surfaceFormatsCore_H
#define surfaceFormatsCore_H

#include "Map.H"
#include "HashSet.H"
#include "labelList.H"
#include "surfZoneList.H"
#include "surfZoneIdentifierList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class ISstream;
class Time;

namespace fileFormats
{

/*---------------------------------------------------------------------------*\
               Class fileFormats::surfaceFormatsCore Declaration
\*---------------------------------------------------------------------------*/

class surfaceFormatsCore
{
protected:

        //- Read non-empty and non-comment line
        static string getLineNoComment(ISstream& is, const char comment='#');


        //- Return a surfZone list with a single entry, the size of which
        //- corresponds to that of the container
        template<class Container>
        static List<surfZone> oneZone
        (
            const Container& container,
            const word& name = "zone0"
        )
        {
            return List<surfZone>(1, surfZone(name, container.size(), 0, 0));
        }

public:

    // Static Data

        //- The file extension corresponding to 'native' surface format
        //  Normally "ofs" (mnemonic: OF = OpenFOAM, S = Surface)
        static word nativeExt;


    // Static Member Functions

        //- Helper function when checking if a file extension is supported.
        static bool checkSupport
        (
            const wordHashSet& available,
            const word& ext,
            const bool verbose,
            const word& functionName
        );

        //- Return the local file name (within time directory)
        //  NEEDS FIXING
        static fileName localMeshFileName
        (
            const word& surfName = word::null
        );

        //- Find instance with surfName
        //  NEEDS FIXING
        static fileName findMeshInstance
        (
            const Time&,
            const word& surfName = word::null
        );

        //- Find mesh file with surfName
        //  NEEDS FIXING
        static fileName findMeshFile
        (
            const Time&,
            const word& surfName = word::null
        );


    // Constructors

        //- Construct null
        surfaceFormatsCore() = default;


    //- Destructor
    virtual ~surfaceFormatsCore() = default;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fileFormats
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
