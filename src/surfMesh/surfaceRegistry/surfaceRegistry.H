/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Released 2009-2011 OpenCFD Ltd.
    Copyright (C) 2011 OpenFOAM Foundation
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

Class
    Foam::surfaceRegistry

Description
    Wraps the normal objectRegistry with a local instance for %surfaces.

SourceFiles
    surfaceRegistry.C

\*---------------------------------------------------------------------------*/

#ifndef surfaceRegistry_H
#define surfaceRegistry_H

#include "objectRegistry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class surfaceRegistry Declaration
\*---------------------------------------------------------------------------*/

class surfaceRegistry
:
    public objectRegistry
{

    // Private Member Functions

        //- No copy construct
        surfaceRegistry(const surfaceRegistry&) = delete;

        //- No copy assignment
        surfaceRegistry& operator=(const surfaceRegistry&) = delete;


public:

    //- Runtime type information
    TypeName("surfaceRegistry");

        //- The prefix to local: %surfaces
        static const word prefix;

        //- The default surface name: %default
        static word defaultName;


    // Constructors

        //- Construct for the given objectRegistry and named surface
        surfaceRegistry
        (
            const objectRegistry& obr,
            const word& surfName = word::null
        );


    //- Destructor
    virtual ~surfaceRegistry() = default;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
