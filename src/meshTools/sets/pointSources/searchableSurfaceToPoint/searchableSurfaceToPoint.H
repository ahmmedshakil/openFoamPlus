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

Class
    Foam::searchableSurfaceToPoint

Description
    A topoSetPointSource to select mesh points within a searchableSurface.

    \heading Dictionary parameters
    \table
        Property    | Description                           | Required | Default
        surfaceType | The searchable surface type           | yes   |
        surfaceName | Name for the IOobject                 | no    | mesh-name
        surface     | Same as 'surfaceType'                 | no    |
    \endtable

SourceFiles
    searchableSurfaceToPoint.C

\*---------------------------------------------------------------------------*/

#ifndef searchableSurfaceToPoint_H
#define searchableSurfaceToPoint_H

#include "topoSetPointSource.H"
#include "searchableSurface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class searchableSurfaceToPoint Declaration
\*---------------------------------------------------------------------------*/

class searchableSurfaceToPoint
:
    public topoSetPointSource
{
    // Private Data

        //- Add usage string
        static addToUsageTable usage_;

        //- The searchableSurface
        autoPtr<searchableSurface> surf_;


    // Private Member Functions

        void combine(topoSet& set, const bool add) const;


protected:

    // Protected Member Functions

        //- Retrieve surface name from dictionary entry
        static word getSurfaceName
        (
            const dictionary& dict,
            const word& defaultName
        );


public:

    //- Runtime type information
    TypeName("searchableSurfaceToPoint");


    // Constructors

        //- Construct surface-type from dictionary
        searchableSurfaceToPoint
        (
            const word& surfaceType,
            const polyMesh& mesh,
            const dictionary& dict
        );

        //- Construct from dictionary
        searchableSurfaceToPoint
        (
            const polyMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~searchableSurfaceToPoint() = default;


    // Member Functions

        virtual void applyToSet
        (
            const topoSetSource::setAction action,
            topoSet& set
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
