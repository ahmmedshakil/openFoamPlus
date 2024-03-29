/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Released 2004-2011 OpenCFD Ltd.
    Copyright (C) 2011 OpenFOAM Foundation
    Modified code Copyright (C) 2018 OpenCFD Ltd.
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
    Foam::cellToCell

Description
    A topoSetCellSource to select the cells from another cellSet.

    \heading Dictionary parameters
    \table
        Property    | Description                          | Required | Default
        sets        | The input cellSet names              | possibly |
        set         | The input cellSet name               | possibly |
    \endtable

Note
    Must specify "sets" or "set" (highest to lowest precedence).

SourceFiles
    cellToCell.C

\*---------------------------------------------------------------------------*/

#ifndef cellToCell_H
#define cellToCell_H

#include "topoSetCellSource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class cellToCell Declaration
\*---------------------------------------------------------------------------*/

class cellToCell
:
    public topoSetCellSource
{
    // Private data

        //- Add usage string
        static addToUsageTable usage_;

        //- Names of sets to use
        wordList names_;


public:

    //- Runtime type information
    TypeName("cellToCell");

    // Constructors

        //- Construct from components
        cellToCell(const polyMesh& mesh, const word& setName);

        //- Construct from dictionary
        cellToCell(const polyMesh& mesh, const dictionary& dict);

        //- Construct from Istream
        cellToCell(const polyMesh& mesh, Istream& is);


    //- Destructor
    virtual ~cellToCell() = default;


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
