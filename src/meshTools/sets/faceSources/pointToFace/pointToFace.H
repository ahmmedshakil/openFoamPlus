/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Released 2004-2011 OpenCFD Ltd.
    Copyright (C) 2011-2012 OpenFOAM Foundation
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
    Foam::pointToFace

Description
    A topoSetFaceSource to select faces based on use of points.

    \heading Dictionary parameters
    \table
        Property    | Description                       | Required | Default
        sets        | The pointSet names to use         | possibly |
        set         | The pointSet name to use          | possibly |
        option      | Selection type (all/any/edge)     | yes   |
    \endtable

Note
    Must specify "sets" or "set" (highest to lowest precedence).

SourceFiles
    pointToFace.C

\*---------------------------------------------------------------------------*/

#ifndef pointToFace_H
#define pointToFace_H

#include "topoSetFaceSource.H"
#include "Enum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class pointToFace Declaration
\*---------------------------------------------------------------------------*/

class pointToFace
:
    public topoSetFaceSource
{

public:

        //- Enumeration defining the valid options
        enum pointAction
        {
            ANY,
            ALL,
            EDGE
        };


private:

        //- Add usage string
        static addToUsageTable usage_;

        static const Enum<pointAction> pointActionNames_;

        //- Names of sets to use
        wordList names_;

        //- Option
        pointAction option_;


    // Private Member Functions

        //- Depending on face to cell option add to or delete from cellSet.
        void combine(topoSet& set, const bool add, const word& setName) const;


public:

    //- Runtime type information
    TypeName("pointToFace");

    // Constructors

        //- Construct from components
        pointToFace
        (
            const polyMesh& mesh,
            const word& setName,
            const pointAction option
        );

        //- Construct from dictionary
        pointToFace(const polyMesh& mesh, const dictionary& dict);

        //- Construct from Istream
        pointToFace(const polyMesh& mesh, Istream& is);


    //- Destructor
    virtual ~pointToFace() = default;


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
