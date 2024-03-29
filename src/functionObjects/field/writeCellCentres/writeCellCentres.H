/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016 OpenFOAM Foundation
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
    Foam::functionObjects::writeCellCentres

Group
    grpFieldFunctionObjects

Description
    Writes the cell-centres volVectorField and the three component fields as
    volScalarFields.

Usage
    Example of function object specification:
    \verbatim
    writeCellCentres
    {
        type        writeCellCentres;
        libs        ("libfieldFunctionObjects.so");
        ...
    }
    \endverbatim

    Where the entries comprise:
    \table
        Property  | Description                 | Required    | Default value
        type      | type name: writeCellCentres | yes       |
    \endtable

See also
    Foam::functionObjects::fvMeshFunctionObject

SourceFiles
    writeCellCentres.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_writeCellCentres_H
#define functionObjects_writeCellCentres_H

#include "fvMeshFunctionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class writeCellCentres Declaration
\*---------------------------------------------------------------------------*/

class writeCellCentres
:
    public fvMeshFunctionObject
{
    // Private member functions

        //- No copy construct
        writeCellCentres(const writeCellCentres&) = delete;

        //- No copy assignment
        void operator=(const writeCellCentres&) = delete;


public:

    //- Runtime type information
    TypeName("writeCellCentres");


    // Constructors

        //- Construct from Time and dictionary
        writeCellCentres
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~writeCellCentres();


    // Member Functions

        //- Read the cell-centre rate data
        virtual bool read(const dictionary&);

        //- Do nothing
        virtual bool execute();

        //- Write the cell-centre fields
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
