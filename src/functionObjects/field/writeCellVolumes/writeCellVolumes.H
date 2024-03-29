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
    Foam::functionObjects::writeCellVolumes

Group
    grpFieldFunctionObjects

Description
    Writes the cell-volumes volScalarField.

Usage
    Example of function object specification:
    \verbatim
    writeCellVolumes
    {
        type        writeCellVolumes;
        libs        ("libfieldFunctionObjects.so");
        ...
    }
    \endverbatim

    Where the entries comprise:
    \table
        Property  | Description                 | Required    | Default value
        type      | type name: writeCellVolumes | yes       |
    \endtable

See also
    Foam::functionObjects::fvMeshFunctionObject

SourceFiles
    writeCellVolumes.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_writeCellVolumes_H
#define functionObjects_writeCellVolumes_H

#include "fvMeshFunctionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class writeCellVolumes Declaration
\*---------------------------------------------------------------------------*/

class writeCellVolumes
:
    public fvMeshFunctionObject
{
    // Private member functions

        //- No copy construct
        writeCellVolumes(const writeCellVolumes&) = delete;

        //- No copy assignment
        void operator=(const writeCellVolumes&) = delete;


public:

    //- Runtime type information
    TypeName("writeCellVolumes");


    // Constructors

        //- Construct from Time and dictionary
        writeCellVolumes
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~writeCellVolumes();


    // Member Functions

        //- Read the cell-volume data
        virtual bool read(const dictionary&);

        //- Do nothing
        virtual bool execute();

        //- Write the cell-volume fields
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
