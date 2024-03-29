/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::functionObjects::writeObjects

Group
    grpUtilitiesFunctionObjects

Description
    Allows specification of different writing frequency of objects registered
    to the database.

    It has similar functionality as the main time database through the
    \c writeControl setting:
      - timeStep
      - writeTime
      - adjustableRunTime
      - runTime
      - clockTime
      - cpuTime

    It also has the ability to write the selected objects that were defined
    with the respective write mode for the requested \c writeOption, namely:
    \vartable
        autoWrite | objects set to write at output time
        noWrite   | objects set to not write by default
        anyWrite  | any option of the previous two
    \endvartable

Usage
    Example of function object specification:
    \verbatim
    writeObjects1
    {
        type        writeObjects;
        libs        ("libutilityFunctionObjects.so");
        ...
        objects     (obj1 obj2);
        writeOption anyWrite;
    }
    \endverbatim

    Where the entries comprise:
    \table
        Property     | Description             | Required     | Default value
        type         | type name: writeObjects | yes          |
        objects      | objects to write        | yes          |
        writeOption  | only those with this write option | no | anyWrite
    \endtable

    Note: Regular expressions can also be used in \c objects.

See also
    Foam::functionObject
    Foam::functionObjects::timeControl

SourceFiles
    writeObjects.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_writeObjects_H
#define functionObjects_writeObjects_H

#include "functionObject.H"
#include "wordRes.H"
#include "Enum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class objectRegistry;

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                Class functionObjects::writeObjects Declaration
\*---------------------------------------------------------------------------*/

class writeObjects
:
    public functionObject
{
public:

    // Public data types

        //- Re-enumeration defining the write options, based on the original
        //  ones at IOobject::writeOption
        enum writeOption
        {
            AUTO_WRITE,
            NO_WRITE,
            ANY_WRITE
        };

        static const Enum<writeOption> writeOptionNames_;

private:

    // Private data

        //- Reference to registry
        const objectRegistry& obr_;

        //- To only write objects of defined writeOption
        writeOption writeOption_;

        //- Names of objects to control
        wordRes objectNames_;


    // Private Member Functions

        //- No copy construct
        writeObjects(const writeObjects&) = delete;

        //- No copy assignment
        void operator=(const writeObjects&) = delete;


public:

    //- Runtime type information
    TypeName("writeObjects");


    // Constructors

        //- Construct from Time and dictionary
        writeObjects
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~writeObjects() = default;


    // Member Functions

        //- Read the writeObjects data
        virtual bool read(const dictionary&);

        //- Do nothing
        virtual bool execute();

        //- Write the registered objects
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
