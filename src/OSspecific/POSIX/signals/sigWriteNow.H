/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Modified code Copyright (C) 2016-2019 OpenCFD Ltd.
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
    Foam::sigWriteNow

Description
    Signal handler to write once and continue.
    The interrupt is defined by OptimisationSwitches::writeNowSignal

SourceFiles
    sigWriteNow.C

\*---------------------------------------------------------------------------*/

#ifndef sigWriteNow_H
#define sigWriteNow_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class Time;

/*---------------------------------------------------------------------------*\
                         Class sigWriteNow Declaration
\*---------------------------------------------------------------------------*/

class sigWriteNow
{
    // Private Data

        //- Signal number to use
        static int signal_;


    // Private Member Functions

        //- Handler for caught signals
        static void sigHandler(int);


public:

    //- Allow setter access to signal_
    friend class addwriteNowSignalToOpt;


    // Constructors

        //- Construct null
        sigWriteNow();

        //- Construct from components
        sigWriteNow(Time& runTime, bool verbose=false);


    //- Destructor
    ~sigWriteNow();


    // Member Functions

        //- Is active?
        static bool active();

        //- The signal number being used
        static int signalNumber();

        //- Set/reset signal handler
        static void set(bool verbose=false);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
