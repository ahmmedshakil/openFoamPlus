/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2019 PCOpt/NTUA
    Modified code Copyright (C) 2013-2019 FOSS GP
    Modified code Copyright (C) 2019 OpenCFD Ltd.
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
    Foam::RASTurbulenceModel

Description
    Solves for a RAS turbulence model, with constant U and phi values

\*---------------------------------------------------------------------------*/

#ifndef RASTurbulenceModel_H
#define RASTurbulenceModel_H

#include "incompressiblePrimalSolver.H"
#include "SIMPLEControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class RASTurbulenceModel Declaration
\*---------------------------------------------------------------------------*/

class RASTurbulenceModel
:
    public incompressiblePrimalSolver
{
private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        RASTurbulenceModel(const RASTurbulenceModel&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const RASTurbulenceModel&) = delete;


protected:

    // Protected data

        //- Solver control
        autoPtr<SIMPLEControl> solverControl_;


public:


    // Static Data Members

        //- Run-time type information
        TypeName("RASTurbulenceModel");


    // Constructors

        //- Construct from mesh and dictionary
        RASTurbulenceModel
        (
            fvMesh& mesh,
            const word& managerType,
            const dictionary& dict
        );


    //- Destructor
    virtual ~RASTurbulenceModel() = default;


    // Member Functions

        // Evolution

            //- Execute one iteration of the solution algorithm
            virtual void solveIter();

            //- Main control loop
            virtual void solve();

            //- Looper (advances iters, time step)
            virtual bool loop();

            //- Read average iteration
            virtual bool writeData(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
