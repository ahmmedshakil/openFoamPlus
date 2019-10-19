/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2019 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2013-2016 OpenFOAM Foundation
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
    Foam::functionObjects::CourantNo

Group
    grpFieldFunctionObjects

Description
    This function object calculates and outputs the Courant number, i.e. \c Co:

    \f[
        \text{Co} =
            \frac{\delta t}{2}
            \frac{\sum_\text{facei} | \phi_\text{facei} |}{V_\text{cell} \rho}
    \f]
    where
    * \c Co the Courant number as a \c volScalarField
    * \f$ \delta t\f$ the current time-step of the simulation
    * \f$ \phi_\text{facei} \f$ the volumetric flux through face \c facei
    * \f$ \text{facei} \f$ the index of the given cell face
    * \f$ V_\text{cell} \f$ the cell volume
    * \f$ \rho \f$ the fluid density

    Reference:
    \verbatim
        Standard expression:
            Courant, R., Friedrichs, K., & Lewy, H. (1967).
            On the partial difference equations of mathematical physics.
            IBM Journal of Research and Development, 11(2), 215-234.
            doi:10.1147/rd.112.0215
    \endverbatim

Usage
    Example of function object specification to calculate the Courant number:
    \verbatim
    Co1
    {
        // Mandatory entries
        type        CourantNo;
        libs        ("libfieldFunctionObjects.so");
    
        // Optional fieldExpression entries
    }
    \endverbatim

    Where the entries comprise:
    \table
        Property     | Description             | Required    | Default value
        type         | Type name: CourantNo    | yes         |
        rho          | Name of density field   | no          | rho
        field        | Name of flux field      | no          | phi
        result       | Name of Courant number field | no     | Co
        log          | Log to standard output  | no          | yes
    \endtable

Note
    The field is stored on the mesh database using the default name \c Co
    so that it can be retrieved and used for other applications.

See also
    Foam::functionObjects::fieldExpression
    Foam::functionObjects::fvMeshFunctionObject

SourceFiles
    CourantNo.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_CourantNo_H
#define functionObjects_CourantNo_H

#include "fieldExpression.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                          Class CourantNo Declaration
\*---------------------------------------------------------------------------*/

class CourantNo
:
    public fieldExpression
{
    // Private Data

        //- Name of density field (optional)
        word rhoName_;


    // Private Member Functions

        //- Divide the Courant number by rho if required
        tmp<volScalarField::Internal> byRho
        (
            const tmp<volScalarField::Internal>& Co
        ) const;

        //- Calculate the Courant number field and return true if successful
        virtual bool calc();


public:

    //- Runtime type information
    TypeName("CourantNo");


    // Constructors

        //- Construct from Time and dictionary
        CourantNo
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~CourantNo() = default;


    // Member Functions

        //- Read the CourantNo data
        virtual bool read(const dictionary&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
