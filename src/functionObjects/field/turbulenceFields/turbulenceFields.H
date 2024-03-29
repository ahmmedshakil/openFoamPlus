/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2016 OpenFOAM Foundation
    Modified code Copyright (C) 2015-2018 OpenCFD Ltd.
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
    Foam::functionObjects::turbulenceFields

Group
    grpFieldFunctionObjects

Description
    Stores turbulence fields on the mesh database for further manipulation.

    Fields are stored as copies of the original, with the prefix
    "tubulenceModel:", e.g.:

    \verbatim
    turbulenceModel:R
    \endverbatim

Usage
    Example of function object specification:
    \verbatim
    turbulenceFields1
    {
        type        turbulenceFields;
        libs        ("libfieldFunctionObjects.so");
        ...
        fields
        (
            R
            devRhoReff
        );
    }
    \endverbatim

    Where the entries comprise:
    \table
        Property     | Description                 | Required | Default value
        type         | type name: turbulenceFields | yes      |
        fields       | fields to store (see below) | yes      |
    \endtable

    Where \c fields can include:
    \plaintable
        k           | turbulence kinetic energy
        epsilon     | turbulence kinetic energy dissipation rate
        omega       | turbulence specific dissipation rate
        nuTilda     | turbulence modified viscosity
        nut         | turbulence viscosity (incompressible)
        nuEff       | effective turbulence viscosity (incompressible)
        mut         | turbulence viscosity (compressible)
        muEff       | effective turbulence viscosity (compressible)
        alphat      | turbulence thermal diffusivity (compressible)
        alphaEff    | effective turbulence thermal diffusivity (compressible)
        R           | Reynolds stress tensor
        devReff     | Deviatoric part of the effective Reynolds stress
        devRhoReff  | Divergence of the Reynolds stress
        L           | turbulence length scale
        I           | turbulence intensity
    \endplaintable

See also
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::timeControl

SourceFiles
    turbulenceFields.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_turbulenceFields_H
#define functionObjects_turbulenceFields_H

#include "fvMeshFunctionObject.H"
#include "HashSet.H"
#include "Enum.H"
#include "volFieldsFwd.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                      Class turbulenceFields Declaration
\*---------------------------------------------------------------------------*/

class turbulenceFields
:
    public fvMeshFunctionObject
{
public:

    enum compressibleField
    {
        cfK,
        cfEpsilon,
        cfOmega,
        cfNuTilda,
        cfMut,
        cfMuEff,
        cfAlphat,
        cfAlphaEff,
        cfR,
        cfDevRhoReff,
        cfL,
        cfI
    };
    static const Enum<compressibleField> compressibleFieldNames_;

    enum incompressibleField
    {
        ifK,
        ifEpsilon,
        ifOmega,
        ifNuTilda,
        ifNut,
        ifNuEff,
        ifR,
        ifDevReff,
        ifL,
        ifI
    };
    static const Enum<incompressibleField> incompressibleFieldNames_;

    static const word modelName;


protected:

    // Protected data

        //- Fields to load
        wordHashSet fieldSet_;


    // Protected Member Functions

        //- Return true if compressible turbulence model is identified
        bool compressible();

        //- Process the turbulence field
        template<class Type>
        void processField
        (
            const word& fieldName,
            const tmp<GeometricField<Type, fvPatchField, volMesh>>& tvalue
        );

        //- Return omega calculated from k and epsilon
        template<class Model>
        tmp<volScalarField> omega(const Model& model) const;

        //- Return nuTilda calculated from k and omega
        template<class Model>
        tmp<volScalarField> nuTilda(const Model& model) const;

        //- Return L calculated from k and epsilon
        template<class Model>
        tmp<volScalarField> L(const Model& model) const;

        //- Return I calculated from k and U
        template<class Model>
        tmp<volScalarField> I(const Model& model) const;


private:

    // Private member functions

        //- No copy construct
        turbulenceFields(const turbulenceFields&) = delete;

        //- No copy assignment
        void operator=(const turbulenceFields&) = delete;


public:

    //- Runtime type information
    TypeName("turbulenceFields");


    // Constructors

        //- Construct from Time and dictionary
        turbulenceFields
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~turbulenceFields();


    // Member Functions

        //- Read the controls
        virtual bool read(const dictionary&);

        //- Calculate turbulence fields
        virtual bool execute();

        //- Do nothing.
        //  The turbulence fields are registered and written automatically
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "turbulenceFieldsTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
