/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015 OpenFOAM Foundation
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
    Foam::LESModels::kOmegaSSTIDDES

Group
    grpDESTurbulence

Description
    k-omega-SST IDDES turbulence model for incompressible and compressible
    flows

    Reference:
    \verbatim
        Gritskevich, M.S., Garbaruk, A.V., Schuetze, J., Menter, F.R. (2011)
        Development of DDES and IDDES Formulations for the k-omega
        Shear Stress Transport Model, Flow, Turbulence and Combustion,
        pp. 1-19
    \endverbatim

SourceFiles
    kOmegaSSTIDDES.C

\*---------------------------------------------------------------------------*/

#ifndef kOmegaSSTIDDES_H
#define kOmegaSSTIDDES_H

#include "kOmegaSSTDES.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
                       class kOmegaSSTIDDES Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class kOmegaSSTIDDES
:
    public kOmegaSSTDES<BasicTurbulenceModel>
{
    // Private Member Functions

        //- Check that the supplied delta is an IDDESDelta
        const IDDESDelta& setDelta() const;

        tmp<volScalarField> alpha() const;
        tmp<volScalarField> ft(const volScalarField& magGradU) const;
        tmp<volScalarField> fl(const volScalarField& magGradU) const;

        tmp<volScalarField> rd
        (
            const volScalarField& nur,
            const volScalarField& magGradU
        ) const;

        //- Delay function
        tmp<volScalarField> fdt(const volScalarField& magGradU) const;

        //- No copy construct
        kOmegaSSTIDDES(const kOmegaSSTIDDES&) = delete;

        //- No copy assignment
        void operator=(const kOmegaSSTIDDES&) = delete;


protected:

    // Protected data

        // Model coefficients

            dimensionedScalar Cdt1_;
            dimensionedScalar Cdt2_;
            dimensionedScalar Cl_;
            dimensionedScalar Ct_;

        // Fields

            const IDDESDelta& IDDESDelta_;

    // Protected Member Functions

        //- Length scale
        virtual tmp<volScalarField> dTilda
        (
            const volScalarField& magGradU,
            const volScalarField& CDES
        ) const;


public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;


    //- Runtime type information
    TypeName("kOmegaSSTIDDES");


    // Constructors

        //- Construct from components
        kOmegaSSTIDDES
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& propertiesName = turbulenceModel::propertiesName,
            const word& type = typeName
        );


    //- Destructor
    virtual ~kOmegaSSTIDDES()
    {}


    // Member Functions

        //- Re-read model coefficients if they have changed
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "kOmegaSSTIDDES.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
