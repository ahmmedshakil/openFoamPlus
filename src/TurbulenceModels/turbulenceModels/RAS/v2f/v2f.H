/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2012-2016 OpenFOAM Foundation
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
    Foam::RASModels::v2f

Group
    grpRASTurbulence

Description
    The v2-f turbulence closure model for incompressible and compressible flows.

    The model is a three-transport-equation linear-eddy-viscosity turbulence
    closure model alongside an elliptic relaxation equation:
      - Turbulent kinetic energy, \c k
      - Turbulent kinetic energy dissipation rate, \c epsilon
      - Wall-normal fluctuating velocity scale, \c v2
      - Elliptic relaxation factor, \c f

    The model employs the variant (i.e. N=6) ensuring zero \c f on walls.

    The set of wall boundary conditions can be:

        k       = kLowReWallFunction
        epsilon = epsilonWallFunction
        v2      = v2WallFunction
        f       = fWallFunction

    These are applicable to both low- and high-Reynolds number flows.

    Inlet values can be approximated by:

        v2      = 2/3 k
        f       = zero-gradient

    References:
    \verbatim
        Standard model (Tag:LK):
            Lien, F. S., & Kalitzin, G. (2001).
            Computations of transonic flow with the v2-f turbulence model.
            International Journal of Heat and Fluid Flow, 22(1), 53-61.
            doi:10.1016/S0142-727X(00)00073-4

        Typo correction for Eq. 4 (Tag:BL):
            Billard, F., & Laurence, D. (2012).
            A robust k−ε−v2/k elliptic blending turbulence model applied to
            near-wall, separated and buoyant flows.
            International Journal of Heat and Fluid Flow, 33(1), 45-58.
            doi:10.1016/j.ijheatfluidflow.2011.11.003

        Limiter on the turbulent viscosity and v2 (Modification-I) (Tag:D):
            Davidson, L., Nielsen, P., & Sveningsson, A. (2003).
            Modifications of the v2-f model for computing the flow in a
            3D wall jet.
            In Hanjalić, K., Nagano, Y., & Tummers, M. J. (Eds.),
            Proceedings of the 4th International Symposium on Turbulence, Heat
            and Mass Transfer, 577-584, October 12-17, 2003, Antalya, Turkey

        For the Rapid Distortion Theory-based (RDT) compression term (E):
            El Tahry, S. H. (1983).
            k-epsilon equation for compressible reciprocating engine flows.
            Journal of Energy, 7(4), 345-353.
            doi:10.2514/3.48086
    \endverbatim

    The default model coefficients are (LK:Eq. 17; BL:Table 5):
    \verbatim
        v2fCoeffs
        {
            Cmu         0.22;
            CmuKEps     0.09;
            C1          1.4;
            C2          0.3;
            Ceta        70.0;
            Ceps1a      1.4;
            Ceps1b      1.0;
            Ceps1c      0.05;
            Ceps2       1.9;
            CRDT       -0.33;
            CT          6.0;
            CL          0.23;
            sigmaK      1.0;
            sigmaEps    1.3;
            sigmaV2     1.0;
        }
    \endverbatim

Note
    If the kLowReWallFunction is employed, a velocity variant of the turbulent
    viscosity wall function should be used, e.g. nutUWallFunction.  Turbulence
    k variants (nutk...) for this case will not behave correctly.

See also
    Foam::RASModels::v2fBase
    Foam::RASModels::kEpsilon
    Foam::kLowReWallFunctionFvPatchScalarField
    Foam::epsilonWallFunctionFvPatchScalarField
    Foam::v2WallFunctionFvPatchScalarField
    Foam::fWallFunctionFvPatchScalarField

SourceFiles
    v2f.C

\*---------------------------------------------------------------------------*/

#ifndef v2f_H
#define v2f_H

#include "v2fBase.H"
#include "RASModel.H"
#include "eddyViscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                             Class v2f Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class v2f
:
    public eddyViscosity<RASModel<BasicTurbulenceModel>>,
    public v2fBase
{
    // Private Member Functions

        //- No copy construct
        v2f(const v2f&) = delete;

        //- No copy assignment
        void operator=(const v2f&) = delete;


protected:

    // Protected Data

        // Model coefficients

            dimensionedScalar Cmu_;
            dimensionedScalar CmuKEps_;
            dimensionedScalar C1_;
            dimensionedScalar C2_;
            dimensionedScalar Ceta_;
            dimensionedScalar Ceps1a_;
            dimensionedScalar Ceps1b_;
            dimensionedScalar Ceps1c_;
            dimensionedScalar Ceps2_;
            dimensionedScalar CRDT_;
            dimensionedScalar CT_;
            dimensionedScalar CL_;
            dimensionedScalar sigmaK_;
            dimensionedScalar sigmaEps_;
            dimensionedScalar sigmaV2_;


        // Fields

            //- Turbulent kinetic energy [m2/s2]
            volScalarField k_;

            //- Turbulent kinetic energy dissipation rate [m2/s3]
            volScalarField epsilon_;

            //- The wall-normal fluctuating velocity scale [m2/s2]
            volScalarField v2_;

            //- Elliptic relaxation factor [1/s]
            volScalarField f_;

            //- Turbulent time scale [s]
            volScalarField T_;


        // Bounding values

            dimensionedScalar v2Min_;
            dimensionedScalar fMin_;
            dimensionedScalar TMin_;
            dimensionedScalar L2Min_;


    // Protected Member Functions

        //- Update nut with the latest available k, epsilon, v2 and T
        virtual void correctNut();

        //- Return the turbulent time scale
        tmp<volScalarField> Ts
        (
            const volScalarField& C
        ) const;

        //- Return the turbulent length scale
        tmp<volScalarField> Ls
        (
            const volScalarField& C
        ) const;


public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;


    //- Runtime type information
    TypeName("v2f");


    // Constructors

        //- Construct from components
        v2f
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
    virtual ~v2f() = default;


    // Member Functions

        //- Re-read model coefficients if they have changed
        virtual bool read();

        //- Return the effective diffusivity for k
        tmp<volScalarField> DkEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    "DkEff",
                    this->nut_/sigmaK_ + this->nu()
                )
            );
        }

        //- Return the effective diffusivity for epsilon
        tmp<volScalarField> DepsilonEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    "DepsilonEff",
                    this->nut_/sigmaEps_ + this->nu()
                )
            );
        }

        //- Return the effective diffusivity for v2
        tmp<volScalarField> Dv2Eff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    "Dv2Eff",
                    this->nut_/sigmaV2_ + this->nu()
                )
            );
        }

        //- Return the turbulent kinetic energy field
        virtual tmp<volScalarField> k() const
        {
            return k_;
        }

        //- Return the turbulent kinetic energy dissipation rate field
        virtual tmp<volScalarField> epsilon() const
        {
            return epsilon_;
        }

        //- Return the wall-normal fluctuating velocity scale field
        virtual tmp<volScalarField> v2() const
        {
            return v2_;
        }

        //- Return the elliptic relaxation factor field
        virtual tmp<volScalarField> f() const
        {
            return f_;
        }

        //- Solve the transport equations and correct the turbulent viscosity
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "v2f.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
