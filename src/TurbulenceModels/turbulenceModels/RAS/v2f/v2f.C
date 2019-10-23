/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2012-2017 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "v2f.H"
#include "fvOptions.H"
#include "bound.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
tmp<volScalarField> v2f<BasicTurbulenceModel>::Ts
(
    const volScalarField& C
) const
{
    // SAF: limiting thermo->nu(). If psiThermo is used rho might be < 0
    // temporarily when p < 0 then nu < 0 which needs limiting
    // (BL:Eq. 8 & Table 5)
    return
        max
        (
            min
            (
                k_/epsilon_,
                0.6*k_/C
            ),
            CT_*sqrt
            (
                max
                (
                    this->nu(),
                    dimensionedScalar(this->nu()().dimensions(), Zero)
                )/epsilon_
            )
        );
}


template<class BasicTurbulenceModel>
tmp<volScalarField> v2f<BasicTurbulenceModel>::Ls
(
    const volScalarField& C
) const
{
    // SAF: limiting thermo->nu(). If psiThermo is used rho might be < 0
    // temporarily when p < 0 then nu < 0 which needs limiting
    // (BL:Eq. 8 & Table 5)
    return
        CL_*max
        (
            min
            (
                pow(k_, 1.5)/epsilon_,
                pow(k_, 1.5)/C
            ),
            Ceta_*pow025
            (
                pow3
                (
                    max
                    (
                        this->nu(),
                        dimensionedScalar(this->nu()().dimensions(), Zero)
                    )
                )/epsilon_
            )
        );
}


template<class BasicTurbulenceModel>
void v2f<BasicTurbulenceModel>::correctNut()
{
    // (D:Eq. 7)
    this->nut_ =
        min
        (
            CmuKEps_*sqr(k_)/epsilon_,
            Cmu_*v2_*T_
        );
    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_).correct(this->nut_);

    BasicTurbulenceModel::correctNut();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
v2f<BasicTurbulenceModel>::v2f
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName,
    const word& type
)
:
    eddyViscosity<RASModel<BasicTurbulenceModel>>
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName
    ),
    v2fBase(),

    Cmu_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Cmu",
            this->coeffDict_,
            0.22
        )
    ),
    CmuKEps_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "CmuKEps",
            this->coeffDict_,
            0.09
        )
    ),
    C1_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "C1",
            this->coeffDict_,
            1.4
        )
    ),
    C2_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "C2",
            this->coeffDict_,
            0.3
        )
    ),
    Ceta_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Ceta",
            this->coeffDict_,
            70.0
        )
    ),
    Ceps1a_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Ceps1a",
            this->coeffDict_,
            1.4
        )
    ),
    Ceps1b_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Ceps1b",
            this->coeffDict_,
            1.0
        )
    ),
    Ceps1c_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Ceps1c",
            this->coeffDict_,
            0.050
        )
    ),
    Ceps2_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Ceps2",
            this->coeffDict_,
            1.9
        )
    ),
    CRDT_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "CRDT",
            this->coeffDict_,
            -0.33
        )
    ),
    CT_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "CT",
            this->coeffDict_,
            6.0
        )
    ),
    CL_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "CL",
            this->coeffDict_,
            0.23
        )
    ),
    sigmaK_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "sigmaK",
            this->coeffDict_,
            1.0
        )
    ),
    sigmaEps_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "sigmaEps",
            this->coeffDict_,
            1.3
        )
    ),
    sigmaV2_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "sigmaV2",
            this->coeffDict_,
            1.0
        )
    ),

    k_
    (
        IOobject
        (
            IOobject::groupName("k", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    epsilon_
    (
        IOobject
        (
            IOobject::groupName("epsilon", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    v2_
    (
        IOobject
        (
            IOobject::groupName("v2", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    f_
    (
        IOobject
        (
            IOobject::groupName("f", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    T_
    (
        IOobject
        (
            "T",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        this->mesh_,
        dimensionedScalar(dimTime, Zero)
    ),

    v2Min_(dimensionedScalar("v2Min", v2_.dimensions(), SMALL)),
    fMin_(dimensionedScalar("fMin", f_.dimensions(), Zero)),
    TMin_(dimensionedScalar("TMin", dimTime, SMALL)),
    L2Min_(dimensionedScalar("L2Min", sqr(dimLength), SMALL))
{
    bound(k_, this->kMin_);
    bound(epsilon_, this->epsilonMin_);
    bound(v2_, v2Min_);
    bound(f_, fMin_);

    if (type == typeName)
    {
        this->printCoeffs(type);
    }

    if
    (
        mag(sigmaK_.value()) < VSMALL
     || mag(sigmaEps_.value()) < VSMALL
     || mag(sigmaV2_.value()) < VSMALL
    )
    {
        FatalErrorInFunction
            << "Tiny magnitude for a model constant:" << nl
            << "sigmaK = " << sigmaK_ << nl
            << "sigmaEps = " << sigmaEps_ << nl
            << "sigmaV2 = " << sigmaV2_ << nl
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool v2f<BasicTurbulenceModel>::read()
{
    if (eddyViscosity<RASModel<BasicTurbulenceModel>>::read())
    {
        Cmu_.readIfPresent(this->coeffDict());
        CmuKEps_.readIfPresent(this->coeffDict());
        C1_.readIfPresent(this->coeffDict());
        C2_.readIfPresent(this->coeffDict());
        Ceta_.readIfPresent(this->coeffDict());
        Ceps1a_.readIfPresent(this->coeffDict());
        Ceps1b_.readIfPresent(this->coeffDict());
        Ceps1c_.readIfPresent(this->coeffDict());
        Ceps2_.readIfPresent(this->coeffDict());
        CRDT_.readIfPresent(this->coeffDict());
        CT_.readIfPresent(this->coeffDict());
        CL_.readIfPresent(this->coeffDict());
        sigmaK_.readIfPresent(this->coeffDict());
        sigmaEps_.readIfPresent(this->coeffDict());
        sigmaV2_.readIfPresent(this->coeffDict());

        return true;
    }

    return false;
}


template<class BasicTurbulenceModel>
void v2f<BasicTurbulenceModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }

    // Construct local convenience references
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    const volVectorField& U = this->U_;
    volScalarField& nut = this->nut_;
    fv::options& fvOptions(fv::options::New(this->mesh_));

    eddyViscosity<RASModel<BasicTurbulenceModel>>::correct();

    const volScalarField::Internal divU
    (
        fvc::div
        (
            fvc::absolute(this->phi(), U)
        )
    );

    // (LK:p. 54; D:Eq. 2)
    tmp<volTensorField> tgradU(fvc::grad(U));

    const volSymmTensorField S("Sij", symm(tgradU()));

    const volScalarField G(this->GName(), nut*(2.0*(dev(S) && tgradU)));

    // (LK:Eqs. 4-5)
    volScalarField C
    (
        "C",
        (sqrt(6.0)*Cmu_)*(sqrt(S && S)*v2_) + this->epsilonMin_
    );

    // (LK:Eq. 4)
    T_ = this->Ts(C);
    bound(T_, TMin_);

    // (LK:Eq. 5)
    const volScalarField L2(type() + ":L2", sqr(Ls(C)) + L2Min_);

    // (LK:Eq. 15; D:Eq. 6)
    const volScalarField::Internal v2fAlpha
    (
        type() + ":alpha",
        1.0/T_()*((C1_ - 6.0)*v2_() - 2.0/3.0*k_()*(C1_ - 1.0))
    );

    // (LK:Eq. 17)
    const volScalarField::Internal Ceps1
    (
        "Ceps1",
        Ceps1a_*(Ceps1b_ + Ceps1c_*sqrt(k_()/v2_()))
    );


    // Update epsilon (and possibly G) at the wall
    epsilon_.boundaryFieldRef().updateCoeffs();

    // Turbulent kinetic energy dissipation rate equation (LK:Eq. 7)
    // k/T ~ epsilon
    tmp<fvScalarMatrix> epsEqn
    (
        fvm::ddt(alpha, rho, epsilon_)
      + fvm::div(alphaRhoPhi, epsilon_)
      - fvm::laplacian(alpha*rho*DepsilonEff(), epsilon_)
     ==
        Ceps1*alpha()*rho()*G()/T_()
      - fvm::SuSp
        (
            (2.0/3.0*Ceps1 + CRDT_)*alpha()*rho()*divU,
            epsilon_
        )
      - fvm::Sp(Ceps2_*alpha()*rho()/T_(), epsilon_)
      + fvOptions(alpha, rho, epsilon_)
    );

    epsEqn.ref().relax();
    fvOptions.constrain(epsEqn.ref());
    epsEqn.ref().boundaryManipulate(epsilon_.boundaryFieldRef());
    solve(epsEqn);
    fvOptions.correct(epsilon_);
    bound(epsilon_, this->epsilonMin_);


    // Turbulent kinetic energy equation (LK:Eq. 6)
    // epsilon/k ~ 1/Ts (LK:p. 55)
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(alpha, rho, k_)
      + fvm::div(alphaRhoPhi, k_)
      - fvm::laplacian(alpha*rho*DkEff(), k_)
     ==
        alpha()*rho()*G()
      - fvm::SuSp((2.0/3.0)*alpha()*rho()*divU, k_)
      - fvm::Sp(alpha()*rho()/T_(), k_)
      + fvOptions(alpha, rho, k_)
    );

    kEqn.ref().relax();
    fvOptions.constrain(kEqn.ref());
    solve(kEqn);
    fvOptions.correct(k_);
    bound(k_, this->kMin_);


    // Elliptic relaxation function equation (LK:Eq. 15)
    tmp<fvScalarMatrix> fEqn
    (
      - fvm::laplacian(f_)
     ==
      - fvm::Sp(1.0/L2(), f_)
      - (v2fAlpha - C2_*G())/(L2()*k_())
      - (2.0/3.0*C2_)*divU/L2()
    );

    fEqn.ref().relax();
    fvOptions.constrain(fEqn.ref());
    solve(fEqn);
    fvOptions.correct(f_);
    bound(f_, fMin_);


    // Wall-normal fluctuating velocity scale equation (LK:Eq. 14; D:Eq. 6)
    tmp<fvScalarMatrix> v2Eqn
    (
        fvm::ddt(alpha, rho, v2_)
      + fvm::div(alphaRhoPhi, v2_)
      - fvm::laplacian(alpha*rho*Dv2Eff(), v2_)
      ==
        alpha()*rho()*min(k_()*f_(), C2_*(G() - 2.0/3.0*k_()*divU) - v2fAlpha)
      - fvm::Sp(6.0*alpha()*rho()*epsilon_()/k_(), v2_)
      + fvOptions(alpha, rho, v2_)
    );

    v2Eqn.ref().relax();
    fvOptions.constrain(v2Eqn.ref());
    solve(v2Eqn);
    fvOptions.correct(v2_);
    bound(v2_, v2Min_);

    // Update the turbulent time scale
    C = (sqrt(6.0)*Cmu_)*(sqrt(S && S)*v2_) + this->epsilonMin_;
    T_ = this->Ts(C);

    correctNut();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
