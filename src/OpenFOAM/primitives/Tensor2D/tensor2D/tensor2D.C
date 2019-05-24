/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2010, 2019 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "tensor2D.H"
#include "quadraticEqn.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<>
const char* const Foam::tensor2D::vsType::typeName = "tensor2D";

template<>
const char* const Foam::tensor2D::vsType::componentNames[] =
{
    "xx", "xy",
    "yx", "yy"
};

template<>
const Foam::tensor2D Foam::tensor2D::vsType::vsType::zero
(
    tensor2D::uniform(0)
);

template<>
const Foam::tensor2D Foam::tensor2D::vsType::one
(
    tensor2D::uniform(1)
);

template<>
const Foam::tensor2D Foam::tensor2D::vsType::max
(
    tensor2D::uniform(VGREAT)
);

template<>
const Foam::tensor2D Foam::tensor2D::vsType::min
(
    tensor2D::uniform(-VGREAT)
);

template<>
const Foam::tensor2D Foam::tensor2D::vsType::rootMax
(
    tensor2D::uniform(ROOTVGREAT)
);

template<>
const Foam::tensor2D Foam::tensor2D::vsType::rootMin
(
    tensor2D::uniform(-ROOTVGREAT)
);

template<>
const Foam::tensor2D Foam::tensor2D::I
(
    1, 0,
    0, 1
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::Vector2D<Foam::complex> Foam::eigenValues
(
    const tensor2D& T
)
{
    const scalar a = T.xx();
    const scalar b = T.xy();
    const scalar c = T.yx();
    const scalar d = T.yy();

    const scalar trace = a + d;

    // (Jeannerod et al. (2013), p. 2246)
    scalar w = b*c;
    scalar e = std::fma(-b, c, w);
    scalar f = std::fma(a, d, -w);
    const scalar determinant = f + e;

    // Square-distance between two eigenvalues
    const scalar gapSqr = std::fma(-4.0, determinant, sqr(trace));

    // (Ford, (2014), Section 8.4.2.)
    // Eigenvalues are effectively real
    if (0 <= gapSqr)
    {
        scalar firstRoot = 0.5*(trace - sign(-trace)*Foam::sqrt(gapSqr));

        #ifdef FULLDEBUG
            if (mag(firstRoot) < SMALL)
            FatalErrorInFunction
                << "Almost-zero-valued root is detected."
                << abort(FatalError);
        #endif

        Vector2D<complex> lambdas
        (
            complex(firstRoot, 0),
            complex(determinant/firstRoot, 0)
        );

        // Sort the eigenvalues into ascending order
        if (mag(lambdas.y()) < mag(lambdas.x()))
        {
            Swap(lambdas.x(), lambdas.y());
        }

        return lambdas;
    }
    // Eigenvalues are complex
    else
    {
        complex lambda(0.5*trace, 0.5*Foam::sqrt(mag(gapSqr)));

        return Vector2D<complex>
        (
            lambda,
            lambda.conjugate()
        );
    }
}


Foam::Vector2D<Foam::complex> Foam::eigenVector
(
    const tensor2D& T,
    const complex& lambda,
    const Vector2D<complex>& direction
)
{
    // (Knill, (2004))
    if (VSMALL < mag(T.yx()))
    {
        Vector2D<complex> eVec(lambda - complex(T.yy()), complex(T.yx()));
        return eVec/mag(eVec);
    }
    else if (VSMALL < mag(T.xy()))
    {
        Vector2D<complex> eVec(complex(T.xy()), lambda - complex(T.xx()));
        return eVec/mag(eVec);
    }
    else
    {
        return direction;
    }
}


Foam::Tensor2D<Foam::complex> Foam::eigenVectors
(
    const tensor2D& T,
    const Vector2D<complex>& lambdas
)
{
    Vector2D<complex> Ux(pTraits<complex>::one, Zero);
    Vector2D<complex> Uy(Zero, pTraits<complex>::one);

    return Tensor2D<complex>
    (
        eigenVector(T, lambdas.x(), Ux),
        eigenVector(T, lambdas.y(), Uy)
    );
}


Foam::Tensor2D<Foam::complex> Foam::eigenVectors
(
    const tensor2D& T
)
{
    const Vector2D<complex> lambdas(eigenValues(T));

    return eigenVectors(T, lambdas);
}


// ************************************************************************* //
