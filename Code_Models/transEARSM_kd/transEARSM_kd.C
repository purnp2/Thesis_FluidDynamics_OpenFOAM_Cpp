/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
     \\/     M anipulation  |
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

#include "transEARSM_kd.H"
#include "fvOptions.H"
#include "bound.H"
#include "wallDist.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //



// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //



template<class BasicTurbulenceModel>
void transEARSM_kd<BasicTurbulenceModel>::correctNut()
{
  correctNonlinearStress(fvc::grad(this->U_));
}


template<class BasicTurbulenceModel>
volScalarField transEARSM_kd<BasicTurbulenceModel>::N
(
    const volScalarField& A3p,
    const volScalarField& P1,
    const volScalarField& P2
) const
{
    volScalarField N = A3p / 3.0;

    forAll(N, i)
    {
        if (P2[i] < 0)
        {
            N[i] += 2*pow(sqr(P1[i]) - P2[i], 1./6.)
                * cos( 1./3.*acos( P1[i]/sqrt(sqr(P1[i]) - P2[i]))) ;
        }
        else
        {
            scalar a = max(P1[i] + sqrt(P2[i]), 0.0);
            scalar b = P1[i] - sqrt(P2[i]);
            N[i] += pow(a, 1./3.) + sign(b) * pow(fabs(b), 1./3.);
        }
    };

    forAll(N.boundaryField(), patchi)
    {
        fvPatchScalarField& pN = N.boundaryFieldRef()[patchi];
        const fvPatchScalarField& pP1 = P1.boundaryField()[patchi];
        const fvPatchScalarField& pP2 = P2.boundaryField()[patchi];

        forAll(pN, i)
        {
            if (pP2[i] < 0)
            {
                pN[i] += 2*pow(sqr(pP1[i]) - pP2[i], 1./6.)
                    * cos( 1./3.*acos( pP1[i]/sqrt(sqr(pP1[i]) - pP2[i]))) ;
            }
            else
            {
                scalar a = max(pP1[i] + sqrt(pP2[i]), 0.0);
                scalar b = pP1[i] - sqrt(pP2[i]);
                pN[i] += pow(a, 1./3.) + sign(b) * pow(fabs(b), 1./3.);
            }
        };
        
    };

    return N;
}


template<class BasicTurbulenceModel>
void transEARSM_kd<BasicTurbulenceModel>::correctNonlinearStress(const volTensorField& gradU)
{
    scalar Ctau = 6.0;
    volScalarField tau(
        max
        (
            1.0 / (this->betaStar_ * this->omega_),
            Ctau * sqrt(this->nu() / (this->betaStar_ * max(this->k_, this->kMin_) * this->omega_))
        ));
    
    volSymmTensorField S("S", tau * dev(symm(gradU)));
    volTensorField     W("W", -tau * skew(gradU));
    // NOTE: Wij = 1/2(dui/dxj - duj/dxi) = - skew(grad(U))

    volScalarField IIS  = tr(S & S);
    volScalarField IIW  = tr(W & W);
    volScalarField IV   = tr(S & W & W);

    scalar Neq = 81.0 / 20.0;
    scalar CDiff = 2.2;
    volScalarField beta1eq = - 6.0/5.0 * Neq / (sqr(Neq) - 2*IIW);
    volScalarField A3p = 9.0/5.0 + 9.0/4.0 * CDiff * max(1 + beta1eq*IIS, 0.0);
    volScalarField P1 = (sqr(A3p)/27 + (9.0/20.0)*IIS - (2.0/3.0)*IIW) * A3p;
    volScalarField P2 = sqr(P1) - pow3(sqr(A3p)/9 + 0.9*IIS + (2.0/3.0)*IIW);
    
    volScalarField N = this->N(A3p, P1, P2);

    volScalarField Q = 5.0/6.0*(sqr(N) - 2*IIW)*(2*sqr(N)-IIW);

    volScalarField beta1 = -N*(2.0*sqr(N) - 7.0*IIW) / Q;
    volScalarField beta3 = -12.0 * IV / (N * Q);
    volScalarField beta4 = -2.0 * (sqr(N)  - 2.0*IIW) / Q;
    volScalarField beta6 = -6.0 * N / Q;
    volScalarField beta9 =  6.0 / Q;

    
    volScalarField Cmu = - 0.5 * (beta1 + IIW * beta6);

    this->nut_ = Cmu * this->k_ * tau;
    this->nut_.correctBoundaryConditions();

    
    this->nonlinearStress_ = this->k_ * symm
    (
        beta3 * ( (W & W) - (1.0/3.0) * IIW * I )
        + beta4 * ( (S & W) - (W & S) )
        + beta6 * ( (S & W & W) + (W & W & S) - IIW * S - (2.0/3.0) * IV * I)
        + beta9 * ( (W & S & W & W) - (W & W & S & W) )
    );

    this->nonlinearStress_.correctBoundaryConditions();

    BasicTurbulenceModel::correctNut();

}

/*  // Dr Furst's formula for f_ss i.e. shearShelterFactor
template<class BasicTurbulenceModel>
volScalarField transEARSM_kd<BasicTurbulenceModel>::shearShelterFactor
(
    const volScalarField& W2
)
{
    volScalarField shearShelterFactor
    (
        exp( -pow( Css_ * this->nu() * (sqrt(W2)/k_) ,2) )
    );
    return shearShelterFactor;

}   */

// Dr Kubacki's formula for f_ss i.e. shearShelterFactor
template<class BasicTurbulenceModel>
volScalarField transEARSM_kd<BasicTurbulenceModel>::shearShelterFactor
(
    const volScalarField& S2,
    const volScalarField& W2
)
{
    tmp<volScalarField> f_k = 1.0 - tanh(k_/(Ck_ * this->nu() * omega_));
    tmp<volScalarField> tmp1
    (
        sqrt(S2) - sqrt(W2)
	    +dimensionedScalar("small",dimensionSet(0,0,-1,0,0), SMALL)
    );
    tmp<volScalarField> chi = tanh( (- sqrt(W2) * tmp1 ) / (Ckh_ * pow( betaStar_* omega_,2.0)));
    tmp<volScalarField> tmp2
    (
        f_k * chi
		+dimensionedScalar("small", dimensionSet(0,0,0,0,0), SMALL)
    );
    tmp<volScalarField> f_ss = exp( -pow((Cs_* (1.0 + Ca_ * tmp2) * this->nu()/ (sqrt(k_) * y_ )),2));
    Info << "Upto here- KD shearShelteringFactor formula used !!!"<< endl;
    return f_ss;
}

/*  //Dr Furst's defination for Intermittency Factor
template<class BasicTurbulenceModel>
volScalarField transEARSM_kd<BasicTurbulenceModel>::IF_gamma
(
    const volScalarField& W2
)
{
    volScalarField IF_gamma
    (
        min
		(
			max
			(
				(1/this->nu()) * (k_/sqrt(W2)) - CT_,
				0.0
			),
			1.0
		)

    );
    return IF_gamma;

}   */

//Dr Furst's defination for Intermittency Factor
template<class BasicTurbulenceModel>
volScalarField transEARSM_kd<BasicTurbulenceModel>::IF_gamma()
{
    volScalarField IF_gamma
    (
        min
		(
			max
			(
				(sqrt(k_) * y_)/(Agamma_ * this->nu()) - 1.0,
				0.0
			),
			1.0
		)

    );
    Info << "Upto here- KD intermittency Factor !!!"<< endl;
    return IF_gamma;

}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
transEARSM_kd<BasicTurbulenceModel>::transEARSM_kd
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
    nonlinearEddyViscosity<RASModel<BasicTurbulenceModel> >
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
    
    Agamma_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Agamma",
            this->coeffDict_,
            12.0
        )
    ),

    Cs_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cs",
            this->coeffDict_,
            21.0
        )
    ),

    Ca_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Ca",
            this->coeffDict_,
            1.0
        )
    ),

    Ckh_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Ckh",
            this->coeffDict_,
            10.0
        )
    ),

    Ck_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Ck",
            this->coeffDict_,
            10
        )
    ),

    betaStar_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "betaStar",
            this->coeffDict_,
            0.09
        )
    ),

    gamma_Coeff_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "gamma_Coeff",
            this->coeffDict_,
            0.553
        )
    ),

    beta_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "beta",
            this->coeffDict_,
            0.075
        )
    ),

    alphaK_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaK",
            this->coeffDict_,
            1.01
        )
    ),

    alphaOmega_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaOmega",
            this->coeffDict_,
            1.0
        )
    ),

    alphaD_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaD",
            this->coeffDict_,
            0.52
        )
    ),

    CT_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CT",
            this->coeffDict_,
            1.8125
        )
    ),

    Css_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Css",
            this->coeffDict_,
            2.75
        )
    ),

    A0_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "A0",
            this->coeffDict_,
            -0.72
        )
    ),

    k_
    (
        IOobject
        (
            IOobject::groupName("k", U.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),

    omega_
    (
        IOobject
        (
            IOobject::groupName("omega", U.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),

    y_(wallDist::New(this->mesh_).y())

{
    bound(k_, this->kMin_);
    bound(omega_, this->omegaMin_);

    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool transEARSM_kd<BasicTurbulenceModel>::read()
{
    if (nonlinearEddyViscosity<RASModel<BasicTurbulenceModel> >::read())
    {    
        Agamma_.readIfPresent(this->coeffDict());
        Cs_.readIfPresent(this->coeffDict());
        Ca_.readIfPresent(this->coeffDict());
        Ckh_.readIfPresent(this->coeffDict());
        Ck_.readIfPresent(this->coeffDict());
        betaStar_.readIfPresent(this->coeffDict());
        gamma_Coeff_.readIfPresent(this->coeffDict());
        beta_.readIfPresent(this->coeffDict());
        alphaK_.readIfPresent(this->coeffDict());
        alphaOmega_.readIfPresent(this->coeffDict());
        alphaD_.readIfPresent(this->coeffDict());;
        CT_.readIfPresent(this->coeffDict());
        Css_.readIfPresent(this->coeffDict());
        A0_.readIfPresent(this->coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


template<class BasicTurbulenceModel>
void transEARSM_kd<BasicTurbulenceModel>::validate()
{
    this->correctNut();
}


template<class BasicTurbulenceModel>
void transEARSM_kd<BasicTurbulenceModel>::correct()
{

    if (!this->turbulence_)
    {
        return;
    }

    nonlinearEddyViscosity<RASModel<BasicTurbulenceModel> >::correct();

    // Local references
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    const volVectorField& U = this->U_;
    volScalarField& nut = this->nut_;
    fv::options& fvOptions(fv::options::New(this->mesh_));
    
    volScalarField divU(fvc::div(fvc::absolute(this->phi(), U)));

    tmp<volTensorField> tgradU = fvc::grad(U);

    volScalarField G
    (
        this->GName(),
        (nut * dev(twoSymm(tgradU())) - this->nonlinearStress_) && tgradU()
    );

    omega_.boundaryFieldRef().updateCoeffs();
    
    volScalarField gradKgradOmegaByOmega
    (
        (fvc::grad(k_) & fvc::grad(omega_)) / omega_
    );

    volScalarField S2(2*magSqr(symm(tgradU()))); //read as S-square
    volScalarField W2(2*magSqr(-skew(tgradU()))); //read as W-square
    //volScalarField f_ss( this->shearShelterFactor(W2) ); //Dr Furst formula
    volScalarField f_ss( this->shearShelterFactor( S2 , W2 )); //Dr Kubacki formula
    //volScalarField IF_gamma( this->IF_gamma(W2) ); //Defination by Dr Furst
    volScalarField IF_gamma(this->IF_gamma()); // Defination by Dr Kubacki


    {

        tmp<volScalarField> CDOmega = this->alphaD_ * alpha * rho *
            max( gradKgradOmegaByOmega, dimensionedScalar("zero", inv(sqr(dimTime)), 0.0));

        tmp<fvScalarMatrix> omegaEqn
        (
            fvm::ddt(alpha, rho, omega_)
          + fvm::div(alphaRhoPhi, omega_)
          - fvm::laplacian(alpha*rho*DomegaEff(), omega_)
         ==
            alpha*rho* this->gamma_Coeff_ * omega_ / max(k_, this->kMin_) * (G*f_ss)
            - fvm::SuSp((2.0/3.0)*alpha*rho*this->gamma_Coeff_ *divU, omega_)
            - fvm::Sp(alpha*rho*this->beta_*omega_, omega_)
            + CDOmega
            + fvOptions(alpha, rho, omega_)            
        );

        omegaEqn.ref().relax();
        fvOptions.constrain(omegaEqn.ref());
        omegaEqn.ref().boundaryManipulate(omega_.boundaryFieldRef());
        solve(omegaEqn);
        fvOptions.correct(omega_);
        bound(omega_, this->omegaMin_);
    }
    
    
    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(alpha, rho, k_)
      + fvm::div(alphaRhoPhi, k_)
      - fvm::laplacian(alpha*rho*DkEff(), k_)
     ==
        alpha*rho*(G*f_ss)*IF_gamma
      - fvm::SuSp((2.0/3.0)*alpha*rho*divU, k_)
        - fvm::Sp(betaStar_*alpha*rho*omega_, k_)
      + fvOptions(alpha, rho, k_)
    );

    kEqn.ref().relax();
    fvOptions.constrain(kEqn.ref());
    solve(kEqn);
    fvOptions.correct(k_);
    bound(k_, this->kMin_);

    correctNonlinearStress(tgradU());
    
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //