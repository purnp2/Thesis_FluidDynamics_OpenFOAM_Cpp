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

Class
    Foam::RASModels::transEARSM_kd

Group
    grpRASTurbulence

Description
   Implementation of the transEARSM_kd turbulence model,
   which is combination of Hellsten2005 EARSM model and 
   Kubacki-Dick2017 transition model.

    This turbulence model is described in:
    \verbatim
        Hellsten, Antti K. (2005). 
        "New Advanced K-W Turbulence Model for High-Lift Aerodynamics." 
        AIAA Journal 43 (9): 1857–69. doi:10.2514/1.13754.
        and
        Kubacki S. and Dick E. (2017) .
        "Transition Models for Turbomachinery Boundary Layer Flows: A Review"
        International Journal of Turbomachinery Propulsion and Power
    \endverbatim

    The default model coefficients are
    \verbatim
        transEARSM_kdCoeffs
        {
            Agamma          0.09
            Cs              0.518
            Ca              0.44
            Ckh             10
            Ck              10  //Changed from 6 to 10 in updated model.
            betaStar        0.09
            gamma_Coeff     0.553
            beta            0.075
            alphaK          1.01
            alphaOmega      1.0
            alphaD          0.52
            CT              1.8125
            Css             2.75    
            A0             -0.72
        }
    \endverbatim

Author:
    Purn Prakash <purn.prakash.stud@pw.edu.pl>
    Credits: EARSM model used published by
    Dr Jiri Furst <Jiri.Furst@gmail.com>

SourceFiles
    EARSM.C (by Dr Furst)

\*---------------------------------------------------------------------------*/

#ifndef transEARSM_kd_H
#define transEARSM_kd_H

#include "RASModel.H"
#include "nonlinearEddyViscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                           Class transEARSM_kd Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class transEARSM_kd
:
    public nonlinearEddyViscosity<RASModel<BasicTurbulenceModel> >
{
    // Private memmber functions
    
    // Disallow default bitwise copy construct and assignment
    transEARSM_kd(const transEARSM_kd&);
    transEARSM_kd& operator=(const transEARSM_kd&);
    
    
protected:
    
    // Protected data
    
    // Model coefficients

    dimensionedScalar Agamma_;
    dimensionedScalar Cs_;
    dimensionedScalar Ca_;
    dimensionedScalar Ckh_;
    dimensionedScalar Ck_;
    dimensionedScalar betaStar_;

    dimensionedScalar gamma_Coeff_;
    dimensionedScalar beta_;
    dimensionedScalar alphaK_;
    dimensionedScalar alphaOmega_;
    dimensionedScalar alphaD_;

    dimensionedScalar CT_;
    dimensionedScalar Css_;
    
    dimensionedScalar A0_;
    
    // Fields
    
    volScalarField k_;
    volScalarField omega_;

    //- Wall distance
    //  Note: different to wall distance in parent RASModel
    //  which is for near-wall cells only
    const volScalarField& y_;


    // Protected Member Functions

    //virtual volScalarField shearShelterFactor(const volScalarField& W2); //Furst
    virtual volScalarField shearShelterFactor(const volScalarField& SSSSS, const volScalarField& OOOOO); //Kubacki
    //virtual volScalarField IF_gamma(const volScalarField& W2); //Furst formula for intermettency Factor
    volScalarField IF_gamma(); // Dr Kubacki original formula for Intm. Factor
    
    virtual void correctNut();


    virtual volScalarField N
    (
        const volScalarField& A3p,
        const volScalarField& P1,
        const volScalarField& P2
    ) const;
    
   
    virtual void correctNonlinearStress(const volTensorField& gradU);

public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;


    //- Runtime type information
    TypeName("transEARSM_kd");


    // Constructors

        //- Construct from components
        transEARSM_kd
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
    virtual ~transEARSM_kd()
    {}


    // Member Functions
    
    //- Re-read model coefficients if they have changed
    virtual bool read();
    
    //- Return the effective diffusivity for k
    tmp<volScalarField> DkEff() const
        {
            return tmp<volScalarField>
                (
                    new volScalarField("DkEff", alphaK_*this->nut_ + this->nu())
                );
        }
    
    //- Return the effective diffusivity for omega
    tmp<volScalarField> DomegaEff() const
        {
            return tmp<volScalarField>
                (
                    new volScalarField("DomegaEff", alphaOmega_*this->nut_ + this->nu())
                );
        }
    
    //- Return the total fluctuation kinetic energy
    virtual tmp<volScalarField> k() const
        {
            return k_;
        }
    
    
    //- Return the turbulence specific dissipation rate
    virtual tmp<volScalarField> omega() const
        {
            return omega_;
        }
    
    //- Return the total fluctuation kinetic energy dissipation rate
    virtual tmp<volScalarField> epsilon() const
        {
            return tmp<volScalarField>
                (
                    new volScalarField
                    (
                    IOobject
                    (
                        "epsilon",
                        this->mesh_.time().timeName(),
                        this->mesh_
                    ),
                    betaStar_*k_*omega_,
                    omega_.boundaryField().types() 
                    )
                );
            
        }
    
    //- Validate the turbulence fields after construction
    //  Update turbulence viscosity and other derived fields as requires
    virtual void validate();
    
    //- Solve the turbulence equations and correct the turbulence viscosity
    virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#ifdef NoRepository
#   include "transEARSM_kd.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
