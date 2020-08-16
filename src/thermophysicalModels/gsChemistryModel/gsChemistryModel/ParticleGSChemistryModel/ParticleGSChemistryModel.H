/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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
    Foam::ParticleGSChemistryModel

Description
    Extends base chemistry model by adding a thermo package, and ODE functions.
    Introduces chemistry equation system and evaluation of chemical source
    terms.

SourceFiles
    ParticleGSChemistryModelI.H
    ParticleGSChemistryModel.C

\*---------------------------------------------------------------------------*/

#ifndef ParticleGSChemistryModel_H
#define ParticleGSChemistryModel_H

#include "ReactionList.H"
#include "SurfaceReactionList.H"
#include "ODESystem.H"
#include "volFields.H"
#include "rhoReactionThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class fvMesh;

/*---------------------------------------------------------------------------*\
                      Class ParticleGSChemistryModel Declaration
\*---------------------------------------------------------------------------*/

template<class ReactionThermo, class GThermoType, class SThermoType>
class ParticleGSChemistryModel
:
    public ReactionThermo,
    public ODESystem
{
    // Private Member Functions

        //- Solve the reaction system for the given time step
        //  of given type and return the characteristic time
        template<class DeltaTType>
        scalar solve(const DeltaTType& deltaT);

        //- Initialize surface coverages
        void initializeSurface();

        //- Solve surface species QSSA in celli
        void solveSurfaceQSSA
        (
            const label celli,
            const scalar treshold,
            const label maxiter,
            const scalar deltaT
        );


protected:

    typedef SThermoType thermoType;


    // Protected data

        //- Reference to the field of specie mass fractions
        PtrList<volScalarField>& Y_;

        //- Reference to the field of solid species coverages
        PtrList<volScalarField>& Ys_;

        //- Thermodynamic data of the gas species
        const PtrList<GThermoType>& gasSpecieThermo_;

        //- Thermodynamic data of the solid species
        const PtrList<SThermoType>& solidSpecieThermo_;

        //- Gas reactions
        const ReactionList<GThermoType> gasReactions_;

        //- Solid reactions
        const SurfaceReactionList<SThermoType> solidReactions_;

        //- Number of species
        label nGasSpecie_;
        label nSolidSpecie_;
        label nSpecie_;

        //- Number of reactions
        label nGasReaction_;
        label nSolidReaction_;
        label nReaction_;

        //- Temperature below which the reaction rates are assumed 0
        scalar Treact_;

        //- Temporary concentration field
        mutable scalarField c_;

        //- Temporary rate-of-change of concentration field
        mutable scalarField dcdt_;

        //- Flag to indicate if initialisation is taking place
        bool initialize_;

        //- Catalyst surface-area-to-volume
        const scalar catarea_;

        //- Catalyst porosity
        const scalar por_;

        //- Gas density in pores [kg/m3]
        scalar rhogs_;

        //- Use instantaneous rate
        const bool instantaneousRates_;

public:

    //- Runtime type information
    TypeName("gsParticle");


    // Constructors

        //- Construct from thermo
        ParticleGSChemistryModel(rhoReactionThermo& thermo);

        //- Disallow default bitwise copy construction
        ParticleGSChemistryModel(ParticleGSChemistryModel&);


    //- Destructor
    virtual ~ParticleGSChemistryModel();


    // Member Functions

        //- The reactions
        inline const PtrList<Reaction<GThermoType>>& gasReactions() const;
        inline const PtrList<SurfaceReaction<SThermoType>>& solidReactions() const;

        //- Thermodynamic data of the species
        inline const PtrList<GThermoType>& gasSpecieThermo() const;
        inline const PtrList<SThermoType>& solidSpecieThermo() const;

        //- The number of species
        virtual inline label nSpecie() const;

        //- The number of reactions
        virtual inline label nReaction() const;

        //- Temperature below which the reaction rates are assumed 0
        inline scalar Treact() const;

        //- Temperature below which the reaction rates are assumed 0
        inline scalar& Treact();

        //- dc/dt = omega, rate of change in concentration, for each species
        virtual void omega
        (
            const scalar p,
            const scalar T,
            const scalarField& c,
            const label li,
            scalarField& dcdt
        ) const;

        //- Calculates the reaction rates
        virtual void calculate();


        // Chemistry model functions (overriding abstract functions in
        // basicChemistryModel.H)

            //- Return const access to the chemical source terms for gas specie, i
            inline const volScalarField::Internal& RRg
            (
                const label i
            ) const;

            //- Return const access to the chemical source terms for solid specie, i
            inline const volScalarField::Internal& RRs
            (
                const label i
            ) const;

            //- Return const access to the chemical source terms for specie, i
            inline const volScalarField::Internal& RR
            (
                const label i
            ) const;

            //- Return non const access to chemical source terms [kg/m^3/s]
            virtual volScalarField::Internal& RR
            (
                const label i
            );

            //- Return reaction rate of the speciei in reactionI
            virtual tmp<volScalarField::Internal> calculateRR
            (
                const label reactionI,
                const label speciei
            ) const;

            //- Solve the reaction system for the given time step
            //  and return the characteristic time
            virtual scalar solve(const scalar deltaT);

            //- Solve the reaction system for the given time step
            //  and return the characteristic time
            virtual scalar solve(const scalarField& deltaT);

            //- Return the chemical time scale
            virtual tmp<volScalarField> tc() const;

            //- Return the heat release rate [kg/m/s^3]
            virtual tmp<volScalarField> Qdot() const;


        // ODE functions (overriding abstract functions in ODE.H)

            //- Number of ODE's to solve
            inline virtual label nEqns() const;

            virtual void derivatives
            (
                const scalar t,
                const scalarField& c,
                const label li,
                scalarField& dcdt
            ) const;

            virtual void jacobian
            (
                const scalar t,
                const scalarField& c,
                const label li,
                scalarField& dcdt,
                scalarSquareMatrix& J
            ) const;

            virtual void solve
            (
                scalar& p,
                scalar& T,
                scalarField& c,
                const label li,
                scalar& deltaT,
                scalar& subDeltaT
            ) const = 0;

        // Other functions

            //- Correct catalytic reaction rates
            virtual void correctRates
            (
                PtrList<volScalarField>& Yg,
                const volScalarField& rhog,
                const volScalarField& Tg,
                const volScalarField& pg
            );

            //- Set gas phase mass fractions (e.g. from other phase)
            virtual void setGasMassFractions(const PtrList<volScalarField>& Yg);

            //- Disable surface species equations
            virtual void disableSurfaceSpeciesEqn();

            //- Correct catalytic wall fluxes
            virtual void correctCatalyticWallFluxes
            (
                PtrList<volScalarField>& Yg,
                const volScalarField& rhog,
                volScalarField& heg,
                volScalarField& Tg,
                const volScalarField& pg,
                const wordList& patches
            );

            //- Return catalytic wall species flux on patchi
            virtual scalarField catalyticWallSpeciesFlux
            (
                const word speciei,
                const label patchi
            ) const;

            //- Return catalytic wall heat flux on patchi
            virtual scalarField catalyticWallHeatFlux(const label patchi) const;

            //- Set volume fraction alpha
            virtual void setAlpha(const volScalarField& alpha);

            //- Set temperatures (e.g. from other phase)
            virtual void setTemperature(const volScalarField& Tg);

            //- Get rate in one cell or for one particle,
            //  given the conditions (Y, Ys, rho, T, p).
            //  This function returns the reaction heat
            //  Optionally, initialize surface coverages
            virtual scalar getRatesQdotI
            (
                const scalar deltaT,
                const label celli,
                scalarField& RR,
                const scalarField& Yg,
                scalarField& Ys,
                const scalar rhog,
                const scalar Tg,
                scalar Ts,
                scalar p,
                const bool initialize = false
            );


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const ParticleGSChemistryModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ParticleGSChemistryModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ParticleGSChemistryModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
