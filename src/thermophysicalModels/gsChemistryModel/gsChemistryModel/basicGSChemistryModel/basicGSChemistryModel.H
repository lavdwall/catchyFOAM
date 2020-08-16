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
    Foam::basicGSChemistryModel

Description
    Base class for chemistry models

SourceFiles
    basicGSChemistryModelI.H
    basicGSChemistryModel.C

\*---------------------------------------------------------------------------*/

#ifndef basicGSChemistryModel_H
#define basicGSChemistryModel_H

#include "IOdictionary.H"
#include "Switch.H"
#include "scalarField.H"
#include "volFields.H"
#include "rhoReactionThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class fvMesh;

/*---------------------------------------------------------------------------*\
                   class basicGSChemistryModel Declaration
\*---------------------------------------------------------------------------*/

class basicGSChemistryModel
:
    public IOdictionary
{
protected:

    // Protected data

        //- Reference to the mesh database
        const fvMesh& mesh_;

        //- Thermo
        rhoReactionThermo& thermo_;

        //- Chemistry activation switch
        Switch chemistry_;

        //- Initial chemical time step
        const scalar deltaTChemIni_;

        //- Maximum chemical time step
        const scalar deltaTChemMax_;

        //- Latest estimation of integration step
        volScalarField::Internal deltaTChem_;


    // Protected Member Functions

        //- Return non-const access to the latest estimation of integration
        //  step, e.g. for multi-chemistry model
        inline volScalarField::Internal& deltaTChem();

        //- Correct function - updates due to mesh changes
        void correct();


public:

    //- Runtime type information
    TypeName("basicGSChemistryModel");

    //- Thermo type
    typedef rhoReactionThermo reactionThermo;

    //- Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        basicGSChemistryModel,
        thermo,
        (rhoReactionThermo& thermo),
        (thermo)
    );


    // Constructors

        //- Construct from thermo
        basicGSChemistryModel(rhoReactionThermo& thermo);

        //- Construct as copy (not implemented)
        basicGSChemistryModel(const basicGSChemistryModel&);


    // Selectors
    static autoPtr<basicGSChemistryModel> New(rhoReactionThermo& thermo);


    //- Destructor
    virtual ~basicGSChemistryModel();


    // Member Functions

        //- Return const access to the mesh database
        inline const fvMesh& mesh() const;

        //- Return const access to the thermo package
        inline const rhoReactionThermo& thermo() const;

        //- Return non-const access to the thermo package
        inline rhoReactionThermo& thermo();

        //- Chemistry activation switch
        inline Switch chemistry() const;

        //- The number of species
        virtual label nSpecie() const = 0;

        //- The number of reactions
        virtual label nReaction() const = 0;

        //- Return the latest estimation of integration step
        inline const volScalarField::Internal& deltaTChem() const;


        // Functions to be derived in derived classes

            // Fields

                //- Return const access to chemical source terms for solids [kg/m^3/s]
                virtual const volScalarField::Internal& RRs
                (
                    const label i
                ) const = 0;

                //- Return const access to chemical source terms for gases [kg/m^3/s]
                virtual const volScalarField::Internal& RRg
                (
                    const label i
                ) const = 0;

                //- Return const access to chemical source terms [kg/m^3/s]
                virtual const volScalarField::Internal& RR
                (
                    const label i
                ) const = 0;

                //- Return access to chemical source terms [kg/m^3/s]
                virtual volScalarField::Internal& RR
                (
                    const label i
                ) = 0;

                //- Return reaction rate of the speciei in reactioni
                virtual tmp<volScalarField::Internal> calculateRR
                (
                    const label reactioni,
                    const label speciei
                ) const = 0;


            // Chemistry solution

                //- Calculates the reaction rates
                virtual void calculate() = 0;

                //- Solve the reaction system for the given time step
                //  and return the characteristic time
                virtual scalar solve(const scalar deltaT) = 0;

                //- Solve the reaction system for the given time step
                //  and return the characteristic time
                virtual scalar solve(const scalarField& deltaT) = 0;

                //- Return the chemical time scale
                virtual tmp<volScalarField> tc() const = 0;

                //- Return the heat release rate [kg/m/s^3]
                virtual tmp<volScalarField> Qdot() const = 0;

            // Other functions

                //- Correct catalytic reaction rates
                virtual void correctRates
                (
                    PtrList<volScalarField>& Yg,
                    const volScalarField& rhog,
                    const volScalarField& Tg,
                    const volScalarField& pg
                ) = 0;

                //- Set gas phase mass fractions (e.g. from other phase)
                virtual void setGasMassFractions(const PtrList<volScalarField>& Yg) = 0;

                //- Disable surface species equations
                virtual void disableSurfaceSpeciesEqn() = 0;

                //- Correct catalytic wall fluxes
                virtual void correctCatalyticWallFluxes
                (
                    PtrList<volScalarField>& Yg,
                    const volScalarField& rhog,
                    volScalarField& heg,
                    volScalarField& Tg,
                    const volScalarField& pg,
                    const wordList& patches
                ) = 0;

                //- Return catalytic wall species flux on patchi
                virtual scalarField catalyticWallSpeciesFlux
                (
                    const word speciei,
                    const label patchi
                ) const = 0;

                //- Return catalytic wall heat flux on patchi
                virtual scalarField catalyticWallHeatFlux(const label patchi) const = 0;

                //- Set volume fraction alpha
                virtual void setAlpha(const volScalarField& alpha) = 0;

                //- Set temperatures (e.g. from other phase)
                virtual void setTemperature(const volScalarField& Tg) = 0;

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
                ) = 0;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const basicGSChemistryModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "basicGSChemistryModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
