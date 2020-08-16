/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2019 OpenFOAM Foundation
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
    Foam::ReactiveGasExchangePhaseSystem

Description
    Class to provide interfacial heat and mass transfer between a number of
    phases according to a interface composition model.

SourceFiles
    ReactiveGasExchangePhaseSystem.C

\*---------------------------------------------------------------------------*/

#ifndef ReactiveGasExchangePhaseSystem_H
#define ReactiveGasExchangePhaseSystem_H

#include "phaseSystem.H"

#include "fvm.H"
#include "phaseModel.H"
#include "ThermoPhaseModel.H"
#include "MovingPhaseModel.H"
#include "StationaryPhaseModel.H"
#include "ReactingGSPhaseModel.H"
#include "rhoReactionThermo.H"
#include "basicGSChemistryModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class ReactiveGasExchangePhaseSystem Declaration
\*---------------------------------------------------------------------------*/

template<class BasePhaseSystem>
class ReactiveGasExchangePhaseSystem
:
    public BasePhaseSystem
{
protected:

    // Protected data
    const word gasPhaseName_;

    const List<word> exchangePhaseModels_;

    const bool multiphaseReactionHeat_;

    // Protected member functionstypedef
    typedef ReactingGSPhaseModel
                <
                    MovingPhaseModel
                    <
                        ThermoPhaseModel<phaseModel, rhoReactionThermo>
                    >,
                    basicGSChemistryModel
                >

        reactingGSPhaseModel;

    typedef ReactingGSPhaseModel
                <
                    StationaryPhaseModel
                    <
                        ThermoPhaseModel<phaseModel, rhoReactionThermo>
                    >,
                    basicGSChemistryModel
                >

        stationaryReactingGSPhaseModel;


public:

    // Constructors

        //- Construct from fvMesh
        ReactiveGasExchangePhaseSystem(const fvMesh&);


    //- Destructor
    virtual ~ReactiveGasExchangePhaseSystem();


    // Member Functions

        //- Return the mass transfer rate for an interface
        virtual tmp<volScalarField> dmdtf(const phasePairKey& key) const;

        //- Return the mass transfer rates for each phase
        virtual PtrList<volScalarField> dmdts() const;

        //- Return the momentum transfer matrices for the cell-based algorithm
        virtual autoPtr<phaseSystem::momentumTransferTable> momentumTransfer();

        //- Return the momentum transfer matrices for the face-based algorithm
        virtual autoPtr<phaseSystem::momentumTransferTable> momentumTransferf();

        //- Return the heat transfer matrices
        virtual autoPtr<phaseSystem::heatTransferTable> heatTransfer() const;

        //- Return the specie transfer matrices
        virtual autoPtr<phaseSystem::specieTransferTable>
            specieTransfer() const;

        //- Correct the fluid properties other than those listed below
        virtual void correct();

        //- Read base phaseProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ReactiveGasExchangePhaseSystem.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
