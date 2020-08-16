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

\*---------------------------------------------------------------------------*/

#include "ReactiveGasExchangePhaseSystem.H"
#include "interfaceCompositionModel.H"
#include "heatTransferModel.H"
#include "diffusiveMassTransferModel.H"
#include "fvmSup.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::ReactiveGasExchangePhaseSystem<BasePhaseSystem>::
ReactiveGasExchangePhaseSystem
(
    const fvMesh& mesh
)
:
    BasePhaseSystem(mesh),
    gasPhaseName_
    (
        this->subDict("reactiveGasExchange").lookup("gasPhase")
    ),
    exchangePhaseModels_(this->subDict("reactiveGasExchange").lookup("exchangePhases")),
    multiphaseReactionHeat_(this->subDict("reactiveGasExchange").lookupOrDefault("multiphaseReactionHeat", false))
{
    Info << nl<< "Multiphase reaction heat is " << (multiphaseReactionHeat_? "on" : "off") << nl << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::ReactiveGasExchangePhaseSystem<BasePhaseSystem>::
~ReactiveGasExchangePhaseSystem()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::tmp<Foam::volScalarField>
Foam::ReactiveGasExchangePhaseSystem<BasePhaseSystem>::dmdtf
(
    const phasePairKey& key
) const
{
    tmp<volScalarField> tDmdtf = BasePhaseSystem::dmdtf(key);
    return tDmdtf;
}


template<class BasePhaseSystem>
Foam::PtrList<Foam::volScalarField>
Foam::ReactiveGasExchangePhaseSystem<BasePhaseSystem>::dmdts() const
{
    PtrList<volScalarField> dmdts(BasePhaseSystem::dmdts());
    return dmdts;
}


template<class BasePhaseSystem>
Foam::autoPtr<Foam::phaseSystem::momentumTransferTable>
Foam::ReactiveGasExchangePhaseSystem<BasePhaseSystem>::
momentumTransfer()
{
    autoPtr<phaseSystem::momentumTransferTable> eqnsPtr =
        BasePhaseSystem::momentumTransfer();
    return eqnsPtr;
}


template<class BasePhaseSystem>
Foam::autoPtr<Foam::phaseSystem::momentumTransferTable>
Foam::ReactiveGasExchangePhaseSystem<BasePhaseSystem>::
momentumTransferf()
{
    autoPtr<phaseSystem::momentumTransferTable> eqnsPtr =
        BasePhaseSystem::momentumTransferf();
    return eqnsPtr;
}


template<class BasePhaseSystem>
Foam::autoPtr<Foam::phaseSystem::heatTransferTable>
Foam::ReactiveGasExchangePhaseSystem<BasePhaseSystem>::
heatTransfer() const
{
    autoPtr<phaseSystem::heatTransferTable> eqnsPtr =
        BasePhaseSystem::heatTransfer();

    if (multiphaseReactionHeat_)
    {
        phaseSystem::heatTransferTable& eqns = eqnsPtr();

        forAll(exchangePhaseModels_, phasei)
        {
            if (!this->phaseModels_[exchangePhaseModels_[phasei]].stationary())
            {
                const reactingGSPhaseModel& phase = dynamic_cast<const reactingGSPhaseModel&>(this->phaseModels_[exchangePhaseModels_[phasei]]);

                Info<< " min/max reactive heat exchange "<< exchangePhaseModels_[phasei] << ": "
                    << min(phase.Qdot()).value() << " / " << max(phase.Qdot()).value() << endl;

                *eqns[exchangePhaseModels_[phasei]] -= phase.Qdot();
                *eqns[gasPhaseName_] += phase.Qdot();
            }
            else
            {
                const stationaryReactingGSPhaseModel& phase = dynamic_cast<const stationaryReactingGSPhaseModel&>(this->phaseModels_[exchangePhaseModels_[phasei]]);

                Info<< " min/max reactive heat exchange "<< exchangePhaseModels_[phasei] << ": "
                    << min(phase.Qdot()).value() << " / " << max(phase.Qdot()).value() << endl;

                *eqns[exchangePhaseModels_[phasei]] -= phase.Qdot();
                *eqns[gasPhaseName_] += phase.Qdot();
            }
        }
    }
    return eqnsPtr;
}


template<class BasePhaseSystem>
Foam::autoPtr<Foam::phaseSystem::specieTransferTable>
Foam::ReactiveGasExchangePhaseSystem<BasePhaseSystem>::
specieTransfer() const
{
    autoPtr<phaseSystem::specieTransferTable> eqnsPtr =
        BasePhaseSystem::specieTransfer();

    phaseSystem::specieTransferTable& eqns = eqnsPtr();

    const phaseModel& gasPhase = this->phaseModels_[gasPhaseName_];
    forAll(exchangePhaseModels_, phasei)
    {
        const phaseModel& phase = this->phaseModels_[exchangePhaseModels_[phasei]];
        const PtrList<volScalarField>& Yg = gasPhase.Y();

        forAll(Yg, i)
        {
            const word name = Yg[i].name();
            const word otherName = (IOobject::groupName(IOobject::member(name),phase.name()));

            *eqns[name] += phase.Rg(Yg[i]);
            *eqns[otherName] -= phase.Rg(phase.Y(IOobject::member(name)));
        }
    }

    return eqnsPtr;
}


template<class BasePhaseSystem>
void Foam::ReactiveGasExchangePhaseSystem<BasePhaseSystem>::
correct()
{
    BasePhaseSystem::correct();
}


template<class BasePhaseSystem>
bool Foam::ReactiveGasExchangePhaseSystem<BasePhaseSystem>::read()
{
    if (BasePhaseSystem::read())
    {
        bool readOK = true;

        // Models ...

        return readOK;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
