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

#include "CatalyticMassTransferPhaseSystem.H"
#include "interfaceCompositionModel.H"
#include "heatTransferModel.H"
#include "diffusiveMassTransferModel.H"
#include "fvmSup.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::CatalyticMassTransferPhaseSystem<BasePhaseSystem>::
CatalyticMassTransferPhaseSystem
(
    const fvMesh& mesh
)
:
    BasePhaseSystem(mesh)
{
    this->generatePairsAndSubModels
    (
        "diffusiveMassTransfer",
        diffusiveMassTransferModels_
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::CatalyticMassTransferPhaseSystem<BasePhaseSystem>::
~CatalyticMassTransferPhaseSystem()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::tmp<Foam::volScalarField>
Foam::CatalyticMassTransferPhaseSystem<BasePhaseSystem>::dmdtf
(
    const phasePairKey& key
) const
{
    tmp<volScalarField> tDmdtf = BasePhaseSystem::dmdtf(key);
    return tDmdtf;
}


template<class BasePhaseSystem>
Foam::PtrList<Foam::volScalarField>
Foam::CatalyticMassTransferPhaseSystem<BasePhaseSystem>::dmdts() const
{
    PtrList<volScalarField> dmdts(BasePhaseSystem::dmdts());
    return dmdts;
}


template<class BasePhaseSystem>
Foam::autoPtr<Foam::phaseSystem::momentumTransferTable>
Foam::CatalyticMassTransferPhaseSystem<BasePhaseSystem>::
momentumTransfer()
{
    autoPtr<phaseSystem::momentumTransferTable> eqnsPtr =
        BasePhaseSystem::momentumTransfer();
    return eqnsPtr;
}


template<class BasePhaseSystem>
Foam::autoPtr<Foam::phaseSystem::momentumTransferTable>
Foam::CatalyticMassTransferPhaseSystem<BasePhaseSystem>::
momentumTransferf()
{
    autoPtr<phaseSystem::momentumTransferTable> eqnsPtr =
        BasePhaseSystem::momentumTransferf();
    return eqnsPtr;
}


template<class BasePhaseSystem>
Foam::autoPtr<Foam::phaseSystem::heatTransferTable>
Foam::CatalyticMassTransferPhaseSystem<BasePhaseSystem>::
heatTransfer() const
{
    autoPtr<phaseSystem::heatTransferTable> eqnsPtr =
        BasePhaseSystem::heatTransfer();
    return eqnsPtr;
}


template<class BasePhaseSystem>
Foam::autoPtr<Foam::phaseSystem::specieTransferTable>
Foam::CatalyticMassTransferPhaseSystem<BasePhaseSystem>::
specieTransfer() const
{
    autoPtr<phaseSystem::specieTransferTable> eqnsPtr =
        BasePhaseSystem::specieTransfer();

    phaseSystem::specieTransferTable& eqns = eqnsPtr();

    forAll(this->phaseModels_, phasei)
    {
        const phaseModel& phase = this->phaseModels_[phasei];

        const PtrList<volScalarField>& Yi = phase.Y();

        forAll(Yi, i)
        {
            eqns.insert
            (
                Yi[i].name(),
                new fvScalarMatrix(Yi[i], dimMass/dimTime)
            );
        }
    }

    // Mass transfer across the interface
    forAllConstIter
    (
        diffusiveMassTransferModelTable,
        diffusiveMassTransferModels_,
        diffusiveMassTransferModelIter
    )
    {
        const phasePair& pair(this->phasePairs_[diffusiveMassTransferModelIter.key()]);

        const phaseModel& gas = pair.continuous();
        const phaseModel& cat = pair.dispersed();

        const PtrList<volScalarField>& Ygi = gas.Y();
        const PtrList<volScalarField>& Yci = cat.Y();

        const volScalarField K(diffusiveMassTransferModelIter()->K());

        forAll(Ygi, i)
        {
            const word gasname
            (
                Ygi[i].name()
            );

            const word catname
            (
                Yci[i].name()
            );

            volScalarField gKD = K*gas.thermo().alpha();

            *eqns[gasname] += gKD*Yci[i] - fvm::Sp(gKD, eqns[gasname]->psi());
            *eqns[catname] += gKD*Ygi[i] - fvm::Sp(gKD, eqns[catname]->psi());
        }
    }

    return eqnsPtr;
}


template<class BasePhaseSystem>
void Foam::CatalyticMassTransferPhaseSystem<BasePhaseSystem>::
correct()
{
    BasePhaseSystem::correct();
}


template<class BasePhaseSystem>
bool Foam::CatalyticMassTransferPhaseSystem<BasePhaseSystem>::read()
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
