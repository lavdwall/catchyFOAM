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

#include "ReactingGSPhaseModel.H"
#include "phaseSystem.H"
#include "fvMatrix.H"
#include "localEulerDdtScheme.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasePhaseModel, class ReactionType>
Foam::ReactingGSPhaseModel<BasePhaseModel, ReactionType>::ReactingGSPhaseModel
(
    const phaseSystem& fluid,
    const word& phaseName,
    const bool referencePhase,
    const label index
)
:
    BasePhaseModel(fluid, phaseName, referencePhase, index),
    reaction_(ReactionType::New(this->thermo_())),
    gasExchange_(false),
    chemistryT_(false)
{
    word phaseSystemType(fluid.lookup("type"));
    if (phaseSystemType.find("reactiveGasExchange")!=std::string::npos)
    {
        Info << "Found reactiveGasExchange phase system" << endl;
        gasExchange_ = true;
        gasPhaseName_ = word(fluid.subDict("reactiveGasExchange").lookup("gasPhase"));
        forAll(this->thermo_->composition().Y(), i)
        {
            this->thermo_->composition().setInactive(i);
        }
        chemistryT_ = fluid.subDict("reactiveGasExchange").
                      lookupOrDefault("useGasTemperatureForChemistry", false);

        if (!fluid.subDict("reactiveGasExchange").lookupOrDefault("coveragesTransport", true))
        {
            forAll(this->thermo_->composition().Ys(), i)
            {
                this->thermo_->composition().setInactive(i + this->thermo_->composition().Y().size());
            }
        }
    }
    if (!gasExchange_ && reaction_->lookupOrDefault("twoPhaseAlpha",false))
    {
        FatalErrorInFunction << "Mass transfer model cannot be combined " <<
          "with coupled chemistry approach: " << nl <<
          "either disable twoPhaseAlpha in chemistryProperties."<< phaseName <<
          "or, select reactiveGasExchange phase system in phaseProperties!"
          << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasePhaseModel, class ReactionType>
Foam::ReactingGSPhaseModel<BasePhaseModel, ReactionType>::~ReactingGSPhaseModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasePhaseModel, class ReactionType>
void Foam::ReactingGSPhaseModel<BasePhaseModel, ReactionType>::correctSpecies()
{
    if (gasExchange_)
    {
        reaction_->setGasMassFractions(this->fluid().phases()[gasPhaseName_].Y());
    }

    BasePhaseModel::correctSpecies();
}


template<class BasePhaseModel, class ReactionType>
void Foam::ReactingGSPhaseModel<BasePhaseModel, ReactionType>::correctReactions()
{
    if (gasExchange_) reaction_->setAlpha(*this);
    if (gasExchange_)
    {
        if (chemistryT_)
        {
            reaction_->setTemperature(this->fluid().phases()[gasPhaseName_].thermo().T());
        }
    }

    if (fv::localEulerDdt::enabled(this->mesh()))
    {
        const scalarField& rDeltaT =
            fv::localEulerDdt::localRDeltaT(this->mesh());

        if (this->fluid().subDict(this->name()).found("maxIntegrationTime"))
        {
            scalar maxIntegrationTime
            (
                readScalar(this->fluid().subDict(this->name())
                    .lookup("maxIntegrationTime"))
            );

            reaction_->solve
            (
                min(1.0/rDeltaT, maxIntegrationTime)()
            );
        }
        else
        {
            reaction_->solve((1.0/rDeltaT)());
        }
    }
    else
    {
        reaction_->solve(this->mesh().time().deltaTValue());
    }

    BasePhaseModel::correctReactions();
}


template<class BasePhaseModel, class ReactionType>
Foam::tmp<Foam::fvScalarMatrix>
Foam::ReactingGSPhaseModel<BasePhaseModel, ReactionType>::R
(
    volScalarField& Yi
) const
{
    tmp<fvScalarMatrix> tSu(new fvScalarMatrix(Yi, dimMass/dimTime));
    fvScalarMatrix& Su = tSu.ref();

    const label specieI = this->thermo_->composition().species()[Yi.member()];


    Su += reaction_->RR(specieI);

    if (gasExchange_)
      return tSu;
    else
      return (*this)*tSu;
}


template<class BasePhaseModel, class ReactionType>
Foam::tmp<Foam::fvScalarMatrix>
Foam::ReactingGSPhaseModel<BasePhaseModel, ReactionType>::Rs
(
    volScalarField& Ysi
) const
{
    tmp<fvScalarMatrix> tSu(new fvScalarMatrix(Ysi, dimVolume/dimTime));
    fvScalarMatrix& Su = tSu.ref();

    const label specieI_S = this->thermo_->composition().solidSpecies()[Ysi.member()];
    const label nG_ = this->thermo_->composition().Y().size();
    Su += reaction_->RR(specieI_S+nG_);

    return tSu;
}


template<class BasePhaseModel, class ReactionType>
Foam::tmp<Foam::fvScalarMatrix>
Foam::ReactingGSPhaseModel<BasePhaseModel, ReactionType>::Rg
(
    const volScalarField& Yi
) const
{
    tmp<fvScalarMatrix> tSu(new fvScalarMatrix(Yi, dimMass/dimTime));
    fvScalarMatrix& Su = tSu.ref();

    const label specieI = this->thermo_->composition().species()[Yi.member()];


    Su += reaction_->RR(specieI);

    return tSu;
}


template<class BasePhaseModel, class ReactionType>
Foam::tmp<Foam::volScalarField>
Foam::ReactingGSPhaseModel<BasePhaseModel, ReactionType>::Qdot() const
{
    if (gasExchange_)
      return reaction_->Qdot();
    else
      return (*this)*reaction_->Qdot();
}


// ************************************************************************* //
