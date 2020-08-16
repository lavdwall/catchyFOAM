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

#include "MultiComponentGSPhaseModel.H"

#include "phaseSystem.H"

#include "fvmDdt.H"
#include "fvmDiv.H"
#include "fvmSup.H"
#include "fvmLaplacian.H"
#include "fvcDdt.H"
#include "fvcDiv.H"
#include "fvMatrix.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasePhaseModel>
Foam::MultiComponentGSPhaseModel<BasePhaseModel>::MultiComponentGSPhaseModel
(
    const phaseSystem& fluid,
    const word& phaseName,
    const bool referencePhase,
    const label index
)
:
    BasePhaseModel(fluid, phaseName, referencePhase, index),
    residualAlpha_
    (
        "residualAlpha",
        dimless,
        fluid.mesh().solverDict("Yi")
    ),
    inertIndex_(-1),
    inertSIndex_(-1)
{
    const word inertSpecie
    (
        this->thermo_->lookupOrDefault("inertSpecie", word::null)
    );

    if (inertSpecie != word::null)
    {
        inertIndex_ = this->thermo_->composition().species()[inertSpecie];
    }

    const word inertSolidSpecie
    (
        this->thermo_->lookupOrDefault("inertSolidSpecie", word::null)
    );

    if (inertSolidSpecie != word::null)
    {
        inertSIndex_ = this->thermo_->composition().solidSpecies()[inertSolidSpecie];
    }

    PtrList<volScalarField>& Y = this->thermo_->composition().Y();

    forAll(Y, i)
    {
        if (i != inertIndex_ && this->thermo_->composition().active(i))
        {
            const label j = YActive_.size();
            YActive_.resize(j + 1);
            YActive_.set(j, &Y[i]);
        }
    }

    PtrList<volScalarField>& Ys = this->thermo_->composition().Ys();

    forAll(Ys, i)
    {
        if (i != inertSIndex_ && this->thermo_->composition().active(i+Y.size()))
        {
            const label j = YsActive_.size();
            YsActive_.resize(j + 1);
            YsActive_.set(j, &Ys[i]);
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasePhaseModel>
Foam::MultiComponentGSPhaseModel<BasePhaseModel>::~MultiComponentGSPhaseModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasePhaseModel>
void Foam::MultiComponentGSPhaseModel<BasePhaseModel>::correctSpecies()
{
    BasePhaseModel::correctSpecies();

    volScalarField Yt
    (
        IOobject
        (
            IOobject::groupName("Yt", this->name()),
            this->fluid().mesh().time().timeName(),
            this->fluid().mesh()
        ),
        this->fluid().mesh(),
        dimensionedScalar(dimless, 0)
    );

    PtrList<volScalarField>& Yi = YRef();

    forAll(Yi, i)
    {
        if (i != inertIndex_)
        {
            Yi[i].max(0);
            Yt += Yi[i];
        }
    }

    if (inertIndex_ != -1)
    {
        Yi[inertIndex_] = scalar(1) - Yt;
        Yi[inertIndex_].max(0);
    }
    else
    {
        forAll(Yi, i)
        {
            Yi[i] /= Yt;
        }
    }

    volScalarField Yst
    (
        IOobject
        (
            IOobject::groupName("Yst", this->name()),
            this->fluid().mesh().time().timeName(),
            this->fluid().mesh()
        ),
        this->fluid().mesh(),
        dimensionedScalar("zero", dimless, 0)
    );

    PtrList<volScalarField>& Ysi = YsRef();

    forAll(Ysi, i)
    {
        if (i != inertSIndex_)
        {
            Ysi[i].max(0);
            Yst += Ysi[i];
        }
    }

    if (inertSIndex_ != -1)
    {
        Ysi[inertSIndex_] = scalar(1) - Yst;
        Ysi[inertSIndex_].max(0);
    }
    else
    {
        forAll(Ysi, i)
        {
            Ysi[i] /= Yst;
        }
    }
}


template<class BasePhaseModel>
bool Foam::MultiComponentGSPhaseModel<BasePhaseModel>::pure() const
{
    return false;
}


template<class BasePhaseModel>
Foam::tmp<Foam::fvScalarMatrix>
Foam::MultiComponentGSPhaseModel<BasePhaseModel>::YiEqn(volScalarField& Yi)
{
    const volScalarField& alpha = *this;
    const volScalarField& rhog =
        this->fluid().mesh().objectRegistry::lookupObject<volScalarField>(IOobject::groupName("thermo:rhog",this->name()));
    const surfaceScalarField alphaRhogPhi
    (
        this->alphaPhi()*fvc::interpolate(rhog)
    );

    return
    (
        fvm::ddt(alpha, rhog, Yi)
      + fvm::div(alphaRhogPhi, Yi, "div(" + alphaRhogPhi.name() + ",Yi)")
     ==
        this->R(Yi)

      + fvc::ddt(residualAlpha_*rhog, Yi)
      - fvm::ddt(residualAlpha_*rhog, Yi)
    );
}


template<class BasePhaseModel>
const Foam::PtrList<Foam::volScalarField>&
Foam::MultiComponentGSPhaseModel<BasePhaseModel>::Y() const
{
    return this->thermo_->composition().Y();
}


template<class BasePhaseModel>
const Foam::volScalarField&
Foam::MultiComponentGSPhaseModel<BasePhaseModel>::Y(const word& name) const
{
    return this->thermo_->composition().Y(name);
}


template<class BasePhaseModel>
Foam::PtrList<Foam::volScalarField>&
Foam::MultiComponentGSPhaseModel<BasePhaseModel>::YRef()
{
    return this->thermo_->composition().Y();
}


template<class BasePhaseModel>
const Foam::UPtrList<Foam::volScalarField>&
Foam::MultiComponentGSPhaseModel<BasePhaseModel>::YActive() const
{
    return YActive_;
}


template<class BasePhaseModel>
Foam::UPtrList<Foam::volScalarField>&
Foam::MultiComponentGSPhaseModel<BasePhaseModel>::YActiveRef()
{
    return YActive_;
}


template<class BasePhaseModel>
Foam::tmp<Foam::fvScalarMatrix>
Foam::MultiComponentGSPhaseModel<BasePhaseModel>::YsiEqn(volScalarField& Ysi)
{
    const volScalarField& alpha = *this;
    const surfaceScalarField alphaRhoPhi(this->alphaRhoPhi());
    const volScalarField& rho = this->thermo().rho();

    return
    (
        fvm::ddt(alpha, rho, Ysi)
      + fvm::div(alphaRhoPhi, Ysi, "div(" + alphaRhoPhi.name() + ",Ysi)")
     ==
        alpha*rho*this->Rs(Ysi)

      + fvc::ddt(residualAlpha_*rho, Ysi)
      - fvm::ddt(residualAlpha_*rho, Ysi)
    );
}


template<class BasePhaseModel>
const Foam::PtrList<Foam::volScalarField>&
Foam::MultiComponentGSPhaseModel<BasePhaseModel>::Ys() const
{
    return this->thermo_->composition().Ys();
}


template<class BasePhaseModel>
const Foam::volScalarField&
Foam::MultiComponentGSPhaseModel<BasePhaseModel>::Ys(const word& name) const
{
    return this->thermo_->composition().Ys(name);
}


template<class BasePhaseModel>
Foam::PtrList<Foam::volScalarField>&
Foam::MultiComponentGSPhaseModel<BasePhaseModel>::YsRef()
{
    return this->thermo_->composition().Ys();
}


template<class BasePhaseModel>
const Foam::UPtrList<Foam::volScalarField>&
Foam::MultiComponentGSPhaseModel<BasePhaseModel>::YsActive() const
{
    return YsActive_;
}


template<class BasePhaseModel>
Foam::UPtrList<Foam::volScalarField>&
Foam::MultiComponentGSPhaseModel<BasePhaseModel>::YsActiveRef()
{
    return YsActive_;
}


// ************************************************************************* //
