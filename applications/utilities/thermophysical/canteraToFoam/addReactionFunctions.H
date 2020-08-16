enum reactionType
{
    irreversible,
    reversible
};

template<class ReactionRateType>
void addReactionType
(
    ReactionList<gasHThermoPhysics>& reactions_,
    speciesTable& speciesTable_,
    HashPtrTable<gasHThermoPhysics>& speciesThermo_,
    const reactionType rType,
    List<specieCoeffs>& lhs,
    List<specieCoeffs>& rhs,
    const ReactionRateType& rr
)
{
    switch (rType)
    {
        case irreversible:
        {
            reactions_.append
            (
                new IrreversibleReaction
                <gasHThermoPhysics, ReactionRateType>
                (
                    ReactionProxy<gasHThermoPhysics>
                    (
                        speciesTable_,
                        lhs,
                        rhs,
                        speciesThermo_
                    ),
                    rr
                )
            );
        }
        break;

        case reversible:
        {
            reactions_.append
            (
                new ReversibleReaction
                <gasHThermoPhysics, ReactionRateType>
                (
                    ReactionProxy<gasHThermoPhysics>
                    (
                        speciesTable_,
                        lhs,
                        rhs,
                        speciesThermo_
                    ),
                    rr
                )
            );
        }
        break;

        default:

            NotImplemented;
    }
}


void addGasReaction
(
    shared_ptr<Cantera::Reaction> ctreaction_,
    ReactionList<gasHThermoPhysics>& reactions_,
    speciesTable& speciesTable_,
    HashPtrTable<gasHThermoPhysics>& speciesThermo_
)
{
    Composition ct_rs = ctreaction_->reactants;
    Composition ct_ps = ctreaction_->products;
    Composition ct_os = ctreaction_->orders;
    List<specieCoeffs> lhs(ct_rs.size());
    List<specieCoeffs> rhs(ct_ps.size());
    int nl = 0, nr = 0;

    forAll(speciesTable_, i)
    {
        if(ct_rs.count(speciesTable_[i])>0)
        {
            lhs[nl].index = i;
            lhs[nl].stoichCoeff = ct_rs[speciesTable_[i]];
            lhs[nl].exponent =
                ct_os.count(speciesTable_[i])>0 ?
                ct_os[speciesTable_[i]] :
                lhs[nl].stoichCoeff;
            nl++;
        }
        if(ct_ps.count(speciesTable_[i])>0)
        {
            rhs[nr].index = i;
            rhs[nr].stoichCoeff = ct_ps[speciesTable_[i]];
            rhs[nr].exponent =
                ct_os.count(speciesTable_[i])>0 ?
                ct_os[speciesTable_[i]] :
                rhs[nr].stoichCoeff;
            nr++;
        }
    }

    switch (ctreaction_->reaction_type)
    {
        case ELEMENTARY_RXN:
        {
            Arrhenius arr = dynamic_cast<ElementaryReaction&>(*ctreaction_).rate;
            addReactionType
            (
                reactions_, speciesTable_, speciesThermo_,
                ctreaction_->reversible? reversible : irreversible,
                lhs, rhs,
                ArrheniusReactionRate
                (
                    arr.preExponentialFactor(),
                    arr.temperatureExponent(),
                    arr.activationEnergy_R()
                )
            );
            break;
        }
        case THREE_BODY_RXN:
        {
            ThreeBodyReaction tbr = dynamic_cast<ThreeBodyReaction&>(*ctreaction_);
            scalarList efficiencies(speciesTable_.size(), 1.0);
            forAll(efficiencies, e)
            {
                if (tbr.third_body.efficiencies.count(speciesTable_[e])>0)
                    efficiencies[e] = tbr.third_body.efficiencies[speciesTable_[e]];
            }
            addReactionType
            (
                reactions_, speciesTable_, speciesThermo_,
                ctreaction_->reversible? reversible : irreversible,
                lhs, rhs,
                thirdBodyArrheniusReactionRate
                (
                    tbr.rate.preExponentialFactor(),
                    tbr.rate.temperatureExponent(),
                    tbr.rate.activationEnergy_R(),
                    thirdBodyEfficiencies(speciesTable_, efficiencies)
                )
            );
            break;
        }
        case PLOG_RXN:
        {
            PlogReaction pr = dynamic_cast<PlogReaction&>(*ctreaction_);
            scalarList ps(pr.rate.rates().size());
            scalarList As(ps.size());
            scalarList betas(ps.size());
            scalarList Tas(ps.size());
            forAll (ps, i)
            {
                ps[i] = pr.rate.rates()[i].first;
                As[i] = pr.rate.rates()[i].second.preExponentialFactor();
                betas[i] = pr.rate.rates()[i].second.temperatureExponent();
                Tas[i] = pr.rate.rates()[i].second.activationEnergy_R();
            }
            addReactionType
            (
                reactions_, speciesTable_, speciesThermo_,
                ctreaction_->reversible? reversible : irreversible,
                lhs, rhs,
                pdepArrheniusReactionRate
                (
                    ps, As, betas, Tas
                )
            );
            break;
        }
        case FALLOFF_RXN:
        {
            FalloffReaction forr = dynamic_cast<FalloffReaction&>(*ctreaction_);
            Arrhenius arrlow = forr.low_rate;
            Arrhenius arrhigh = forr.high_rate;
            scalarList efficiencies(speciesTable_.size(), 1.0);
            forAll(efficiencies, e)
            {
                if (forr.third_body.efficiencies.count(speciesTable_[e])>0)
                    efficiencies[e] = forr.third_body.efficiencies[speciesTable_[e]];
            }

            shared_ptr<Falloff> fof = forr.falloff;
            string fofType = fof->type();

            if (fofType.compare("Lindemann") == 0)
            {
                addReactionType
                (
                    reactions_, speciesTable_, speciesThermo_,
                    ctreaction_->reversible? reversible : irreversible,
                    lhs, rhs,
                    FallOffReactionRate
                        <ArrheniusReactionRate, LindemannFallOffFunction>
                    (
                        ArrheniusReactionRate
                        (
                            arrlow.preExponentialFactor(),
                            arrlow.temperatureExponent(),
                            arrlow.activationEnergy_R()
                        ),
                        ArrheniusReactionRate
                        (
                            arrhigh.preExponentialFactor(),
                            arrhigh.temperatureExponent(),
                            arrhigh.activationEnergy_R()
                        ),
                        LindemannFallOffFunction(),
                        thirdBodyEfficiencies(speciesTable_, efficiencies)
                    )
                );
            }
            else if (fofType.compare("Troe") == 0)
            {
                scalarList TroeCoeffs(4);
                fof->getParameters(TroeCoeffs.begin());

                if (TroeCoeffs[3] == 0.0)
                {
                    TroeCoeffs[3] = great;
                }

                addReactionType
                (
                    reactions_, speciesTable_, speciesThermo_,
                    ctreaction_->reversible? reversible : irreversible,
                    lhs, rhs,
                    FallOffReactionRate
                        <ArrheniusReactionRate, TroeFallOffFunction>
                    (
                        ArrheniusReactionRate
                        (
                            arrlow.preExponentialFactor(),
                            arrlow.temperatureExponent(),
                            arrlow.activationEnergy_R()
                        ),
                        ArrheniusReactionRate
                        (
                            arrhigh.preExponentialFactor(),
                            arrhigh.temperatureExponent(),
                            arrhigh.activationEnergy_R()
                        ),
                        TroeFallOffFunction
                        (
                            TroeCoeffs[0],
                            TroeCoeffs[1],
                            TroeCoeffs[2],
                            TroeCoeffs[3]
                        ),
                        thirdBodyEfficiencies(speciesTable_, efficiencies)
                    )
                );
            }
            else if (fofType.compare("SRI") == 0)
            {
                scalarList SRICoeffs(5);
                fof->getParameters(SRICoeffs.begin());

                addReactionType
                (
                    reactions_, speciesTable_, speciesThermo_,
                    ctreaction_->reversible? reversible : irreversible,
                    lhs, rhs,
                    FallOffReactionRate
                        <ArrheniusReactionRate, SRIFallOffFunction>
                    (
                        ArrheniusReactionRate
                        (
                            arrlow.preExponentialFactor(),
                            arrlow.temperatureExponent(),
                            arrlow.activationEnergy_R()
                        ),
                        ArrheniusReactionRate
                        (
                            arrhigh.preExponentialFactor(),
                            arrhigh.temperatureExponent(),
                            arrhigh.activationEnergy_R()
                        ),
                        SRIFallOffFunction
                        (
                            SRICoeffs[0],
                            SRICoeffs[1],
                            SRICoeffs[2],
                            SRICoeffs[3],
                            SRICoeffs[4]
                        ),
                        thirdBodyEfficiencies(speciesTable_, efficiencies)
                    )
                );
            }
            else
            {
                FatalErrorInFunction
                    << "Fall-off function type not implemented"
                    << exit(FatalError);
            }
            break;
        }
        default:
        {
            Info //FatalErrorInFunction
                << "Reaction rate type not implemented"
                << endl;//exit(FatalError);
        }
    }
}



void addSurfaceReaction
(
    shared_ptr<Cantera::SurfPhase> surf_,
    shared_ptr<Cantera::InterfaceKinetics> surfkin_,
    shared_ptr<Cantera::Reaction> ctreaction_,
    SurfaceReactionList<constHThermoPhysics>& reactions_,
    speciesTable& gSpeciesTable_,
    speciesTable& sSpeciesTable_,
    HashPtrTable<constHThermoPhysics>& sSpeciesThermo_
)
{
    if (ctreaction_->reversible)
    {
            FatalErrorInFunction
                << "Reversible surface reactions not supported yet."
                << exit(FatalError);
    }
    Composition ct_rs = ctreaction_->reactants;
    Composition ct_ps = ctreaction_->products;
    Composition ct_os = ctreaction_->orders;
    List<specieCoeffs> lhs, rhs, glhs, grhs;
    int nl = 0, nr = 0, ngl = 0, ngr = 0;

    forAll(sSpeciesTable_, i)
    {
        if(ct_rs.count(sSpeciesTable_[i])>0)
        {
            lhs.append(specieCoeffs());
            lhs[nl].index = i;
            lhs[nl].stoichCoeff = ct_rs[sSpeciesTable_[i]];
            lhs[nl].exponent =
                ct_os.count(sSpeciesTable_[i])>0 ?
                ct_os[sSpeciesTable_[i]] :
                lhs[nl].stoichCoeff;
            nl++;
        }
        if(ct_ps.count(sSpeciesTable_[i])>0)
        {
            rhs.append(specieCoeffs());
            rhs[nr].index = i;
            rhs[nr].stoichCoeff = ct_ps[sSpeciesTable_[i]];
            rhs[nr].exponent =
                ct_os.count(sSpeciesTable_[i])>0 ?
                ct_os[sSpeciesTable_[i]] :
                rhs[nr].stoichCoeff;
            nr++;
        }
    }

    forAll(gSpeciesTable_, i)
    {
        if(ct_rs.count(gSpeciesTable_[i])>0)
        {
            glhs.append(specieCoeffs());
            glhs[ngl].index = i;
            glhs[ngl].stoichCoeff = ct_rs[gSpeciesTable_[i]];
            glhs[ngl].exponent =
                ct_os.count(gSpeciesTable_[i])>0 ?
                ct_os[gSpeciesTable_[i]] :
                glhs[ngl].stoichCoeff;
            ngl++;
        }
        if(ct_ps.count(gSpeciesTable_[i])>0)
        {
            grhs.append(specieCoeffs());
            grhs[ngr].index = i;
            grhs[ngr].stoichCoeff = ct_ps[gSpeciesTable_[i]];
            grhs[ngr].exponent =
                ct_os.count(gSpeciesTable_[i])>0 ?
                ct_os[gSpeciesTable_[i]] :
                grhs[ngr].stoichCoeff;
            ngr++;
        }
    }

    switch (ctreaction_->reaction_type)
    {
        case INTERFACE_RXN:
        {
            Arrhenius arr = dynamic_cast<InterfaceReaction&>(*ctreaction_).rate;
            if (dynamic_cast<InterfaceReaction&>(*ctreaction_).is_sticking_coefficient)
            {
                if (dynamic_cast<InterfaceReaction&>(*ctreaction_).use_motz_wise_correction)
                {
                    FatalErrorInFunction
                        << "Sticking coefficient reaction with motz-wise correction currently not supported "
                        << exit(FatalError);
                }
                // Identify the interface phase
                size_t iInterface = npos;
                size_t min_dim = 4;
                for (size_t n = 0; n < surfkin_->nPhases(); n++) {
                    if (surfkin_->thermo(n).nDim() < min_dim) {
                        iInterface = n;
                        min_dim = surfkin_->thermo(n).nDim();
                    }
                }
                InterfaceReaction r = dynamic_cast<InterfaceReaction&>(*ctreaction_);
                std::string sticking_species = r.sticking_species;
                if (sticking_species == "") {
                    // Identify the sticking species if not explicitly given
                    bool foundStick = false;
                    for (const auto& sp : r.reactants) {
                        size_t iPhase = surfkin_->speciesPhaseIndex(surfkin_->kineticsSpeciesIndex(sp.first));
                        if (iPhase != iInterface) {
                            // Non-interface species. There should be exactly one of these
                            if (foundStick) {
                                FatalErrorInFunction
                                    << "Multiple non-interface species found"
                                    << " in sticking reaction: '" + r.equation() + "'"
                                    << exit(FatalError);
                            }
                            foundStick = true;
                            sticking_species = sp.first;
                        }
                    }
                    if (!foundStick) {
                        FatalErrorInFunction
                            << "No non-interface species found"
                            << " in sticking reaction: '" + r.equation() + "'"
                            << exit(FatalError);
                    }
                }
                double surface_order = 0.0;
                double multiplier = 1.0;
                // Adjust the A-factor
                for (const auto& sp : r.reactants) {
                    size_t iPhase = surfkin_->speciesPhaseIndex(surfkin_->kineticsSpeciesIndex(sp.first));
                    const ThermoPhase& p = surfkin_->thermo(iPhase);
                    size_t k = p.speciesIndex(sp.first);
                    if (sp.first == sticking_species) {
                        multiplier *= Foam::sqrt(GasConstant/(2*Pi*p.molecularWeight(k)));
                    } else {
                        // Non-sticking species. Convert from coverages used in the
                        // sticking probability expression to the concentration units
                        // used in the mass action rate expression. For surface phases,
                        // the dependence on the site density is incorporated when the
                        // rate constant is evaluated, since we don't assume that the
                        // site density is known at this time.
                        double order = getValue(r.orders, sp.first, sp.second);
                        if (&p == std::dynamic_pointer_cast<ThermoPhase>(surf_).get()) {
                            multiplier *= Foam::pow(surf_->size(k), order);
                            surface_order += order;
                        }
                    }
                }
                multiplier /= Foam::pow(surf_->siteDensity(), surface_order);

                scalar preExp = arr.preExponentialFactor()*multiplier;
                scalar tbeta = arr.temperatureExponent() + 0.5;
                scalar EaR = arr.activationEnergy_R();
                arr = Arrhenius(preExp, tbeta, EaR);
            }
            if (!dynamic_cast<InterfaceReaction&>(*ctreaction_).coverage_deps.empty())
            {
                std::map<std::string, CoverageDependency> coverage_deps
                    = dynamic_cast<InterfaceReaction&>(*ctreaction_).coverage_deps;
                List<Tuple2<word, scalarField>> cov;
                forAll(sSpeciesTable_, i)
                {
                    if (coverage_deps.count(sSpeciesTable_[i])>0)
                    {
                        word covspecie = sSpeciesTable_[i];
                        scalarField covcoeffs(4);
                        covcoeffs[0] = coverage_deps[covspecie].a;
                        covcoeffs[1] = coverage_deps[covspecie].m;
                        covcoeffs[2] = coverage_deps[covspecie].E;
                        size_t k = surf_->speciesIndex(covspecie);
                        covcoeffs[3] = surf_->size(k);
                        cov.append(Tuple2<word, scalarField>(covspecie, covcoeffs));
                    }
                }
                Foam::ReactionProxy<constHThermoPhysics> rgas
                (
                    sSpeciesTable_, lhs, rhs, sSpeciesThermo_
                );
                reactions_.append
                (
                    new IrreversibleSurfaceReaction
                    <constHThermoPhysics, coverageDependentSurfaceArrheniusReactionRate>
                    (
                        SurfaceReactionProxy<constHThermoPhysics>
                        (
                            rgas,
                            gSpeciesTable_,
                            glhs,
                            grhs
                        ),
                        coverageDependentSurfaceArrheniusReactionRate
                        (
                            arr.preExponentialFactor(),
                            arr.temperatureExponent(),
                            arr.activationEnergy_R(),
                            sSpeciesTable_,
                            surf_->siteDensity(),
                            gSpeciesTable_.size(),
                            cov
                        )
                    )
                );
            }
            else
            {
                Foam::ReactionProxy<constHThermoPhysics> rgas
                (
                    sSpeciesTable_, lhs, rhs, sSpeciesThermo_
                );
                reactions_.append
                (
                    new IrreversibleSurfaceReaction
                    <constHThermoPhysics, catSurfaceArrheniusReactionRate>
                    (
                        SurfaceReactionProxy<constHThermoPhysics>
                        (
                            rgas,
                            gSpeciesTable_,
                            glhs,
                            grhs
                        ),
                        catSurfaceArrheniusReactionRate
                        (
                            arr.preExponentialFactor(),
                            arr.temperatureExponent(),
                            arr.activationEnergy_R()
                        )
                    )
                );
            }
            break;
        }
        default:
        {
            Info //FatalErrorInFunction
                << "Reaction rate type not implemented"
                << endl;//exit(FatalError);
        }
    }
}
