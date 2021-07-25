/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------

Application
    canteraToFoam

    Example usage:
        canteraToFoam mechanism.cti gas transportProperties chem thermo -surface surface1

Description
    Converts Cantera thermodynamics and reaction data files into
    OpenFOAM format.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "OFstream.H"
#include "OStringStream.H"
#include "IStringStream.H"

#include <cantera/base/config.h>
#if CT_USE_SYSTEM_FMT
    #include <fmt/format.cc>
#endif

#include <cantera/base/Solution.h>
#include <cantera/base/utilities.h>
#include <cantera/thermo/Species.h>
#include <cantera/thermo/IdealGasPhase.h>
#include <cantera/thermo/SurfPhase.h>
#include <cantera/kinetics/GasKinetics.h>
#include <cantera/kinetics/InterfaceKinetics.h>
#include <cantera/transport/MixTransport.h>

#include "ReactionList.H"
#include "ReactionProxy.H"
#include "IrreversibleReaction.H"
#include "ReversibleReaction.H"
#include "ArrheniusReactionRate.H"
#include "pdepArrheniusReactionRate.H"
#include "thirdBodyArrheniusReactionRate.H"
#include "FallOffReactionRate.H"
#include "LindemannFallOffFunction.H"
#include "TroeFallOffFunction.H"
#include "SRIFallOffFunction.H"

#include "SurfaceReactionList.H"
#include "SurfaceReactionProxy.H"
#include "IrreversibleSurfaceReaction.H"
#include "catSurfaceArrheniusReactionRate.H"
#include "coverageDependentSurfaceArrheniusReactionRate.H"
#include "Polynomial.H"

using namespace Foam;
using namespace Cantera;

#include "thermodynamicConstants.H"
using namespace Foam::constant::thermodynamic;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "specie.H"
#include "perfectGas.H"
#include "rhoConst.H"
#include "hConstThermo.H"
#include "janafThermo.H"
#include "sensibleEnthalpy.H"
#include "thermo.H"
#include "sutherlandTransport.H"
#include "constTransport.H"

typedef
constTransport
<
    species::thermo
    <
        hConstThermo
        <
            rhoConst<specie>
        >,
        sensibleEnthalpy
    >
>
constHThermoPhysics;

typedef
sutherlandTransport
<
    species::thermo
    <
        janafThermo
        <
            perfectGas<specie>
        >,
        sensibleEnthalpy
    >
> gasHThermoPhysics;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "addReactionFunctions.H"

int main(int argc, char *argv[])
{
    unsetenv("FOAM_SIGFPE");
    // Increase the precision of the output for JANAF coefficients
    Ostream::defaultPrecision(10);

    argList::validArgs.append("Cantera cti file (input)");
    argList::validArgs.append("Cantera gas phase name (input)");
    argList::validArgs.append("OpenFOAM transport file (input)");
    argList::validArgs.append("OpenFOAM gas chemistry file (output)");
    argList::validArgs.append("OpenFOAM thermodynamics file (output)");

    argList::addBoolOption
    (
        "transport",
        "read Cantera transport data and write polynomial information"
    );

    argList::addOption
    (
        "surface",
        "word",
        "read and convert surface chemistry and thermo"
    );

    argList::addOption
    (
        "transportMode",
        "word",
        "transport model type (CK = Polynomial<4>, other = Polynomial<5>)"
    );

    argList args(argc, argv);

    bool surf_ = args.optionFound("surface");
    bool fm_mode = args.optionFound("transportMode");
    int m_mode = 10;
    if (fm_mode)
    {
        word m_mode_name;
        args.optionReadIfPresent("transportMode", m_mode_name);
        if (m_mode_name == "CK")
        {
            m_mode = 10;
        }
        else
        {
            m_mode = 0;
        }
    }

    // Gas phase

        OFstream gasReactionsFile("" + args[4] + ".gas");
        OFstream gasThermoFile("" + args[5] + ".gas");

        // Create Cantera IdealGasMix object
        Info<< "Reading gas phase mechanism" << endl;
        shared_ptr<Solution> ctgassol =
                newSolution
                (
                    args[1],
                    args[2],
                    args.optionFound("transport") ? "CK_Mix" : ""
                );
        shared_ptr<IdealGasPhase> ctgas_thermo = std::dynamic_pointer_cast<IdealGasPhase>(ctgassol->thermo());
        shared_ptr<GasKinetics> ctgas_kin = std::dynamic_pointer_cast<GasKinetics>(ctgassol->kinetics());

        // Create gas species table
        Info<< indent << "- Reading gas species names" << endl;
        List<word> gasnames(ctgas_thermo->nSpecies());
        forAll(gasnames, gi)
        {
            gasnames[gi] = ctgas_thermo->speciesName(gi);
        }
        List<word> elementnames(ctgas_thermo->nElements());
        forAll(elementnames, ei)
        {
            elementnames[ei] = ctgas_thermo->elementName(ei);
        }
        speciesTable gasSpecies(gasnames);
        speciesTable elements(elementnames);
        gasReactionsFile.writeKeyword("elements")
            << elements << token::END_STATEMENT << nl << nl;
        gasThermoFile.writeKeyword("species")
            << gasSpecies << token::END_STATEMENT << nl << nl;

        // Assign thermodynamic properties
        Info<< indent << "- Reading transport / thermophysical properties" << endl;
        fileName transportFile(args[3]);
        dictionary transportDict_;
        transportDict_.read(IFstream(transportFile)());
        HashPtrTable<gasHThermoPhysics> gasSpeciesThermo;

        word currentSpecieName;
        scalar currentMolecularWeight;
        scalar currentLowT = 0, currentHighT = 0, currentCommonT = 0;
        gasHThermoPhysics::coeffArray highCpCoeffs(scalarList(7));
        gasHThermoPhysics::coeffArray lowCpCoeffs(scalarList(7));
        forAll(gasSpecies, i)
        {
            size_t is; label tType; scalar pRef_;
            doublereal coeffs_ [15];

            currentSpecieName = gasSpecies[i];
            currentMolecularWeight = ctgas_thermo->molecularWeight(i);
	        ctgas_thermo->species(i)->thermo->reportParameters(is, tType, currentLowT, currentHighT, pRef_, &coeffs_[0]);
            currentCommonT = coeffs_[0];

            for (label coefLabel=0; coefLabel<gasHThermoPhysics::nCoeffs_; coefLabel++)
            {
                highCpCoeffs[coefLabel] = coeffs_[1+coefLabel];
                lowCpCoeffs[coefLabel] = coeffs_[8+coefLabel];
            }

            gasSpeciesThermo.insert
            (
                currentSpecieName,
                new gasHThermoPhysics
                (
                    janafThermo<perfectGas<specie>>
                    (
                        specie
                        (
                            currentSpecieName,
                            1.0,
                            currentMolecularWeight
                        ),
                        currentLowT,
                        currentHighT,
                        currentCommonT,
                        highCpCoeffs,
                        lowCpCoeffs,
                        true
                    ),
                    transportDict_.subDict(currentSpecieName)
                )
            );
        }
        //gasSpeciesThermo.write(gasThermoFile);
        OStringStream osg;
        gasSpeciesThermo.write(osg);
        dictionary gtempDict(IStringStream(osg.str())());

        wordList gsList(gtempDict.toc());

        // Add elements
        forAll(gasSpecies, gi)
        {
            dictionary elementsDict("elements");
            scalarField nei(elements.size());
            ctgas_thermo->getAtoms(gi, nei.begin());
            forAll(nei, ei)
            {
                if (nei[ei] != 0)
                {
                    elementsDict.add
                    (
                        elements[ei],
                        nei[ei]
                    );
                }
            }

            gtempDict.subDict(gasSpecies[gi]).add("elements", elementsDict);

            if (transportDict_.subDict(gasSpecies[gi]).subDict("transport").found("mu"))
                gtempDict.subDict(gasSpecies[gi]).subDict("transport")
                    .add("mu", transportDict_.subDict(gasSpecies[gi]).subDict("transport").lookup("mu"));
            if (transportDict_.subDict(gasSpecies[gi]).subDict("transport").found("Pr"))
                gtempDict.subDict(gasSpecies[gi]).subDict("transport")
                    .add("Pr", transportDict_.subDict(gasSpecies[gi]).subDict("transport").lookup("Pr"));
        }

        // Optionally read and write transport data
        if (args.optionFound("transport"))
        {            
            shared_ptr<MixTransport> ctgas_trans = std::dynamic_pointer_cast<MixTransport>(ctgassol->transport());
            scalarField coeffs(8, 0.0);

            forAll(gasSpecies, i)
            {
                ctgas_trans->getViscosityPolynomials(i, coeffs.begin());
                Polynomial<8> muCoeffs_(coeffs);
                ctgas_trans->getConductivityPolynomials(i, coeffs.begin());
                Polynomial<8> kappaCoeffs_(coeffs);

                List<Polynomial<8>> DkjCoeffs_(gasSpecies.size());
                forAll(gasSpecies, j)
                {
                    ctgas_trans->getBinDiffusivityPolynomials(i, j, coeffs.begin());
                    DkjCoeffs_[j] = Polynomial<8>(coeffs);
                }

                gtempDict.subDict(gasSpecies[i]).subDict("transport")
                        .add("muLogCoeffs<8>", muCoeffs_);
                gtempDict.subDict(gasSpecies[i]).subDict("transport")
                        .add("kappaLogCoeffs<8>", kappaCoeffs_);
                gtempDict.subDict(gasSpecies[i]).subDict("transport")
                        .add("DkjLogCoeffs<8>", DkjCoeffs_);
            }
        }

        gtempDict.write(gasThermoFile, false);

        // Create list of gas reactions
        Info<< indent << "- Reading gas phase reactions" << endl;
        ReactionList<gasHThermoPhysics> gasReactions_;
        for(size_t i=0; i<ctgas_kin->nReactions(); i++)
        {
            addGasReaction(ctgas_kin->reaction(i), gasReactions_, gasSpecies, gasSpeciesThermo);
        }
        gasReactions_.write(gasReactionsFile);

    // Solid phase
    if (surf_)
    {
        word surfName;
        args.optionReadIfPresent("surface", surfName);
        OFstream surfReactionsFile("" + args[4] + ".surf");
        OFstream surfThermoFile("" + args[5] + ".surf");

        // Create Cantera Interface object
        Info<< "Reading surface mechanism" << endl;
        shared_ptr<Solution> ctsurfsol = newSolution(args[1], surfName, "", {ctgassol});
        shared_ptr<SurfPhase> ctsurf_thermo = std::dynamic_pointer_cast<SurfPhase>(ctsurfsol->thermo());
        shared_ptr<InterfaceKinetics> ctsurf_kin = std::dynamic_pointer_cast<InterfaceKinetics>(ctsurfsol->kinetics());


        // Create surface species table
        Info<< indent << "- Reading surface species" << endl;
        List<word> surfnames(ctsurf_thermo->nSpecies());
        scalarList surfsizes(ctsurf_thermo->nSpecies());
        for(size_t i=0; i<ctsurf_thermo->nSpecies(); i++)
        {
            surfnames[i] = ctsurf_thermo->speciesName(i);
            surfsizes[i] = ctsurf_thermo->size(i);
        }
        speciesTable solidSpecies(surfnames);
        surfThermoFile.writeKeyword("species")
            << solidSpecies << token::END_STATEMENT << nl << nl;

        // Assign thermodynamic properties
        Info<< indent << "- Reading transport / thermophysical properties" << endl;
        HashPtrTable<constHThermoPhysics> solidSpeciesThermo;

        forAll(solidSpecies, i)
        {
            currentSpecieName = solidSpecies[i];

            solidSpeciesThermo.insert
            (
                currentSpecieName,
                new constHThermoPhysics
                (
                    currentSpecieName,
                    constHThermoPhysics
                    (
                        transportDict_.subDict("solidSpecies")
                    )
                )
            );
        }

        //solidSpeciesThermo.write(surfThermoFile);
        OStringStream oss;
        solidSpeciesThermo.write(oss);
        dictionary stempDict(IStringStream(oss.str())());

        wordList ssList(stempDict.toc());

        forAll(solidSpecies, si)
        {
            dictionary elementsDict("elements");
            scalarField nei(elements.size());
            ctsurf_thermo->getAtoms(si, nei.begin());
            forAll(nei, ei)
            {
                if (nei[ei] != 0)
                {
                    elementsDict.add
                    (
                        elements[ei],
                        nei[ei]
                    );
                }
            }

            stempDict.subDict(solidSpecies[si]).add("elements", elementsDict);
        }

        // Add elements
        forAll(ssList, si)
        {
            stempDict.subDict(ssList[si]).subDict("specie")
                .add("size", surfsizes[solidSpecies[ssList[si]]]);
            stempDict.subDict(ssList[si]).subDict("specie")
                .add("surfaceDensity",  ctsurf_thermo->siteDensity());

            if (transportDict_.subDict("solidSpecies").subDict("transport").found("kappa"))
                stempDict.subDict(ssList[si]).subDict("transport")
                    .add("kappa", transportDict_.subDict("solidSpecies").subDict("transport").lookup("kappa"));
            if (transportDict_.subDict("solidSpecies").subDict("thermodynamics").found("Cv"))
                stempDict.subDict(ssList[si]).subDict("thermodynamics")
                    .add("Cv", transportDict_.subDict("solidSpecies").subDict("thermodynamics").lookup("Cv"));
        }

        stempDict.write(surfThermoFile, false);

        // Create list of surface reactions
        Info<< indent << "- Reading surface reactions" << endl;
        SurfaceReactionList<constHThermoPhysics> surfReactions_;
        for(size_t i=0; i<ctsurf_kin->nReactions(); i++)
        {
            addSurfaceReaction
            (
                ctsurf_thermo,
                ctsurf_kin,
                ctsurf_kin->reaction(i),
                surfReactions_,
                gasSpecies,
                solidSpecies,
                solidSpeciesThermo
            );
        }
        surfReactions_.write(surfReactionsFile);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
