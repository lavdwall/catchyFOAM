tmp<fv::convectionScheme<scalar>> mvConvection
(
    fv::convectionScheme<scalar>::New
    (
        mesh,
        fields,
        phi,
        mesh.divScheme("div(phi,Yi_h)")
    )
);

{
    reaction->correct();
    catReaction->correctRates
    (
        Y, thermo.rho(), thermo.T(), thermo.p()
    );
    Qdot = reaction->Qdot() + catReaction->Qdot();
    volScalarField Yt(0.0*Y[0]);

    forAll(Y, i)
    {
        if (i != inertIndex && composition.active(i))
        {
            volScalarField& Yi = Y[i];

            fvScalarMatrix YiEqn
            (
                fvm::ddt(rho, Yi)
              + mvConvection->fvmDiv(phi, Yi)
              + thermophysicalTransport->divj(Yi)
             ==
                reaction->R(Yi)
              + catReaction->RR(i)
              + fvOptions(rho, Yi)
            );

            YiEqn.relax();

            fvOptions.constrain(YiEqn);

            YiEqn.solve("Yi");

            fvOptions.correct(Yi);

            Yi.max(0.0);
            Yt += Yi;
        }
    }

    if (inertIndex != -1)
    {
        Y[inertIndex] = scalar(1) - Yt;
        Y[inertIndex].max(0.0);
    }
    else
    {
        forAll(Y, i)
        {
            Y[i] /= Yt;
        }
    }
}
