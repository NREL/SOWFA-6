{
    volScalarField nu = laminarTransport.nu();

    forAll(U,ii)
    {
      //U[ii].x() = Uinf;
        U[ii].x() = 0.0;
        U[ii].y() = 0.0;
        U[ii].z() = 0.0;
    }

    for (int i = 0; i < nVortices; i++)
    {
        tC[i] += runTime.deltaT().value();
        xC[i] = xMin + (tC[i] * Uinf);
        if (xC[i] > xMax)
        {
            xC[i] = xMin + (xC[i] - xMax);
            tC[i] = (xC[i] - xMin) / Uinf;
        }


        forAll(mesh.C(),ii)
        {
            scalar xD = mesh.C()[ii].x() - xC[i];
            scalar zD = mesh.C()[ii].z() - zC[i];
            scalar r = Foam::sqrt(xD*xD + zD*zD);
            scalar theta = Foam::atan2(zD,xD);
            scalar Vtheta = (gamma/(r*Foam::constant::mathematical::twoPi))*(1.0 - Foam::exp(-(r * r) / (4.0 * nu[ii] * tC[i])));
            U[ii].x() += Foam::sin(theta) * Vtheta;
            U[ii].z() += Foam::cos(theta) * Vtheta;
        }
    }


    U.correctBoundaryConditions();
 

    phi = fvc::flux(U);

}
