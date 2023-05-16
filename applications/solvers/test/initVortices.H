scalar distBetweenVortices = 100.0;
scalar Uinf = 10.0;
scalar gamma = 100.0;

scalar a=readScalar("TRef", laminarTransport);


// Get mesh bounding information.
boundBox meshBounds = mesh.bounds();
scalar xMin = meshBounds.min().x() - 10.0*distBetweenVortices;
scalar xMax = meshBounds.max().x() + 10.0*distBetweenVortices;
scalar zHalf = 0.5 * (meshBounds.min().z() + meshBounds.max().z());

Info << "xMin/xMax = " << xMin << "/" << xMax << endl;


// Figure out number of and distance between vortices.
label nVortices = round((xMax - xMin) / distBetweenVortices);
distBetweenVortices = (xMax - xMin) / nVortices;

Info << "distance between vortices = " << distBetweenVortices << endl;
Info << "number of vortices = " << nVortices << endl;


// Initialize vortex locations and timings.
List<scalar> tC(nVortices,0.0);
List<scalar> xC(nVortices,0.0);
List<scalar> zC(nVortices,zHalf);

for (int i = 0; i < nVortices; i++)
{
    xC[i] = xMin + i*distBetweenVortices;
    tC[i] = (xC[i] - xMin) / Uinf;
}
