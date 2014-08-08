trilat
======

Solution to the trilateration problem, based on the Gauss-Newton method and implemented in Octave.

Trilateration, the problem of determining a point in space, given the distances (known as radii) to the point from a number of reference points (called beacons). The most famous application is ofcourse the Global Positioning System (GPS).

The package is intended to provide functionality for determining the points, given the coordinates of beacons and the radii corresponding to the respective beacons.



SHORT INTRO TO TRILATERATION
Trilateretion can be visualized as the problem of finding the intersection of spheres (in 3D space). Generally if the space has m dimensions and there are n beacons, trilateration will only have a unique solution if n > m. In practice the sphereres don't intersect at the exact same point. In that case the package will find the least-square approximation of the intersection. If m = n, trilateration will have two solutions. In practice one of these solutions can typically be ignored. The package tackles this problem by letting the user to optionally specify an initial guess x0. If m = n and x0 is sufficiently close to the real solution x, then trilat will find the correct of the two solutions.
