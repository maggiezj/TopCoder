/**
 * Copyright (C) 2014 TopCoder Inc., All Rights Reserved.
 */
package asteroidtracker;

interface IAsteroidTracker {
    public void initialize(double antennaPositions[],
                           double peakGain,
                           double minDistanceGain[],
                           double maxDistanceGain[]);

    public void asteroidAppearance(int asteroidIndex,
                                   double scienceScoreMultiplier,
                                   double reflectivityMultiplier,
                                   double initialImageInformation,
                                   double initialTrajectoryInformation,
                                   double trajectory[]);

    public String nextCommand(double currentTime);
}

