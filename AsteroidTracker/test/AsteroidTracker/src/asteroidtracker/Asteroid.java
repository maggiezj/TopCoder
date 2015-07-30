/**
 * Copyright (C) 2014 TopCoder Inc., All Rights Reserved.
 */
package asteroidtracker;

import java.io.Serializable;

class Asteroid implements Serializable {
    public AsteroidSighting trajectory[];
    public double scienceScoreMultiplier;   // multiplier for image/tracking knowledge score
    public double reflectivityMultiplier;   // multiplier including area, reflectivity, etc.
    public double initialImageInformation;
    public double initialTrajectoryInformation;
}

class AsteroidSighting implements Serializable {
    public double time;
    public Vector3d pos;
}

