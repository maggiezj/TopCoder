/**
 * Copyright (C) 2014 TopCoder Inc., All Rights Reserved.
 */
package asteroidtracker;

import java.awt.Point;
import java.awt.geom.*;
import java.io.Serializable;
import java.io.File;
import java.util.ArrayList;
import java.util.Scanner;
import java.io.FileNotFoundException;
import java.util.List;

class TestCase implements Serializable {
    static public final double SPEED_OF_LIGHT = 299792458.0;
    static public final double MAX_TRANSMITTING_POWER = 40000.0;
    static public final double BACKGROUND_NOISE = 1.0e-30;
    static public final double SHIELDING_MULTIPLIER = 1.0e-27;
    static public final double RELOCATE_SPEED = Math.PI / (10 * 60);
    static public final double RELOCATE_POWER = 5000.0;
    static public final double SIMULATION_TIME = 7 * 24 * 60 * 60;
    static public final double T_MIN = 0.1;
    static public final double CRITICAL_TRACKING_ANGLE = 1.0 * Math.PI / (180.0 * 60);
    static public final double ANTENNA_SAFE_RADIUS = 11.0;

    static public final double Q_LOST = (24 * 60 * 60) / Math.log(2.0);
    static public final double Q_TRAJECTORY = 6.0e3;
    static public final double Q_IMAGE = 10.0e5;
    static public final double IMAGE_SCORE_MULTIPLIER = 30.0;
    static public final double TRAJECTORY_SCORE_MULTIPLIER = 60.0 / SIMULATION_TIME;

    public Point2D.Double antennaPositions[];
    public Asteroid asteroids[];
    public double antennaMinDistance;
    public double antennaMaxDistance;
    public double peakGain;
    public double minDistanceGain[];
    public double maxDistanceGain[];

    public void generate(String path) throws FileNotFoundException {
        Scanner sc = new Scanner(new File(path));

        final int numberOfAntennas = sc.nextInt();

        this.antennaPositions = new Point2D.Double[numberOfAntennas];

        for (int antennaIndex = 0; antennaIndex < numberOfAntennas; ++antennaIndex) {
            final double x = sc.nextDouble();
            final double y = sc.nextDouble();
            this.antennaPositions[antennaIndex] = new Point2D.Double(x, y);
        }

        this.peakGain = sc.nextDouble();

        final int gainLength = sc.nextInt();
        this.minDistanceGain = new double[gainLength];
        this.maxDistanceGain = new double[gainLength];

        {
            for (int i = 0; i < gainLength; ++i) {
                this.minDistanceGain[i] = sc.nextDouble();
            }
        }

        {
            for (int i = 0; i < gainLength; ++i) {
                this.maxDistanceGain[i] = sc.nextDouble();
            }
        }

        final int numberOfAsteroids = sc.nextInt();
        this.asteroids = new Asteroid[numberOfAsteroids];

        for (int asteroidIndex = 0; asteroidIndex < numberOfAsteroids; ++asteroidIndex) {
            Asteroid asteroid = new Asteroid();
            asteroid.scienceScoreMultiplier = sc.nextDouble();
            asteroid.reflectivityMultiplier = sc.nextDouble();
            asteroid.initialImageInformation = sc.nextDouble();
            asteroid.initialTrajectoryInformation = sc.nextDouble();

            final int len = sc.nextInt();

            asteroid.trajectory = new AsteroidSighting[len];

            for (int i = 0; i < len; ++i) {
                AsteroidSighting s = new AsteroidSighting();
                s.pos = new Vector3d();

                s.time = sc.nextDouble();
                s.pos.x = sc.nextDouble();
                s.pos.y = sc.nextDouble();
                s.pos.z = sc.nextDouble();

                asteroid.trajectory[i] = s;
            }

            this.asteroids[asteroidIndex] = asteroid;
        }

        this.antennaMinDistance = Double.MAX_VALUE;
        this.antennaMaxDistance = 0.0;

        for (int i = 0; i < numberOfAntennas; ++i) {
            for (int j = 0; j < numberOfAntennas; ++j) {
                if (i != j) {
                    double dx = antennaPositions[i].x - antennaPositions[j].x;
                    double dy = antennaPositions[i].y - antennaPositions[j].y;
                    double dist = Math.sqrt(dx * dx + dy * dy);

                    this.antennaMinDistance = Math.min(this.antennaMinDistance, dist);
                    this.antennaMaxDistance = Math.max(this.antennaMaxDistance, dist);
                }
            }
        }
    }
}
