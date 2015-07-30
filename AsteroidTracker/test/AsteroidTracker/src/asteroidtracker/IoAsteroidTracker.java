/**
 * Copyright (C) 2014 TopCoder Inc., All Rights Reserved.
 */
package asteroidtracker;

import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.lang.Thread;
import java.awt.geom.*;

/*
 *  Takes a process and implements the IAsteroidTracker interface by communicating via IO
 */
class IoAsteroidTracker implements IAsteroidTracker {
    private BufferedReader reader;
    private PrintWriter writer;

    IoAsteroidTracker(final Process solution) {
        this.reader = new BufferedReader(new InputStreamReader(solution.getInputStream()));
        this.writer = new PrintWriter(solution.getOutputStream());

        /*
         *  The error output stream from the solution should be printed (in standard out)
         */
        final Thread errorStreamRedirector = new Thread()
        {
            private BufferedReader errorReader = new BufferedReader(new InputStreamReader( solution.getErrorStream() ));

            public void run() {
                while (true) {
                    String s;
                    try {
                        s = errorReader.readLine();
                    } catch (Exception e) {
                        // e.printStackTrace();
                        return;
                    }
                    if (s == null) {
                        break;
                    }
                    System.out.println(s);
                }
            }
        };

        errorStreamRedirector.start();
    }

    public void asteroidAppearance(int asteroidIndex,
                                   double scienceScoreMultiplier,
                                   double reflectivityMultiplier,
                                   double initialImageInformation,
                                   double initialTrajectoryInformation,
                                   double trajectory[]) {
        writer.print("A");
        writer.print(" " + asteroidIndex);
        writer.print(" " + scienceScoreMultiplier);
        writer.print(" " + reflectivityMultiplier);
        writer.print(" " + initialImageInformation);
        writer.print(" " + initialTrajectoryInformation);

        int numberOfSightings = trajectory.length / 4;
        writer.print(" " + numberOfSightings);
        writer.println();

        for (int trajectoryIndex = 0; trajectoryIndex < numberOfSightings; ++trajectoryIndex) {
            writer.println(trajectory[trajectoryIndex * 4 + 0] + " "
                           + trajectory[trajectoryIndex * 4 + 1] + " "
                           + trajectory[trajectoryIndex * 4 + 2] + " "
                           + trajectory[trajectoryIndex * 4 + 3]);
        }

        writer.flush();
    }

    public String nextCommand(double currentTime) {
        writer.print("C");
        writer.print(" " + currentTime);
        writer.println();
        writer.flush();

        try {
            return reader.readLine();
        }
        catch (Exception e) {
            System.err.println("ERROR: Unable to read nextCommand from your soulution. Exception: " + e.getMessage());
            System.exit(-1);
            return "";
        }
    }

    /*
     *  Provide all information and data in this testcase
     */
    public void initialize(double antennaPositions[],
                           double peakGain,
                           double minDistanceGain[],
                           double maxDistanceGain[]) {
        int numberOfAntennas = antennaPositions.length / 2;
        writer.println(numberOfAntennas);

        for (int antennaIndex = 0; antennaIndex < numberOfAntennas; ++antennaIndex) {
            writer.println(antennaPositions[antennaIndex * 2 + 0] + " " + antennaPositions[antennaIndex * 2 + 1]);
        }

        writer.println(peakGain);

        writer.println(minDistanceGain.length);
        for (double d : minDistanceGain) {
            writer.println(d);
        }

        for (double d : maxDistanceGain) {
            writer.println(d);
        }

        writer.flush();
    }
}
