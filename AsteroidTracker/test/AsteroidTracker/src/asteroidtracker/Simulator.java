/**
 * Copyright (C) 2014 TopCoder Inc., All Rights Reserved.
 */
package asteroidtracker;

import java.util.ArrayList;
import java.util.List;

import java.util.ArrayDeque;
import java.util.Arrays;
import java.awt.geom.*;

class Simulator {
    public Simulator(TestCase testCase) {
        this.testCase = testCase;
        this.numberOfAntennas = testCase.antennaPositions.length;
        this.numberOfAsteroids = testCase.asteroids.length;

        this.antennaAssignments = new int[numberOfAntennas];
        this.isRelocating = new boolean[numberOfAntennas];
        this.transmittingPowers = new double[numberOfAntennas];
        this.antennaDirections = new Vector3d[numberOfAntennas];

        this.imageInformations = new double[numberOfAsteroids];
        this.trajectoryInformations = new double[numberOfAsteroids];
        this.trajectoryScores = new double[numberOfAsteroids];
        this.receivedSignalStrength = new double[numberOfAsteroids];
        this.outputPowerQueue = new ArrayDeque[numberOfAsteroids];
        this.trajectoryIndices = new int[numberOfAsteroids];

        this.energySpent = 0.0;
        this.time = 0.0;

        for (int antennaIndex = 0; antennaIndex < numberOfAntennas; ++antennaIndex) {
            antennaAssignments[antennaIndex] = -1;
            isRelocating[antennaIndex] = false;
            transmittingPowers[antennaIndex] = 0.0;
            antennaDirections[antennaIndex] = new Vector3d(0.0, 0.0, 1.0); // straight up
        }

        for (int asteroidIndex = 0; asteroidIndex < numberOfAsteroids; ++asteroidIndex) {
            imageInformations[asteroidIndex] = 0.0;
            trajectoryInformations[asteroidIndex] = 0.0;
            trajectoryScores[0] = 0.0;
            receivedSignalStrength[asteroidIndex] = 0.0;
            outputPowerQueue[asteroidIndex] = new ArrayDeque<ReceiveEvent>();
            trajectoryIndices[asteroidIndex] = 0;
        }
    }

    /*
     * For each asteroid, generate the list of antennas tracking it
     */
    List<Integer>[] getSubarrays() {
        List<Integer> subarrays[] = new List[numberOfAsteroids];

        for (int asteroidIndex = 0; asteroidIndex < numberOfAsteroids; ++asteroidIndex) {
            subarrays[asteroidIndex] = new ArrayList<Integer>();
        }

        for (int antennaIndex = 0; antennaIndex < numberOfAntennas; ++antennaIndex) {
            final int asteroidIndex = antennaAssignments[antennaIndex];
            if (asteroidIndex != -1
                    && !isRelocating[antennaIndex]
                    && getAsteroidPosition(asteroidIndex).z >= 0.0) {
                subarrays[asteroidIndex].add(antennaIndex);
            }
        }

        return subarrays;
    }

    /*
     * Get current asteroid position, given an asteroid index
     */
    Vector3d getAsteroidPosition(int asteroidIndex) {
        int trajectoryIndex = trajectoryIndices[asteroidIndex] - 1;
        return trajectoryIndex == -1 ? null : testCase.asteroids[asteroidIndex].trajectory[trajectoryIndex].pos;
    }

    boolean isAntennaTransmitting(int antennaIndex) {
        final int asteroidIndex = antennaAssignments[antennaIndex];
        return asteroidIndex != -1
               && transmittingPowers[antennaIndex] > 0.0
               && !isRelocating[antennaIndex]
               && getAsteroidPosition(asteroidIndex).z >= 0.0;
    }

    /*
     * Get current image score (before multipliers) for an asteroid, given an asteroid index
     */
    double getImageScore(int asteroidIndex) {
        return Math.tanh(imageInformations[asteroidIndex] / testCase.Q_IMAGE);
    }

    /*
     * Return current total score
     */
    double getScore() {
        double score = 0.0;

        for (int asteroidIndex = 0; asteroidIndex < numberOfAsteroids; ++asteroidIndex) {
            score += testCase.asteroids[asteroidIndex].scienceScoreMultiplier *
                     (testCase.IMAGE_SCORE_MULTIPLIER * getImageScore(asteroidIndex) + testCase.TRAJECTORY_SCORE_MULTIPLIER * trajectoryScores[asteroidIndex]);
        }

        score -= energySpent * 1.0e-9;

        return score;
    }

    /*
     *  Interpolate the gain with bilinear interpolation
     */
    public double interpolateGain(double distanceBetweenAntennas, double angle) {
        // Number of elements in the amplitude arrays
        final int len = testCase.minDistanceGain.length;

        // Real_index = index + u, where 0 <= u <= 1
        double u = (len - 1.0) * angle / Math.PI;
        int index = (int)Math.floor(u);
        u -= index;

        // if the angle happens to be exactly 180 degrees:
        if (index == len - 1) {
            index -= 1;
            u += 1.0;
        }

        // 0 <= v <= 1
        final double v = (distanceBetweenAntennas - testCase.antennaMinDistance) / (testCase.antennaMaxDistance - testCase.antennaMinDistance);

        final double p0 = (1.0 - v) * testCase.minDistanceGain[index + 0] + v * testCase.maxDistanceGain[index + 0];
        final double p1 = (1.0 - v) * testCase.minDistanceGain[index + 1] + v * testCase.maxDistanceGain[index + 1];

        return (1.0 - u) * p0 + u * p1;
    }

    /*
     * The induced noise for a given antenna index
     */
    double inducedNoise(int receivingAntennaIndex) {
        double inducedNoiseSum = 0.0;

        final Point2D.Double receivingAntennaPosition = testCase.antennaPositions[receivingAntennaIndex];

        for (int transmittingAntennaIndex = 0; transmittingAntennaIndex < numberOfAntennas; ++transmittingAntennaIndex) {
            if (receivingAntennaIndex == transmittingAntennaIndex
                    || !isAntennaTransmitting(transmittingAntennaIndex)) {
                continue;
            }

            final Point2D.Double transmittingAntennaPosition = testCase.antennaPositions[transmittingAntennaIndex];
            final Vector3d antennaDirection = antennaDirections[transmittingAntennaIndex];

            /*
             * Relative position and distance between antennas
             */
            final double antennaDX = receivingAntennaPosition.x - transmittingAntennaPosition.x;
            final double antennaDY = receivingAntennaPosition.y - transmittingAntennaPosition.y;
            final double distanceBetweenAntennas = Math.sqrt(antennaDX * antennaDX + antennaDY * antennaDY);

            /*
             * Calculate angle between receivingAntenna and antennaDirection, from transmittingAntenna's point of view
             */
            final double scalarProduct = antennaDX * antennaDirection.x + antennaDY * antennaDirection.y;
            final double angle = Math.acos(scalarProduct / distanceBetweenAntennas);

            /*
             * Calculate induced noise for this antenna pair, and update the sum
             */
            final double transmittingPower = transmittingPowers[transmittingAntennaIndex];
            final double inducedNoise = transmittingPower * interpolateGain(distanceBetweenAntennas, angle) / (distanceBetweenAntennas * distanceBetweenAntennas);
            inducedNoiseSum += inducedNoise;
        }

        return inducedNoiseSum;
    }

    /*
     * Get the current signal return for an asteroid, given the number of antennas tracking it
     */
    double effectiveSignalPowerReturn(int asteroidIndex, int N_receive) {
        /*
         * Calculate distance (R) to asteroid and R^4 (R4)
         */
        final double R = getAsteroidPosition(asteroidIndex).length();
        final double R2 = R * R;
        final double R4 = R2 * R2;

        /*
         * Trajectory accuracy penalty
         */
        final double T_multiplier = Math.max(Math.tanh(trajectoryInformations[asteroidIndex] / testCase.Q_TRAJECTORY), TestCase.T_MIN);

        return testCase.peakGain
               * testCase.asteroids[asteroidIndex].reflectivityMultiplier
               * T_multiplier
               * receivedSignalStrength[asteroidIndex]
               * (N_receive * N_receive)
               / R4;
    }

    /*
     *  Check if the beam from one antenna intersects any other antenna
     */
    public void checkAntennaBeamIntersection() throws ConstraintException {
        Vector3d receiverPosition = new Vector3d();
        Vector3d closestBeamPosition = new Vector3d();
        Vector3d relativeBeamPosition = new Vector3d();

        for (int transmittingAntennaIndex = 0; transmittingAntennaIndex < numberOfAntennas; ++transmittingAntennaIndex) {
            if (!isAntennaTransmitting(transmittingAntennaIndex)) {
                continue;
            }

            final Point2D.Double transmitterPosition = testCase.antennaPositions[transmittingAntennaIndex];
            final Vector3d beamDirection = antennaDirections[transmittingAntennaIndex];

            for (int receivingAntennaIndex = 0; receivingAntennaIndex < numberOfAntennas; ++receivingAntennaIndex) {
                if (receivingAntennaIndex == transmittingAntennaIndex) {
                    continue;
                }

                receiverPosition.x = testCase.antennaPositions[receivingAntennaIndex].x - transmitterPosition.x;
                receiverPosition.y = testCase.antennaPositions[receivingAntennaIndex].y - transmitterPosition.y;
                receiverPosition.z = 0.0;

                double closestBeamDistance = beamDirection.dot(receiverPosition);

                // We can only have an intersection if the other antenna is in front of the beam
                if (closestBeamDistance > 0.0) {
                    closestBeamPosition.scale(closestBeamDistance, beamDirection);
                    relativeBeamPosition.sub(receiverPosition, closestBeamPosition);

                    double distanceFromBeam = relativeBeamPosition.length();
                    if (distanceFromBeam < this.testCase.ANTENNA_SAFE_RADIUS) {
                        throw new ConstraintException(transmittingAntennaIndex,
                                                      receivingAntennaIndex,
                                                      "Your solution is violating near-field constraints. Antenna " + transmittingAntennaIndex
                                                      + " is damaging antenna " + receivingAntennaIndex + " with its beam.");
                    }
                }
            }
        }
    }

    /*
     *  Check target proximity
     */
    public void checkTargetProximity() throws ConstraintException {
        List<Integer> subarrays[] = getSubarrays();

        for (int asteroidIndexA = 0; asteroidIndexA < numberOfAsteroids; ++asteroidIndexA) {
            if (subarrays[asteroidIndexA].isEmpty()) {
                continue;
            }

            final Vector3d posA = getAsteroidPosition(asteroidIndexA);

            for (int asteroidIndexB = 0; asteroidIndexB < numberOfAsteroids; ++asteroidIndexB) {
                if (subarrays[asteroidIndexB].isEmpty() || asteroidIndexA == asteroidIndexB) {
                    continue;
                }

                final Vector3d posB = getAsteroidPosition(asteroidIndexB);
                final double angle = posA.angle(posB);
                final double minAngle = TestCase.CRITICAL_TRACKING_ANGLE / Math.min(subarrays[asteroidIndexA].size(), subarrays[asteroidIndexB].size());

                if (angle < minAngle) {
                    throw new ConstraintException("Tracked asteroids " + asteroidIndexA + " and " + asteroidIndexB + " are too close. Angle is only " + angle + ", minimum allowed for this pair is " + minAngle);
                }
            }
        }
    }

    /*
     * Update tracking antennas directions
     */
    void updateTrackingAntennasDirections() {
        for (int antennaIndex = 0; antennaIndex < numberOfAntennas; ++antennaIndex) {
            final int asteroidIndex = antennaAssignments[antennaIndex];

            // If the antenna is not tracking anything, ignore
            if (asteroidIndex == -1 || isRelocating[antennaIndex]) {
                continue;
            }

            // We are following this asteroid already, update direction to follow
            final Vector3d asteroidPosition = getAsteroidPosition(asteroidIndex);
            Vector3d antennaDirection = antennaDirections[antennaIndex];
            antennaDirection.x = asteroidPosition.x;
            antennaDirection.y = asteroidPosition.y;
            antennaDirection.z = asteroidPosition.z;
            antennaDirection.normalize();
        }
    }

    /*
     * Update relocating antennas directions given elapsed time
     */
    void updateRelocatingAntennasDirections(double deltaTime) {
        // If no time has passed (simultaneous events), no change
        if (deltaTime == 0.0) {
            return;
        }

        final double angle = deltaTime * testCase.RELOCATE_SPEED;  // how far we can turn an antenna given the time
        final double cosA = Math.cos(angle);
        final double oneMinusCosA = 1.0 - cosA;
        final double sinA = Math.sin(angle);
        Vector3d rotationAxis = new Vector3d();
        for (int antennaIndex = 0; antennaIndex < numberOfAntennas; ++antennaIndex) {
            final Vector3d antennaDirection = antennaDirections[antennaIndex];
            final int asteroidIndex = antennaAssignments[antennaIndex];

            // If the antenna is not relocating, ignore
            if (asteroidIndex == -1 || !isRelocating[antennaIndex]) {
                continue;
            }

            final Vector3d asteroidPosition = getAsteroidPosition(asteroidIndex);

            /*
             * Slew the antenna towards the asteroid position
             */
            rotationAxis.cross(antennaDirection, asteroidPosition);
            rotationAxis.normalize();
            final double ux = rotationAxis.x;
            final double uy = rotationAxis.y;
            final double uz = rotationAxis.z;

            // rotation matrix for rotation around rotationAxis with specified angle
            final double m00 = cosA + ux * ux * oneMinusCosA;
            final double m01 = ux * uy * oneMinusCosA - uz * sinA;
            final double m02 = ux * uz * oneMinusCosA + uy * sinA;
            final double m10 = ux * uy * oneMinusCosA + uz * sinA;
            final double m11 = cosA + uy * uy * oneMinusCosA;
            final double m12 = uy * uz * oneMinusCosA - ux * sinA;
            final double m20 = ux * uz * oneMinusCosA - uy * sinA;
            final double m21 = uy * uz * oneMinusCosA + ux * sinA;
            final double m22 = cosA + uz * uz * oneMinusCosA;

            // matrix multiplication
            final double newX = m00 * antennaDirection.x + m01 * antennaDirection.y + m02 * antennaDirection.z;
            final double newY = m10 * antennaDirection.x + m11 * antennaDirection.y + m12 * antennaDirection.z;
            final double newZ = m20 * antennaDirection.x + m21 * antennaDirection.y + m22 * antennaDirection.z;

            antennaDirection.x = newX;
            antennaDirection.y = newY;
            antennaDirection.z = newZ;

            antennaDirection.normalize(); // should be close to 1.0 already
        }
    }


    /*
     * Update current transmitted signal output power
     */
    void updateOutputSignalPower() {
        List<Integer> subarrays[] = getSubarrays();

        for (int asteroidIndex = 0; asteroidIndex < numberOfAsteroids; ++asteroidIndex) {
            // Ignore if asteroid not visible yet
            if (trajectoryIndices[asteroidIndex] == 0) {
                continue;
            }

            /*
             * Calculated the combined output power from all antennas
             */
            double amplitudeSum = 0.0;
            for (int antennaIndex : subarrays[asteroidIndex]) {
                amplitudeSum += Math.sqrt(transmittingPowers[antennaIndex]);
            }
            final double combinedOutputPower = amplitudeSum * amplitudeSum;

            /*
             * Calculate distance (R) to asteroid
             */
            final double R = getAsteroidPosition(asteroidIndex).length();

            // Find the power for last queued event
            final ArrayDeque<ReceiveEvent> Q = outputPowerQueue[asteroidIndex];
            final double lastOutputPower = Q.isEmpty() ? receivedSignalStrength[asteroidIndex] : Q.getLast().power;

            /*
             * Check if the received power differs from the last queued event, and add it to the queue in that case
             */
            if (combinedOutputPower != lastOutputPower) {
                final double travelTime = 2.0 * R / testCase.SPEED_OF_LIGHT;
                Q.add(new ReceiveEvent(this.time + travelTime, asteroidIndex, combinedOutputPower));
            }
        }
    }

    /*
     * Update image/trajectory information, energy spent and scores, given elapsed time
     */
    void updateScore(double deltaTime) {
        /*
         * If no time has passed (simultaneous events), no change
         */
        if (deltaTime == 0.0) {
            return;
        }

        /*
         * Update energy spent
         */
        for (int antennaIndex = 0; antennaIndex < numberOfAntennas; ++antennaIndex) {
            // If the antenna is not tracking anything, zero power consumption
            if (antennaAssignments[antennaIndex] == -1) {
                continue;
            }

            if (isRelocating[antennaIndex]) {
                this.energySpent += deltaTime * testCase.RELOCATE_POWER;
            }
            else {
                this.energySpent += deltaTime * transmittingPowers[antennaIndex];
            }
        }

        /*
         * Update accumulated information and score
         */
        final List<Integer>[] subarrays = getSubarrays();
        final double decayMultiplier = Math.exp(-deltaTime / TestCase.Q_LOST);
        for (int asteroidIndex = 0; asteroidIndex < numberOfAsteroids; ++asteroidIndex) {
            final List<Integer> antennaIndices = subarrays[asteroidIndex];
            final int N_receive = antennaIndices.size();

            double informationRate = 0.0;
            if (N_receive > 0) {
                final double powerReturn = effectiveSignalPowerReturn(asteroidIndex, N_receive);

                double inducedNoiseSum = 0.0;
                for (int receivingAntennaIndex : antennaIndices) {
                    inducedNoiseSum += inducedNoise(receivingAntennaIndex);
                }

                final double backgroundNoiseSum = N_receive * TestCase.BACKGROUND_NOISE;

                final double noise = TestCase.SHIELDING_MULTIPLIER * inducedNoiseSum + backgroundNoiseSum;

                final double ratio = powerReturn / noise;
                if (ratio > 1.0) {
                    informationRate = 10.0 * Math.log10(ratio);
                }
            }

            imageInformations[asteroidIndex] += deltaTime * informationRate;

            final double previousTrajectoryInformation = trajectoryInformations[asteroidIndex];
            final double newTrajectoryInformation = decayMultiplier * previousTrajectoryInformation + (1.0 - decayMultiplier) * testCase.Q_LOST * informationRate;
            trajectoryInformations[asteroidIndex] = newTrajectoryInformation;
            trajectoryScores[asteroidIndex] += deltaTime * (Math.tanh(previousTrajectoryInformation / testCase.Q_TRAJECTORY) + Math.tanh(newTrajectoryInformation / testCase.Q_TRAJECTORY)) / 2.0;
        }
    }

    final TestCase testCase;
    final int numberOfAntennas;
    final int numberOfAsteroids;

    int antennaAssignments[];          // which asteroid each antenna is assigned to, or -1 if not assigned
    boolean isRelocating[];            // true if antenna is relocating
    double transmittingPowers[];       // the transmitting power for each antenna
    Vector3d antennaDirections[];      // the main target direction for each antenna

    double imageInformations[];        // image information for each asteroid
    double trajectoryInformations[];   // trajectory information for each asteroid
    double trajectoryScores[];         // trajectory score for each asteroid
    double receivedSignalStrength[];   // the received signal from each asteroid
    ArrayDeque<ReceiveEvent> outputPowerQueue[];  // a queue with new signals that will return in the future from each asteroid
    int trajectoryIndices[];           // which trajectory index we are currently at for each asteroid

    double energySpent;                // total energy spent
    double time;                       // current time in seconds of the simulation
}
