/**
 * Copyright (C) 2014 TopCoder Inc., All Rights Reserved.
 */
package asteroidtracker;

/*
 * Base class for Events: time and comparison
 */
abstract class Event implements Comparable<Event> {
    public Event(double time) {
        this.time = time;
    }

    public int compareTo(Event that) {
        return Double.compare(this.time, that.time);
    }

    public double time;

    /*
     * This function updates the state when the event is executed
     */
    abstract public void execute(Simulator simulator, IAsteroidTracker asteroidTracker) throws ConstraintException;
}


class ClockEvent extends Event {
    public ClockEvent(double time) {
        super(time);
    }

    public void execute(Simulator simulator, IAsteroidTracker asteroidTracker) {
        // do nothing
    }
}


class ChangeTransmittingPowerEvent extends Event {
    private int antennaIndex;
    private double power;

    ChangeTransmittingPowerEvent(double time, int antennaIndex, double power) {
        super(time);
        this.antennaIndex = antennaIndex;
        this.power = power;
    }

    public void execute(Simulator simulator, IAsteroidTracker asteroidTracker) {
        simulator.transmittingPowers[antennaIndex] = power;
    }
}


class StartRelocatingAntennaEvent extends Event {
    int antennaIndex;
    int asteroidIndex;

    StartRelocatingAntennaEvent(double time, int antennaIndex, int asteroidIndex) {
        super(time);
        this.antennaIndex = antennaIndex;
        this.asteroidIndex = asteroidIndex;
    }

    public void execute(Simulator simulator, IAsteroidTracker asteroidTracker) {
        simulator.antennaAssignments[antennaIndex] = asteroidIndex;
        simulator.isRelocating[antennaIndex] = (asteroidIndex != -1);
    }
}


class DoneRelocatingAntennaEvent extends Event {
    private int antennaIndex;

    DoneRelocatingAntennaEvent(double time, int antennaIndex) {
        super(time);
        this.antennaIndex = antennaIndex;
    }

    public void execute(Simulator simulator, IAsteroidTracker asteroidTracker) throws ConstraintException {
        final int asteroidIndex = simulator.antennaAssignments[antennaIndex];

        if (simulator.getAsteroidPosition(asteroidIndex).z < 0.0) {
            throw new ConstraintException(antennaIndex, "Your solution tried to relocate to asteroid " + asteroidIndex + " that is below the horizon, with antenna " + antennaIndex);
        }

        simulator.isRelocating[antennaIndex] = false;
    }
}


class ReceiveEvent extends Event {
    private int asteroidIndex;
    public double power;

    ReceiveEvent(double time, int asteroidIndex, double power) {
        super(time);
        this.asteroidIndex = asteroidIndex;
        this.power = power;
    }

    public void execute(Simulator simulator, IAsteroidTracker asteroidTracker) {
        simulator.receivedSignalStrength[asteroidIndex] = power;
        simulator.outputPowerQueue[asteroidIndex].remove();
    }
}

class AsteroidAppearanceEvent extends Event {
    private int asteroidIndex;

    public AsteroidAppearanceEvent(double time, int asteroidIndex) {
        super(time);
        this.asteroidIndex = asteroidIndex;
    }

    public void execute(Simulator simulator, IAsteroidTracker asteroidTracker) {
        Asteroid asteroid = simulator.testCase.asteroids[asteroidIndex];

        // Set the intial information
        simulator.trajectoryIndices[asteroidIndex] = 1;
        simulator.imageInformations[asteroidIndex] = TestCase.Q_IMAGE * asteroid.initialImageInformation;
        simulator.trajectoryInformations[asteroidIndex] = TestCase.Q_TRAJECTORY * asteroid.initialTrajectoryInformation;

        // Reveal information for asteroid tracker
        int len = asteroid.trajectory.length;
        double trajectory[] = new double[len * 4];

        for (int i = 0; i < len; ++i) {
            trajectory[i * 4 + 0] = asteroid.trajectory[i].time;
            trajectory[i * 4 + 1] = asteroid.trajectory[i].pos.x;
            trajectory[i * 4 + 2] = asteroid.trajectory[i].pos.y;
            trajectory[i * 4 + 3] = asteroid.trajectory[i].pos.z;
        }

        asteroidTracker.asteroidAppearance(asteroidIndex,
                                           asteroid.scienceScoreMultiplier,
                                           asteroid.reflectivityMultiplier,
                                           asteroid.initialImageInformation,
                                           asteroid.initialTrajectoryInformation,
                                           trajectory);
    }
}


class AsteroidUpdatePositionEvent extends Event {
    private int asteroidIndex;

    public AsteroidUpdatePositionEvent(double time, int asteroidIndex) {
        super(time);
        this.asteroidIndex = asteroidIndex;
    }

    public void execute(Simulator simulator, IAsteroidTracker asteroidTracker) {
        simulator.trajectoryIndices[asteroidIndex] += 1;
    }
}
