/**
 * Copyright (C) 2014 TopCoder Inc., All Rights Reserved.
 */
package asteroidtracker;

import java.util.Scanner;
import java.awt.Point;
import java.awt.event.KeyAdapter;
import java.awt.event.KeyEvent;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.awt.geom.Point2D;

import org.apache.commons.cli.*;
import java.util.List;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Queue;
import java.io.IOException;
import java.io.FileNotFoundException;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;

class VisualizerMain {
    static final int DEFAULT_SIMULATION_SPEED = 1;
    static final int DEFAULT_SEED = 1;

    static public void main(String[] args) throws FileNotFoundException, IOException, ClassNotFoundException {
        Options options = new Options();
        options.addOption("help", "h", false, "Print this message.");
        options.addOption("exec", "e", true, "Sets the command you would use to execute your solution. "
                          + "If your compiled solution is an executable file, the command will just be the full path to it.");
        options.addOption("seed", "s", true, "Sets the seed used for test case generation. The default value is " + DEFAULT_SEED + ".");
        options.addOption("pause", "p", false, "Starts visualizer in paused mode.");
        options.addOption("simulation_speed", "S", true, "Sets number of events between visualized frames. The default value is " + DEFAULT_SIMULATION_SPEED + ".");
        options.addOption("no_visualization", "n", false, "Switches the visualization off, leaving only text output.");

        CommandLine cmd = null;
        boolean doPrintHelp = true;
        Process solution = null;

        class CommandLineException extends Exception {
            CommandLineException(String msg) {
                super(msg);
            }
        }

        int seed = DEFAULT_SEED;
        int simulationSpeed = DEFAULT_SIMULATION_SPEED;
        boolean isVisualizing = true;
        boolean isPaused = false;
        try {
            CommandLineParser parser = new BasicParser();
            try {
                cmd = parser.parse(options, args);
                seed = Integer.parseInt(cmd.getOptionValue("seed", Integer.toString(DEFAULT_SEED)));
                simulationSpeed = Integer.parseInt(cmd.getOptionValue("simulation_speed", Integer.toString(DEFAULT_SIMULATION_SPEED)));
                isVisualizing = !cmd.hasOption("no_visualization");
                isPaused = cmd.hasOption("pause");
            }
            catch (ParseException e) {
                throw new CommandLineException("ERROR: Could not parse your arguments. " + e.getMessage());
            }
            catch (NumberFormatException e) {
                throw new CommandLineException("ERROR: Could not parse your arguments. " + e.getMessage());
            }

            if (cmd.hasOption("help")) {
                throw new CommandLineException("help");
            }

            if (!cmd.hasOption("exec")) {
                throw new CommandLineException("ERROR: You did not provide the command to execute your solution. Please use -exec <command> for this.");
            }

            doPrintHelp = false;

            String execCommand = cmd.getOptionValue("exec");
            try {
                solution = Runtime.getRuntime().exec(execCommand);
            }
            catch (IOException e) {
                throw new CommandLineException("ERROR: Unable to execute your solution using the provided command: " + execCommand + ". " + e.getMessage());
            }
        }
        catch (CommandLineException e) {
            System.err.println(e.getMessage());

            if (doPrintHelp) {
                HelpFormatter formatter = new HelpFormatter();
                formatter.printHelp("AsteroidTrackerVisualizer", options);
            }

            System.exit(0);
        }

        final TestCase testCase = new TestCase();
        testCase.generate("testcase" + seed + ".txt");

        IAsteroidTracker asteroidTracker = new IoAsteroidTracker(solution);

        VisualizerMain vis = new VisualizerMain(testCase, asteroidTracker, isVisualizing, isPaused, simulationSpeed);
        vis.run();
    }

    VisualizerMain(TestCase testCase,
                   IAsteroidTracker asteroidTracker,
                   boolean isVisualizing,
                   boolean isPaused,
                   int simulationSpeed) {
        this.testCase = testCase;
        this.asteroidTracker = asteroidTracker;
        this.simulationSpeed = simulationSpeed;
        this.isVisualizing = isVisualizing;

        this.keyMutex = new Object();
        this.pauseMode = isPaused;
        this.keyPressed = false;
    }

    private final TestCase testCase;
    private  IAsteroidTracker asteroidTracker;
    private int simulationSpeed;
    private final boolean isVisualizing;

    private final Object keyMutex;
    private boolean pauseMode;
    private boolean keyPressed;

    public void run() {
        /*
         *  Initialize simulator and asteroid tracker
         */
        Simulator simulator = new Simulator(testCase);
        Drawer drawer = isVisualizing ? makeDrawer(15, 20, simulator) : null;

        /*
         *  Provide all information and data in this testcase
         */
        {
            double antennaPositionsAsArray[] = new double[2 * simulator.numberOfAntennas];
            for (int antennaIndex = 0; antennaIndex < simulator.numberOfAntennas; ++antennaIndex) {
                antennaPositionsAsArray[antennaIndex * 2 + 0] = testCase.antennaPositions[antennaIndex].x;
                antennaPositionsAsArray[antennaIndex * 2 + 1] = testCase.antennaPositions[antennaIndex].y;
            }

            asteroidTracker.initialize(antennaPositionsAsArray,
                                       testCase.peakGain,
                                       testCase.minDistanceGain,
                                       testCase.maxDistanceGain);
        }

        /*
         * Simulation loop
         */
        int frameCount = 0;
        Event nextCommand = null;
        while (simulator.time < testCase.SIMULATION_TIME) {
            List<Event> deterministicEvents = new ArrayList<Event>();
            List<AsteroidAppearanceEvent> appearanceEvents = new ArrayList<AsteroidAppearanceEvent>();

            // Don't exceed end of simulation
            deterministicEvents.add(new ClockEvent(testCase.SIMULATION_TIME));

            /*
             * Add events for when antennas are done relocating
             */
            for (int antennaIndex = 0; antennaIndex < simulator.numberOfAntennas; ++antennaIndex) {
                if (simulator.isRelocating[antennaIndex]) {
                    Vector3d antennaDirection = simulator.antennaDirections[antennaIndex];

                    int asteroidIndex = simulator.antennaAssignments[antennaIndex];
                    Vector3d asteroidPosition = simulator.getAsteroidPosition(asteroidIndex);
                    double angle = antennaDirection.angle(asteroidPosition);

                    double timeLeft = angle / testCase.RELOCATE_SPEED;

                    deterministicEvents.add(new DoneRelocatingAntennaEvent(simulator.time + timeLeft, antennaIndex));
                }
            }

            /*
             * Add events concerning asteroids
             */
            for (int asteroidIndex = 0; asteroidIndex < simulator.numberOfAsteroids; ++asteroidIndex) {
                final int nextTrajectoryIndex = simulator.trajectoryIndices[asteroidIndex];
                AsteroidSighting trajectory[] = testCase.asteroids[asteroidIndex].trajectory;

                if (nextTrajectoryIndex < trajectory.length) {
                    if (nextTrajectoryIndex == 0) {
                        /*
                         * Add event for asteroid appearance
                         */
                        appearanceEvents.add(new AsteroidAppearanceEvent(trajectory[0].time, asteroidIndex));
                    }
                    else {
                        /*
                         * Add events for when asteroids change position
                         */
                        deterministicEvents.add(new AsteroidUpdatePositionEvent(trajectory[nextTrajectoryIndex].time, asteroidIndex));
                    }
                }

                /*
                 * Add events for when signal reception from asteroids change
                 */
                Queue<ReceiveEvent> Q = simulator.outputPowerQueue[asteroidIndex];
                if (!Q.isEmpty()) {
                    deterministicEvents.add(Q.element());
                }
            }

            /*
             * Add user event
             */
            if (nextCommand == null) {
                final String command = asteroidTracker.nextCommand(simulator.time);
                try {
                    nextCommand = parseCommand(simulator, command);
                }
                catch (ConstraintException e) {
                    if (drawer != null) {
                        drawer.setConstraintException(e);
                    }
                    System.err.println("Error at time: " + simulator.time + " while parsing your command: " + command);
                    System.err.println("Error message: " + e.getMessage());
                    break;
                }
            }

            /*
             * Find the event that will happen next
             */
            final Event nextDeterministicEvent = Collections.min(deterministicEvents);
            final Event nextAsteroidAppearance = appearanceEvents.isEmpty() ? null : Collections.min(appearanceEvents);
            Event nextEvent;

            if (nextAsteroidAppearance != null && nextAsteroidAppearance.time <= nextCommand.time && nextAsteroidAppearance.time <= nextDeterministicEvent.time) {
                nextEvent = nextAsteroidAppearance;
                nextCommand = null; // ignore this command
            }
            else if (nextCommand.time < nextDeterministicEvent.time) {
                nextEvent = nextCommand;
                nextCommand = null; // ask for a new command next time
            }
            else {
                nextEvent = nextDeterministicEvent;
                // don't ask for a new command next time
            }

            final double deltaTime = nextEvent.time - simulator.time;

            /*
             * Update state for this time period, update score, then execute next event and update future changes in output power
             */
            try {
                if (deltaTime > 0.0) {
                    simulator.checkAntennaBeamIntersection();
                    simulator.checkTargetProximity();
                    simulator.updateTrackingAntennasDirections();
                    simulator.updateRelocatingAntennasDirections(deltaTime);
                    simulator.updateScore(deltaTime);
                    simulator.time = nextEvent.time;
                }

                nextEvent.execute(simulator, asteroidTracker);
                simulator.updateOutputSignalPower();
            }
            catch (ConstraintException e) {
                System.err.println("Error at time: " + simulator.time);
                System.err.println("Error message: " + e.getMessage());
                if (drawer != null) {
                    drawer.setConstraintException(e);
                }
                break;
            }

            /*
             * Update visualization
             */
            if (drawer != null && deltaTime > 0.0) {
                if (frameCount++ % simulationSpeed == 0 || this.pauseMode) {
                    this.processPause();
                    drawer.repaint();

                    /*
                     *	Limit frame rate
                     */
                    try {
                        Thread.sleep(10);
                    }
                    catch (InterruptedException e) {
                        // do nothing
                    }
                }
            }
        }

        if (drawer != null) {
            drawer.repaint();
        }
        final double finalScore = simulator.time < testCase.SIMULATION_TIME ? 0.0 : simulator.getScore();

        System.out.println("Score: " + finalScore);

        if (drawer == null) {
            System.exit(0);
        }
    }

    static Event parseCommand(Simulator simulator, String command) throws ConstraintException {
        Scanner scanner = new Scanner(command);

        double executeTime = scanner.nextDouble();
        if (executeTime < simulator.time) {
            throw new ConstraintException("Your solution returned a time = " + executeTime + " which is before the current time = " + simulator.time);
        }

        int antennaIndex = scanner.nextInt();
        if (antennaIndex < 0) {
            throw new ConstraintException("Your solution returned a negative antennaIndex = " + antennaIndex);
        }
        else if (antennaIndex >= simulator.numberOfAntennas) {
            throw new ConstraintException("Your solution returned an out of bounds antennaIndex = " + antennaIndex + " (0-based). Number of antennas = " +  simulator.numberOfAntennas);
        }

        String commandType = scanner.next();

        if (commandType.equals("R")) {
            int asteroidIndex = scanner.nextInt();

            if (asteroidIndex < -1 || asteroidIndex >= simulator.numberOfAsteroids) {
                throw new ConstraintException(antennaIndex, "Your solution has returned an invalid asteroid index: " + asteroidIndex
                                              + ". It should be between 0 and " + (simulator.numberOfAsteroids - 1) + " inclusive, or -1 if turning off antenna");
            }

            if (asteroidIndex != -1 && simulator.trajectoryIndices[asteroidIndex] == 0) {
                throw new ConstraintException(antennaIndex, "Your solution tried to relocate to an asteroid that has not yet appeared");
            }

            return new StartRelocatingAntennaEvent(executeTime, antennaIndex, asteroidIndex);
        }
        else if (commandType.equals("P")) {
            double transmittingPower = scanner.nextDouble();

            /*
             *	Check that antenna power is valid
             */
            if (transmittingPower < 0) {
                throw new ConstraintException(antennaIndex, "Your solution is trying to use a negative power: " + transmittingPower + ".");
            }
            else if (transmittingPower > TestCase.MAX_TRANSMITTING_POWER) {
                throw new ConstraintException(antennaIndex, "Your solution is trying to use too much power: " + transmittingPower + ", maximum is: " + TestCase.MAX_TRANSMITTING_POWER + ".");
            }

            return new ChangeTransmittingPowerEvent(executeTime, antennaIndex, transmittingPower);
        }
        else {
            throw new ConstraintException(antennaIndex, "Your solution returned an invalid command type: " + commandType + ". It should be 'R' or 'P'");
        }
    }

    private Drawer makeDrawer(int numberOfCircels, int numberOfLines, Simulator simulator) {
        double[] latitudeCircles = new double[numberOfCircels];
        double[] longitudeLines = new double[numberOfLines];

        double antennaMinX = Double.MAX_VALUE;
        double antennaMaxX = Double.MIN_VALUE;
        double antennaMinY = Double.MAX_VALUE;
        double antennaMaxY = Double.MIN_VALUE;

        for (int antennaIndex = 0; antennaIndex < simulator.numberOfAntennas; ++antennaIndex) {
            double x = testCase.antennaPositions[antennaIndex].x;
            double y = testCase.antennaPositions[antennaIndex].y;

            antennaMinX = Math.min(antennaMinX, x);
            antennaMaxX = Math.max(antennaMaxX, x);
            antennaMinY = Math.min(antennaMinY, y);
            antennaMaxY = Math.max(antennaMaxY, y);
        }

        Point2D.Double normalizedAntennaPositions[] = new Point2D.Double[simulator.numberOfAntennas];
        for (int antennaIndex = 0; antennaIndex < simulator.numberOfAntennas; ++antennaIndex) {
            double rx = testCase.antennaPositions[antennaIndex].x;
            double ry = testCase.antennaPositions[antennaIndex].y;

            double x = 0.5f;
            double y = 0.5f;

            double renderScale = 0.4f;

            if (antennaMaxX > antennaMinX) {
                x += renderScale * (-0.5f + (rx - antennaMinX) / (antennaMaxX - antennaMinX));
            }

            if (antennaMaxY > antennaMinY) {
                y -= renderScale * (-0.5f + (ry - antennaMinY) / (antennaMaxY - antennaMinY));
            }

            normalizedAntennaPositions[antennaIndex] = new Point2D.Double(x, y);
        }

        for (int i = 0; i < numberOfCircels; ++i) {
            double r = (double)i / (double)(numberOfCircels - 1);
            double elevation = (double)Math.PI / 2.0f * (1.0f - r);
            double radius = 0.5f * (double)Math.cos(elevation);
            latitudeCircles[i] = radius;
        }

        for (int i = 0; i < numberOfLines; ++i) {
            double azimuth = 2.0f * (double)Math.PI * (1.0f - (double)i / (double)(numberOfLines - 1));
            longitudeLines[i] = azimuth;
        }

        Drawer drawer = new Drawer();
        drawer.initialize("Asteroid Tracker Visualizer (space to pause, q to quit, +/- to change simulation speed)",
                          300,
                          new VisualizerKeyListener(),
                          new VisualizerWindowListener(),
                          latitudeCircles,
                          longitudeLines,
                          normalizedAntennaPositions,
                          simulator);
        return drawer;
    }

    private class VisualizerWindowListener extends WindowAdapter {
        public void windowClosing(WindowEvent event) {
            System.exit(0);
        }
    }

    private class VisualizerKeyListener extends KeyAdapter {
        public void keyPressed(KeyEvent e) {
            synchronized (keyMutex) {
                switch (e.getKeyChar()) {
                case ' ':
                    pauseMode = !pauseMode;
                    break;

                case 'q':
                    System.exit(0);
                    break;

                case '+':
                    simulationSpeed += 1;
                    break;

                case '-':
                    simulationSpeed = Math.max(1, simulationSpeed - 1);
                    break;
                };

                keyPressed = true;
                keyMutex.notifyAll();
            }
        }
    }

    private void processPause() {
        synchronized (this.keyMutex) {
            if (!this.pauseMode) {
                return;
            }
            this.keyPressed = false;
            while (!this.keyPressed) {
                try {
                    this.keyMutex.wait();
                } catch (InterruptedException e) {
                    // do nothing
                }
            }
        }
    }
}
