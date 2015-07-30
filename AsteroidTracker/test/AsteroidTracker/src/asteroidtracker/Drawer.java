/**
 * Copyright (C) 2014 TopCoder Inc., All Rights Reserved.
 */
package asteroidtracker;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Image;
import java.awt.Insets;
import java.awt.event.KeyListener;
import java.awt.event.WindowListener;
import java.awt.geom.Point2D;
import java.awt.image.BufferedImage;
import java.awt.Point;
import java.util.List;
import javax.swing.JFrame;


class Drawer extends Component {
    static private int antennaRenderRadius = 15;

    private JFrame frame;
    private double[] latitudeCircles;
    private double[] longitudeLines;
    private Point2D.Double[] antennaPositions;
    private Simulator simulator = null;

    private int hemisphereRadius;
    private Image backgroundImage;
    private boolean backgroundDirty;
    private ConstraintException constraintException;

    /*
     *  Provide all information about the antennas and the asteroids in this testcase
     */
    public void initialize(String title,                   // title of the window
                           int hemisphereRadius,           // initial radius of the hemisphere in pixels, but the window should be resizable
                           KeyListener keyListener,        // key listener for this GUI
                           WindowListener windowListener,  // window listener for this GUI
                           double[] latitudeCircles,        // list of normalized radius (between [0, 1]) for the latitude circles
                           double[] longitudeLines,         // list of angels in radians for the longitude lines
                           Point2D.Double[] antennaPositions,  // list of normalized antenna positions (x and y between [0, 1])
                           Simulator simulator)
    {
        this.hemisphereRadius = hemisphereRadius;
        this.backgroundDirty = true;
        this.latitudeCircles = latitudeCircles;
        this.longitudeLines = longitudeLines;
        this.antennaPositions = antennaPositions;
        this.simulator = simulator;
        this.backgroundDirty = true;
        this.constraintException = null;

        frame = new JFrame(title);
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.getContentPane().add(this);
        frame.pack();
        Insets insets = frame.getInsets();
        frame.setSize(new Dimension(insets.left + insets.right+hemisphereRadius * 2 + 200,
                                    insets.top + insets.bottom+hemisphereRadius * 2));
        frame.setResizable(true);
        frame.addKeyListener(keyListener);
        frame.addWindowListener(windowListener);
        frame.setVisible(true);
    }

    public void setConstraintException(ConstraintException constraintException) {
        this.constraintException = constraintException;
    }

    private int normalizedToPixel(double d) {
        int id = (int)(hemisphereRadius * 2 * d);
        return Math.min(Math.max(0, id), hemisphereRadius * 2 - 1);
    }

    private Point normalizedToPixel(Point2D.Double point) {
        return new Point(normalizedToPixel(point.x), normalizedToPixel(point.y));
    }

    public void paint(Graphics g) {
        /*
         * Exit if we haven't been initialized
         */
        if (this.simulator == null) {
            return;
        }

        Graphics2D graphics = (Graphics2D)g;
        Font font = new Font(Font.MONOSPACED, Font.PLAIN, 12);
        g.setFont(font);

        /*
         * Draw background
         */
        {
            Dimension size = this.getSize();
            int radius = Math.min(size.width,size.height) / 2;
            if (radius != this.hemisphereRadius) {
                backgroundDirty = true;
                this.hemisphereRadius = radius;
            }
            if (backgroundDirty) {
                updateBackgroundImage();
                backgroundDirty = false;
            }
            graphics.drawImage(backgroundImage, 0, 0, null);
        }

        /*
         * Project asteroid positions onto circle
         */
        Point[] asteroidPositions = new Point[simulator.numberOfAsteroids];
        for (int asteroidIndex = 0; asteroidIndex < simulator.numberOfAsteroids; ++asteroidIndex) {
            Vector3d asteroidPositionVec = simulator.getAsteroidPosition(asteroidIndex);
            if (asteroidPositionVec == null) {
                continue;
            }

            final double dist = asteroidPositionVec.length();

            Point2D.Double projectedPoint = new Point2D.Double();
            projectedPoint.x = 0.5 + 0.5 * asteroidPositionVec.x / dist;
            projectedPoint.y = 0.5 - 0.5 * asteroidPositionVec.y / dist;

            asteroidPositions[asteroidIndex] = asteroidPositionVec.z > 0.0 ? normalizedToPixel(projectedPoint) : null;
        }


        Point[] antennaAimPositions = new Point[simulator.numberOfAntennas];
        for (int antennaIndex = 0; antennaIndex < simulator.numberOfAntennas; ++antennaIndex) {
            final Point antennaPixel = normalizedToPixel(antennaPositions[antennaIndex]);
            final Vector3d dir = simulator.antennaDirections[antennaIndex];
            antennaAimPositions[antennaIndex] = new Point(antennaPixel.x + (int)Math.round(antennaRenderRadius * dir.x),
                    antennaPixel.y - (int)Math.round(antennaRenderRadius * dir.y));
        }

        /*
         *  Render assignment
         */
        List<Integer>[] subarrays = simulator.getSubarrays();
        graphics.setStroke(new BasicStroke(1));
        for (int antennaIndex = 0; antennaIndex < simulator.numberOfAntennas; ++antennaIndex) {
            int asteroidIndex = simulator.antennaAssignments[antennaIndex];

            if (asteroidIndex != -1 && simulator.transmittingPowers[antennaIndex] > 0.0 && !simulator.isRelocating[antennaIndex]) {
                graphics.setColor(Color.black);

                Point asteroidPixel = asteroidPositions[asteroidIndex];
                if (asteroidPixel != null) {
                    Point antennaAim = antennaAimPositions[antennaIndex];

                    if (constraintException != null && antennaIndex == constraintException.antennaTransmitterIndex) {
                        graphics.setColor(Color.red);
                    }

                    graphics.drawLine(asteroidPixel.x, asteroidPixel.y, antennaAim.x, antennaAim.y);
                }
            }
        }

        /*
         *  Render antenna direction
         */
        for (int antennaIndex = 0; antennaIndex < simulator.numberOfAntennas; ++antennaIndex) {
            int aimPointRadius = 1;
            graphics.setColor(Color.white);
            Point antennaAim = antennaAimPositions[antennaIndex];

            graphics.fillOval(antennaAim.x - aimPointRadius,
                              antennaAim.y - aimPointRadius,
                              aimPointRadius * 2, aimPointRadius * 2);
        }

        /*
         *  Render asteroids
         */
        graphics.setStroke(new BasicStroke(2));
        int asteroidRenderRadius = 5;
        int yOffset = 15;
        int textHeight = graphics.getFontMetrics().getHeight();
        for (int asteroidIndex = 0; asteroidIndex < simulator.numberOfAsteroids; ++asteroidIndex) {
            Point asteroidPos = asteroidPositions[asteroidIndex];
            if (asteroidPos != null) {
                List<Integer> antennaIndices = subarrays[asteroidIndex];
                final int N_receive = antennaIndices.size();

                graphics.setColor(N_receive == 0 ? Color.white : Color.green);

                graphics.drawOval(asteroidPos.x - asteroidRenderRadius, asteroidPos.y - asteroidRenderRadius,
                                  asteroidRenderRadius * 2, asteroidRenderRadius * 2);

                if (N_receive > 0) {
                    final double powerReturn = simulator.effectiveSignalPowerReturn(asteroidIndex, N_receive);

                    double inducedNoiseSum = 0.0;
                    for (int receivingAntennaIndex : antennaIndices) {
                        inducedNoiseSum += simulator.inducedNoise(receivingAntennaIndex);
                    }

                    final double backgroundNoiseSum = N_receive * TestCase.BACKGROUND_NOISE;

                    final double noise = TestCase.SHIELDING_MULTIPLIER * inducedNoiseSum + backgroundNoiseSum;

                    graphics.drawString(String.format("Signal: %.2f", powerReturn * 1.0e30), asteroidPos.x + asteroidRenderRadius + 2, asteroidPos.y - 2 * textHeight);
                    graphics.drawString(String.format("Noise:  %.2f", noise * 1.0e30), asteroidPos.x + asteroidRenderRadius + 2, asteroidPos.y - 1 * textHeight);
                }

                final TestCase testCase = simulator.testCase;
                final double scienceScoreMultiplier = testCase.asteroids[asteroidIndex].scienceScoreMultiplier;
                final double trajectoryInformation = simulator.trajectoryInformations[asteroidIndex] / testCase.Q_TRAJECTORY;
                final double imageScore = scienceScoreMultiplier * testCase.IMAGE_SCORE_MULTIPLIER * simulator.getImageScore(asteroidIndex);

                final double trajectoryScore = scienceScoreMultiplier * testCase.TRAJECTORY_SCORE_MULTIPLIER * simulator.trajectoryScores[asteroidIndex];

                graphics.drawString(String.format("Tr.Inf: %.2f", trajectoryInformation), asteroidPos.x + asteroidRenderRadius + 2, asteroidPos.y + 0 * textHeight);
                graphics.drawString(String.format("Img.Sc: %.2f", imageScore), asteroidPos.x + asteroidRenderRadius + 2, asteroidPos.y + 1 * textHeight);
                graphics.drawString(String.format("Tr.Sc:  %.2f", trajectoryScore), asteroidPos.x + asteroidRenderRadius + 2, asteroidPos.y + 2 * textHeight);
            }
        }

        /*
         * Render antenna transmitting power information
         */
        for (int antennaIndex = 0; antennaIndex < simulator.numberOfAntennas; ++antennaIndex) {
            double transmittingPower = simulator.transmittingPowers[antennaIndex];

            Point2D.Double pos = this.antennaPositions[antennaIndex];
            Point pixel = normalizedToPixel(pos);
            graphics.setColor(Color.yellow);

            if (transmittingPower > 0) {
                graphics.drawString(String.format("Output: %d", (int)Math.round(transmittingPower)), pixel.x + antennaRenderRadius + 4, pixel.y);
            }
        }

        /*
         * Print score and energy information
         */
        paintInformation(g);
    }

    private void paintInformation(Graphics g) {
        Font font = new Font("Sans", Font.BOLD, 16);

        g.setFont(font);
        g.setColor(Color.black);
        int rightOffset = 250;
        int lineHeight = font.getSize() + 4;
        int yPos = lineHeight;
        int leftPos = getSize().width - rightOffset;
        g.drawString(textOfTime((int)Math.round(simulator.time)), leftPos, yPos);

        yPos += lineHeight;
        g.drawString(String.format("Energy used: %5.1f 10^9", simulator.energySpent * 1.0e-9), leftPos, yPos);
        yPos += lineHeight;
        g.drawString(String.format("Score: %3.4f", simulator.getScore()), leftPos, yPos);
        yPos += lineHeight;
    }

    private String textOfTime(int time) {
        int s = time % 60;
        time /= 60;
        int m = time % 60;
        time /= 60;
        int h = time;
        return String.format("Time: %02d:%02d:%02d",h,m,s);
    }

    public void updateBackgroundImage() {
        this.backgroundImage = new BufferedImage(hemisphereRadius * 2, hemisphereRadius * 2, BufferedImage.TYPE_INT_ARGB);

        Graphics2D graphics = (Graphics2D)this.backgroundImage.getGraphics();
        graphics.setColor(Color.blue);
        graphics.fillOval(0, 0, hemisphereRadius * 2, hemisphereRadius * 2);

        /*
         * Latitude circles
         */
        for (double normRadius : this.latitudeCircles) {
            Point pixel = normalizedToPixel(new Point2D.Double(normRadius,normRadius));
            int radius = pixel.x;
            graphics.setColor(Color.white);
            graphics.drawOval(hemisphereRadius - radius, hemisphereRadius - radius,
                              2 * radius, 2 * radius);
        }

        /*
         * Longitude lines
         */
        for (double longitude : this.longitudeLines) {
            Point2D.Double target = new Point2D.Double(0.5f + 0.5f * (double)Math.cos(longitude),
                    0.5f + 0.5f * (double)Math.sin(longitude));
            Point targetPixel = normalizedToPixel(target);
            Point sourcePixel = normalizedToPixel(new Point2D.Double(0.5f, 0.5f));
            graphics.drawLine(sourcePixel.x, sourcePixel.y, targetPixel.x, targetPixel.y);
        }

        /*
         *  Draw antennas
         */
        graphics.setColor(Color.red);
        graphics.setStroke(new BasicStroke(2));
        int yOffset = 15;
        int xOffset = -30;
        for (int antennaIndex = 0; antennaIndex < antennaPositions.length; ++antennaIndex) {
            Point2D.Double point = antennaPositions[antennaIndex];
            Point pixel = normalizedToPixel(point);

            graphics.setColor(Color.red);
            graphics.drawOval(pixel.x - antennaRenderRadius,
                              pixel.y - antennaRenderRadius,
                              antennaRenderRadius * 2, antennaRenderRadius * 2);

            graphics.setColor(Color.yellow);
        }
    }
}
