/**
 * Copyright (C) 2014 TopCoder Inc., All Rights Reserved.
 */
package asteroidtracker;

class ConstraintException extends Exception {
    int antennaTransmitterIndex;
    int antennaReceiverIndex;

    ConstraintException(String msg) {
        super(msg);
        this.antennaTransmitterIndex = -1;
        this.antennaReceiverIndex = -1;
    }

    ConstraintException(int antennaTransmitterIndex, String msg) {
        super(msg);
        this.antennaTransmitterIndex = antennaTransmitterIndex;
        this.antennaReceiverIndex = -1;
    }

    ConstraintException(int antennaTransmitterIndex, int antennaReceiverIndex, String msg) {
        super(msg);
        this.antennaTransmitterIndex = antennaTransmitterIndex;
        this.antennaReceiverIndex = antennaReceiverIndex;
    }
}
