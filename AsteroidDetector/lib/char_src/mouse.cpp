/*****************************************************************************
 *                   ibmcom.c                                                *
 *****************************************************************************
 * DESCRIPTION:    This file contains a set of routines for doing low-level  *
 *        serial communications on the IBM PC.  It was translated            *
 *        directly from Wayne Conrad's IBMCOM.PAS version 3.1, with          *
 *        the goal of near-perfect functional correspondence between         *
 *        the Pascal and C versions.                                         *
 *                                                                           *
 * REVISIONS:    18 OCT 89 - RAC - Original translation from IBMCOM.PAS,     *
 *                  with liberal plagiarism of comments from the             *
 *                  Pascal.                                                  *
 *                                                                           *
 *               4 Jan 94 - Bill J. Gray - added #def's for Watcom and       *
 *                  MSC compilers;  fixed some idiosyncratic braces;         *
 *                  added a #ifdef around the definition for old_vector      *
 *                  because it causes fatal errors in Watcom and MSC.        *
 *                                                                           *
 *               5 Jan 94 - BJG - changed rval for com_rx to -1 to           *
 *                 indicate "no data",  to allow xmission/reception          *
 *                 of binary data.  Also changed to allow setting of         *
 *                 buffer sizes at initialisation time.                      *
 *****************************************************************************/

/* Microsoft mouse data protocol (from January 1998 _Electronics Now_):     */
/* Each data packet contains three bytes.  The high bit in each byte may    */
/* or may not be set;  it's the remaining seven bits that matter.  The      */
/* link is 1200 baud, one stop bit,  no parity.                             */
/*                                                                          */
/* 1st byte    1 LB RB Y7 Y6 X7 X6                                          */
/* 2nd byte    0 X5 X4 X3 X2 X1 X0                                          */
/* 3rd byte    0 Y5 Y4 Y3 Y2 Y1 Y0                                          */
/*                                                                          */
/* As you can see,  the positional data is spread over all three bytes.     */
/* Bit 7 gives you a way to 'resynch' and account for lost bytes.           */

#include    <stdio.h>
#include    <stdlib.h>
#include    <conio.h>
#include    <dos.h>
#define MOUSE_DOT_C
#include    "mouse.h"

#ifdef MICROSOFT_C
#define inportb     _inp
#define outportb    _outp
#define outport     _outpw
#define enable      _enable
#define disable     _disable
#define getvect     _dos_getvect
#define setvect     _dos_setvect
#endif

#ifdef __WATCOMC__
#define inportb     inp
#define outportb    outp
#define outport     outpw
#define enable      _enable
#define disable     _disable
#define getvect     _dos_getvect
#define setvect     _dos_setvect
#endif

/*****************************************************************************
 *                   8250 Definitions                                        *
 *****************************************************************************/

/*      Offsets to various 8250 registers.  Taken from IBM Technical         */
/*      Reference Manual, p. 1-225                                           */

#define TXBUFF  0                       /* Transmit buffer register */
#define RXBUFF  0                       /* Receive buffer register */
#define DLLSB   0                       /* Divisor latch LS byte */
#define DLMSB   1                       /* Divisor latch MS byte */
#define IER     1                       /* Interrupt enable register */
#define IIR     2                       /* Interrupt ID register */
#define LCR     3                       /* Line control register */
#define MCR     4                       /* Modem control register */
#define LSR     5                       /* Line status register */
#define MSR     6                       /* Modem status register */

/*      Modem control register bits                                          */

#define DTR     0x01                    /* Data terminal ready */
#define RTS     0x02                    /* Request to send */
#define OUT1    0x04                    /* Output #1 */
#define OUT2    0x08                    /* Output #2 */
#define LPBK    0x10                    /* Loopback mode bit */

/*      Modem status register bits                                           */

#define DCTS    0x01                    /* Delta clear to send */
#define DDSR    0x02                    /* Delta data set ready */
#define TERI    0x04                    /* Trailing edge ring indicator */
#define DRLSD   0x08                    /* Delta Rx line signal detect */
#define CTS     0x10                    /* Clear to send */
#define DSR     0x20                    /* Data set ready */
#define RI      0x40                    /* Ring indicator */
#define RLSD    0x80                    /* Receive line signal detect */

/*      Line control register bits                                           */

#define DATA5   0x00                    /* 5 Data bits */
#define DATA6   0x01                    /* 6 Data bits */
#define DATA7   0x02                    /* 7 Data bits */
#define DATA8   0x03                    /* 8 Data bits */

#define STOP1   0x00                    /* 1 Stop bit */
#define STOP2   0x04                    /* 2 Stop bits */

#define NOPAR   0x00                    /* No parity */
#define ODDPAR  0x08                    /* Odd parity */
#define EVNPAR  0x18                    /* Even parity */
#define STKPAR  0x28                    /* Stick parity */
#define ZROPAR    0x38            /* Zero parity */

/*      Line status register bits                                            */

#define RDR     0x01                    /* Receive data ready */
#define ERRS    0x1E                    /* All the error bits */
#define TXR     0x20                    /* Transmitter ready */

/*      Interrupt enable register bits                                       */

#define DR      0x01                    /* Data ready */
#define THRE    0x02                    /* Tx buffer empty */
#define RLS     0x04                    /* Receive line status */

/*****************************************************************************
 *                   Names for Numbers                 *
 *****************************************************************************/

#define MAX_PORT    4

#define TRUE        1
#define FALSE        0

/*****************************************************************************
 *                  Global Data                     *
 *****************************************************************************/

/*  UART i/o addresses.  Values depend upon which COMM port is selected  */

static int    uart_data;        /* Data register */
static int    uart_ier;        /* Interrupt enable register */
static int    uart_iir;        /* Interrupt identification register */
static int    uart_lcr;        /* Line control register */
static int    uart_mcr;        /* Modem control register */
static int    uart_lsr;        /* Line status register */
static int    uart_msr;        /* Modem status register */

static unsigned char    com_installed;        /* Flag: Communications routines installed */
static int    intnum;            /* Interrupt vector number for chosen port */
static unsigned char    i8259bit;        /* 8259 bit mask */
static unsigned char    old_i8259_mask;        /* Copy as it was when we were called */
static unsigned char    old_ier;        /* Modem register contents saved for */
static unsigned char    old_mcr;        /*  restoring when we're done */
#if 0
static void (_interrupt _far *)(old_vector)( );    /* Place to save COM1 vector */
#else
static void (interrupt *old_vector)( void);
#endif

#ifndef interrupt
#define interrupt   _interrupt _far
#endif
void interrupt    mouse_interrupt_driver( void);

int global_mouse_x, global_mouse_y, global_mouse_button;

static unsigned n_received_chars;

/*****************************************************************************
 *                 com_install( )                     *
 *****************************************************************************
 * DESCRIPTION:    Installs the communications drivers.                   *
 *                                         *                              *
 * SYNOPSIS:    status = com_install(int portnum);                        *
 *        int    portnum;    Desired port number                          *
 *        int    status;        0 = Successful installation               *
 *                    1 = Invalid port number                             *
 *                    2 = No UART for specified port                      *
 *                    3 = Drivers already installed                       *
 *                                                                        *
 * REVISIONS:    18 OCT 89 - RAC - Translated from IBMCOM.PAS             *
 *****************************************************************************/

static const int    uart_base[] =    { 0x3F8, 0x2F8, 0x3E8, 0x2E8 };
static const unsigned char    intnums[] =    { 0x0C,  0x0B,  0x0C,  0x0B };
static const unsigned char    i8259levels[] =    { 4,     3,     4,     3 };

#define DIVISOR ((unsigned)(115200L / 1200L))

int mouse_install( int portnum)
{
   if (com_installed)                /* Drivers already installed */
      return 3;
   if( (portnum < 1) || (portnum > MAX_PORT)) /* Port number out of bounds */
      return 1;

   uart_data = uart_base[portnum-1];        /* Set UART I/O addresses */
   uart_ier  = uart_data + IER;        /*  for the selected comm */
   uart_iir  = uart_data + IIR;        /*  port */
   uart_lcr  = uart_data + LCR;
   uart_mcr  = uart_data + MCR;
   uart_lsr  = uart_data + LSR;
   uart_msr  = uart_data + MSR;
   intnum    = intnums[portnum-1];        /* Ditto for interrupt */
   i8259bit  = (unsigned char)( 1 << i8259levels[portnum-1]);
                             /*  vector and 8259 bit mask */

   old_ier = (unsigned char)inportb( uart_ier);
                                      /* Return an error if we */
   outportb( uart_ier, 0);            /*  can't access the UART */
   if( inportb(uart_ier) != 0)
      return 2;

   disable( );                    /* Save the original 8259 */
   old_i8259_mask = (unsigned char)inportb( 0x21);
                                  /*  mask, then disable the */
   outportb(0x21, old_i8259_mask | i8259bit);    /*  8259 for this interrupt */
   enable( );

   n_received_chars = 0;
   global_mouse_x = global_mouse_y  = global_mouse_button = 0;

   old_vector = getvect( intnum);        /* Save old COMM vector, */
   setvect( intnum, &mouse_interrupt_driver);    /*  then install a new one, */
   com_installed = TRUE;            /*  and note that we did */

   outportb(uart_lcr, DATA8 + NOPAR + STOP1);  /* 8 data, no parity, 1 stop */

   disable( );                    /* Save MCR, then enable */
   old_mcr = (unsigned char)inportb(uart_mcr);
                                     /*  interrupts onto the bus, */
   outportb( uart_mcr,               /*  activate RTS and leave */
        (old_mcr & DTR) | (OUT2 + RTS));    /*  DTR the way it was */
   enable( );

   outportb(uart_ier, DR);            /* Enable receive interrupts */

   disable( );                    /* Now enable the 8259 for */
   outportb(0x21, inportb(0x21) & ~i8259bit);    /*  this interrupt */
   enable( );

                                  /* Now raise DTR */
   disable( );
   outportb(uart_mcr, inportb(uart_mcr) | DTR);
   enable( );

   disable( );                               /* Interrupts off */
   outportb( uart_lcr,                       /* Set up to load baud rate */
            inportb( uart_lcr) | 0x80);      /*  divisor into UART */
   outport( uart_data, DIVISOR);             /* Do so */
   outportb( uart_lcr,                       /* Back to normal UART ops */
            inportb( uart_lcr) & ~0x80);
   enable( );                /* Interrupts back on */

   return 0;                    /* Successful installation */
}                        /* End com_install( ) */

/*****************************************************************************
 *                 com_install( )                     *
 *****************************************************************************
 * DESCRIPTION:    Deinstalls the communications drivers completely, without *
 *        changing the baud rate or DTR.  It tries to leave the              *
 *        interrupt vectors and enables and everything else as they          *
 *        were when the driver was installed.                                *
 *                                                                           *
 * NOTE:    This function MUST be called before returning to DOS, so the     *
 *        interrupt vector won't point to our driver anymore, since it       *
 *        will surely get overwritten by some other transient program        *
 *        eventually.                                                        *
 *                                                                           *
 * REVISIONS:    18 OCT 89 - RAC - Translated from IBMCOM.PAS                *
 *****************************************************************************/

void mouse_deinstall( void)
{
   if (com_installed)              /* Don't de-install twice! */
      {
      outportb( uart_mcr, old_mcr);        /* Restore the UART */
      outportb( uart_ier, old_ier);        /*  registers ... */
      disable( );
      outportb( 0x21,                /*  ... the 8259 interrupt */
          (inportb( 0x21)  & ~i8259bit) | /*  mask ... */
          (old_i8259_mask &  i8259bit));
      enable( );
      setvect( intnum, old_vector);
                                     /*  ... and the comm */
      com_installed = FALSE;            /*  interrupt vector */
      }                    /* End com_installed */
}                        /* End com_deinstall( ) */

/*****************************************************************************
 *                com_interrupt_driver( )                                     *
 *****************************************************************************
 * DESCRIPTION:    Handles communications interrupts.  The UART will         *
 *        interrupt whenever a character has been received or when it is     *
 *        ready to transmit another character.  This routine responds by     *
 *        sticking received characters into the receive queue and            *
 *        yanking characters to be transmitted from the transmit queue       *
 *                                                                           *
 * REVISIOSN:    18 OCT 89 - RAC - Translated from the Pascal.               *
 *****************************************************************************/

void interrupt mouse_interrupt_driver( void)
{
   unsigned char    iir;     /* Local copy if IIR */
   int c;                    /* Local character variable */
   static unsigned char  received_chars[20];

/*  While bit 0 of the IIR is 0, there remains an interrupt to process  */

   while( !((iir = (unsigned char)inportb(uart_iir)) & 1))
      {    /* While there is an int ... */
      switch (iir)
         {                /* Branch on interrupt type */
         case 0:                /* Modem status interrupt */
            inportb(uart_msr);        /* Just clear the interrupt */
            break;

         case 2:                /* Transmit register empty */
            outportb(uart_ier,        /*  Turn off transmit interrupts */
                     inportb(uart_ier) & ~2);
                                /* End 'tx buffer not empty */
            break;

         case 4:                /* Received data interrupt */
            c = inportb( uart_data) & 127;    /* Grab received character */
            if( c & 64)    /* First byte in 3-byte pack has this bit set */
               n_received_chars = 0;
            received_chars[n_received_chars++] = (unsigned char)c;
            if( n_received_chars > 4)
               n_received_chars = 4;
            if( n_received_chars == 3)    /* we have a "for real" one here */
               {
               char dx, dy;

               dx = (char)( received_chars[1] |
                                       (received_chars[0] << 6));
               dy = (char)( received_chars[2] |
                                       ((received_chars[0] & 0xfc) << 4));
               global_mouse_x += dx;
               global_mouse_y += dy;
               global_mouse_button = (received_chars[0] >> 4) & 3;
               }

            break;

         case 6:                /* Line status interrupt */
            inportb(uart_lsr);        /* Just clear the interrupt */
            break;
         }                    /* End switch */
      }                    /* End 'is an interrupt' */
   outportb( 0x20, 0x20);            /* Send EOI to 8259 */
}                        /* End mouse_interrupt_driver( ) */
