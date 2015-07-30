/*
Module : AAMOONMAXDECLINATIONS.H
Purpose: Implementation for the algorithms which obtain the dates and values for Maximum declination of the Moon
Created: PJN / 13-01-2004

Copyright (c) 2004 - 2013 by PJ Naughter (Web: www.naughter.com, Email: pjna@naughter.com)

All rights reserved.

Copyright / Usage Details:

You are allowed to include the source code in any product (commercial, shareware, freeware or otherwise) 
when your product is released in binary form. You are allowed to modify the source code in any way you want 
except you cannot modify the copyright details at the top of each module. If you want to distribute source 
code with your application, then you are only allowed to distribute versions released by the author. This is 
to maintain a single distribution point for the source code. 

*/


/////////////////////// Macros / Defines //////////////////////////////////////

#if _MSC_VER > 1000
#pragma once
#endif

#ifndef __AAMOONMAXDECLINATIONS_H__
#define __AAMOONMAXDECLINATIONS_H__

#ifndef AAPLUS_EXT_CLASS
#define AAPLUS_EXT_CLASS
#endif


////////////////////// Classes ////////////////////////////////////////////////

class AAPLUS_EXT_CLASS CAAMoonMaxDeclinations
{
public:
//Static methods
  static double K(double Year);
  static double MeanGreatestDeclination(double k, bool bNortherly);
  static double MeanGreatestDeclinationValue(double k);
  static double TrueGreatestDeclination(double k, bool bNortherly);
  static double TrueGreatestDeclinationValue(double k, bool bNortherly);
};

#endif //__AAMOONMAXDECLINATIONS_H__
