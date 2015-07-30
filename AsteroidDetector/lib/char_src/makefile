charon.exe: charon.obj a10_extr.obj blink.obj colors.obj findstar.obj &
   gsc_act.obj isophote.obj load_img.obj matcher.obj miscell.obj &
   mouse.obj grab_gsc.obj regions3.obj uncrunch.obj untycho.obj &
   lsquare.obj find_obj.obj settings.obj ucac2.obj ucac3.obj ucac4.obj &
   dispimg.obj \checkmem\watcheck.obj \utils\wbmouse.obj
   wlink name charon.exe option map=charon.map @charon.lnk >> err
   type err

CFLAGS=-j -w4 -ox -dDEBUG_MEM -I\myincl -I\gsc13 -dPIXEL_32 -4r -s -zq

.cpp.obj:
   wcc386 $(CFLAGS) $< >> err
   type err

a10_extr.obj: a10_extr.cpp

charon.obj: charon.cpp

settings.obj: settings.cpp

ucac2.obj: \ucac2\ucac2.cpp
   wcc386 $(CFLAGS) -I\ucac2 \ucac2\ucac2.cpp >> err

ucac3.obj: \ucac3\ucac3.cpp
   wcc386 $(CFLAGS) -I\ucac3 \ucac3\ucac3.cpp >> err

ucac4.obj: \ucac4\ucac4.c
   wcc386 $(CFLAGS) -I\ucac4 \ucac4\ucac4.c   >> err

blink.obj: blink.cpp

colors.obj: \lunar\posttest\colors.cpp
   wcc386 $(CFLAGS) \lunar\posttest\colors.cpp >> err

findstar.obj: findstar.cpp

find_obj.obj: find_obj.cpp

gsc_act.obj: \gsc13\gsc_act.cpp
   wcc386 $(CFLAGS) \gsc13\gsc_act.cpp >> err

matcher.obj: matcher.cpp findstar.h \myincl\lsquare.h
   wcc386 $(CFLAGS) matcher.cpp >> err
   type err

isophote.obj: isophote.cpp

miscell.obj: miscell.cpp

load_img.obj: load_img.cpp

lsquare.obj: lsquare.cpp

dispimg.obj: dispimg.cpp

grab_gsc.obj: grab_gsc.cpp findstar.h
   wcc386 $(CFLAGS) -I\gsc\code grab_gsc.cpp >> err
   type err

regions3.obj: \gsc\code\regions3.cpp
   wcc386 $(CFLAGS) -I\gsc\code \gsc\code\regions3.cpp >> err
   type err

uncrunch.obj: \gsc\code\uncrunch.cpp
   wcc386 $(CFLAGS) -I\gsc\code \gsc\code\uncrunch.cpp >> err
   type err

untycho.obj: \gsc\code\untycho.cpp
   wcc386 $(CFLAGS) -I\gsc\code \gsc\code\untycho.cpp >> err
   type err

mouse.obj: \ibmcom\mouse.cpp
   wcc386 $(CFLAGS) \ibmcom\mouse.cpp >> err

\checkmem\watcheck.obj: \checkmem\checkmem.cpp
   wcc386 -DDOS -DDEBUG_MEM -ox -4r -s -w4 -I\myincl -Fo\checkmem\watcheck.obj \checkmem\checkmem.cpp

\utils\wbmouse.obj: \utils\bmouse.cpp
   wcc386 -DDOS -DDEBUG_MEM -ox -4r -s -w4 -I\myincl -Fo\utils\wbmouse.obj \utils\bmouse.cpp
