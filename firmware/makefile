# ChipWhisperer firmware makefile for the OSIDH-LD glitch target.
# Place this folder under chipwhisperer/firmware/mcu/ (e.g. .../mcu/osidh-glitch/)
# then build with:
#   make PLATFORM=CWLITEARM CRYPTO_TARGET=NONE SS_VER=SS_VER_1_1
TARGET = target_osidh
SRC += target_osidh.c

# no crypto primitive is needed on the target
CRYPTO_TARGET = NONE

# SimpleSerial v1.1 (target_osidh.c uses the (buf,len) handler signature)
SS_VER ?= SS_VER_1_1

include ../simpleserial/Makefile.simpleserial
FIRMWAREPATH = ../.
include $(FIRMWAREPATH)/Makefile.inc
