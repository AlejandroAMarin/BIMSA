###############################################################################
# Definitions
###############################################################################
CC=gcc
CC_FLAGS=-Wall -g -O3

###############################################################################
# Rules
###############################################################################
all: wfa.edit.alignment wfa.edit.distance wfa.edit.distance.nopeeling

# Building rules
wfa.edit.alignment: wfa.edit.alignment.c
	$(CC) $(CC_FLAGS) $< -o $@
wfa.edit.distance: wfa.edit.distance.c
	$(CC) $(CC_FLAGS) $< -o $@
wfa.edit.distance.nopeeling: wfa.edit.distance.nopeeling.c
	$(CC) $(CC_FLAGS) $< -o $@
