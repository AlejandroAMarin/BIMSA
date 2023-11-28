###############################################################################
# Rules
###############################################################################
all: generate_dataset

clean:
	rm generate_dataset

generate_dataset: generate_dataset.c
	$(CC) $(CC_FLAGS) generate_dataset.c -o generate_dataset $(LD_FLAGS) -lm
