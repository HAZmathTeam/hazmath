#
TARGET   = ingrid
CCOMPILE = gcc
CFLAGS   = -g -I. -I../../../include -L../../../lib
COBJS = main_input.o

$(TARGET): $(COBJS)
	$(CCOMPILE) -C $(CFLAGS) -o $(TARGET) $(COBJS) -lhazmath -lm
%.o:	%.c
	$(CCOMPILE) $(CFLAGS) -o $@ -c $<
clean:
	+rm -rf $(TARGET) *.o *~ *.mod ./*.dSYM 

