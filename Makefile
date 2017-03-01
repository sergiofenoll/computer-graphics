CC = g++
CFLAGS = -std=c++11 -Wall -Wextra -g3
LFLAGS =
EXEC = engine
CPP_FILES = $(wildcard src/*.cc)
HED_FILES = $(OBJ_FILES:.o=.d)
OBJ_FILES = $(addprefix obj/,$(notdir $(CPP_FILES:.cc=.o)))

$(EXEC) : $(OBJ_FILES)
	$(CC) $(LFLAGS) -o $@ $^

-include $(HED_FILES)

obj/%.o : src/%.cc
	$(CC) $(CFLAGS) -MM -MT $@ -MF $(patsubst %.o,%.d,$@)  $<
	$(CC) $(CFLAGS) -c -o $@ $<

clean:
	find obj/ -name "*.o" -delete
	find obj/ -name "*.d" -delete
	find . -name $(EXEC) -delete

targz:
	tar -czf $(EXEC).tar.gz src/ obj/ Makefile
zip:
	zip $(EXEC).zip src/ obj/ Makefile
