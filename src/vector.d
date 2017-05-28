vector.o: src/vector.cc src/vector.hh
	$(CC) $(CXXFLAGS) -c $< -o $@
