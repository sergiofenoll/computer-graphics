l_parser.o: src/l_parser.cc src/l_parser.hh
	$(CC) $(CXXFLAGS) -c $< -o $@
