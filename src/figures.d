figures.o: src/figures.cc src/figures.hh src/easy_image.hh src/vector.hh \
 src/l_parser.hh src/color.hh src/ini_configuration.hh
	$(CC) $(CXXFLAGS) -c $< -o $@
