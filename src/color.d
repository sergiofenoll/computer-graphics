color.o: src/color.cc src/color.hh src/easy_image.hh src/vector.hh \
 src/figures.hh src/l_parser.hh src/ini_configuration.hh
	$(CC) $(CXXFLAGS) -c $< -o $@
