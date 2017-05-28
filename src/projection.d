projection.o: src/projection.cc src/projection.hh src/vector.hh \
 src/figures.hh src/easy_image.hh src/l_parser.hh src/color.hh \
 src/ini_configuration.hh
	$(CC) $(CXXFLAGS) -c $< -o $@
