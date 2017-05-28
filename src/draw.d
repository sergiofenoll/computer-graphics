draw.o: src/draw.cc src/draw.hh src/easy_image.hh src/vector.hh \
 src/color.hh src/figures.hh src/l_parser.hh src/ini_configuration.hh \
 src/projection.hh
	$(CC) $(CXXFLAGS) -c $< -o $@
