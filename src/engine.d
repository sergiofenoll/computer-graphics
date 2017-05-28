engine.o: src/engine.cc src/easy_image.hh src/vector.hh \
 src/ini_configuration.hh src/color.hh src/l_parser.hh src/figures.hh \
 src/draw.hh src/projection.hh
	$(CC) $(CXXFLAGS) -c $< -o $@
