ini_configuration.o: src/ini_configuration.cc src/ini_configuration.hh \
 src/color.hh src/easy_image.hh src/vector.hh
	$(CC) $(CXXFLAGS) -c $< -o $@
