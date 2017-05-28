easy_image.o: src/easy_image.cc src/easy_image.hh src/vector.hh
	$(CC) $(CXXFLAGS) -c $< -o $@
