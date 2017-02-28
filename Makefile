TARGET  := prospproj pinhole raster3d projmatrix glprojmatrix glorthoprojmatrix camerarays

SRC_FILE := main.cpp

.PHONY: all

all : $(TARGET)

%:%.cpp
	g++ -o $@ $^ -std=c++11

clean:
	rm -f $(TARGET) *.o *.svg *.ppm .*.swp

    
