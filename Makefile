TARGET  := prospproj pinhole raster3d projmatrix glprojmatrix glorthoprojmatrix camerarays simpleshapes raybox raytri \
           raytracepolymesh \
           raytracetransform \
           shading \
           phong

SRC_FILE := main.cpp

.PHONY: all

all : $(TARGET)

%:%.cpp
	g++ -o $@ $^ -std=c++11 -O3

clean:
	rm -f $(TARGET) *.o *.svg *.ppm .*.swp

    
