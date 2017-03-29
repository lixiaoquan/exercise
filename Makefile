TARGET  := prospproj pinhole raster3d projmatrix glprojmatrix glorthoprojmatrix camerarays simpleshapes raybox raytri \
           raytracepolymesh \
           raytracetransform \
           shading \
           phong \
           mcsim \
           mcintegration

SRC_FILE := main.cpp

.PHONY: all

all : $(TARGET)

%:%.cpp
	g++ -o $@ $^ -std=gnu++11 -O3

clean:
	rm -f $(TARGET) *.o *.svg *.ppm .*.swp

    
