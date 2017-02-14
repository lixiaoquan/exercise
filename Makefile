TARGET   := render

SRC_FILE := main.cpp

.PHONY: all

all : $(TARGET) prospproj

$(TARGET) : $(SRC_FILE)
	g++ -o $@ $^

prospproj : prospproj.o 
	g++ $^ -o $@ 

%.o:%.cpp
	g++ -c -o $@ $^ -std=c++11

clean:
	rm -f $(TARGET) *.o prospproj

    
