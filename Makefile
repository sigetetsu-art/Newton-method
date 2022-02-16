CC = g++
CFLAG = -O4 -Wall

all: newton_method

newton_method: Newton_method.cpp
						$(CC) $(CFLAG) Newton_method.cpp -o newton_method

clean:
		rm -f newton_method *.o*~