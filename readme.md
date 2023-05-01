# Assignment 2
## Libraries Required
* `OpenMP`
* `OpenGL`
## This folder Contains
*  [Source File: many-body-sim-program.cpp](many-body-sim-program.cpp)
*  [Helper File: first.py](first.py)

*  [Source File: graphics-program.cpp](graphics-program.cpp)
*  [Helper File: glad.c](glad.c)
*  [Result: scalability.txt](scalability.txt)
*  [Result: simulation_log.txt](scalability.txt)
## Execution Process
>Make sure that the input file is in same directory as the source file 

Firstly run the helper python script to get numerical values out of the input file
```
$ python first.py  
```
Then you get the numerical data in `inp.txt`

Now, compile & run the many-body-sim.cpp to generate the `trajectory.txt` file using:
```
$ g++ many-body-sim-program.cpp -o a -fopenmp && ./a
```
## Example Result
Time taken to exceute is shown in the terminal like this:
```
Number of threads : 8   
Number of steps : 2400
Time taken (parallel part) for 8 threads : 24.37128
Average Time taken for each step threads : 0.0101547
```
## For Graphic Visualisation
>Make sure OpenGL is installed and glad.c is in the same directory as the source file.

>Note: For guidelines for installing OpenGL click [here](https://medium.com/geekculture/a-beginners-guide-to-setup-opengl-in-linux-debian-2bfe02ccd1e).


Now, compile the `graphics-program.cpp` to visualise the `trajectory.txt` using the following command
```
$ g++ graphics-program.cpp glad.c -lGL -lGLU -lglut -ldl -lglfw
```
To run the program simply execute the compiled file by
```
$ ./a.out
```
