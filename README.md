# double-pendulum
This C++ script models the behaviour of a double-pendulum which is a chaotic system.  

The script takes a text file as input called parameters.txt  
This text file contains the following values:
* The mass (kg) of the first pendulum
* The mass (kg) of the second pendulum
* The length (m) of the first pendulum's rod
* The length (m) of the second pendulum's rod
* The angle (rad) of the first pendulum  
* The angle (rad) of the second pendulum
* The magnitude of the time step (s) that will be used
* The time (s) at which the simulation should stop

The system is integrated using a fourth-order Runge-Kutta numerical scheme.  
For each time step the current time as well as the x, y co-ordinates of each of the masses are written to a text file called output.txt