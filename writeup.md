**Path Planning**
---

**Path Planning Project**

The goals of this project are the following:

* Design a path planner that is able to create smooth, safe paths for the car to follow along a 3 lane highway with traffic
* Run the implemented algorithm on the simulator
* Summarize the results with a written report



## Project Basics
In this project, I used C++ to write a path planning algorithm keep the vehicle inside its lane, avoid hitting other cars, and pass slower moving traffic all by using localization, sensor fusion, and map data.

See files in the `src` folder for the primary C++ files making up this project.
---

### The code in steps
All the code in one file **main.cpp**. I haven't added new code files, however I add 2 functions:
* **successor_states:** this function gets all possible successors for current states *"line 170"*
* **get_cost:** this function returns the cost of any successor state passed *"line 205"*

Here is the algorithm code in steps:

1. I initialized the values of the starting  velocity, lane and state *"lane 297:301"*

2. I get all the successors of all current state using `successor_states` function and iterate on them and use the counter to return the same, "KL" or "PLC", state for at least 20 counts to make the transition between states smoother  *"lane 352:362"*

3. I get the lane and velocity of the successors *"lane 368:383"*

4. I generate 5 points to represent my path, first 2 points are the last 2 points of the previous path and second 3 points are generated according to lane and velocity, then I use `spline library` to curve fit these 5 points *"lane 397:447"*

5. I generated the path points. the first part of points I take them from the remaining points from the previous path and the other point I generated them using the spline function *"lane 452:473"*

6. I get the cost of the successors using `get_cost` function and I divided the cost to 3 parts; lane change cost, velocity from reference cost and collision cost *"lane 485"*

7. I get the successors with the minimum cost and change the
lane and velocity according to it and pass its path to the simulator  *"lane 489:528"*

### Results
The path planner successfully is able to keep inside its lane, avoid hitting other cars, and pass slower moving traffic.

And here is a link to a video showing the algorithm working:
[Path planning project video](https://youtu.be/6C9_n_T0MBU)
