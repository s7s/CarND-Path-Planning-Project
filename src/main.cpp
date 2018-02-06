#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "json.hpp"
#include "spline.h"

using namespace std;

// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.find_first_of("}");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

double distance(double x1, double y1, double x2, double y2)
{
	return sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
}
int ClosestWaypoint(double x, double y, const vector<double> &maps_x, const vector<double> &maps_y)
{

	double closestLen = 100000; //large number
	int closestWaypoint = 0;

	for(int i = 0; i < maps_x.size(); i++)
	{
		double map_x = maps_x[i];
		double map_y = maps_y[i];
		double dist = distance(x,y,map_x,map_y);
		if(dist < closestLen)
		{
			closestLen = dist;
			closestWaypoint = i;
		}

	}

	return closestWaypoint;

}

int NextWaypoint(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{

	int closestWaypoint = ClosestWaypoint(x,y,maps_x,maps_y);

	double map_x = maps_x[closestWaypoint];
	double map_y = maps_y[closestWaypoint];

	double heading = atan2((map_y-y),(map_x-x));

	double angle = fabs(theta-heading);
  angle = min(2*pi() - angle, angle);

  if(angle > pi()/4)
  {
    closestWaypoint++;
  if (closestWaypoint == maps_x.size())
  {
    closestWaypoint = 0;
  }
  }

  return closestWaypoint;
}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> getFrenet(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{
	int next_wp = NextWaypoint(x,y, theta, maps_x,maps_y);

	int prev_wp;
	prev_wp = next_wp-1;
	if(next_wp == 0)
	{
		prev_wp  = maps_x.size()-1;
	}

	double n_x = maps_x[next_wp]-maps_x[prev_wp];
	double n_y = maps_y[next_wp]-maps_y[prev_wp];
	double x_x = x - maps_x[prev_wp];
	double x_y = y - maps_y[prev_wp];

	// find the projection of x onto n
	double proj_norm = (x_x*n_x+x_y*n_y)/(n_x*n_x+n_y*n_y);
	double proj_x = proj_norm*n_x;
	double proj_y = proj_norm*n_y;

	double frenet_d = distance(x_x,x_y,proj_x,proj_y);

	//see if d value is positive or negative by comparing it to a center point

	double center_x = 1000-maps_x[prev_wp];
	double center_y = 2000-maps_y[prev_wp];
	double centerToPos = distance(center_x,center_y,x_x,x_y);
	double centerToRef = distance(center_x,center_y,proj_x,proj_y);

	if(centerToPos <= centerToRef)
	{
		frenet_d *= -1;
	}

	// calculate s value
	double frenet_s = 0;
	for(int i = 0; i < prev_wp; i++)
	{
		frenet_s += distance(maps_x[i],maps_y[i],maps_x[i+1],maps_y[i+1]);
	}

	frenet_s += distance(0,0,proj_x,proj_y);

	return {frenet_s,frenet_d};

}

// Transform from Frenet s,d coordinates to Cartesian x,y
vector<double> getXY(double s, double d, const vector<double> &maps_s, const vector<double> &maps_x, const vector<double> &maps_y)
{
	int prev_wp = -1;

	while(s > maps_s[prev_wp+1] && (prev_wp < (int)(maps_s.size()-1) ))
	{
		prev_wp++;
	}

	int wp2 = (prev_wp+1)%maps_x.size();

	double heading = atan2((maps_y[wp2]-maps_y[prev_wp]),(maps_x[wp2]-maps_x[prev_wp]));
	// the x,y,s along the segment
	double seg_s = (s-maps_s[prev_wp]);

	double seg_x = maps_x[prev_wp]+seg_s*cos(heading);
	double seg_y = maps_y[prev_wp]+seg_s*sin(heading);

	double perp_heading = heading-pi()/2;

	double x = seg_x + d*cos(perp_heading);
	double y = seg_y + d*sin(perp_heading);

	return {x,y};

}


// this function gets all possible successors for current states
// KL:keep lane, PLC: prepare lane change, LCR: lane change right, LCL: lane change left
vector<string> successor_states(string state,int lane, int &counter){
	vector<string> successors;

	// counter to make the state at least continue to 100 counts (PLC&KL states)
	int count_to_change_state=20;
	if ((counter<count_to_change_state)&&(state=="KL"||state=="PLC")){
		successors.push_back(state);
		return successors;
	}

	if (state=="KL"){
		successors.push_back("KL");
		successors.push_back("PLC");			
	}
	if (state=="PLC"){
		successors.push_back("KL");
		successors.push_back("PLC");
		if(lane!=2){
			successors.push_back("LCR");
		}
		if(lane!=0){
			successors.push_back("LCL");
		}
	}
	if (state=="LCR"){
		successors.push_back("KL");
	}
	if (state=="LCL"){
		successors.push_back("KL");
	}
	counter=0; // reset the counter after the 100 counts 
	return successors;
}

// this function returns the cost of any successor state passed
double get_cost (int lane, int next_state_lane, double ref_vel, double next_state_vel, vector<vector<double>> sensor_fusion, double start_path_s, double end_path_s, int path_num_points){
	double cost=0;
	double lane_cost_weight=100;
	double velocity_cost_weight=10;
	double collision_cost_weight=10000000000;
	double epscilon=0.0000001; // to ensute not to devide by zero
	//change lane cost
	double lane_cost=lane_cost_weight*+abs(lane-next_state_lane);
	// velocity diffrence from the refrence cost
	double velocity_cost=velocity_cost_weight*(abs(next_state_vel- ref_vel)/ref_vel);
	// collision cost
	double collisions=0;
	double danger_dist=35;
	double start_danger_dist=20;
	//loop over all the cars
	for(int i=0;i<sensor_fusion.size();i++){
		double d= sensor_fusion[i][6];


		// check if this car in the state lane
		if(d<((next_state_lane*4)+4)&&d>((next_state_lane*4))){
			double vx=sensor_fusion[i][3];
			double vy=sensor_fusion[i][4];
			double v=sqrt(pow(vx,2)+pow(vy,2)); //calculate the velocity
			double start_s=sensor_fusion[i][5];
			double end_s=start_s+path_num_points*0.02*v; // calculate the s value at the end of my path
			// if the lane will change ,check the the cars near of my car at the start of the path
			if(next_state_lane!=lane){
				if(abs(start_s- start_path_s)<start_danger_dist){
					collisions+=1/(abs(start_s- start_path_s)+epscilon);
					collisions+=1000*(next_state_vel)/ref_vel;
					collisions+=1/(abs(end_s- end_path_s)+epscilon);
					collisions+=1000*(next_state_vel)/ref_vel;
				}
				if((start_s<start_path_s)&&(end_s>end_path_s)){
					collisions+=100;
				}
				if((start_s>start_path_s)&&(end_s<end_path_s)){
					collisions+=100;
				}
			}
			if((abs(end_s- end_path_s)<danger_dist)&&(start_s>start_path_s)){
				collisions+=1/(abs(end_s- end_path_s)+epscilon);
				collisions+=1000*(next_state_vel)/ref_vel;
				break; // if the collision will happen once ,it's enogh
			}
		}		
	}
	double collision_cost=collision_cost_weight*collisions;

	cost=lane_cost+velocity_cost+collision_cost;
	return cost;
}


int main() {
  uWS::Hub h;

  // Load up map values for waypoint's x,y,s and d normalized normal vectors
  vector<double> map_waypoints_x;
  vector<double> map_waypoints_y;
  vector<double> map_waypoints_s;
  vector<double> map_waypoints_dx;
  vector<double> map_waypoints_dy;

  // Waypoint map to read from
  string map_file_ = "../data/highway_map.csv";
  // The max s value before wrapping around the track back to 0
  double max_s = 6945.554;

  ifstream in_map_(map_file_.c_str(), ifstream::in);

  string line;
  while (getline(in_map_, line)) {
  	istringstream iss(line);
  	double x;
  	double y;
  	float s;
  	float d_x;
  	float d_y;
  	iss >> x;
  	iss >> y;
  	iss >> s;
  	iss >> d_x;
  	iss >> d_y;
  	map_waypoints_x.push_back(x);
  	map_waypoints_y.push_back(y);
  	map_waypoints_s.push_back(s);
  	map_waypoints_dx.push_back(d_x);
  	map_waypoints_dy.push_back(d_y);
  }

  int lane =1;
  double scan_velocity=0;
  //KL:keep lane, PLC: prepare lane change, LCR: lane change right, LCL: lane change left
  vector<string> states={"KL","PLC","LCR","LCL"};
  string state=states[0];
  double ref_vel=49;
  double low_vel=30;
  int counter=0;
  h.onMessage([&counter,&scan_velocity,&low_vel,&states,&state,&lane,&ref_vel,&map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                     uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    //auto sdata = string(data).substr(0, length);
    //cout << sdata << endl;
    if (length && length > 2 && data[0] == '4' && data[1] == '2') {

      auto s = hasData(data);

      if (s != "") {
        auto j = json::parse(s);
        
        string event = j[0].get<string>();
        
        if (event == "telemetry") {
          // j[1] is the data JSON object
          
        	// Main car's localization Data
          	double car_x = j[1]["x"];
          	double car_y = j[1]["y"];
          	double car_s = j[1]["s"];
          	double car_d = j[1]["d"];
          	double car_yaw = j[1]["yaw"];
          	double car_speed = j[1]["speed"];

          	// Previous path data given to the Planner
          	auto previous_path_x = j[1]["previous_path_x"];
          	auto previous_path_y = j[1]["previous_path_y"];
          	// Previous path's end s and d values 
          	double end_path_s = j[1]["end_path_s"];
          	double end_path_d = j[1]["end_path_d"];

          	// Sensor Fusion Data, a list of all other cars on the same side of the road.
          	auto sensor_fusion = j[1]["sensor_fusion"];

          	json msgJson;

          	vector<double> next_x_vals;
          	vector<double> next_y_vals;


          	//define a path made up of (x,y) points that the car will visit sequentially every .02 seconds
        	  double prev_path_size=previous_path_x.size();

        	  //get successors of current state
        	  vector<string>successors= successor_states(state, lane, counter);
        	  int next_state_lane;
        	  double next_state_vel;
        	  vector<vector<double>> x_paths;
        	  vector<vector<double>> y_paths;
        	  vector<double> costs;
						string next_state;
						int path_num_points=50;

						// loop on all successors to check which has the minimum cost
						for(int j=0;j<successors.size();j++){
							vector<double> current_x_path;
							vector<double> current_y_path;
							double velocity =scan_velocity;
							next_state=successors[j];

							if (next_state=="KL"){
								next_state_lane=lane;
								next_state_vel=ref_vel;
							}
							if (next_state=="PLC"){
								next_state_lane=lane;
								next_state_vel=low_vel;
							}
							if (next_state=="LCR"){
								next_state_lane=lane+1;
								next_state_vel=ref_vel;
							}
							if (next_state=="LCL"){
								next_state_lane=lane-1;
								next_state_vel=ref_vel;
							}

							// to change the velocity gradually
							if(velocity<next_state_vel){
								velocity+=0.35;
							}
							else if (velocity>next_state_vel)
							{
								velocity-=0.35;
							}

	        	  //get the five spline points
	        	  //first 2 points are the last points
	        	  //second 3 points are generated
	        	  vector<double> x_pts;
	        	  vector<double> y_pts;
	        	  double ref_x;
	        	  double ref_y;
	        	  double ref_yaw;
	        	  double prev_ref_x;
	        	  double prev_ref_y;

	        	  if(prev_path_size<2){
	        	  	ref_x=car_x;
	        	  	ref_y=car_y;
	        	  	ref_yaw=car_yaw;
	        	  	prev_ref_x=car_x-cos(car_yaw);
	        	  	prev_ref_y=car_y-sin(car_yaw);
	        	  	x_pts.push_back(prev_ref_x);
	        	  	x_pts.push_back(ref_x);
	        	  	y_pts.push_back(prev_ref_y);
	        	  	y_pts.push_back(ref_y);
	        	  }
	        	  else{
	        	  	ref_x=previous_path_x[prev_path_size-1];
	        	  	ref_y=previous_path_y[prev_path_size-1];
	        	  	prev_ref_x=previous_path_x[prev_path_size-2];;
	        	  	prev_ref_y=previous_path_y[prev_path_size-2];;
	        	  	ref_yaw=atan2(ref_y - prev_ref_y,ref_x - prev_ref_x);
	        	  	x_pts.push_back(prev_ref_x);
	        	  	x_pts.push_back(ref_x);
	        	  	y_pts.push_back(prev_ref_y);
	        	  	y_pts.push_back(ref_y);
	        	  }

	        	  double dist_inc=45;
						  for(int i = 0; i < 3; i++)
						    {
						    	double next_s= car_s+(i+1)*dist_inc;
						    	vector<double> xy=getXY(next_s, next_state_lane*4+2, map_waypoints_s, map_waypoints_x, map_waypoints_y);
				          x_pts.push_back(xy[0]);
				          y_pts.push_back(xy[1]);
						    }


						 	for(int i=0; i<x_pts.size();i++){
						 		double shift_x=x_pts[i]-ref_x;
						 		double shift_y=y_pts[i]-ref_y;
						 		x_pts[i]=shift_x*cos(-ref_yaw)-shift_y*sin(-ref_yaw);
						 		y_pts[i]=shift_x*sin(-ref_yaw)+shift_y*cos(-ref_yaw);
						 	}

						 	//generate spline
						 	tk::spline s;
						 	s.set_points(x_pts,y_pts);

						 	//generate path points
						 	//first pack of points are the previous path points
						 	//second pack of points are the spline points
						 	for(int i=0; i<prev_path_size;i++){
						 		current_x_path.push_back(previous_path_x[i]);
						 		current_y_path.push_back(previous_path_y[i]);
						 	}

						 	double target_x=30;
						 	double target_y=s(target_x);
						 	double target_dist=sqrt(pow(target_x,2)+pow(target_y,2));
						 	double N=target_dist/(0.02*velocity/2.24);

						 	for(int i=0; i<path_num_points-prev_path_size;i++){
						 		double next_x= (i+1)*(target_x/N);
						 		double next_y=s(next_x);
						 		double point_x=next_x;
						 		double point_y=next_y;
						 		next_x=point_x*cos(ref_yaw)-point_y*sin(ref_yaw)+ref_x;
						 		next_y=point_x*sin(ref_yaw)+point_y*cos(ref_yaw)+ref_y;
						 		current_x_path.push_back(next_x);
						 		current_y_path.push_back(next_y);
						 	}
	        	  x_paths.push_back(current_x_path);
	        	  y_paths.push_back(current_y_path);
	        	  
	        	  double start_yaw = atan2(current_y_path[1] - current_y_path[0], current_x_path[1] - current_x_path[0]);
	        	  double end_yaw = atan2(current_y_path[path_num_points-1] - current_y_path[path_num_points-2], current_x_path[path_num_points-1] - current_x_path[path_num_points-2]);

							vector<double>sd_start= getFrenet(current_x_path[0], current_y_path[0], start_yaw, map_waypoints_x, map_waypoints_y);
							vector<double>sd_end= getFrenet(current_x_path[path_num_points-1], current_y_path[path_num_points-1], end_yaw, map_waypoints_x, map_waypoints_y);

						 	double start_path_s=sd_start[0];
						 	double end_path_s=sd_end[0];

						 	// append the cost of the current successor
						 	costs.push_back(get_cost(lane, next_state_lane, ref_vel, next_state_vel, sensor_fusion, start_path_s, end_path_s,path_num_points));
						}	

						// calculate which successor has the minimum cost
						double min_cost=costs[0];
						int min_cost_idx=0;
						for(int i=0; i<costs.size();i++){
				 			//cout<<i<<"\t"<<costs[i]<<endl;
							if (costs[i]<min_cost){
								min_cost_idx=i;
							}
						}
						//get the path of the minimum cost successor
				 		next_x_vals=x_paths[min_cost_idx];
				 		next_y_vals=y_paths[min_cost_idx];

						next_state=successors[min_cost_idx];
						state=next_state;

						//apply the action of the minimum cost successor
						if (next_state=="KL"){
							lane=lane;
							next_state_vel=ref_vel;
						}
						if (next_state=="PLC"){
							lane=lane;
							next_state_vel=low_vel;
						}
						if (next_state=="LCR"){
							lane++;
							next_state_vel=ref_vel;
						}
						if (next_state=="LCL"){
							lane--;
							next_state_vel=ref_vel;
						}

						if(scan_velocity<next_state_vel){
							scan_velocity+=0.35;
						}
						else if (scan_velocity>next_state_vel)
						{
							scan_velocity-=0.35;
						}

						//increase the change state counter
						counter++;

          	//end

          	msgJson["next_x"] = next_x_vals;
          	msgJson["next_y"] = next_y_vals;

          	auto msg = "42[\"control\","+ msgJson.dump()+"]";

          	//this_thread::sleep_for(chrono::milliseconds(1000));
          	ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
          
        }
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}
