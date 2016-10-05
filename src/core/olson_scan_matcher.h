#include "sensor_data.h"
#include "state_data.h"
#include "geometry_utils.h"
#include "maps/grid_cell_strategy.h"
#include "maps/grid_map.h"
//#include "tiny_grid_cells.h"

#include <iostream>

struct Point{
  double x,y;
  Point(double x, double y) {
    this->x = x; this->y = y;
  }
};

class Olson_scan_matcher : public GridScanMatcher{
public:
  double radius; // TODO replace radius & sigma in one param
  const double sigma;
  double meters_per_cell;
  std::shared_ptr<GridMap> table;
  TransformedLaserScan prev_scan;
  const std::shared_ptr<GridCellFactory> factory = std::shared_ptr<GridCellFactory>(new TinyBaseCellFactory());
  RobotState prev_pose;
  RobotState prev_init_pose;

public:
  Olson_scan_matcher(double _radius, double _sigma, double _meters_per_cell) :
                     GridScanMatcher(nullptr), radius(_radius), sigma(_sigma),
                     meters_per_cell(_meters_per_cell*5.0), table(nullptr) {
    radius = meters_per_cell*2.01;
  }

  double normal_variation(double dx, double dy, double dz = 0, double sig = 1.0) {
    return std::exp(-0.5 * (QUAD(dx) + QUAD(dy) + QUAD(dz))/QUAD(sig));
  }

  void blur_probability(Point point) {
    if(table == nullptr) {
      return;
    }

    double step = table->scale()*0.1;
    for (double i = -radius; i <= radius; i += step) {
      double range_j = std::sqrt(QUAD(radius) - QUAD(i));
      for (double j = -range_j; j <= range_j; j += step) {
        double prev_value = table->cell_value(DiscretePoint2D(i + point.x, j+point.y));
        double value = prev_value + normal_variation(i,j,0.0, sigma)/4.0;
        value = value > 1.0 ? 1.0 : value;
        table->update_cell(table->world_to_cell(i + point.x, j + point.y), Occupancy(value, 1.0), 1.0);
      }
    }
  }



  std::vector<Point> scan_to_cartezian(const RobotState& pose,
                                        const TransformedLaserScan& scan,
                                        double &x_min, double &x_max,
                                        double &y_min, double &y_max) {
    std::vector<Point> result;
    x_min = 10000;
    y_min = 10000;
    x_max = -10000;
    y_max = -10000;
    int i = 0;
    for (auto point : scan.points) {
      i++;
      if(!point.is_occupied) {
        continue;
      }
      if(i%20) continue;
      double x_world = pose.x + point.range * std::cos(point.angle+pose.theta);
      x_max = x_max < x_world ? x_world : x_max;
      x_min = x_world < x_min ? x_world : x_min;

      double y_world = pose.y + point.range * std::sin(point.angle+pose.theta);
      y_max = y_max < y_world ? y_world : y_max;
      y_min = y_world < y_min ? y_world : y_min;

      result.push_back(Point(x_world, y_world));
    }
    return result;
  }

  void generate_table(const RobotState& pose,
                      const TransformedLaserScan& scan) {
    GridMapParams table_params;
    table_params.meters_per_cell = meters_per_cell;
    double x_min, y_min, x_max, y_max;
    std::vector<Point> current_cartezian_scan =
                    scan_to_cartezian(pose, scan, x_min, x_max, y_min, y_max);

    table_params.height = y_max - y_min + 2.0 * radius;
    table_params.width  = x_max - x_min + 2.0 * radius;
    ////std::cout << "table height: "<<table_params.height << "table width: "<<table_params.width << std::endl;

    table = std::shared_ptr<GridMap>(new GridMap(factory, table_params));

    for (auto point : current_cartezian_scan) {
      blur_probability(point);
    }
    //std::cout << "Table (" << table->width() << "x" << table->height() << ") was generated" << std::endl;
  }

  void set_scale(double meters_per_cell) {
    this->meters_per_cell = meters_per_cell;
  }

  double esimate_scan_prob_old(const TransformedLaserScan& scan, double x, double y, double theta) {
    double a,b,c,d;
    std::vector<Point> current_cartezian_scan =
                        scan_to_cartezian(RobotState(x,y,theta), scan,a,b,c,d);
    double result = 0;
    for (auto point : current_cartezian_scan) {
      DiscretePoint2D cell_coord = table->world_to_cell(point.x, point.y);
      double point_prob = table->cell_value(cell_coord);
      result += point_prob;
    }
    return result/current_cartezian_scan.size();
  }

  double esimate_scan_prob(const TransformedLaserScan& scan, const RobotState &pose, const GridMap& map){

    int i =0;
    double result = 0;
    for (auto point : scan.points) {
         i++;
         if(!point.is_occupied) {
           continue;
         }
         if(i%20) continue;
         double x_world = pose.x + point.range * std::cos(point.angle+pose.theta);
         //x_max = x_max < x_world ? x_world : x_max;
         //x_min = x_world < x_min ? x_world : x_min;

         double y_world = pose.y + point.range * std::sin(point.angle+pose.theta);
         //y_max = y_max < y_world ? y_world : y_max;
         //y_min = y_world < y_min ? y_world : y_min;

         DiscretePoint2D cell_coord = map.world_to_cell(x_world, y_world);
         double point_prob = map.cell_value(cell_coord);
         result += std::log(point_prob);
       }
    result /= i;
    return result;

  }

  double process_scan(const RobotState &init_pose,
                      const TransformedLaserScan &scan,
                      const GridMap &map,
                      RobotState &pose_delta) override {
    // if it is a first call of this function;
    //table = std::shared_ptr<GridMap>(&map);
    if(prev_scan.points.empty()) {
      prev_scan = scan;
      //pose_delta = init_pose;
      prev_pose = init_pose;
      prev_init_pose = init_pose;
      return 0;
    }
    //std::cout << "!" << std::endl;
    double max_prob = -1.0;
    double mdx = 0,mdy = 0,mdt = 0;
    GridMapParams params;

    double odom_x =  scan.d_x;
    double odom_y = scan.d_y;
    double odom_theta = scan.d_yaw;

    RobotState prob_pose(prev_pose.x + odom_x, prev_pose.y + odom_y, prev_pose.theta + odom_theta);
    generate_table(prev_init_pose, prev_scan);

    //std::cout << "odom_x: " << odom_x << " odom_y: " << odom_y << "odom_t: " << odom_theta << std::endl;
    int i = 0;
    for (double dx = -std::abs(odom_x); dx < std::abs(odom_x); dx += std::abs(odom_x)/10.0) {
      for (double dy = -std::abs(odom_y); dy < std::abs(odom_y); dy += std::abs(odom_y)/10.0) {
        for (double dt = -std::abs(odom_theta)*1.0; dt < std::abs(odom_theta)*1.0; dt += std::abs(odom_theta)/10.0) {
          //double scan_probability = do_bla_bla(scan, RobotState(prob_pose.x + dx,prob_pose.y + dy,prob_pose.theta + dt), map);
          //std::cout << scan_probability << std::endl << std::endl;
          //double scan_probability = esimate_scan_prob(scan, init_pose.x + dx,init_pose.y + dy,init_pose.theta + dt);
          double scan_probability = std::exp(esimate_scan_prob(scan, RobotState(init_pose.x + dx,init_pose.y + dy,init_pose.theta + dt), map));
          double pos_prob = scan_probability * normal_variation(dx, dy, dt,0.5);
          //std::cout << pos_prob << " ";
          ////if(std::abs(dx + 0.06666) < 0.001 && std::abs(dy) < 0.001 && std::abs(dt+0.13333) < 0.001)
          ////  std::cout << "he thinks, that the best probability is: " <<pos_prob << std::endl;
          ////if(std::abs(dx) < 0.00001 && std::abs(dy) < 0.00001 && std::abs(dt) < 0.00001)
          ////  std::cout << i << "  " << pos_prob << std::endl;
          if (pos_prob > max_prob) {
            max_prob = pos_prob;
            mdx = dx;
            mdy = dy;
            mdt = dt;
            //std::cout << "delta_x: " << mdx << " delta_y: " << mdy << "delta_t: " << mdt << std::endl;
            //std::cout << "odom_theta: " << odom_theta<< " delta_t: " << mdt << std::endl;

          }
        }
      }
      //if(!i % 10) { std::cout<<"1" << std::endl;}
      //i++;
    }
    //std::cout << "odom_x: " << odom_x << " odom_y: " << odom_y << " odom_t: " << odom_theta << std::endl;
    //std::cout << "  my_x: " << odom_x + mdx << "   my_y: " << odom_y + mdy << "   my_t: " << odom_theta + mdt << std::endl;
    //std::cout << "delta_x: " << mdx << " delta_y: " << mdy << "delta_t: " << mdt << std::endl;
    /*pose_delta.x = mdx;
    pose_delta.y = mdy;
    pose_delta.theta = mdt;*/

    pose_delta.x = mdx;
    pose_delta.y = mdy;
    pose_delta.theta = mdt;

    //pose_delta.x = mdx;
    //pose_delta.y = mdy;
    //pose_delta.theta = mdt;d
    //prev_pose.update(odom_x +mdx, odom_y +mdy, odom_theta +mdt);

    prev_scan = scan;
    prev_init_pose = init_pose;
    prev_init_pose.update(mdx,mdy,mdt);
    return -max_prob;
  }

  void reset_state() {
    table.reset();
  };
};



/*class
private:

};*/
