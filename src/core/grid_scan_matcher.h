/**
 * \file
 * \brief Defines some interfaces for Scan Matchers
 * There are class GridScanMatcherObserver, ScanCostEstimator, GridScanMatcher
 */

#ifndef __GRID_SCAN_MATCHER_H
#define __GRID_SCAN_MATCHER_H

#include <vector>
#include <memory>
#include <algorithm>

#include "state_data.h"
#include "sensor_data.h"

/**
 * \brief Interface of observer on Scan Matcher
 */
class GridScanMatcherObserver {
public:
  /**
   * A callback that is invoked on scan matching start
   * \param RobotState Pose of a robot
   * \param TransformedLaserScan A laser scan with a transformation
   * \param GridMap A grid map that is used by the matcher
   */
  virtual void on_matching_start(const RobotState &,           /*pose*/
                                 const TransformedLaserScan &, /*scan*/
                                 const GridMap &) {}    /*map*/

  /**
   * A callback that is invoked on scan testing
   * \param RobotState Pose of a robot
   * \param TransformedLaserScan A laser scan with a transformation
   */
  virtual void on_scan_test(const RobotState &,           /*pose*/
                            const TransformedLaserScan &, /*scan*/
                            double) {};                  /*score*/

  /**
   * A callback that is invoked on updating a pose of robot
   * \param RobotState Pose of a robot
   * \param TransformedLaserScan A laser scan with a transformation
   */
  virtual void on_pose_update(const RobotState &,            /*pose*/
                              const TransformedLaserScan &,  /*scan*/
                              double) {};                    /*score*/

  /**
   * A callback that is invoked on scan matching stops
   * \param RobotState Pose of a robot
   */
  virtual void on_matching_end(const RobotState &, /*delta*/
                               double) {};         /*best_score*/
};

/**
 * \brief Interface of Estimator of Scan Cost.
 * Cost - is a number that complies to a scan. As greater cost than worser scan
 */
class ScanCostEstimator {
public:

  /**
   * A callback that is invoked on estimating the cost of scan
   * \param pose Pose of a robot
   * \param scan A laser scan with a transformation
   * \param map A grid map that is used by the matcher
   */
  virtual double estimate_scan_cost(const RobotState &pose,
                                    const TransformedLaserScan &scan,
                                    const GridMap &map,
                                    double min_cost) = 0;
};

/**
 * \brief  Class that matches scans
 * Performes scan adjustment by altering robot pose in order to maximize
 * correspondence between the scan and a grid map. The rule of correspondence
 * computation is defined in ScanCostEstimator subclasses
 */
class GridScanMatcher {
public:
  GridScanMatcher(std::shared_ptr<ScanCostEstimator> estimator) :
    _cost_estimator(estimator) {}
  /**
   * A callback that is invoked on processing of scan
   * \param init_pose Pose of a robot
   * \param scan A laser scan with a transformation
   * \param map A grid map that is used by the matcher
   * \pose_delta Difference between real robot pose and estimation of pose
   */
  virtual double process_scan(const RobotState &init_pose,
                              const TransformedLaserScan &scan,
                              const GridMap &map,
                              RobotState &pose_delta) = 0;

  /// A callback that is invoked on reset of scan matcher state
  virtual void reset_state() {};

  /// Adds an observer
  void subscribe(std::shared_ptr<GridScanMatcherObserver> obs) {
    _observers.push_back(obs);
  }

  /// Removes an observer
  void unsubscribe(std::shared_ptr<GridScanMatcherObserver> obs) {
    // TODO: replace with more ideomatic C++
    std::vector<std::weak_ptr<GridScanMatcherObserver>> new_observers;
    for (auto &raw_obs : GridScanMatcher::observers()) {
      auto obs_ptr = raw_obs.lock();
      if (obs_ptr && obs_ptr != obs) {
        new_observers.push_back(raw_obs);
      }
    }
    _observers = new_observers;
  }

protected:
  /// Returns pointer to the cost estimator
  std::shared_ptr<ScanCostEstimator> cost_estimator() {
    return _cost_estimator;
  }

  /// Returns a reference to vector of pointers on observers
  std::vector<std::weak_ptr<GridScanMatcherObserver>> & observers() {
    return _observers;
  }
private:
  std::vector<std::weak_ptr<GridScanMatcherObserver>> _observers;
  std::shared_ptr<ScanCostEstimator> _cost_estimator;
};

#endif
