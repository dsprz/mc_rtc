/*
 * Copyright 2015-2022 CNRS-UM LIRMM, CNRS-AIST JRL
 */

#pragma once

#include <mc_tasks/SplineTrajectoryTask.h>

#include <mc_trajectory/BSpline.h>




namespace mc_tasks
{

/*! \brief Track a bezier curve with a robot surface.
 *
 * SplineTrajectoryTask takes care of much of the logic (task target updates,
 * orientation waypoint handling, etc.), and brings in all the functionalities
 * from TrajectoryTaskGeneric. This BSplineTrajectoryTask only handles specific
 * aspect of the Bezier curve.
 */
struct MC_TASKS_DLLAPI BSplineTrajectoryTask : public SplineTrajectoryTask<BSplineTrajectoryTask>
{
  friend struct SplineTrajectoryTask<BSplineTrajectoryTask>;
  using SplineTrajectoryBase = SplineTrajectoryTask<BSplineTrajectoryTask>;
  using waypoints_t = mc_trajectory::BSpline::waypoints_t;

public:
  using SplineTrajectoryBase::target;
  /**
   * \brief Creates a trajectory that follows a bspline curve
   *
   * \param frame Control frame
   * \param duration Duration of motion (eg time it takes to go from the current surface position to the curve's final
   * point)
   * \param stiffness Stiffness of the underlying TrajectoryTask (position
   * and orientation)
   * \param weight Task weight
   * \param target Final world pose to reach
   * \param posWp Waypoints in position
   * \param oriWp Waypoints in orientation specified as pairs of (time, orientation).
   */
  BSplineTrajectoryTask(const mc_rbdyn::RobotFrame & frame,
                        double duration,
                        double stiffness,
                        double weight,
                        const sva::PTransformd & target,
                        const waypoints_t & posWp = {},
                        const std::vector<std::pair<double, Eigen::Matrix3d>> & oriWp = {},
                      const bool &verbose_active = false);

  /**
   * \brief Creates a trajectory that follows a bspline curve
   *
   * \param robots Robots controlled by the task
   * \param robotIndex Which robot is controlled
   * \param surfaceName Surface controlled by the task, should belong to the controlled robot
   * \param duration Duration of motion (eg time it takes to go from the current surface position to the curve's final
   * point)
   * \param stiffness Stiffness of the underlying TrajectoryTask (position
   * and orientation)
   * \param weight Task weight
   * \param target Final world pose to reach
   * \param posWp Waypoints in position
   * \param oriWp Waypoints in orientation specified as pairs of (time, orientation).
   */
  BSplineTrajectoryTask(const mc_rbdyn::Robots & robots,
                        unsigned int robotIndex,
                        const std::string & surfaceName,
                        double duration,
                        double stiffness,
                        double weight,
                        const sva::PTransformd & target,
                        const waypoints_t & posWp = {},
                        const std::vector<std::pair<double, Eigen::Matrix3d>> & oriWp = {},
                      const bool &verbose_active = false);


  typedef Eigen::Vector3d Point;
  typedef Point point_t;
  typedef ndcurves::curve_constraints<point_t> curve_constraints_t;
  
  BSplineTrajectoryTask(const mc_rbdyn::RobotFrame & frame,
                      double duration,
                      double stiffness,
                      double weight,
                      const sva::PTransformd & target,
                      const curve_constraints_t & constr = {},
                      const waypoints_t & posWp = {},
                      const std::vector<std::pair<double, Eigen::Matrix3d>> & oriWp = {},
                    const bool &verbose_active = false);


  /*! \brief const accessor to the underlying spline (used by SplineTrajectoryTask)
   *
   * \returns The spline
   */
  const mc_trajectory::BSpline & spline() const { return bspline; };
  /*! \brief accessor to the underlying spline (used by SplineTrajectoryTask)
   *
   * \returns The spline
   */
  mc_trajectory::BSpline & spline() { return bspline; };

  /*! \brief Add interactive GUI elements to control the curve waypoints
   */
  void addToGUI(mc_rtc::gui::StateBuilder & gui);

  /** \brief Control points for the bezier curve (position)
   *
   * \param posWp Vector of position control points for the bezier curve.
   * Shouldn't include the starting and target position (use target() instead).
   */
  void posWaypoints(const waypoints_t & posWp);
  bool verbose_active = false;

protected:
  /*! \brief Sets the curve target pose
   * \param target Target pose for the curve
   */
  void targetPos(const Eigen::Vector3d & target);

  /** \brief Returns the curves target position */
  const Eigen::Vector3d & targetPos() const;

protected:
  mc_trajectory::BSpline bspline;
};

} // namespace mc_tasks
