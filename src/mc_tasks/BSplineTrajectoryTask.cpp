/*
 * Copyright 2015-2022 CNRS-UM LIRMM, CNRS-AIST JRL
 */

#include <mc_tasks/BSplineTrajectoryTask.h>

#include <mc_tasks/MetaTaskLoader.h>

#include <mc_trajectory/BSpline.h>
#include <mc_trajectory/InterpolatedRotation.h>

#include <mc_rtc/deprecated.h>

namespace mc_tasks
{
  //  Helene et Jimmy
  typedef Eigen::Vector3d Point;
  typedef Point point_t;
  typedef ndcurves::curve_constraints<point_t> curve_constraints_t;
using BSpline = mc_trajectory::BSpline;

BSplineTrajectoryTask::BSplineTrajectoryTask(const mc_rbdyn::Robots & robots,
                                             unsigned int robotIndex,
                                             const std::string & surfaceName,
                                             double duration,
                                             double stiffness,
                                             double weight,
                                             const sva::PTransformd & target,
                                             const waypoints_t & posWp,
                                             const std::vector<std::pair<double, Eigen::Matrix3d>> & oriWp,
                                            const bool &verbose_active)
                                            
: BSplineTrajectoryTask(robots.robot(robotIndex).frame(surfaceName), duration, stiffness, weight, target, posWp, oriWp)
{
  if(verbose_active)
  {
    std::cout << "Using the BSplineTrajectoryTask constructor WITHOUT constraints" << std::endl;
  }
}

BSplineTrajectoryTask::BSplineTrajectoryTask(const mc_rbdyn::RobotFrame & frame,
                                             double duration,
                                             double stiffness,
                                             double weight,
                                             const sva::PTransformd & target,
                                             const waypoints_t & posWp,
                                             const std::vector<std::pair<double, Eigen::Matrix3d>> & oriWp,
                                            const bool &verbose_active)
: SplineTrajectoryTask<BSplineTrajectoryTask>(frame, duration, stiffness, weight, target.rotation(), oriWp),
  bspline(duration, frame.position().translation(), target.translation(), posWp)
{
  type_ = "bspline_trajectory";
  name_ = "bspline_trajectory_" + frame.robot().name() + "_" + frame.name();
  if(verbose_active)
  {
    std::cout << "Using the BSplineTrajectoryTask constructor WITHOUT constraints" << std::endl;
  }

}

// Add the declaration in /home/jimmyvu/Documents/mc_rtc/src/mc_rtc/include/mc_tasks/BSplineTrajectoryTask.h
BSplineTrajectoryTask::BSplineTrajectoryTask(const mc_rbdyn::RobotFrame & frame,
                      double duration,
                      double stiffness,
                      double weight,
                      const sva::PTransformd & target,
                      const curve_constraints_t & constr,
                      const waypoints_t & posWp,
                      const std::vector<std::pair<double, Eigen::Matrix3d>> & oriWp,
                    const bool &verbose_active)
: SplineTrajectoryTask<BSplineTrajectoryTask>(frame, 
                                              duration, 
                                              stiffness, 
                                              weight,
                                               target.rotation(),
                                               oriWp),
  bspline(duration,
     frame.position().translation(), 
     target.translation(), 
     constr, 
     posWp,
    verbose_active)
{
  type_ = "bspline_trajectory";
  name_ = "bspline_trajectory_" + frame.robot().name() + "_" + frame.name();
  if(verbose_active)
  {
    std::cout << "Using the BSplineTrajectoryTask constructor with constraints" << std::endl;
    std::cout << "In BSplineTrajectoryTask.cpp, printing constr.end_vel.z() : " << constr.end_vel.z() << std::endl;
    // std::cout << "In BSplineTrajectoryTask.cpp, printing bspline.get_bezier().constr.end_vel.z() : " << bspline.get_bezier()->constr_.end_vel.z() << std::endl;
  }
}

void BSplineTrajectoryTask::posWaypoints(const BSpline::waypoints_t & posWp)
{
  bspline.waypoints(posWp);
}

void BSplineTrajectoryTask::targetPos(const Eigen::Vector3d & target)
{
  bspline.target(target);
}

const Eigen::Vector3d & BSplineTrajectoryTask::targetPos() const
{
  return bspline.target();
}

void BSplineTrajectoryTask::addToGUI(mc_rtc::gui::StateBuilder & gui)
{
  SplineTrajectoryBase::addToGUI(gui);
  bspline.addToGUI(gui, {"Tasks", name_});
}

} // namespace mc_tasks

namespace
{
static auto registered = mc_tasks::MetaTaskLoader::register_load_function(
    "bspline_trajectory",
    [](mc_solver::QPSolver & solver, const mc_rtc::Configuration & config)
    {
      sva::PTransformd finalTarget_;
      mc_tasks::BSplineTrajectoryTask::waypoints_t waypoints;
      std::vector<std::pair<double, Eigen::Matrix3d>> oriWp;
      const auto robotIndex = robotIndexFromConfig(config, solver.robots(), "bspline_trajectory");

      bool has_targetSurface = config.has("targetSurface");
      bool has_targetFrame = config.has("targetFrame");

      if(has_targetSurface || has_targetFrame)
      { // Target defined from a target surface, with an offset defined
        // in the surface coordinates
        if(has_targetSurface)
        {
          mc_rtc::log::deprecated("ExactCubicTrajectoryTaskLoading", "targetSurface", "targetFrame");
        }
        const auto & c = config(has_targetSurface ? "targetSurface" : "targetFrame");
        const auto & targetSurfaceName = c(has_targetSurface ? "surface" : "frame");
        const auto & robot = robotFromConfig(c, solver.robots(), "bspline_trajectory::targetFrame");

        const sva::PTransformd & targetFrame = robot.frame(targetSurfaceName).position();
        const Eigen::Vector3d trans = c("translation", Eigen::Vector3d::Zero().eval());
        const Eigen::Matrix3d rot = c("rotation", Eigen::Matrix3d::Identity().eval());
        sva::PTransformd offset(rot, trans);
        finalTarget_ = offset * targetFrame;

        if(c.has("controlPoints"))
        {
          // Control points offsets defined wrt to the target surface frame
          const auto & controlPoints = c("controlPoints");
          waypoints.resize(controlPoints.size());
          for(unsigned int i = 0; i < controlPoints.size(); ++i)
          {
            const Eigen::Vector3d wp = controlPoints[i];
            sva::PTransformd X_offset(wp);
            waypoints[i] = (X_offset * targetFrame).translation();
          }
        }

        if(c.has("oriWaypoints"))
        {
          std::vector<std::pair<double, Eigen::Matrix3d>> oriWaypoints = c("oriWaypoints");
          for(const auto & wp : oriWaypoints)
          {
            const sva::PTransformd offset{wp.second};
            const sva::PTransformd ori = offset * targetFrame;
            oriWp.push_back(std::make_pair(wp.first, ori.rotation()));
          }
        }
      }
      else
      { // Absolute target pose
        finalTarget_ = config("target");

        if(config.has("controlPoints"))
        {
          // Control points defined in world coordinates
          const auto & controlPoints = config("controlPoints");
          waypoints.resize(controlPoints.size());
          for(unsigned int i = 0; i < controlPoints.size(); ++i)
          {
            const Eigen::Vector3d wp = controlPoints[i];
            waypoints[i] = wp;
          }
        }

        oriWp = config("oriWaypoints", std::vector<std::pair<double, Eigen::Matrix3d>>{});
      }

      const auto & frame = [&]() -> const mc_rbdyn::RobotFrame &
      {
        if(config.has("surface"))
        {
          mc_rtc::log::deprecated("ExactCubicTrajectoryTask", "surface", "frame");
          return solver.robots().robot(robotIndex).frame(config("surface"));
        }
        return solver.robots().robot(robotIndex).frame(config("frame"));
      }();

      auto t =
          std::make_shared<mc_tasks::BSplineTrajectoryTask>(frame, config("duration", 10.), config("stiffness", 100.),
                                                            config("weight", 500.), finalTarget_, waypoints, oriWp);
      t->load(solver, config);
      const auto displaySamples = config("displaySamples", t->displaySamples());
      t->displaySamples(displaySamples);
      t->pause(config("paused", false));
      return t;
    });
} // namespace
