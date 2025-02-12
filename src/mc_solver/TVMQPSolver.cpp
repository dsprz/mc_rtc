/*
 * Copyright 2015-2022 CNRS-UM LIRMM, CNRS-AIST JRL
 */

#include <mc_solver/TVMQPSolver.h>

#include <mc_solver/DynamicsConstraint.h>

#include <mc_tasks/MetaTask.h>

#include <mc_tvm/ContactFunction.h>
#include <mc_tvm/Robot.h>

#include <mc_rtc/gui/Force.h>

#include <tvm/solver/defaultLeastSquareSolver.h>
#include <tvm/task_dynamics/ProportionalDerivative.h>

#include <iostream>
#include <fstream>

#include <iomanip>  // For std::setw
#include <string>
namespace mc_solver
{

inline static Eigen::MatrixXd discretizedFrictionCone(double muI)
{
  Eigen::MatrixXd C(4, 3);
  double mu = muI / std::sqrt(2);
  C << Eigen::Matrix2d::Identity(), Eigen::Vector2d::Constant(mu), -Eigen::Matrix2d::Identity(),
      Eigen::Vector2d::Constant(mu);
  return C;
}

TVMQPSolver::TVMQPSolver(mc_rbdyn::RobotsPtr robots, double dt)
: QPSolver(robots, dt, Backend::TVM), solver_(tvm::solver::DefaultLSSolverOptions{})
{
}

TVMQPSolver::TVMQPSolver(double dt) : QPSolver(dt, Backend::TVM), solver_(tvm::solver::DefaultLSSolverOptions{}) {}

size_t TVMQPSolver::getContactIdx(const mc_rbdyn::Contact & contact)
{
  for(size_t i = 0; i < contacts_.size(); ++i)
  {
    if(contacts_[i] == contact) { return i; }
  }
  return contacts_.size();
}

void TVMQPSolver::setContacts(ControllerToken, const std::vector<mc_rbdyn::Contact> & contacts)
{
  for(const auto & c : contacts) { addContact(c); }
  size_t i = 0;
  for(auto it = contacts_.begin(); it != contacts_.end();)
  {
    const auto & c = *it;
    if(std::find(contacts.begin(), contacts.end(), c) == contacts.end())
    {
      const std::string & r1 = robots().robot(c.r1Index()).name();
      const std::string & r1S = c.r1Surface()->name();
      const std::string & r2 = robots().robot(c.r2Index()).name();
      const std::string & r2S = c.r2Surface()->name();
      logger_->removeLogEntry("contact_" + r1 + "::" + r1S + "_" + r2 + "::" + r2S);
      if(gui_) { gui_->removeElement({"Contacts", "Forces"}, fmt::format("{}::{}/{}::{}", r1, r1S, r2, r2S)); }
      it = removeContact(i);
    }
    else
    {
      ++i;
      ++it;
    }
  }
}

const sva::ForceVecd TVMQPSolver::desiredContactForce(const mc_rbdyn::Contact & id) const
{
  const auto & r1 = robot(id.r1Index());
  auto it1 = dynamics_.find(r1.name());
  if(it1 != dynamics_.end()) { return it1->second->dynamicFunction().contactForce(r1.frame(id.r1Surface()->name())); }
  const auto & r2 = robot(id.r2Index());
  auto it2 = dynamics_.find(r2.name());
  if(it2 != dynamics_.end()) { return it2->second->dynamicFunction().contactForce(r2.frame(id.r2Surface()->name())); }
  return sva::ForceVecd::Zero();
}

double TVMQPSolver::solveTime()
{
  return solve_dt_.count();
}

double TVMQPSolver::solveAndBuildTime()
{
  return solve_dt_.count();
}

bool TVMQPSolver::run_impl(FeedbackType fType)
{
  switch(fType)
  {
    case FeedbackType::None:
      return runOpenLoop();
    case FeedbackType::Joints:
      return runJointsFeedback(false);
    case FeedbackType::JointsWVelocity:
      return runJointsFeedback(true);
    case FeedbackType::ObservedRobots:
      return runClosedLoop(true);
    case FeedbackType::ClosedLoopIntegrateReal:
      return runClosedLoop(false);
    default:
      mc_rtc::log::error("FeedbackType set to unknown value");
      return false;
  }
}

bool TVMQPSolver::runCommon()
{
  for(auto & c : constraints_) { c->update(*this); }
  for(auto & t : metaTasks_)
  {
    t->update(*this);
    t->incrementIterInSolver();
  }

  auto start_t = mc_rtc::clock::now();
  auto r = solver_.solve(problem_);
  solve_dt_ = mc_rtc::clock::now() - start_t;
  return r;
}

bool TVMQPSolver::runOpenLoop()
{
  int i = 0;
  for(auto & robot : *robots_p)
  {
     std::cout << "i: " << i << std::endl;
     std::cout << robot.name() << std ::endl; 

    if(i == 0)
    {
      std :: cout  << robot.name() << std::endl;
 
// Prints tvm q 
      std::ofstream pos("tvmrobot-q.csv");

        if(pos.is_open())
        {
          pos << std::setw(15) << robot.tvmRobot().q()->value() << " ";
          pos.close();
          std::cout << "Matrix saved to tvmrobot-q.csv" << std::endl;
        }else{
          std::cerr << "Error opening tvmpos" << std::endl;
        } 

// Prints tvm alpha
      std::ofstream speed("tvmrobot-alpha.csv");

      if(speed.is_open())
      {
        speed << std::setw(15) << robot.tvmRobot().alpha()->value() << " ";
        speed.close();
        std::cout << "Matrix saved to tvmrobot-alpha.csv" << std::endl;
      }else{
        std::cerr << "Error opening tvmspeed" << std::endl;
      } 
// robot q
        std::ofstream posR("robot-q.csv");

        if(posR.is_open())
        {
          for(const auto& row: robot.mbc().q)
          {
            for(const auto& col: row)
            {
              posR << std::setw(15) << col << " ";
            }
            posR << "\n"; 
          }
          posR.close();
          std::cout << "Matrix saved to robot-q.csv" << std::endl;
        }else{
          std::cerr << "Error opening pos" << std::endl;
        } 
// Prints alpha
      std::ofstream speedR("robot-alpha.csv");

      if(speedR.is_open())
      {
          for(const auto& row: robot.mbc().alpha)
          {
            for(const auto& col: row)
            {
              speedR << std::setw(15) << col << " ";
            }
            speedR << "\n"; 
          }
        speedR.close();
        std::cout << "Matrix saved to robot-alpha.csv" << std::endl;
      }else{
        std::cerr << "Error opening speed" << std::endl;
      } 
    }
    i++;

  }

  if(runCommon())
  {
    std::cout << "(*robots_p)[0].name(): " << (*robots_p)[0].name() << std::endl;
    int j = 0;
    for(auto & robot : *robots_p)
    {
      auto & mb = robot.mb();
      if(mb.nrDof() > 0) { updateRobot(robot); }

      std::cout << robot.name() << std ::endl;

    if(j==0)
    {
    //H
        std::ofstream mass("matrix-tvm-H.csv");

        if(mass.is_open())
        {
          mass << robot.tvmRobot().H().matrix();
          mass.close();
          std::cout << "Matrix saved to matrix-tvm-H.csv" << std::endl;
        }else{
          std::cerr << "Error opening matrix-tvm-H" << std::endl;
        }  

        std::ofstream cor("matrix-tvm-C.csv");
    // C
        if(cor.is_open())
        {
          cor << robot.tvmRobot().C().matrix();
          cor.close();
          std::cout << "Matrix saved to matrix-tvm-C.csv" << std::endl;
        }else{
          std::cerr << "Error opening matrix-tvm-C" << std::endl;
        }

        // Prints robot qdotdot
        std::ofstream accR("robot-alphaD.csv");

        if(accR.is_open())
        {
            for(const auto& row: robot.mbc().alphaD)
            {
              for(const auto& col: row)
              {
                accR << std::setw(15) << col << " ";
              }
              accR << "\n"; 
            }
          accR.close();
          std::cout << "Matrix saved to robot-alphaD.csv" << std::endl;
        }else{
          std::cerr << "Error opening speed" << std::endl;
        } 

//print robot tau

        std::ofstream tau("robot-tau.csv");

        if(tau.is_open())
        {
            for(const auto& row: robot.mbc().jointTorque)
            {
              for(const auto& col: row)
              {
                tau << std::setw(15) << col << " ";
              }
              tau << "\n"; 
            }
          tau.close();
          std::cout << "Matrix saved to robot-tau.csv" << std::endl;
        }else{
          std::cerr << "Error opening tau" << std::endl;
        } 

//print forces

        std::ofstream force("robot-force.csv");

        if(force.is_open())
        {
          for(const auto& f: robot.mbc().force)
          {

            force << f;
            force << "\n";

          }
          force.close();
          std::cout << "Matrix saved to robot-force.csv" << std::endl;
        }else{
          std::cerr << "Error opening force" << std::endl;
        } 

        // Prints tvm alphaD
    std::ofstream acc("tvmrobot-alphaD.csv");

    if(acc.is_open())
    {
      acc << std::setw(15) << robot.tvmRobot().alphaD()->value() << " ";
      acc.close();
      std::cout << "Matrix saved to tvmrobot-alphaD.csv" << std::endl;
    }else{
      std::cerr << "Error opening acc" << std::endl;
    } 

    // Prints tvm torque
      std::ofstream torque("tvmrobot-torque.csv");

      if(torque.is_open())
      {
        torque << std::setw(15) << robot.tvmRobot().tau()->value() << " ";
        torque.close();
        std::cout << "Matrix saved to tvmrobot-torque.csv" << std::endl;
      }else{
        std::cerr << "Error opening torque" << std::endl;
      } 

    // Prints tvm contact forces
      std::ofstream forces("tvmrobot-tauexternal.csv");

      if(forces.is_open())
      {
        forces << std::setw(15) << robot.tvmRobot().tauExternal().value() << " ";
        forces.close();
        std::cout << "Matrix saved to tvm-tauexternal.csv" << std::endl;
      }else{
        std::cerr << "Error opening tau external" << std::endl;
      } 


    }  
        j++; 
    }
    return true;
  }
  return false;
}

bool TVMQPSolver::runJointsFeedback(bool wVelocity)
{
  if(control_q_.size() < robots().size())
  {
    prev_encoders_.resize(robots().size());
    encoders_alpha_.resize(robots().size());
    control_q_.resize(robots().size());
    control_alpha_.resize(robots().size());
  }
  for(size_t i = 0; i < robots().size(); ++i)
  {
    auto & robot = robots_p->robot(i);
    control_q_[i] = robot.q();
    control_alpha_[i] = robot.alpha();
    const auto & encoders = robot.encoderValues();
    if(encoders.size())
    {
      // FIXME Not correct for every joint types
      if(prev_encoders_[i].size() == 0)
      {
        prev_encoders_[i] = robot.encoderValues();
        encoders_alpha_[i].resize(prev_encoders_[i].size());
      }
      for(size_t j = 0; j < encoders.size(); ++j)
      {
        encoders_alpha_[i][j] = (encoders[j] - prev_encoders_[i][j]) / timeStep;
        prev_encoders_[i][j] = encoders[j];
      }
      const auto & rjo = robot.module().ref_joint_order();
      for(size_t j = 0; j < rjo.size(); ++j)
      {
        auto jI = robot.jointIndexInMBC(j);
        if(jI == -1) { continue; }
        robot.q()[static_cast<size_t>(jI)][0] = encoders[j];
        if(wVelocity) { robot.alpha()[static_cast<size_t>(jI)][0] = encoders_alpha_[i][j]; }
      }
      robot.forwardKinematics();
      robot.forwardVelocity();
      robot.forwardAcceleration();
    }
  }
  if(runCommon())
  {
    for(size_t i = 0; i < robots_p->size(); ++i)
    {
      auto & robot = robots_p->robot(i);
      if(robot.mb().nrDof() == 0) { continue; }
      robot.q() = control_q_[i];
      robot.alpha() = control_alpha_[i];
      updateRobot(robot);
    }
    return true;
  }
  return false;
}

bool TVMQPSolver::runClosedLoop(bool integrateControlState)
{
  if(control_q_.size() < robots().size())
  {
    control_q_.resize(robots().size());
    control_alpha_.resize(robots().size());
  }

  for(size_t i = 0; i < robots().size(); ++i)
  {
    auto & robot = robots().robot(i);
    const auto & realRobot = realRobots().robot(i);

    // Save old integrator state
    if(integrateControlState)
    {
      control_q_[i] = robot.mbc().q;
      control_alpha_[i] = robot.mbc().alpha;
    }

    // Set robot state from estimator
    robot.mbc().q = realRobot.mbc().q;
    robot.mbc().alpha = realRobot.mbc().alpha;
    robot.forwardKinematics();
    robot.forwardVelocity();
    robot.forwardAcceleration();

    std::cout << "i: " << i << std::endl; 

    if(i == 0)
    {
      std :: cout  << robot.name() << std::endl;
 
// Prints tvm q 
      std::ofstream pos("tvmrobot-q.csv");

        if(pos.is_open())
        {
          pos << std::setw(15) << robot.tvmRobot().q()->value() << " ";
          pos.close();
          std::cout << "Matrix saved to tvmrobot-q.csv" << std::endl;
        }else{
          std::cerr << "Error opening tvmpos" << std::endl;
        } 

// Prints tvm alpha
      std::ofstream speed("tvmrobot-alpha.csv");

      if(speed.is_open())
      {
        speed << std::setw(15) << robot.tvmRobot().alpha()->value() << " ";
        speed.close();
        std::cout << "Matrix saved to tvmrobot-alpha.csv" << std::endl;
      }else{
        std::cerr << "Error opening tvmspeed" << std::endl;
      } 
// robot q
        std::ofstream posR("robot-q.csv");

        if(posR.is_open())
        {
          for(const auto& row: robot.mbc().q)
          {
            for(const auto& col: row)
            {
              posR << std::setw(15) << col << " ";
            }
            posR << "\n"; 
          }
          posR.close();
          std::cout << "Matrix saved to robot-q.csv" << std::endl;
        }else{
          std::cerr << "Error opening pos" << std::endl;
        } 
// Prints alpha
      std::ofstream speedR("robot-alpha.csv");

      if(speedR.is_open())
      {
          for(const auto& row: robot.mbc().alpha)
          {
            for(const auto& col: row)
            {
              speedR << std::setw(15) << col << " ";
            }
            speedR << "\n"; 
          }
        speedR.close();
        std::cout << "Matrix saved to robot-alpha.csv" << std::endl;
      }else{
        std::cerr << "Error opening speed" << std::endl;
      } 
    }

  }

  std::cout << "runCommon" << std::endl;
  // Solve QP and integrate
  if(runCommon())
  {
    // exit(0);
    for(size_t i = 0; i < robots_p->size(); ++i)
    {
      auto & robot = robots_p->robot(i);
      if(robot.mb().nrDof() == 0) { continue; }
      if(integrateControlState)
      {
        robot.q() = control_q_[i];
        robot.alpha() = control_alpha_[i];
      }
      updateRobot(robot);

// Prints robot qdotdot
        std::ofstream accR("robot-alphaD.csv");

        if(accR.is_open())
        {
            for(const auto& row: robot.mbc().alphaD)
            {
              for(const auto& col: row)
              {
                accR << std::setw(15) << col << " ";
              }
              accR << "\n"; 
            }
          accR.close();
          std::cout << "Matrix saved to robot-alphaD.csv" << std::endl;
        }else{
          std::cerr << "Error opening speed" << std::endl;
        } 

//print robot tau

        std::ofstream tau("robot-tau.csv");

        if(tau.is_open())
        {
            for(const auto& row: robot.mbc().jointTorque)
            {
              for(const auto& col: row)
              {
                tau << std::setw(15) << col << " ";
              }
              tau << "\n"; 
            }
          tau.close();
          std::cout << "Matrix saved to robot-tau.csv" << std::endl;
        }else{
          std::cerr << "Error opening tau" << std::endl;
        } 

//print forces

        std::ofstream force("robot-force.csv");

        if(force.is_open())
        {
          for(const auto& f: robot.mbc().force)
          {

            force << f;
            force << "\n";

          }
          force.close();
          std::cout << "Matrix saved to robot-force.csv" << std::endl;
        }else{
          std::cerr << "Error opening force" << std::endl;
        } 

  
// H
    std::ofstream mass("matrix-tvm-H.csv");

    if(mass.is_open())
    {
      mass << robot.tvmRobot().H().matrix();
      mass.close();
      std::cout << "Matrix saved to matrix-tvm-H.csv" << std::endl;
    }else{
      std::cerr << "Error opening file" << std::endl;
    }  

    std::ofstream cor("matrix-tvm-C.csv");
// C
    if(cor.is_open())
    {
      cor << robot.tvmRobot().C().matrix();
      cor.close();
      std::cout << "Matrix saved to matrix-tvm-C.csv" << std::endl;
    }else{
      std::cerr << "Error opening file" << std::endl;
    }  

// Prints tvm alphaD
    std::ofstream acc("tvmrobot-alphaD.csv");

    if(acc.is_open())
    {
      acc << std::setw(15) << robot.tvmRobot().alphaD()->value() << " ";
      acc.close();
      std::cout << "Matrix saved to tvmrobot-alphaD.csv" << std::endl;
    }else{
      std::cerr << "Error opening acc" << std::endl;
    } 

// Prints tvm torque
  std::ofstream torque("tvmrobot-torque.csv");

  if(torque.is_open())
  {
    torque << std::setw(15) << robot.tvmRobot().tau()->value() << " ";
    torque.close();
    std::cout << "Matrix saved to tvmrobot-torque.csv" << std::endl;
  }else{
    std::cerr << "Error opening torque" << std::endl;
  } 

// Prints tvm contact forces
  std::ofstream forces("tvmrobot-tauexternal.csv");

  if(forces.is_open())
  {
    forces << std::setw(15) << robot.tvmRobot().tauExternal().value() << " ";
    forces.close();
    std::cout << "Matrix saved to tvm-tauexternal.csv" << std::endl;
  }else{
    std::cerr << "Error opening tau external" << std::endl;
  } 

 }
    return true;
  }
  return false;
}

void TVMQPSolver::updateRobot(mc_rbdyn::Robot & robot)
{
  auto & tvm_robot = robot.tvmRobot();
  rbd::vectorToParam(tvm_robot.tau()->value(), robot.controlTorque());
  rbd::vectorToParam(tvm_robot.alphaD()->value(), robot.alphaD());
  robot.eulerIntegration(timeStep);
  robot.forwardKinematics();
  robot.forwardVelocity();
  robot.forwardAcceleration();
}

void TVMQPSolver::addDynamicsConstraint(mc_solver::DynamicsConstraint * dyn)
{
  const auto & r = robot(dyn->robotIndex());
  if(dynamics_.count(r.name()))
  {
    mc_rtc::log::error_and_throw("Only one dynamic constraint can be added for a given robot and {} already has one",
                                 r.name());
  }
  dynamics_[r.name()] = dyn;
  for(size_t i = 0; i < contacts_.size(); ++i)
  {
    const auto & contact = contacts_[i];
    auto & data = contactsData_[i];
    bool isR1 = contact.r1Index() == dyn->robotIndex();
    bool isR2 = contact.r2Index() == dyn->robotIndex();
    if(isR1 || isR2)
    {
      const auto & r1 = robot(contact.r1Index());
      const auto & r2 = robot(contact.r2Index());
      const auto & s1 = *contact.r1Surface();
      const auto & s2 = *contact.r2Surface();
      const auto & f1 = r1.frame(s1.name());
      const auto & f2 = r2.frame(s2.name());
      const auto C = discretizedFrictionCone(contact.friction());
      // FIXME Debug mc_rbdyn::intersection
      // auto s1Points = mc_rbdyn::intersection(s1, s2);
      const auto & s1Points = s1.points();
      if(isR1) { addContactToDynamics(r1.name(), f1, s1Points, data.f1_, data.f1Constraints_, C, 1.0); }
      if(isR2)
      {
        std::vector<sva::PTransformd> s2Points;
        s2Points.reserve(s1Points.size());
        auto X_b2_b1 =
            r1.mbc().bodyPosW[r1.bodyIndexByName(f1.body())] * r2.mbc().bodyPosW[r2.bodyIndexByName(f2.body())].inv();
        for(const auto & X_b1_p : s1Points) { s2Points.push_back(X_b1_p * X_b2_b1); }
        addContactToDynamics(r2.name(), f2, s2Points, data.f2_, data.f2Constraints_, C, -1.0);
      }
    }
  }

  logger_->addLogEntry("TVM matrix H " + (*robots_p)[0].name() + std::to_string(dynamics_.size()), [&, this](){return (*robots_p)[0].tvmRobot().H().norm();});
  logger_->addLogEntry("TVM matrix C " + (*robots_p)[0].name() + std::to_string(dynamics_.size()), [&, this](){return (*robots_p)[0].tvmRobot().C().norm();});
}

void TVMQPSolver::removeDynamicsConstraint(mc_solver::ConstraintSet * maybe_dyn)
{
  for(const auto & [r, dyn] : dynamics_)
  {
    if(static_cast<const mc_solver::ConstraintSet *>(dyn) == maybe_dyn)
    {
      return removeDynamicsConstraint(static_cast<mc_solver::DynamicsConstraint *>(maybe_dyn));
    }
  }
}

void TVMQPSolver::removeDynamicsConstraint(mc_solver::DynamicsConstraint * dyn)
{
  const auto & r = robot(dyn->robotIndex());
  dynamics_.erase(r.name());
  for(size_t i = 0; i < contacts_.size(); ++i)
  {
    const auto & contact = contacts_[i];
    auto & data = contactsData_[i];
    auto clearContacts = [&](const std::string & robot, tvm::VariableVector & forces,
                             std::vector<tvm::TaskWithRequirementsPtr> & constraints)
    {
      if(robot != r.name()) { return; }
      for(auto & c : constraints) { problem_.remove(*c); }
      constraints.clear();
      forces = tvm::VariableVector();
    };
    const auto & r1 = robot(contact.r1Index());
    clearContacts(r1.name(), data.f1_, data.f1Constraints_);
    const auto & r2 = robot(contact.r2Index());
    clearContacts(r2.name(), data.f2_, data.f2Constraints_);
  }
}

void TVMQPSolver::addContactToDynamics(const std::string & robot,
                                       const mc_rbdyn::RobotFrame & frame,
                                       const std::vector<sva::PTransformd> & points,
                                       tvm::VariableVector & forces,
                                       std::vector<tvm::TaskWithRequirementsPtr> & constraints,
                                       const Eigen::MatrixXd & frictionCone,
                                       double dir)
{
  auto it = dynamics_.find(robot);
  if(it == dynamics_.end()) { return; }
  if(constraints.size())
  {
    // FIXME Instead of this we should be able to change C
    for(const auto & c : constraints) { problem_.remove(*c); }
    constraints.clear();
  }
  else
  {
    it->second->removeFromSolverImpl(*this);
    auto & dyn = it->second->dynamicFunction();
    forces = dyn.addContact(frame, points, dir);
    it->second->addToSolverImpl(*this);
  }
  for(int i = 0; i < forces.numberOfVariables(); ++i)
  {
    auto & f = forces[i];
    constraints.push_back(problem_.add(dir * frictionCone * f >= 0.0, {tvm::requirements::PriorityLevel(0)}));
  }
}

auto TVMQPSolver::addVirtualContactImpl(const mc_rbdyn::Contact & contact) -> std::tuple<size_t, bool>
{
  bool hasWork = false;
  auto idx = getContactIdx(contact);
  if(idx < contacts_.size())
  {
    const auto & oldContact = contacts_[idx];
    if(oldContact.dof() == contact.dof() && oldContact.friction() == contact.friction())
    {
      return std::make_tuple(idx, hasWork);
    }
    hasWork = contact.friction() != oldContact.friction();
    contacts_[idx] = contact;
  }
  else
  {
    hasWork = true;
    contacts_.push_back(contact);
  }
  auto & data = idx < contactsData_.size() ? contactsData_[idx] : contactsData_.emplace_back();
  const auto & r1 = robot(contact.r1Index());
  const auto & r2 = robot(contact.r2Index());
  const auto & f1 = r1.frame(contact.r1Surface()->name());
  const auto & f2 = r2.frame(contact.r2Surface()->name());
  if(!data.contactConstraint_) // New contact
  {
    auto contact_fn = std::make_shared<mc_tvm::ContactFunction>(f1, f2, contact.dof());
    data.contactConstraint_ = problem_.add(contact_fn == 0., tvm::task_dynamics::PD(1.0 / dt(), 1.0 / dt()),
                                           {tvm::requirements::PriorityLevel(0)});
    logger_->addLogEntry(fmt::format("contact_{}::{}_{}::{}", r1.name(), f1.name(), r2.name(), f2.name()),
                         [this, contact]() { return desiredContactForce(contact); });
    gui_->addElement({"Contacts", "Forces"},
                     mc_rtc::gui::Force(
                         fmt::format("{}::{}/{}::{}", r1.name(), f1.name(), r2.name(), f2.name()), [this, contact]()
                         { return desiredContactForce(contact); }, [&f1]() { return f1.position(); }));
  }
  else
  {
    auto contact_fn = std::static_pointer_cast<mc_tvm::ContactFunction>(data.contactConstraint_->task.function());
    contact_fn->dof(contact.dof());
  }
  return std::make_tuple(idx, hasWork);
}

void TVMQPSolver::addContact(const mc_rbdyn::Contact & contact)
{
  size_t idx = contacts_.size();
  bool hasWork = false;
  std::tie(idx, hasWork) = addVirtualContactImpl(contact);
  if(!hasWork) { return; }
  auto & data = contactsData_[idx];
  const auto & r1 = robot(contact.r1Index());
  const auto & r2 = robot(contact.r2Index());
  const auto & s1 = *contact.r1Surface();
  const auto & s2 = *contact.r2Surface();
  const auto & f1 = r1.frame(s1.name());
  const auto & f2 = r2.frame(s2.name());
  // FIXME Let the user decide how much the friction cone should be discretized
  auto C = discretizedFrictionCone(contact.friction());
  auto addContactForce = [&](const std::string & robot, const mc_rbdyn::RobotFrame & frame,
                             const std::vector<sva::PTransformd> & points, tvm::VariableVector & forces,
                             std::vector<tvm::TaskWithRequirementsPtr> & constraints, double dir)
  { addContactToDynamics(robot, frame, points, forces, constraints, C, dir); };
  // FIXME These points computation are a waste of time if they are not needed
  // FIXME Debug mc_rbdyn::intersection
  // auto s1Points = mc_rbdyn::intersection(s1, s2);
  auto s1Points = s1.points();
  addContactForce(r1.name(), f1, s1Points, data.f1_, data.f1Constraints_, 1.0);
  std::vector<sva::PTransformd> s2Points;
  s2Points.reserve(s1Points.size());
  auto X_b2_b1 =
      r1.mbc().bodyPosW[r1.bodyIndexByName(f1.body())] * r2.mbc().bodyPosW[r2.bodyIndexByName(f2.body())].inv();
  std::transform(s1Points.begin(), s1Points.end(), std::back_inserter(s2Points),
                 [&](const auto & X_b1_p) { return X_b1_p * X_b2_b1; });
  addContactForce(r2.name(), f2, s2Points, data.f2_, data.f2Constraints_, -1.0);
}

auto TVMQPSolver::removeContact(size_t idx) -> ContactIterator
{
  auto & contact = contacts_[idx];
  auto & data = contactsData_[idx];
  const auto & r1 = robot(contact.r1Index());
  auto r1DynamicsIt = dynamics_.find(r1.name());
  if(r1DynamicsIt != dynamics_.end())
  {
    r1DynamicsIt->second->removeFromSolverImpl(*this);
    r1DynamicsIt->second->dynamicFunction().removeContact(r1.frame(contact.r1Surface()->name()));
    r1DynamicsIt->second->addToSolverImpl(*this);
  }
  const auto & r2 = robot(contact.r2Index());
  auto r2DynamicsIt = dynamics_.find(r2.name());
  if(r2DynamicsIt != dynamics_.end())
  {
    r2DynamicsIt->second->removeFromSolverImpl(*this);
    r2DynamicsIt->second->dynamicFunction().removeContact(r2.frame(contact.r2Surface()->name()));
    r2DynamicsIt->second->addToSolverImpl(*this);
  }
  for(const auto & c : data.f1Constraints_) { problem_.remove(*c); }
  for(const auto & c : data.f2Constraints_) { problem_.remove(*c); }
  if(data.contactConstraint_)
  {
    problem_.remove(*data.contactConstraint_);
    data.contactConstraint_.reset();
  }
  contactsData_.erase(contactsData_.begin() + static_cast<decltype(contacts_)::difference_type>(idx));
  return contacts_.erase(contacts_.begin() + static_cast<decltype(contacts_)::difference_type>(idx));
}

} // namespace mc_solver
