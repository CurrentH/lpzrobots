/***************************************************************************
 *   Copyright (C) 2005 by Robot Group Leipzig                             *
 *    martius@informatik.uni-leipzig.de                                    *
 *    fhesse@informatik.uni-leipzig.de                                     *
 *    der@informatik.uni-leipzig.de                                        *
 *    frankguettler@gmx.de                                                 *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 *                                                                         *
 *   $Log$
 *   Revision 1.2  2007-02-12 13:30:40  martius
 *   dog looks allready nicer
 *
 *   Revision 1.1  2007/01/26 12:07:09  martius
 *   orientationsensor added
 *
 *   Revision 1.2  2006/07/14 12:23:55  martius
 *   selforg becomes HEAD
 *
 *   Revision 1.1.2.2  2006/06/25 17:01:57  martius
 *   remove old simulations
 *   robots get names
 *
 *   Revision 1.1.2.1  2006/06/10 20:12:45  martius
 *   simulation for uwo (unknown walking object)
 *
 *
 ***************************************************************************/
#include <stdio.h>

// include ode library
#include <ode/ode.h>

// include noisegenerator (used for adding noise to sensorvalues)
#include <selforg/noisegenerator.h>

// include simulation environment stuff
#include "simulation.h"

// include agent (class for holding a robot, a controller and a wiring)
#include "odeagent.h"

// used wiring
#include <selforg/one2onewiring.h>

// used robot
#include "vierbeiner.h"

// used arena
#include "playground.h"
// used passive spheres
#include "passivesphere.h"
#include "joint.h"

// used controller
//#include <selforg/invertnchannelcontroller.h>
#include <selforg/invertmotornstep.h>
#include <selforg/sinecontroller.h>

// fetch all the stuff of lpzrobots into scope
using namespace lpzrobots;


class ThisSim : public Simulation {
public:


  Joint* fixator;

  // starting function (executed once at the beginning of the simulation loop)
  void start(const OdeHandle& odeHandle, const OsgHandle& osgHandle, GlobalData& global) 
  {
    setCameraHomePos(Pos(1.53837, 4.73003, 1.27411),  Pos(154.844, -9.01605, 0));
    // initialization
    // - set noise to 0.0
    // - register file chess.ppm as a texture called chessTexture (used for the wheels)
    global.odeConfig.setParam("controlinterval",4);
    global.odeConfig.setParam("noise",0.05);
    global.odeConfig.setParam("realtimefactor",4);
    //    global.odeConfig.setParam("gravity", 0);
    //    global.odeConfig.setParam("cameraspeed", 250);
    //  int chessTexture = dsRegisterTexture("chess.ppm");

    // use Playground as boundary:
    // - create pointer to playground (odeHandle contains things like world and space the 
    //   playground should be created in; odeHandle is generated in simulation.cpp)
    // - setting geometry for each wall of playground: 
    //   setGeometry(double length, double width, double	height)
    // - setting initial position of the playground: setPosition(double x, double y, double z)
    // - push playground in the global list of obstacles(globla list comes from simulation.cpp)
    Playground* playground = new Playground(odeHandle, osgHandle, osg::Vec3(30, 0.2, 0.5));
    playground->setPosition(osg::Vec3(0,0,0)); // playground positionieren und generieren
    global.obstacles.push_back(playground);

    // add passive spheres as obstacles
    // - create pointer to sphere (with odehandle, osghandle and 
    //   optional parameters radius and mass,where the latter is not used here) )
    // - set Pose(Position) of sphere 
    // - set a texture for the sphere
    // - add sphere to list of obstacles
    for (int i=0; i<= 1/*2*/; i+=2){
      PassiveSphere* s1 = new PassiveSphere(odeHandle, osgHandle, 0.5);
      s1->setPosition(osg::Vec3(-4.5+i*4.5,0,0));
      s1->setTexture("Images/dusty.rgb");
      global.obstacles.push_back(s1);
    }

    VierBeinerConf conf = VierBeiner::getDefaultConf();
    //    conf.frictionGround = 1;
    conf.jointLimit = M_PI/8; 
    conf.legNumber = 4;
    VierBeiner* dog = new VierBeiner(odeHandle, osgHandle,conf, "Dog");
    dog->place(osg::Matrix::translate(0,0,0.4));
    global.configs.push_back(dog);

    Primitive* trunk = dog->getMainPrimitive();
    fixator = new FixedJoint(trunk, global.environment);
    fixator->init(odeHandle, osgHandle);

    // use Nimm4 vehicle as robot:
    // - create pointer to nimm4 (with odeHandle and osg Handle and possible other settings, see nimm4.h)
    // - place robot
    //OdeRobot* vehiInvertMotorSpacecle = new Nimm4(odeHandle, osgHandle);
    //vehicle->place(Pos(0,2,0));

    // create pointer to controller
    // push controller in global list of configurables
    // AbstractController *controller = new SineController();
    InvertMotorNStepConf cc = InvertMotorNStep::getDefaultConf();
    cc.useS=true;
    AbstractController *controller = new InvertMotorNStep(cc);
    controller->setParam("sinerate",50);
    controller->setParam("phaseshift",1);
    controller->setParam("adaptrate",0);
    controller->setParam("epsC",0.05);
    controller->setParam("epsA",0.01);
    controller->setParam("steps",1);
    controller->setParam("s4avg",2);
    controller->setParam("kwta",4);
    controller->setParam("inhibition",0.01);
    
    global.configs.push_back(controller);
  
    // create pointer to one2onewiring
    One2OneWiring* wiring = new One2OneWiring(new ColorUniformNoise(0.1));

    // create pointer to agent
    // initialize pointer with controller, robot and wiring
    // push agent in globel list of agents
    OdeAgent* agent = new OdeAgent(plotoptions);
    agent->init(controller, dog, wiring);
    global.agents.push_back(agent);
  
    showParams(global.configs);
  }

  // add own key handling stuff here, just insert some case values
  virtual bool command(const OdeHandle&, const OsgHandle&, GlobalData& globalData, int key, bool down)
  {
    if (down) { // only when key is pressed, not when released
      switch ( (char) key )
	{
	case 'x': 
	  if(fixator) delete fixator;
	  fixator=0;	 
	  return true;
	  break;
	default:
	  return false;
	  break;
	}
    }
    return false;
  }
};


int main (int argc, char **argv)
{ 
  ThisSim sim;
  return sim.run(argc, argv) ? 0 : 1;

}
 