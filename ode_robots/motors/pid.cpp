/***************************************************************************
 *   Copyright (C) 2005-2011 LpzRobots development team                    *
 *    Georg Martius  <georg dot martius at web dot de>                     *
 *    Frank Guettler <guettler at informatik dot uni-leipzig dot de        *
 *    Frank Hesse    <frank at nld dot ds dot mpg dot de>                  *
 *    Ralf Der       <ralfder at mis dot mpg dot de>                       *
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
 ***************************************************************************/
#include <ode-dbl/ode.h>
#include <iostream>


#include "pid.h"
using namespace std;

namespace lpzrobots {
/*
  PID::PID ( double KP , double KI , double KD)
    : KP(KP), KI(KI), KD(KD)
  {
    P=D=integrator=I=0;

    targetposition = 0;
    derivative     = 0;
    position       = 0;
    lastposition   = 0;
    error          = 0;
    tau            = 1000;
    lasttime       = -1;
    last2position  = 0;
    lasterror      = 0;
    force          = 0;
  }
*/
  PID::PID ( double KPpos, double KIpos, double KDpos, double KPvel, double KIvel, double KDvel )
      : KPpos(KPpos), KIpos(KIpos), KDpos(KDpos), KPvel(KPvel), KIvel(KIvel), KDvel(KDvel)
    {
	  targetposition_pos = 0;
	  targetposition_vel = 0;
      derivative_inner = 0;
      derivative_outer = 0;
      integrator_inner = 0;
      integrator_outer = 0;
      vel = 0;
      pos = 0;
      lastvel = 0;
      lastpos = 0;
      last2vel = 0;
      last2pos = 0;
      error_pos = 0;
      error_vel = 0;
      lasterror_pos = 0;
      lasterror_vel = 0;

      tau            = 1000;
      lasttime       = -1;
      force          = 0;
    }

  void PID::setKP(double KP){
	  this->KP = KP;
	}
  void PID::setKI(double KI){
      this->KI = KI;
    }
  void PID::setKD(double KD){
	  this->KD = KD;
	}

  void PID::setTargetPosition ( double newpos )
  {
    targetposition = newpos;
  }

  double PID::getTargetPosition ( )
  {
    return targetposition;
  }

  /* This is the old implementation. Do not change in order not to brake all
      the simulations
   */
  double PID::step ( double newsensorval, double time)
  {
    if(lasttime != -1 && time - lasttime > 0 ){

    	lastposition = position;
    	position = newsensorval;
    	double stepsize=time-lasttime;

		error = targetposition - position;
		derivative = (error - lasterror) / stepsize;
		integrator += stepsize * error;

		if( integrator > 10 )
		{
			integrator = 10;
		}
		else if( integrator < -10 )
		{
			integrator = -10;
		}

		force = (error*KP) + (derivative*KD) + (integrator*KI);

		/*int limits = 10;

	    if( force > limits )
	    	force = limits;
	    else if( force < -limits )
	    	force = -limits;*/

    } else
    {
    	force=0;
    }

    lasttime=time;
    lasterror = error;
    //cout << "Kp " << KP << endl;
    //cout << "Ki " << KI << endl;
    //cout << "Kd " << KD << endl;
    return force;
  }

  // This is the new implementation used by the center and velocity servos
  double PID::stepNoCutoff ( double newsensorval, double time)
  {
    if(lasttime != -1 && time - lasttime > 0 ){
      lastposition = position;
      position = newsensorval;
      double stepsize=time-lasttime;

      lasterror = error;
      error = targetposition - position;
      derivative += ((lasterror - error) / stepsize - derivative)*0.2; // Georg: Who put the 0.2 here!?

      P = error;
      I *= (1-1/tau);
      I += stepsize * error * KI;
      D = -derivative * KD;
      force = KP*(P + I + D);
    } else {
      force=0;
    }
    lasttime=time;
    return force;
  }

  // This is the new implementation used for the velocity servos (velocity control)
  // no I term and velocity is bound such that we cannot overshoot in one step
  double PID::stepVelocity ( double newsensorval, double time)
  {
    // force is here a nominal velocity

    if(lasttime != -1 && time - lasttime > 0 )
    {
      lastposition = position;
      position = newsensorval;
      double stepsize=time-lasttime;

      lasterror = error;
      error = targetposition - position;

      P = error;
      if(KD!=0.0){
        derivative += ((lasterror - error) / stepsize - derivative);
        D = -derivative * KD;
        force = KP*(P + D);
      } else
        force = KP*P;
      // limit the velocity
      if(stepsize*fabs(force) > fabs(error))
      {
        force = error/stepsize;
      }
    } else
    {
      force=0;
    }

    lasttime=time;

    return force;
  }
  double PID::stepPositionVelocity ( double newsensorval, double motorvel, double time)
  {
	double stepsize = time - lasttime;

	//	The force is here a nominal velocity
	if( lasttime != -1 && time - lasttime > 0 )
	{
		//	The outer loop.
		lastpos = pos;
		pos = newsensorval;

		error_pos = targetposition_pos - pos;
		derivative_outer = (error_pos - lasterror_pos) / stepsize;
		integrator_outer += stepsize * error_pos;

		if( integrator_outer > 10 )
		{
			integrator_outer = 10;
		}
		else if( integrator_outer < -10 )
		{
			integrator_outer = -10;
		}

		double targetvelocity = (error_pos*KPpos) + (derivative_outer*KDpos) + (integrator_outer*KIpos);

		//	The inner loop
		lastvel = vel;
		vel = motorvel;

		error_vel = targetvelocity - vel;
		derivative_inner = (error_vel - lasterror_vel) / stepsize;
		integrator_inner += stepsize * error_vel;

		if( integrator_inner > 10 )
		{
			integrator_inner = 10;
		}
		else if( integrator_inner < -10 )
		{
			integrator_inner = -10;
		}

		force = (error_vel*KPvel) + (derivative_inner*KDvel) + (integrator_inner*KIvel);
	}
	else
	{
		force = 0;
	}

	lasttime = time;
	lasterror_pos = error_pos;
	lasterror_vel = error_vel;

	return force;
  }

}
