/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "helper_functions.h"

using std::string;
using std::vector;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  num_particles = 100;  // TODO: Set the number of particles
  std::default_random_engine gen; //Random number engine
  // Normal Gaussian distribution for x, y & theta
  std::normal_distribution<double> dist_x(x, std[0]);
  std::normal_distribution<double> dist_y(y, std[1]);
  std::normal_distribution<double> dist_theta(theta, std[2]);
  // Initialize particles
  for(int i = 0; i < num_particles; ++i)
  {
    Particle p;
    p.id = i+1;
    p.x = dist_x(gen);
    p.y = dist_y(gen);
    p.theta = dist_theta(gen);
    p.weight = 1.0;
    this->weights.push_back(p.weight);
    this->particles.push_back(p);
  }
  this->is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */
  // Predict position and heading angle of each particle using bicycle model
  double pred_x; // Prediction for Position in x coord.
  double pred_y; // Prediction for Position in y coord.
  double pred_theta; // Prediction for Heading Angle
  std::default_random_engine gen; //Random number engine
  for(int i = 0; i < num_particles; ++i)
  {
    pred_x     = this->particles[i].x + velocity/yaw_rate *
                 (sin(this->particles[i].theta + yaw_rate*delta_t)
                 -sin(this->particles[i].theta));
    pred_y     = this->particles[i].y + velocity/yaw_rate *
                 (cos(this->particles[i].theta)
                 -cos(this->particles[i].theta + yaw_rate*delta_t));
    pred_theta = this->particles[i].theta + yaw_rate * delta_t;
    // Normal Gaussian distribution for x, y & theta predictions
    std::normal_distribution<double> dist_x(pred_x, std_pos[0]);
    std::normal_distribution<double> dist_y(pred_y, std_pos[1]);
    std::normal_distribution<double> dist_theta(pred_theta, std_pos[2]);
    // Add random Gaussian noise to predictions
    this->particles[i].x     = dist_x(gen);
    this->particles[i].y     = dist_y(gen);
    this->particles[i].theta = dist_theta(gen);
  }
}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */

  // Perform Nearest Neighbor method for Data Association
  double distance;
  int id;
  for(int i = 0; i < observations.size(); ++i)
  {
    distance = std::numeric_limits<double>::infinity();
    for(int j = 0; j < predicted.size(); j++)
    {
      double new_dist = dist(observations[i].x, observations[i].y,
                             predicted[j].x, predicted[j].y);
      if(new_dist < distance)
      {
        id = predicted[j].id;
        distance = new_dist; 
      }
    }
    observations[i].id = id;
  }
 }

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */
  // For landmarks (in map frame)
  vector<LandmarkObs> landmarks_vec;
  LandmarkObs landmark;
  // For observations (in map frame)
  vector<LandmarkObs> observations_vec;
  LandmarkObs observation;

  // For debug
  vector<int> associations;
  vector<double> sense_x;
  vector<double> sense_y;
  this->weights.clear();
  for(int i = 0; i < this->num_particles; ++i)
  {
    for(int j = 0; j < map_landmarks.landmark_list.size(); ++j)
    {
      // Check sensor range to eleminate unrelated landmarks
      if(std::abs(particles[i].x - map_landmarks.landmark_list[j].x_f) <= sensor_range && 
      std::abs(particles[i].y - map_landmarks.landmark_list[j].y_f) <= sensor_range)
      {
        // Related landmarks in map frame
        landmark.id = map_landmarks.landmark_list[j].id_i;
        landmark.x  = map_landmarks.landmark_list[j].x_f;
        landmark.y  = map_landmarks.landmark_list[j].y_f;
        landmarks_vec.push_back(landmark);
      }
    }
    for(int k = 0; k < observations.size(); ++k)
    {
      // Coordinate transformation from vehicle frame to map frame for observations
      observation.x = particles[i].x + (cos(particles[i].theta) * observations[k].x) - 
                      (sin(particles[i].theta) * observations[k].y);
      observation.y = particles[i].y + (sin(particles[i].theta) * observations[k].x) + 
                      (cos(particles[i].theta) * observations[k].y);
      observations_vec.push_back(observation);
    }
    dataAssociation(landmarks_vec, observations_vec);
    for(int l = 0; l < observations_vec.size(); ++l)
    {
      for(int m = 0; m < landmarks_vec.size(); ++m)
      {
        if(observations_vec[l].id == landmarks_vec[m].id)
        {
          float term1          = 1 / (2*M_PI*std_landmark[0]*std_landmark[1]);
          float term2_1        = pow(observations_vec[l].x - landmarks_vec[m].x, 2) 
                                 / (2*std_landmark[0]*std_landmark[0]);
          float term2_2        = pow(observations_vec[l].y - landmarks_vec[m].y, 2) 
                                 / (2*std_landmark[1]*std_landmark[1]);
          particles[i].weight *= term1 * exp(-(term2_1 + term2_2));
          
          associations.push_back(landmarks_vec[m].id);
          sense_x.push_back(observations_vec[l].x);
          sense_y.push_back(observations_vec[l].y);
          break;
        }
      }
    }
    SetAssociations(particles[i], associations, sense_x, sense_y);
    
    associations.clear();
    sense_x.clear();
    sense_y.clear();
    landmarks_vec.clear();
    observations_vec.clear();
    this->weights.push_back(particles[i].weight);
  }  
}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
  // Create the distribution with particle weights
  std::discrete_distribution<> d(weights.begin(), weights.end());
  std::default_random_engine gen; //Random number engine

  // Create resampled particles
  vector<Particle> particles_resampled;

  // Resampling
  for(int i = 0; i < this->num_particles; ++i)
  {
    int id = d(gen);
    particles_resampled.push_back(this->particles[id]);
  }

  this->particles = particles_resampled;

  particles_resampled.clear();

}

void ParticleFilter::SetAssociations(Particle& particle, 
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}