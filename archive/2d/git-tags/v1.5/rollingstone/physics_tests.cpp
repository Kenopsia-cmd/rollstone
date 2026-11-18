#include "simulation.hpp"
#include <iostream>
#include <cstdlib>
using namespace rock;
void require(bool test,const char* message){if(!test){std::cerr<<"FAIL: "<<message<<'\n';std::exit(1);}}
Config flat(){Config c;c.terrain.points={{0,0},{120,0}};return c;}
int main(){
    Config c;auto base=simulate(c);auto repeated=simulate(c);
    require(base.reached,"base run reaches toe collector");
    require(base.runoutX>c.terrain.toeX()+1,"toe is not an artificial stopping boundary");
    require(base.detectorHeight==repeated.detectorHeight,"deterministic repeated run");
    Config fine=c;fine.dt=.001;auto refined=simulate(fine);
    require(std::abs(base.detectorHeight-refined.detectorHeight)<.15,"collector height converges on dt halving");
    require(std::abs(base.detectorSpeed-refined.detectorSpeed)/base.detectorSpeed<.03,"collector speed converges on dt halving");
    double initialE=base.trajectory.front().energy+base.finalState.mass*G*base.trajectory.front().z;
    double largestE=initialE;
    for(const auto& point:base.trajectory)largestE=std::max(largestE,point.energy+base.finalState.mass*G*point.z);
    require(largestE<=initialE*1.005,"passive terrain does not materially create mechanical energy");
    c=flat();c.initialSpeed=2.;c.material.rollResistance=.1;
    auto stopping=simulate(c);require(stopping.status=="stopped","flat terrain with resistance comes to rest");
    require(stopping.finalState.v.x==0,"stationary endpoint has zero velocity");
    c=flat();c.material={0,0,0,0};State s=initialState(c);s.p.z=30;s.v={0,0};s.omega=0;
    for(int i=0;i<100;++i)step(s,c);
    double t=100*c.dt;require(std::abs(s.v.z+G*t)<1e-10,"free-flight gravity integration");
    require(std::abs(s.p.z-(30-.5*G*t*t))<.003,"free-flight position integration");
    c=Config();c.shape=Shape::Polygon;c.sides=5;
    auto polygon=simulate(c);State p=initialState(c);
    require(std::abs(p.inertia/(p.mass*c.radius*c.radius)-(2+std::cos(2*PI/5))/6)<1e-12,"actual polygon polar inertia");
    double bottom0=supportRadius(c,0,{0,1}),bottom1=supportRadius(c,.2,{0,1});
    require(std::abs(bottom0-bottom1)>.001,"polygon contact support rotates with shape");
    require(std::abs(polygon.runoutX-base.runoutX)>.01||std::abs(polygon.detectorSpeed-base.detectorSpeed)>.01,"polygon dynamics differ from decorative sphere");
    std::cout<<"PASS sphere and polygon physics tests\n";
}
