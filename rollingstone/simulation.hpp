#pragma once
#include "model.hpp"
#include <limits>

namespace rock {
inline Terrain::Terrain():points{{0,55.6},{10,51.9},{11.25,50.6},{30,35.6},
    {40,34.4},{50,23.1},{60,20.6},{70,15.6},{80,10},{90,8.75},{100,2.88},{120,0}}{}
inline void checkMaterial(const Material& m){
    if(!std::isfinite(m.rollResistance)||!std::isfinite(m.friction)||
       !std::isfinite(m.restitution)||!std::isfinite(m.drag)||m.rollResistance<0||
       m.rollResistance>1||m.friction<0||m.friction>2||m.restitution<0||m.restitution>1||m.drag<0||m.drag>10)
        throw std::invalid_argument("Invalid surface material (finite coefficients required)");
}
inline void Terrain::validate()const{
    if(points.size()<2)throw std::invalid_argument("Terrain needs at least two vertices");
    for(size_t i=0;i<points.size();++i){
        if(!std::isfinite(points[i].x)||!std::isfinite(points[i].z))throw std::invalid_argument("Nonfinite terrain coordinate");
        if(i&&points[i].x<=points[i-1].x)throw std::invalid_argument("Terrain x coordinates must increase strictly");
    }
    if(!surfaces.empty()&&surfaces.size()!=points.size()-1)throw std::invalid_argument("One surface material required per terrain segment");
    for(const auto& m:surfaces)checkMaterial(m);
}
inline size_t Terrain::segment(double x)const{
    if(x<=points.front().x)return 0;
    for(size_t i=0;i+1<points.size();++i)if(x<points[i+1].x)return i;
    return points.size()-2;
}
inline double Terrain::height(double x)const{
    if(x<=points.front().x)return points.front().z;
    if(x>=points.back().x)return points.back().z;
    size_t i=segment(x);return points[i].z+(points[i+1].z-points[i].z)*(x-points[i].x)/(points[i+1].x-points[i].x);
}
inline double Terrain::angle(double x)const{
    if(x>=points.back().x)return 0;
    size_t i=segment(x);return -std::atan2(points[i+1].z-points[i].z,points[i+1].x-points[i].x);
}
inline void validate(const Config& c){
    c.terrain.validate();checkMaterial(c.material);
    auto finite=[](double x){return std::isfinite(x);};
    if(!finite(c.radius)||c.radius<.05||c.radius>10||!finite(c.initialSpeed)||c.initialSpeed<0||c.initialSpeed>200||
       !finite(c.density)||c.density<=0||c.density>100000||!finite(c.dt)||c.dt<.00005||c.dt>.01||
       !finite(c.maxTime)||c.maxTime<=0||c.maxTime>600||!finite(c.cv)||c.cv<0||c.cv>1||
       !finite(c.binHeight)||c.binHeight<=0||c.binHeight>10||!finite(c.gridWidth)||c.gridWidth<.05||
       !finite(c.sampleInterval)||c.sampleInterval<c.dt||!finite(c.detectorX)||
       !finite(c.shapePhase))throw std::invalid_argument("Invalid simulation configuration or integration step");
    if(c.detectorX>=0&&c.detectorX<=c.terrain.points.front().x)throw std::invalid_argument("Collector must be downstream of release");
    if(c.sides<3||c.sides>12)throw std::invalid_argument("Convex polygon needs 3 to 12 sides");
    if(c.barrier.kind!=BarrierKind::None)throw std::invalid_argument("This version does not yet support barriers");
}
inline std::vector<Vec2> bodyVertices(const Config& c,double angle){
    std::vector<Vec2> vertices;
    int count=c.shape==Shape::Sphere?48:c.sides;
    for(int i=0;i<count;++i){double a=angle+2*PI*i/count;vertices.push_back({c.radius*std::cos(a),c.radius*std::sin(a)});}
    return vertices;
}
inline double supportRadius(const Config& c,double angle,Vec2 normal){
    if(c.shape==Shape::Sphere)return c.radius;
    double support=0;for(Vec2 vertex:bodyVertices(c,angle))support=std::max(support,-dot(vertex,normal));
    return support;
}
inline double kineticEnergy(const State& s){return .5*s.mass*dot(s.v,s.v)+.5*s.inertia*s.omega*s.omega;}
inline State initialState(const Config& c){
    validate(c);State s;
    s.mass=c.density*4./3.*PI*c.radius*c.radius*c.radius;
    s.inertia=.4*s.mass*c.radius*c.radius;
    // Polygon is a uniform 2D lamina with the same mass as the reference sphere.
    // Its true polar inertia depends on side count and circumradius.
    if(c.shape==Shape::Polygon)s.inertia=s.mass*c.radius*c.radius*(2+std::cos(2*PI/c.sides))/6;
    s.angle=c.shapePhase;
    double x=c.terrain.points.front().x,th=c.terrain.angle(x);
    s.p={x,c.terrain.height(x)+c.radius/std::cos(th)};
    if(c.shape==Shape::Polygon){
        s.p.z=-std::numeric_limits<double>::infinity();
        for(Vec2 vertex:bodyVertices(c,s.angle))s.p.z=std::max(s.p.z,c.terrain.height(x+vertex.x)-vertex.z);
    }
    s.v={c.initialSpeed*std::cos(th),-c.initialSpeed*std::sin(th)};
    s.omega=-c.initialSpeed/c.radius;
    return s;
}
inline Vec2 rotationalVelocity(double omega,Vec2 arm){return {-omega*arm.z,omega*arm.x};}
inline void impulse(State& s,Vec2 arm,Vec2 normal,double magnitude){
    s.v=s.v+normal*(magnitude/s.mass);s.omega+=cross(arm,normal)*magnitude/s.inertia;
}
inline double effectiveInverseMass(const State& s,Vec2 arm,Vec2 direction){
    double lever=cross(arm,direction);return 1/s.mass+lever*lever/s.inertia;
}
inline void step(State& s,const Config& c){
    if(!s.moving)return;
    const double dt=c.dt;
    const Material& m=c.terrain.surfaces.empty()?c.material:c.terrain.surfaces[c.terrain.segment(s.p.x)];
    double kd=.5*1.225*m.drag*PI*c.radius*c.radius/s.mass;
    double speed=length(s.v);s.v=s.v*(1/(1+kd*speed*dt));
    s.v.z-=G*dt;s.p=s.p+s.v*dt;s.angle+=s.omega*dt;s.time+=dt;
    double th=c.terrain.angle(s.p.x);Vec2 n{std::sin(th),std::cos(th)};
    struct Contact{Vec2 arm,normal;double target=0,jn=0,jt=0;Material mat;};
    std::vector<Contact> contacts;
    if(c.shape==Shape::Sphere){
        double gap=(s.p.z-c.terrain.height(s.p.x))*n.z-c.radius;
        if(gap<=1e-6){
            s.p=s.p+n*std::max(0.,-gap);
            contacts.push_back({n*(-c.radius),n,0,0,0,m});
        }
    }else{
        double correction=0;
        for(Vec2 arm:bodyVertices(c,s.angle)){
            Vec2 p=s.p+arm;double angle=c.terrain.angle(p.x);Vec2 normal{std::sin(angle),std::cos(angle)};
            double verticalGap=p.z-c.terrain.height(p.x);
            if(verticalGap<=1e-5){
                Material material=c.terrain.surfaces.empty()?c.material:c.terrain.surfaces[c.terrain.segment(p.x)];
                contacts.push_back({arm,normal,0,0,0,material});
                correction=std::max(correction,-verticalGap);
            }
        }
        s.p.z+=correction;
    }
    for(auto& contact:contacts){
        double vn=dot(s.v+rotationalVelocity(s.omega,contact.arm),contact.normal);
        contact.target=vn<-.25?-contact.mat.restitution*vn:0;
    }
    // Accumulated projected impulses resolve all contacting polygon vertices.
    // Restitution target is frozen per impact, so iterations do not add bounce energy.
    for(int iteration=0;iteration<8;++iteration)for(auto& contact:contacts){
        Vec2 normal=contact.normal,tangent{normal.z,-normal.x};
        double vn=dot(s.v+rotationalVelocity(s.omega,contact.arm),normal);
        double oldJn=contact.jn;
        contact.jn=std::max(0.,oldJn+(contact.target-vn)/effectiveInverseMass(s,contact.arm,normal));
        impulse(s,contact.arm,normal,contact.jn-oldJn);
        double slip=dot(s.v+rotationalVelocity(s.omega,contact.arm),tangent);
        double oldJt=contact.jt;
        contact.jt=std::clamp(oldJt-slip/effectiveInverseMass(s,contact.arm,tangent),-contact.mat.friction*contact.jn,contact.mat.friction*contact.jn);
        impulse(s,contact.arm,tangent,contact.jt-oldJt);
    }
    if(!contacts.empty()){
        double dw=m.rollResistance*s.mass*G*n.z*c.radius*dt/s.inertia;
        s.omega=std::copysign(std::max(0.,std::abs(s.omega)-dw),s.omega);
        bool canRest=c.shape==Shape::Polygon||std::abs(std::tan(th))<=m.rollResistance+1e-8;
        if(length(s.v)<.05&&std::abs(s.omega)*c.radius<.05&&canRest)s.restTime+=dt;else s.restTime=0;
        if(s.restTime>.5){s.v={0,0};s.omega=0;s.moving=false;s.stopReason="stopped";}
    }else s.restTime=0;
    if(s.p.x<c.terrain.points.front().x-20||s.p.x>std::max(c.terrain.toeX()+80.,c.detectorX+30.)){
        s.moving=false;s.stopReason="domain_exit";
    }
    if(!std::isfinite(s.p.x)||!std::isfinite(s.p.z)||!std::isfinite(s.omega))throw std::runtime_error("Nonfinite physics state");
}
inline Sample sample(const State& s){return {s.time,s.p.x,s.p.z,s.v.x,s.v.z,s.angle,s.omega,kineticEnergy(s),s.barrierDeflection};}
inline RunResult simulate(const Config& c){
    validate(c);RunResult result;result.config=c;State s=initialState(c);
    result.trajectory.push_back(sample(s));double nextSample=c.sampleInterval;
    double detector=c.detectorX<0?c.terrain.toeX():c.detectorX;
    result.runoutX=s.p.x;result.maxSpeed=length(s.v);result.maxEnergy=kineticEnergy(s);
    const int steps=static_cast<int>(std::ceil(c.maxTime/c.dt));
    for(int i=0;i<steps&&s.moving;++i){
        State previous=s;step(s,c);
        result.runoutX=std::max(result.runoutX,s.p.x);
        result.maxSpeed=std::max(result.maxSpeed,length(s.v));
        result.maxEnergy=std::max(result.maxEnergy,kineticEnergy(s));
        result.maxBounceHeight=std::max(result.maxBounceHeight,std::max(0.,s.p.z-c.terrain.height(s.p.x)-supportRadius(c,s.angle,{0,1})));
        if(!result.reached&&previous.p.x<detector&&s.p.x>=detector){
            double f=(detector-previous.p.x)/(s.p.x-previous.p.x);
            result.reached=true;result.detectorHeight=previous.p.z+(s.p.z-previous.p.z)*f-c.terrain.toeZ();
            Vec2 crossingVelocity=previous.v+(s.v-previous.v)*f;
            double crossingOmega=previous.omega+(s.omega-previous.omega)*f;
            result.detectorSpeed=length(crossingVelocity);
            result.detectorEnergy=.5*s.mass*dot(crossingVelocity,crossingVelocity)+.5*s.inertia*crossingOmega*crossingOmega;
            result.detectorTime=previous.time+c.dt*f;
            result.inBand=result.detectorHeight>=0&&result.detectorHeight<=10;
        }
        if(s.time+1e-9>=nextSample){result.trajectory.push_back(sample(s));nextSample+=c.sampleInterval;}
    }
    if(result.trajectory.back().t<s.time-1e-9)result.trajectory.push_back(sample(s));
    if(s.moving){s.stopReason="time_limit";s.moving=false;}
    result.status=s.stopReason;result.finalState=s;return result;
}
}
