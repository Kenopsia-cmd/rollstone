// 山区滚石模拟：单文件C++，支持Visual Studio C++17/C++20 + EasyX。
// 无参数运行直接打开画布；所有项目自定义代码已内联。
// B 批量模拟，E 导出CSV，A 参数研究，G 地形研究，M 多颗回放，F 能量曲线，O 输出目录。
// 左侧实验/地表/地形/形状/拦挡页调参数，调整后按B重算。Space暂停，R重播。
// 命令行 --help 查看本版本支持的批量、研究和参数选项。
// 图形、回放、概率和CSV使用同一计算数据。只需把本文件作为工程唯一主.cpp。
// 这是二维教学近似模型，未作工程标定。
// 字体与布局清晰度修正；保留字体复用、轨迹抽稀和动态研究表校验。
#ifndef NOMINMAX
#define NOMINMAX
#endif
#ifndef UNICODE
#define UNICODE
#endif
#ifndef _UNICODE
#define _UNICODE
#endif
#ifdef _MSC_VER
#pragma execution_character_set("utf-8")
#pragma comment(lib,"user32.lib")
#pragma comment(lib,"gdi32.lib")
#pragma comment(lib,"shell32.lib")
#endif
#define ROCK_STAGE 8
#define ROCK_VERSION "combined"

// ==================== model.hpp ====================

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <map>
#include <random>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace rock {
constexpr double PI=3.14159265358979323846;
constexpr double G=9.81;
struct Vec2 { double x=0,z=0; };
inline Vec2 operator+(Vec2 a,Vec2 b){return {a.x+b.x,a.z+b.z};}
inline Vec2 operator-(Vec2 a,Vec2 b){return {a.x-b.x,a.z-b.z};}
inline Vec2 operator*(Vec2 a,double b){return {a.x*b,a.z*b};}
inline Vec2 operator/(Vec2 a,double b){return {a.x/b,a.z/b};}
inline double dot(Vec2 a,Vec2 b){return a.x*b.x+a.z*b.z;}
inline double cross(Vec2 a,Vec2 b){return a.x*b.z-a.z*b.x;}
inline double length(Vec2 a){return std::hypot(a.x,a.z);}
struct Material { double rollResistance=.02,friction=.45,restitution=.45,drag=1.; };
inline Material materialPreset(int preset){
    if(preset==2)return {.08,.60,.25,1.};
    if(preset==3)return {.12,.75,.20,1.};
    return {.02,.45,.45,1.};
}
enum class Shape { Sphere, Polygon };
enum class BarrierKind { None, Rigid, Elastic };
struct Barrier {
    BarrierKind kind=BarrierKind::None;
    double x=112.,height=6.,restitution=.15,stiffness=2e6,damping=2e5,maxDeflection=3.;
    double mass=1000.; // equivalent translating barrier mass, kg
};
struct Terrain {
    std::vector<Vec2> points;
    std::vector<Material> surfaces; // empty: Config::material; otherwise one material per segment
    Terrain();
    void validate() const;
    size_t segment(double x) const;
    double height(double x) const;
    double angle(double x) const;
    double toeX() const {return points.back().x;}
    double toeZ() const {return points.back().z;}
};
struct Config {
    Terrain terrain;
    Material material;
    double radius=2.,initialSpeed=1.5,density=2600.,dt=.002,maxTime=60.,cv=.1;
    double detectorX=-1.,binHeight=.5,gridWidth=1.,sampleInterval=.02;
    double terrainWidthScale=1.,terrainHeightScale=1.;
    bool randomSurface=false;
    Shape shape=Shape::Sphere;
    int sides=5;
    double shapePhase=0.;
    Barrier barrier;
};
struct State {
    Vec2 p,v;
    double angle=0,omega=0,mass=0,inertia=0,time=0;
    bool moving=true;
    double barrierDeflection=0,barrierDissipation=0,barrierVelocity=0,restTime=0;
    double maxBarrierDeflection=0,peakBarrierForce=0;
    double barrierForce=0,contactImpulse=0; // spring/damper force N; cumulative barrier normal impulse N s
    bool barrierHit=false,barrierFailed=false;
    bool barrierPassed=false,barrierOvertopped=false,barrierIntercepted=false;
    uint64_t integrationSteps=0; // actual legacy steps or adaptive collision substeps
    std::string stopReason;
};
struct Sample {
    double t=0,x=0,z=0,vx=0,vz=0,angle=0,omega=0,energy=0,deflection=0;
    bool barrierHit=false,barrierFailed=false;
    double barrierVelocity=0,barrierForce=0,barrierDissipation=0,contactImpulse=0;
};
struct RunResult {
    Config config;
    std::vector<Sample> trajectory;
    State finalState;
    bool reached=false,inBand=false;
    double detectorHeight=-1,detectorSpeed=0,detectorEnergy=0,detectorTime=0;
    double runoutX=0,maxSpeed=0,maxEnergy=0,maxBounceHeight=0;
    std::string status;
};
struct BatchResult {
    Config base;
    uint32_t seed=42;
    std::vector<RunResult> runs;
    std::vector<int> heightCounts;
    std::map<std::pair<int,int>,int> routeCells;
    size_t reached=0,inBand=0,aboveBand=0,belowBand=0,notReached=0;
};
void validate(const Config& c);
std::vector<Vec2> bodyVertices(const Config& c,double angle);
double supportRadius(const Config& c,double angle,Vec2 normal);
double kineticEnergy(const State& s);
State initialState(const Config& c);
void step(State& s,const Config& c);
RunResult simulate(const Config& c);
Config randomize(const Config& base,std::mt19937& rng);
BatchResult runBatch(const Config& base,int count,uint32_t seed);
void exportBatch(const BatchResult& batch,const std::string& outputDir);
std::string batchSummary(const BatchResult& batch);
int runParameterStudy(const Config& base,int count,uint32_t seed,const std::string& outputDir);
int runTerrainStudy(const Config& base,int count,uint32_t seed,const std::string& outputDir);
}


// ==================== simulation.hpp ====================

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
       !finite(c.binHeight)||c.binHeight<.05||c.binHeight>10||!finite(c.gridWidth)||c.gridWidth<.05||c.gridWidth>1000||
       !finite(c.sampleInterval)||c.sampleInterval<c.dt||!finite(c.detectorX)||
       !finite(c.shapePhase))throw std::invalid_argument("Invalid simulation configuration or integration step");
    if(c.detectorX>10000 || (c.detectorX<0 && c.detectorX!=-1))throw std::invalid_argument("Collector x must be -1 (automatic) or within the domain limit of 10000m");
    if(c.detectorX>=0&&c.detectorX<=c.terrain.points.front().x)throw std::invalid_argument("Collector must be downstream of release");
    if(c.shape!=Shape::Sphere&&c.shape!=Shape::Polygon)throw std::invalid_argument("Unknown rock shape");
#if ROCK_STAGE < 5
    if(c.shape!=Shape::Sphere)throw std::invalid_argument("Polygon support starts at task 5");
#endif
    if(c.sides<3||c.sides>12)throw std::invalid_argument("Regular convex polygon needs 3 to 12 sides");
    if(c.barrier.kind!=BarrierKind::None&&c.barrier.kind!=BarrierKind::Rigid&&c.barrier.kind!=BarrierKind::Elastic)
        throw std::invalid_argument("Unknown barrier kind");
#if ROCK_STAGE < 6
    if(c.barrier.kind!=BarrierKind::None)throw std::invalid_argument("Barrier support starts at task 6");
#elif ROCK_STAGE < 7
    if(c.barrier.kind==BarrierKind::Elastic)throw std::invalid_argument("Elastic barrier support starts at task 7");
#endif
    const auto& b=c.barrier;
    if(b.kind!=BarrierKind::None&&(!finite(b.x)||b.x<=c.terrain.points.front().x||b.x>10000||
       !finite(b.height)||b.height<=0||b.height>100||!finite(b.restitution)||b.restitution<0||b.restitution>1))
        throw std::invalid_argument("Invalid finite-height barrier");
    if(b.kind==BarrierKind::Elastic&&(!finite(b.mass)||b.mass<1||b.mass>1e8||!finite(b.stiffness)||
       b.stiffness<=0||b.stiffness>1e10||!finite(b.damping)||b.damping<0||b.damping>1e10||
       !finite(b.maxDeflection)||b.maxDeflection<=0||b.maxDeflection>100))
        throw std::invalid_argument("Invalid elastic panel mass, stiffness, damping or displacement capacity");
}
inline std::vector<Vec2> bodyVertices(const Config& c,double angle){
    std::vector<Vec2> vertices;
    const int count=c.shape==Shape::Sphere?48:c.sides;
    for(int i=0;i<count;++i){double a=angle+2*PI*i/count;vertices.push_back({c.radius*std::cos(a),c.radius*std::sin(a)});}
    return vertices;
}
inline double supportRadius(const Config& c,double angle,Vec2 normal){
    if(c.shape==Shape::Sphere)return c.radius;
    double support=0;for(const auto& vertex:bodyVertices(c,angle))support=std::max(support,-dot(vertex,normal));
    return support;
}
inline double kineticEnergy(const State& s){return .5*s.mass*dot(s.v,s.v)+.5*s.inertia*s.omega*s.omega;}
inline State initialState(const Config& c){
    validate(c);State s;
    s.mass=c.density*4./3.*PI*c.radius*c.radius*c.radius;
    s.inertia=.4*s.mass*c.radius*c.radius;
    // A regular polygon is a uniform lamina with the reference sphere's mass.
    if(c.shape==Shape::Polygon)s.inertia=s.mass*c.radius*c.radius*(2+std::cos(2*PI/c.sides))/6;
    s.angle=c.shapePhase;
    double x=c.terrain.points.front().x,th=c.terrain.angle(x);
    s.p={x,c.terrain.height(x)+c.radius/std::cos(th)};
    if(c.shape==Shape::Polygon){
        s.p.z=-std::numeric_limits<double>::infinity();
        const auto vertices=bodyVertices(c,s.angle);
        for(const auto& a:vertices)s.p.z=std::max(s.p.z,c.terrain.height(x+a.x)-a.z);
        // Include terrain vertices under lower polygon edges, not only rock vertices.
        for(size_t i=0;i<vertices.size();++i){const auto a=vertices[i],b=vertices[(i+1)%vertices.size()];
            if(b.x<=a.x+1e-12)continue;
            for(const auto& p:c.terrain.points)if(p.x>=x+a.x&&p.x<=x+b.x){
                const double f=(p.x-x-a.x)/(b.x-a.x);
                s.p.z=std::max(s.p.z,p.z-(a.z+(b.z-a.z)*f));
            }
        }
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
inline double simulationDomainRight(const Config& c){
    double right=std::max(c.terrain.toeX()+80.,c.detectorX+30.);
    if(c.barrier.kind!=BarrierKind::None)right=std::max(right,c.barrier.x+
        (c.barrier.kind==BarrierKind::Elastic?c.barrier.maxDeflection:0)+2*c.radius+30.);
    return right;
}
inline void stepSphereLegacy(State& s,const Config& c){
    if(!s.moving)return;
    ++s.integrationSteps;
    const double dt=c.dt;
    const Material& m=c.terrain.surfaces.empty()?c.material:c.terrain.surfaces[c.terrain.segment(s.p.x)];
    double kd=.5*1.225*m.drag*PI*c.radius*c.radius/s.mass;
    double speed=length(s.v);s.v=s.v*(1/(1+kd*speed*dt));
    s.v.z-=G*dt;s.p=s.p+s.v*dt;s.angle+=s.omega*dt;s.time+=dt;
    double th=c.terrain.angle(s.p.x);Vec2 n{std::sin(th),std::cos(th)},t{n.z,-n.x};
    double gap=(s.p.z-c.terrain.height(s.p.x))*n.z-c.radius;
    bool contact=gap<=1e-6;
    if(contact){
        s.p=s.p+n*std::max(0.,-gap);
        Vec2 arm=n*(-c.radius);
        double vn=dot(s.v+rotationalVelocity(s.omega,arm),n);
        double normalImpulse=0;
        if(vn<0){
            double e=vn<-.25?m.restitution:0.;
            normalImpulse=-(1+e)*vn/effectiveInverseMass(s,arm,n);
            impulse(s,arm,n,normalImpulse);
        }
        double slip=dot(s.v+rotationalVelocity(s.omega,arm),t);
        double jt=-slip/effectiveInverseMass(s,arm,t);
        jt=std::clamp(jt,-m.friction*normalImpulse,m.friction*normalImpulse);
        impulse(s,arm,t,jt);
        // Rolling resistance is a torque opposing rotation, never a forward force.
        double dw=m.rollResistance*s.mass*G*n.z*c.radius*dt/s.inertia;
        s.omega=std::copysign(std::max(0.,std::abs(s.omega)-dw),s.omega);
        if(length(s.v)<.05&&std::abs(s.omega)*c.radius<.05&&std::abs(std::tan(th))<=m.rollResistance+1e-8){
            s.v={0,0};s.omega=0;s.moving=false;s.stopReason="stopped";
        }
    }
    if(s.p.x<c.terrain.points.front().x-20||s.p.x>std::max(c.terrain.toeX()+80.,c.detectorX+30.)){
        s.moving=false;s.stopReason="domain_exit";
    }
    if(!std::isfinite(s.p.x)||!std::isfinite(s.p.z)||!std::isfinite(s.omega))throw std::runtime_error("Nonfinite physics state");
}
namespace detail {
constexpr double contactSlop=1e-8;
struct Contact {Vec2 arm,normal;double depth=0,jn=0,jt=0,target=0;Material material;};
inline Vec2 unit(Vec2 a,Vec2 fallback={0,1}){const double d=length(a);return d>1e-12?a/d:fallback;}
inline Material surface(const Config& c,double x){return c.terrain.surfaces.empty()?c.material:c.terrain.surfaces[c.terrain.segment(x)];}
inline Vec2 groundNormal(const Terrain& t,double x){
    if(x<t.points.front().x||x>=t.points.back().x)return {0,1};
    const double a=t.angle(x);return {std::sin(a),std::cos(a)};
}
inline void append(std::vector<Contact>& out,Vec2 arm,Vec2 normal,double depth,const Material& material){
    if(depth < -contactSlop)return;
    for(const auto& q:out)if(length(q.arm-arm)<1e-7&&dot(q.normal,normal)>.999999)return;
    out.push_back({arm,normal,std::max(0.,depth),0,0,0,material});
}
inline std::vector<Contact> groundContacts(const State& s,const Config& c,double& verticalCorrection){
    std::vector<Contact> out;verticalCorrection=0;
    if(c.shape==Shape::Polygon){
        const auto vertices=bodyVertices(c,s.angle);
        for(const auto& arm:vertices){const auto p=s.p+arm;const auto normal=groundNormal(c.terrain,p.x);
            const double overlap=c.terrain.height(p.x)-p.z;
            if(overlap>=-contactSlop){append(out,arm,normal,overlap*normal.z,surface(c,p.x));verticalCorrection=std::max(verticalCorrection,overlap);}
        }
        // A heightfield corner can enter a polygon edge while every rock vertex is clear.
        // Testing all breakpoints of the lower hull gives the exact vertical separation.
        for(size_t i=0;i<vertices.size();++i){const Vec2 a=vertices[i],b=vertices[(i+1)%vertices.size()],edge=b-a;
            if(edge.x<=1e-12)continue;const Vec2 normal=unit({-edge.z,edge.x});
            for(const auto& p:c.terrain.points){
                const double localX=p.x-s.p.x;if(localX<a.x||localX>b.x)continue;
                const double f=(localX-a.x)/edge.x;const Vec2 arm=a+edge*f;
                const double overlap=p.z-s.p.z-arm.z;
                if(overlap>=-contactSlop){append(out,arm,normal,overlap*normal.z,surface(c,p.x));verticalCorrection=std::max(verticalCorrection,overlap);}
            }
        }
    }else{
        // Segment faces plus their endpoints resolve corners for sphere/barrier runs.
        // The historical sphere-without-barrier path remains unchanged above.
        const auto& points=c.terrain.points;
        for(size_t i=0;i+1<points.size();++i){const Vec2 a=points[i],edge=points[i+1]-a;
            const double u=dot(s.p-a,edge)/dot(edge,edge);
            if(u<0||u>1)continue;const Vec2 normal=unit({-edge.z,edge.x});
            append(out,normal*(-c.radius),normal,c.radius-dot(s.p-a,normal),surface(c,.5*(a.x+points[i+1].x)));
        }
        if(s.p.x<=points.front().x)append(out,{0,-c.radius},{0,1},points.front().z+c.radius-s.p.z,surface(c,points.front().x));
        if(s.p.x>=points.back().x)append(out,{0,-c.radius},{0,1},points.back().z+c.radius-s.p.z,surface(c,points.back().x));
        for(const auto& p:points){const Vec2 d=s.p-p;const double distance=length(d);
            if(distance<c.radius+contactSlop&&d.z>=0){const Vec2 normal=unit(d);append(out,normal*(-c.radius),normal,c.radius-distance,surface(c,p.x));}}
    }
    return out;
}
inline void groundImpulse(State& s,std::vector<Contact>& contacts){
    const double before=kineticEnergy(s);
    for(auto& contact:contacts){const double vn=dot(s.v+rotationalVelocity(s.omega,contact.arm),contact.normal);
        contact.target=vn<-.25?-contact.material.restitution*vn:0;}
    for(int iteration=0;iteration<16;++iteration)for(auto& contact:contacts){
        const Vec2 n=contact.normal,t{n.z,-n.x};
        const double vn=dot(s.v+rotationalVelocity(s.omega,contact.arm),n),old=contact.jn;
        contact.jn=std::max(0.,old+(contact.target-vn)/effectiveInverseMass(s,contact.arm,n));
        impulse(s,contact.arm,n,contact.jn-old);
        const double slip=dot(s.v+rotationalVelocity(s.omega,contact.arm),t),oldT=contact.jt;
        contact.jt=std::clamp(oldT-slip/effectiveInverseMass(s,contact.arm,t),-contact.material.friction*contact.jn,contact.material.friction*contact.jn);
        impulse(s,contact.arm,t,contact.jt-oldT);
    }
    // A finite-iteration simultaneous restitution solve must never create kinetic energy.
    // Scaling preserves separating normal velocities and bounds solver error passively.
    const double after=kineticEnergy(s);
    if(after>before*(1+1e-12)&&after>0){const double f=std::sqrt(before/after);s.v=s.v*f;s.omega*=f;}
}
inline bool pointInPolygon(Vec2 p,const std::vector<Vec2>& vertices){
    for(size_t i=0;i<vertices.size();++i)if(cross(vertices[(i+1)%vertices.size()]-vertices[i],p-vertices[i]) < -contactSlop)return false;
    return true;
}
inline std::vector<Contact> barrierContacts(const State& s,const Config& c,double approachSide){
    std::vector<Contact> out;if(c.barrier.kind==BarrierKind::None||s.barrierFailed)return out;
    const double wall=c.barrier.x+(c.barrier.kind==BarrierKind::Elastic?s.barrierDeflection:0);
    const double base=c.terrain.height(c.barrier.x),top=base+c.barrier.height;
    const Material contactMaterial{0,0,c.barrier.restitution,0};
    if(c.shape==Shape::Sphere){
        const Vec2 point{wall,std::clamp(s.p.z,base,top)},delta=s.p-point;
        const double d=length(delta);const Vec2 n=unit(delta,{approachSide,0});
        append(out,n*(-c.radius),n,c.radius-d,contactMaterial);return out;
    }
    const auto vertices=bodyVertices(c,s.angle);
    // Actual vertices against the wall face; clipped edge intersections at the caps
    // are handled separately using the polygon edge normal, never a fake horizontal one.
    for(const auto& arm:vertices){const Vec2 p=s.p+arm;
        if(p.z>base+contactSlop&&p.z<top-contactSlop&&std::abs(s.p.x-wall)<=c.radius+contactSlop)
            append(out,arm,{approachSide,0},-approachSide*(p.x-wall),contactMaterial);
    }
    for(const Vec2 point:std::vector<Vec2>{{wall,base},{wall,top}}){const Vec2 local=point-s.p;
        if(!pointInPolygon(local,vertices))continue;
        double closest=std::numeric_limits<double>::infinity();Vec2 arm,normal;
        for(size_t i=0;i<vertices.size();++i){const Vec2 a=vertices[i],edge=vertices[(i+1)%vertices.size()]-a;
            const Vec2 n=unit({-edge.z,edge.x});const double d=dot(local-a,n);
            if(d<closest){closest=d;normal=n;arm=local-n*d;}}
        append(out,arm,normal,closest,contactMaterial);
    }
    return out;
}
inline double panelStored(const State& s,const Config& c){
    return c.barrier.kind==BarrierKind::Elastic&&!s.barrierFailed?
        .5*c.barrier.stiffness*s.barrierDeflection*s.barrierDeflection+.5*c.barrier.mass*s.barrierVelocity*s.barrierVelocity:0;
}
inline void advancePanel(State& s,const Config& c,double h){
    if(c.barrier.kind!=BarrierKind::Elastic||s.barrierFailed)return;
    const auto& b=c.barrier;const double denominator=1+h*b.damping/b.mass+h*h*b.stiffness/b.mass;
    s.barrierVelocity=(s.barrierVelocity-h*b.stiffness*s.barrierDeflection/b.mass)/denominator;
    s.barrierDeflection+=h*s.barrierVelocity;
    s.barrierDissipation+=h*b.damping*s.barrierVelocity*s.barrierVelocity;
}
inline void finishPanel(State& s,const Config& c){
    if(c.barrier.kind!=BarrierKind::Elastic||s.barrierFailed){s.barrierForce=0;return;}
    const auto& b=c.barrier;s.barrierForce=b.stiffness*s.barrierDeflection+b.damping*s.barrierVelocity;
    s.peakBarrierForce=std::max(s.peakBarrierForce,std::abs(s.barrierForce));
    s.maxBarrierDeflection=std::max(s.maxBarrierDeflection,std::abs(s.barrierDeflection));
    if(std::abs(s.barrierDeflection)>b.maxDeflection){
        s.barrierDissipation+=panelStored(s,c);s.barrierFailed=true;s.barrierVelocity=0;s.barrierForce=0;
    }
}
inline void barrierImpulse(State& s,const Config& c,double h,double approachSide){
    auto contacts=barrierContacts(s,c,approachSide);if(contacts.empty())return;
    s.barrierHit=true;const bool elastic=c.barrier.kind==BarrierKind::Elastic;
    const auto& b=c.barrier;const double denominator=elastic?1+h*b.damping/b.mass+h*h*b.stiffness/b.mass:1;
    const double panelBefore=s.barrierVelocity;
    for(auto& q:contacts){const double vn=dot(s.v+rotationalVelocity(s.omega,q.arm)-Vec2{s.barrierVelocity,0},q.normal);
        q.target=!elastic&&vn<-.25?-b.restitution*vn:0;}
    for(int iteration=0;iteration<16;++iteration)for(auto& q:contacts){
        const double vn=dot(s.v+rotationalVelocity(s.omega,q.arm)-Vec2{elastic?s.barrierVelocity:0,0},q.normal);
        const double inverse=effectiveInverseMass(s,q.arm,q.normal)+(elastic?q.normal.x*q.normal.x/(b.mass*denominator):0);
        const double old=q.jn;q.jn=std::max(0.,old+(q.target-vn)/inverse);const double delta=q.jn-old;
        impulse(s,q.arm,q.normal,delta);
        if(elastic){
            const double dv=-q.normal.x*delta/(b.mass*denominator);s.barrierVelocity+=dv;s.barrierDeflection+=h*dv;
            s.p=s.p+q.normal*(delta/s.mass*h);s.angle+=cross(q.arm,q.normal)*delta/s.inertia*h;
        }
    }
    for(const auto& q:contacts)s.contactImpulse+=q.jn;
    if(elastic)s.barrierDissipation+=h*b.damping*(s.barrierVelocity*s.barrierVelocity-panelBefore*panelBefore);
    // Resolve remaining geometry against the current panel position. The panel itself
    // is never displaced by projection, so correction cannot invent spring energy.
    for(int i=0;i<8;++i){auto remaining=barrierContacts(s,c,approachSide);if(remaining.empty())break;
        const auto worst=std::max_element(remaining.begin(),remaining.end(),[](const Contact& a,const Contact& b){return a.depth<b.depth;});
        if(worst->depth<1e-9)break;s.p=s.p+worst->normal*(worst->depth+1e-9);
    }
}
inline double advancedStepLimit(const State& s,const Config& c){
    const double boundarySpeed=length(s.v)+std::abs(s.omega)*c.radius+std::abs(s.barrierVelocity)+.1;
    double h=std::min(.0005,.04*c.radius/boundarySpeed);
    h=std::min(h,.01/(std::abs(s.omega)+1e-9));
    if(c.barrier.kind==BarrierKind::Elastic&&!s.barrierFailed)h=std::min(h,.05*std::sqrt(c.barrier.mass/c.barrier.stiffness));
    return h;
}
inline void advancedSubstep(State& s,const Config& c,double h){
    if(++s.integrationSteps>5000000)throw std::runtime_error("One trajectory exceeded five million collision substeps; reduce its duration or stiffness");
    const State previous=s;
    const double previousEnergy=kineticEnergy(s)+s.mass*G*s.p.z+panelStored(s,c)+s.barrierDissipation;
    const double previousWall=c.barrier.x+(c.barrier.kind==BarrierKind::Elastic?s.barrierDeflection:0);
    const double approachSide=s.p.x<=previousWall?-1.:1.;
    advancePanel(s,c,h);
    const auto material=surface(c,s.p.x);
    const double kd=.5*1.225*material.drag*PI*c.radius*c.radius/s.mass;
    s.v=s.v*(1/(1+kd*length(s.v)*h));s.v.z-=G*h;s.p=s.p+s.v*h;s.angle+=s.omega*h;s.time+=h;
    double correction=0;auto contacts=groundContacts(s,c,correction);
    if(c.shape==Shape::Polygon)s.p.z+=std::max(0.,correction);
    else for(int k=0;k<4;++k){auto current=groundContacts(s,c,correction);if(current.empty())break;
        const auto deepest=std::max_element(current.begin(),current.end(),[](const Contact& a,const Contact& b){return a.depth<b.depth;});
        if(deepest->depth<1e-9)break;s.p=s.p+deepest->normal*deepest->depth;}
    if(!contacts.empty()){
        groundImpulse(s,contacts);
        const auto normal=groundNormal(c.terrain,s.p.x);
        const double dw=material.rollResistance*s.mass*G*normal.z*c.radius*h/s.inertia;
        s.omega=std::copysign(std::max(0.,std::abs(s.omega)-dw),s.omega);
    }
    barrierImpulse(s,c,h,approachSide);
    // Passive energy guard for discretization/projection error. Initial overlaps are
    // invalid starting states; the normal initialized configurations have no overlap.
    const double panelK=c.barrier.kind==BarrierKind::Elastic&&!s.barrierFailed?.5*c.barrier.mass*s.barrierVelocity*s.barrierVelocity:0;
    const double currentK=kineticEnergy(s)+panelK;
    const double fixedEnergy=s.mass*G*s.p.z+panelStored(s,c)-panelK+s.barrierDissipation;
    const double allowed=std::max(0.,previousEnergy-fixedEnergy);
    if(currentK>allowed+1e-8&&currentK>0){const double f=std::sqrt(allowed/currentK);s.v=s.v*f;s.omega*=f;s.barrierVelocity*=f;}
    finishPanel(s,c);
    if(c.barrier.kind!=BarrierKind::None){
        const double wall=c.barrier.x+(c.barrier.kind==BarrierKind::Elastic&&!s.barrierFailed?s.barrierDeflection:0);
        const double left=s.p.x-supportRadius(c,s.angle,{1,0});
        if(left>wall+1e-6){s.barrierPassed=true;if(!s.barrierFailed)s.barrierOvertopped=true;}
    }
    const double slope=c.terrain.angle(s.p.x);
    const bool canRest=c.shape==Shape::Polygon||std::abs(std::tan(slope))<=material.rollResistance+1e-8||s.barrierHit;
    if(!contacts.empty()&&canRest&&length(s.v)<.025&&std::abs(s.omega)*c.radius<.025&&std::abs(s.barrierVelocity)<.025)s.restTime=previous.restTime+h;
    else s.restTime=0;
    if(s.restTime>.5){
        s.v={0,0};s.omega=0;s.moving=false;s.stopReason="stopped";
        if(s.barrierHit&&!s.barrierFailed&&!s.barrierPassed){
            // Stable near-wall contact, not merely failure to reach the collector.
            State near=s;const double wall=c.barrier.x+(c.barrier.kind==BarrierKind::Elastic?s.barrierDeflection:0);
            near.p.x+=(near.p.x<wall?.01:-.01);
            if(!barrierContacts(near,c,approachSide).empty()){s.barrierIntercepted=true;s.stopReason="intercepted";}
        }
    }
    if(s.p.x<c.terrain.points.front().x-20||s.p.x>simulationDomainRight(c)){s.moving=false;s.stopReason="domain_exit";}
    if(!std::isfinite(s.p.x)||!std::isfinite(s.p.z)||!std::isfinite(s.v.x)||!std::isfinite(s.v.z)||!std::isfinite(s.omega)||!std::isfinite(s.barrierDeflection))
        throw std::runtime_error("Nonfinite extended physics state");
}
} // namespace detail
template<class Observer> inline void stepObserved(State& s,const Config& c,Observer&& observe){
    if(!s.moving)return;
    if(c.shape==Shape::Sphere&&c.barrier.kind==BarrierKind::None){const State previous=s;stepSphereLegacy(s,c);observe(previous,s);return;}
    double remaining=c.dt;int substeps=0;
    while(remaining>1e-12&&s.moving){
        const double h=std::min(remaining,detail::advancedStepLimit(s,c));
        if(++substeps>20000||h<1e-10)throw std::runtime_error("Extended collision substep limit exceeded");
        const State previous=s;detail::advancedSubstep(s,c,h);observe(previous,s);remaining-=h;
    }
}
inline void step(State& s,const Config& c){stepObserved(s,c,[](const State&,const State&){});}
inline Sample sample(const State& s){
    Sample p{s.time,s.p.x,s.p.z,s.v.x,s.v.z,s.angle,s.omega,kineticEnergy(s),s.barrierDeflection};
    p.barrierHit=s.barrierHit;p.barrierFailed=s.barrierFailed;p.barrierVelocity=s.barrierVelocity;
    p.barrierForce=s.barrierForce;p.barrierDissipation=s.barrierDissipation;p.contactImpulse=s.contactImpulse;return p;
}
inline RunResult simulate(const Config& c){
    validate(c);RunResult result;result.config=c;State s=initialState(c);
    result.trajectory.push_back(sample(s));double nextSample=c.sampleInterval;
    double detector=c.detectorX<0?c.terrain.toeX():c.detectorX;
    result.runoutX=s.p.x;result.maxSpeed=length(s.v);result.maxEnergy=kineticEnergy(s);
    const int steps=static_cast<int>(std::ceil(c.maxTime/c.dt));
    for(int i=0;i<steps&&s.moving;++i){
        stepObserved(s,c,[&](const State& previous,const State&){
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
            const double crossingStep=c.shape==Shape::Sphere&&c.barrier.kind==BarrierKind::None?c.dt:s.time-previous.time;
            result.detectorTime=previous.time+crossingStep*f;
            result.inBand=result.detectorHeight>=0&&result.detectorHeight<=10;
        }
        if(s.time+1e-9>=nextSample){result.trajectory.push_back(sample(s));nextSample+=c.sampleInterval;}
        });
    }
    if(result.trajectory.back().t<s.time-1e-9)result.trajectory.push_back(sample(s));
    if(s.moving){s.stopReason="time_limit";s.moving=false;}
    result.status=s.stopReason;result.finalState=s;return result;
}
}


// ==================== batch.hpp ====================

// ==================== batch.hpp ====================

#include <iomanip>
#include <set>
#include <sstream>

namespace rock {
inline double normalBounded(double mean,double cv,double lo,double hi,std::mt19937& rng) {
    if(cv==0 || mean==0) return mean;
    std::normal_distribution<double> distribution(mean,std::abs(mean)*cv);
    for(int attempt=0;attempt<10000;++attempt){
        const double v=distribution(rng);
        if(std::isfinite(v) && v>=lo && v<=hi)return v;
    }
    throw std::runtime_error("Unable to sample physical parameter bounds");
}
inline Config randomize(const Config& base,std::mt19937& rng) {
    Config c=base;
    c.initialSpeed=normalBounded(base.initialSpeed,base.cv,0,200,rng);
    if(base.randomSurface){
        auto sample=[&](Material& m){
            m.rollResistance=normalBounded(m.rollResistance,base.cv,0,1,rng);
            m.friction=normalBounded(m.friction,base.cv,0,2,rng);
            m.restitution=normalBounded(m.restitution,base.cv,0,1,rng);
            m.drag=normalBounded(m.drag,base.cv,0,10,rng);
        };
        sample(c.material);
        for(auto& m:c.terrain.surfaces)sample(m);
    }
    return c;
}
inline double nominalIntegrationStep(const Config& c){
    return c.shape==Shape::Sphere&&c.barrier.kind==BarrierKind::None?c.dt:std::min(c.dt,.0005);
}
inline BatchResult runBatch(const Config& base,int count,uint32_t seed) {
    validate(base);
    if(count<1 || count>5000)throw std::invalid_argument("runs must be in [1,5000]");
    if(count*std::ceil(base.maxTime/nominalIntegrationStep(base))>2e8)
        throw std::invalid_argument("Requested integration work exceeds 200 million steps; reduce runs or time horizon");
    if(count*base.maxTime/base.sampleInterval>4e6)
        throw std::invalid_argument("Too many saved trajectory samples; reduce runs or increase --sample-interval");
    BatchResult b;b.base=base;b.seed=seed;
    b.heightCounts.assign(static_cast<size_t>(std::ceil(10/base.binHeight)),0);
    uint64_t actualSteps=0;
    for(int i=0;i<count;++i){
        std::seed_seq trialSeed{seed,static_cast<uint32_t>(i)};
        std::mt19937 rng(trialSeed);
        auto r=simulate(randomize(base,rng));
        actualSteps+=r.finalState.integrationSteps;
        if(actualSteps>200000000ULL)throw std::runtime_error("Actual collision work exceeded 200 million substeps; reduce runs");
        if(r.reached){
            ++b.reached;
            if(r.inBand){++b.inBand;const size_t k=std::min(b.heightCounts.size()-1,
                    static_cast<size_t>(std::max(0.,r.detectorHeight)/base.binHeight));++b.heightCounts[k];}
            else if(r.detectorHeight>10)++b.aboveBand;
            else ++b.belowBand;
        }else ++b.notReached;
        std::set<std::pair<int,int>> visited;
        // Each trial contributes at most once to a cell. Interpolate saved segments,
        // so fast motion cannot skip a cell merely because output sampling is sparse.
        auto mark=[&](double x,double z){
            const double h=z-base.terrain.toeZ();
            if(h<0 || h>10 || x<base.terrain.points.front().x || x>simulationDomainRight(base))return;
            const int ix=static_cast<int>(std::floor(x/base.gridWidth));
            const int iz=std::min(static_cast<int>(b.heightCounts.size())-1,static_cast<int>(h/base.binHeight));
            visited.insert({ix,iz});
        };
        for(size_t j=0;j<r.trajectory.size();++j){
            const auto& a=r.trajectory[j];
            if(j==0){mark(a.x,a.z);continue;}
            const auto& p=r.trajectory[j-1];
            int steps=std::max(1,static_cast<int>(std::ceil(std::max(std::abs(a.x-p.x)/base.gridWidth,
                  std::abs(a.z-p.z)/base.binHeight)*2)));
            for(int k=0;k<=steps;++k){const double f=double(k)/steps;mark(p.x+(a.x-p.x)*f,p.z+(a.z-p.z)*f);}
        }
        for(const auto& cell:visited)++b.routeCells[cell];
        b.runs.push_back(std::move(r));
    }
    return b;
}
inline std::pair<double,double> wilson(size_t k,size_t n){
    if(!n)return {0,0};
    const double z=1.95996398454,p=double(k)/n,d=1+z*z/n;
    const double center=(p+z*z/(2*n))/d;
    const double half=z*std::sqrt(p*(1-p)/n+z*z/(4.*n*n))/d;
    return {std::max(0.,center-half),std::min(1.,center+half)};
}
inline std::string batchSummary(const BatchResult& b){
    std::ostringstream s;s<<std::fixed<<std::setprecision(4);
    const size_t n=b.runs.size();const auto ci=wilson(b.inBand,n);
    s<<"runs="<<n<<" seed="<<b.seed<<"\nreached="<<b.reached<<" in_0_10m="<<b.inBand
     <<" above_10m="<<b.aboveBand<<" below_0m="<<b.belowBand<<" not_reached="<<b.notReached
     <<"\nP(in_0_10m)="<<(n?double(b.inBand)/n:0)<<" Wilson95%=["<<ci.first<<","<<ci.second<<"]\n";
    for(size_t i=0;i<b.heightCounts.size();++i)
        s<<i*b.base.binHeight<<".."<<std::min(10.,(i+1)*b.base.binHeight)<<" m: "
         <<b.heightCounts[i]<<" P="<<(n?double(b.heightCounts[i])/n:0)<<"\n";
#if ROCK_STAGE >= 6
    if(b.base.barrier.kind!=BarrierKind::None){
        size_t hits=0,intercepted=0,failed=0,overtopped=0;
        for(const auto& r:b.runs){hits+=r.finalState.barrierHit;intercepted+=r.finalState.barrierIntercepted;
            failed+=r.finalState.barrierFailed;overtopped+=r.finalState.barrierOvertopped;}
        s<<"barrier_hit="<<hits<<" stable_interception="<<intercepted<<" failed="<<failed<<" overtopped="<<overtopped<<'\n';
        s<<"P(stable_interception)="<<(n?double(intercepted)/n:0)<<"; not_reached is not interception.\n";
    }
#endif
    return s.str();
}
}


// ==================== export.hpp ====================

#if ROCK_STAGE >= 2
// ==================== export.hpp ====================



#include <filesystem>
#include <fstream>
#include <iomanip>
#include <locale>

namespace rock {
// C++20 uses char8_t for UTF-8 strings. Keep filesystem paths encoded correctly
// without relying on u8path(), which is deprecated in C++20.
inline std::filesystem::path pathFromUtf8(const std::string& value) {
#ifdef __cpp_char8_t
    return std::filesystem::path(std::u8string(value.begin(), value.end()));
#else
    return std::filesystem::u8path(value);
#endif
}
inline std::string pathToUtf8(const std::filesystem::path& value) {
    const auto bytes = value.u8string();
    return std::string(reinterpret_cast<const char*>(bytes.data()), bytes.size());
}
inline std::ofstream csvFile(const std::filesystem::path& p){
    std::ofstream f(p,std::ios::binary);
    if(!f)throw std::runtime_error("Cannot write "+pathToUtf8(p));
    f.imbue(std::locale::classic());f<<std::setprecision(12);return f;
}
inline const char* shapeName(Shape s){return s==Shape::Sphere?"sphere":"polygon";}
inline const char* barrierName(BarrierKind k){return k==BarrierKind::None?"none":k==BarrierKind::Rigid?"rigid":"elastic";}
#if ROCK_STAGE >= 6
inline const char* barrierOutcome(const RunResult& r){
    if(r.config.barrier.kind==BarrierKind::None)return "disabled";
    const auto& s=r.finalState;
    if(s.barrierIntercepted)return "intercepted";
    if(s.barrierFailed&&s.barrierPassed)return "failed_passage";
    if(s.barrierOvertopped&&s.barrierPassed)return "overtopped";
    if(s.barrierFailed)return "failed_unresolved";
    if(s.barrierHit)return "hit_unresolved";
    return "not_encountered";
}
#endif
inline void exportBatch(const BatchResult& b,const std::string& outputDir){
    namespace fs=std::filesystem;const auto out=pathFromUtf8(outputDir);
    if(fs::exists(out) && !fs::is_empty(out))throw std::runtime_error("Output folder must be new or empty: "+outputDir);
    fs::create_directories(out);
    const auto& c=b.base;const size_t n=b.runs.size();
    auto bins=csvFile(out/"height_probability.csv");
    bins<<"height_lower_m,height_upper_m,trial_count,total_trials,probability,probability_given_in_band,wilson95_low,wilson95_high\n";
    for(size_t i=0;i<b.heightCounts.size();++i){const auto ci=wilson(b.heightCounts[i],n);
        bins<<i*c.binHeight<<','<<std::min(10.,(i+1)*c.binHeight)<<','<<b.heightCounts[i]<<','<<n<<','
            <<double(b.heightCounts[i])/n<<','<<(b.inBand?double(b.heightCounts[i])/b.inBand:0)<<','<<ci.first<<','<<ci.second<<'\n';}
    auto grid=csvFile(out/"route_probability.csv");
    grid<<"x_lower_m,x_upper_m,height_lower_m,height_upper_m,visited_trials,total_trials,probability\n";
    for(const auto& item:b.routeCells)grid<<item.first.first*c.gridWidth<<','<<(item.first.first+1)*c.gridWidth<<','
        <<item.first.second*c.binHeight<<','<<std::min(10.,(item.first.second+1)*c.binHeight)<<','<<item.second<<','<<n<<','<<double(item.second)/n<<'\n';
    auto runs=csvFile(out/"runs.csv");
    runs<<"trial_id,seed,initial_speed_m_s,radius_m,roll_resistance,friction,restitution,drag,shape,sides,reached,in_0_10m,detector_height_m,detector_speed_m_s,detector_energy_J,detector_time_s,runout_x_m,max_speed_m_s,max_energy_J,max_bounce_height_m,duration_s,status,barrier_hit,barrier_failed,max_deflection_m,peak_support_force_N,dissipated_energy_J";
#if ROCK_STAGE >= 6
    runs<<",barrier_contact_impulse_N_s,barrier_passed,barrier_overtopped,barrier_intercepted,barrier_outcome,integration_steps";
#endif
    runs<<'\n';
    auto trajectory=csvFile(out/"trajectories.csv");
    trajectory<<"trial_id,t_s,x_m,z_m,height_above_toe_m,vx_m_s,vz_m_s,angle_rad,omega_rad_s,kinetic_energy_J,barrier_deflection_m";
#if ROCK_STAGE >= 6
    trajectory<<",barrier_hit,barrier_failed,barrier_velocity_m_s,barrier_support_force_N,barrier_dissipated_energy_J,barrier_contact_impulse_N_s";
#endif
#if ROCK_STAGE >= 8
    trajectory<<",rock_potential_energy_J,rock_mechanical_energy_J,barrier_elastic_energy_J,barrier_kinetic_energy_J,system_energy_plus_dissipation_J";
#endif
    trajectory<<'\n';
    auto surfaces=csvFile(out/"sampled_surfaces.csv");
    surfaces<<"trial_id,segment_index,roll_resistance,friction,restitution,drag\n";
    for(size_t i=0;i<n;++i){const auto& r=b.runs[i];const auto& m=r.config.material;const auto& s=r.finalState;
        runs<<i+1<<','<<b.seed<<','<<r.config.initialSpeed<<','<<r.config.radius<<','<<m.rollResistance<<','<<m.friction<<','<<m.restitution<<','<<m.drag<<','<<shapeName(r.config.shape)<<','<<r.config.sides<<','<<r.reached<<','<<r.inBand<<',';
        if(r.reached)runs<<r.detectorHeight<<','<<r.detectorSpeed<<','<<r.detectorEnergy<<','<<r.detectorTime;
        else runs<<",,,";
        runs<<','<<r.runoutX<<','<<r.maxSpeed<<','<<r.maxEnergy<<','<<r.maxBounceHeight<<','<<s.time<<','<<r.status<<','<<s.barrierHit<<','<<s.barrierFailed<<','<<s.maxBarrierDeflection<<','<<s.peakBarrierForce<<','<<s.barrierDissipation;
#if ROCK_STAGE >= 6
        runs<<','<<s.contactImpulse<<','<<s.barrierPassed<<','<<s.barrierOvertopped<<','<<s.barrierIntercepted<<','<<barrierOutcome(r)<<','<<s.integrationSteps;
#endif
        runs<<'\n';
        for(const auto& t:r.trajectory){
            trajectory<<i+1<<','<<t.t<<','<<t.x<<','<<t.z<<','<<t.z-c.terrain.toeZ()<<','<<t.vx<<','<<t.vz<<','<<t.angle<<','<<t.omega<<','<<t.energy<<','<<t.deflection;
#if ROCK_STAGE >= 6
            trajectory<<','<<t.barrierHit<<','<<t.barrierFailed<<','<<t.barrierVelocity<<','<<t.barrierForce<<','<<t.barrierDissipation<<','<<t.contactImpulse;
#endif
#if ROCK_STAGE >= 8
            const double potential=s.mass*G*(t.z-c.terrain.toeZ());
            const bool panel=r.config.barrier.kind==BarrierKind::Elastic&&!t.barrierFailed;
            const double stored=panel?.5*r.config.barrier.stiffness*t.deflection*t.deflection:0.;
            const double panelK=panel?.5*r.config.barrier.mass*t.barrierVelocity*t.barrierVelocity:0.;
            trajectory<<','<<potential<<','<<t.energy+potential<<','<<stored<<','<<panelK<<','<<t.energy+potential+stored+panelK+t.barrierDissipation;
#endif
            trajectory<<'\n';
        }
        for(size_t j=0;j+1<r.config.terrain.points.size();++j){
            const auto& sm=r.config.terrain.surfaces.empty()?m:r.config.terrain.surfaces[j];
            surfaces<<i+1<<','<<j<<','<<sm.rollResistance<<','<<sm.friction<<','<<sm.restitution<<','<<sm.drag<<'\n';
        }
    }
    auto terrain=csvFile(out/"terrain.csv");terrain<<"point_index,x_m,z_m\n";
    for(size_t i=0;i<c.terrain.points.size();++i)terrain<<i<<','<<c.terrain.points[i].x<<','<<c.terrain.points[i].z<<'\n';
    auto metadata=csvFile(out/"config.csv");metadata<<"parameter,value\n"
        <<"schema_version,2\nprogram_version,"<<ROCK_VERSION<<"\nseed,"<<b.seed<<"\nruns,"<<n<<"\ndt_s,"<<c.dt<<"\nmax_time_s,"<<c.maxTime
        <<"\nsample_interval_s,"<<c.sampleInterval<<"\ncoefficient_of_variation,"<<c.cv<<"\nrandom_surface,"<<c.randomSurface
        <<"\nradius_m,"<<c.radius<<"\ninitial_speed_m_s,"<<c.initialSpeed<<"\ndensity_kg_m3,"<<c.density
        <<"\nroll_resistance,"<<c.material.rollResistance<<"\nfriction,"<<c.material.friction<<"\nrestitution,"<<c.material.restitution<<"\ndrag,"<<c.material.drag
        <<"\ndetector_x_m,"<<(c.detectorX<0?c.terrain.toeX():c.detectorX)<<"\ndetector_reference_z_m,"<<c.terrain.toeZ()
        <<"\nbin_height_m,"<<c.binHeight<<"\ngrid_width_m,"<<c.gridWidth<<"\nshape,"<<shapeName(c.shape)<<"\nsides,"<<c.sides<<"\nshape_phase_rad,"<<c.shapePhase
        <<"\nbarrier,"<<barrierName(c.barrier.kind)<<"\nbarrier_x_m,"<<c.barrier.x<<"\nbarrier_height_m,"<<c.barrier.height
        <<"\nbarrier_restitution,"<<c.barrier.restitution<<"\nstiffness_N_m,"<<c.barrier.stiffness<<"\ndamping_N_s_m,"<<c.barrier.damping<<"\nmax_deflection_m,"<<c.barrier.maxDeflection<<"\nbarrier_mass_kg,"<<c.barrier.mass
        <<"\nterrain_width_scale,"<<c.terrainWidthScale<<"\nterrain_height_scale,"<<c.terrainHeightScale<<'\n';
    auto summary=csvFile(out/"summary.txt");summary<<batchSummary(b);
    auto readme=csvFile(out/"READ_ME.txt");
    readme<<"UTF-8 CSV; dot decimal; SI units (m,s,kg,N,J).\n"
       <<"Detector: first left-to-right centre crossing; height measured from toe elevation.\n"
       <<"All unconditional probabilities use ALL trials (including not reached) as denominator.\n"
       <<"Bins are [low,high), except the last bin includes 10m. Missing detector fields mean no crossing.\n"
       <<"Conditional bin probability uses trials in 0..10m as denominator (zero when none).\n"
       <<"Route grid = per-trial centre path visits / N; each trial counted once per cell. Cells do not sum to one.\n"
       <<"Trajectory segments are interpolated for the grid; finer sample intervals improve spatial accuracy.\n"
       <<"Zero-hit cells are omitted; seed reproducibility assumes the same executable/standard library.\n"
       <<"Gaussian sigma=abs(mean)*CV; physically invalid draws are rejected (truncated normal).\n"
       <<"Zero-mean variables stay zero. Stored mean/CV are parent-normal parameters, not post-truncation moments.\n"
       <<"Polygon geometry is a regular convex lamina with circumscribed radius and equivalent sphere mass.\n"
       <<"Extended collision dynamics use adaptive substeps; near-start/stop polygon paths may be timestep-sensitive.\n"
       <<"peak_support_force_N is the elastic spring/damper reaction, not impulse divided by dt; rigid value is 0.\n"
       <<"Barrier contact impulse is cumulative normal impulse (N s). Stable interception is distinct from not reaching the detector.\n"
       <<"failed_unresolved and hit_unresolved are censored outcomes, not successful interception.\n"
       <<"Energy is relative to toe elevation; passive contact correction and implicit integration can dissipate additional numerical energy.\n"
       <<"This is a 2-D educational model, not a calibrated engineering design or barrier certification tool.\n";
#if ROCK_STAGE >= 6
    auto outcomes=csvFile(out/"barrier_outcomes.csv");
    outcomes<<"outcome,trial_count,total_trials,probability,wilson95_low,wilson95_high\n";
    std::map<std::string,size_t> counts;
    for(const char* key:{"disabled","intercepted","failed_passage","overtopped","failed_unresolved","hit_unresolved","not_encountered"})counts[key]=0;
    for(const auto& r:b.runs)++counts[barrierOutcome(r)];
    for(const auto& item:counts){const auto ci=wilson(item.second,n);outcomes<<item.first<<','<<item.second<<','<<n<<','<<double(item.second)/n<<','<<ci.first<<','<<ci.second<<'\n';}
    outcomes.flush();if(!outcomes)throw std::runtime_error("Failed to write barrier outcomes");
#endif
    for(auto* f:{&bins,&grid,&runs,&trajectory,&surfaces,&terrain,&metadata,&summary,&readme}){
        f->flush();if(!*f)throw std::runtime_error("CSV output failed (disk full or I/O error)");
    }
}
}

#endif

// ==================== studies.hpp ====================

#if ROCK_STAGE >= 3
// ==================== studies.hpp ====================


#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <limits>
#include <locale>
#include <sstream>

namespace rock {
namespace study {
inline double missing(){return std::numeric_limits<double>::quiet_NaN();}
struct Moments {
    size_t n=0; double mean=0,m2=0;
    void add(double x){if(!std::isfinite(x))return; ++n; double d=x-mean; mean+=d/double(n); m2+=d*(x-mean);}
    double value()const{return n?mean:missing();}
    double sd()const{return n>1?std::sqrt(std::max(0.,m2/double(n-1))):missing();}
    double se()const{return n>1?sd()/std::sqrt(double(n)):missing();}
};
inline std::pair<double,double> interval(const Moments& m,bool probability){
    if(!m.n)return {missing(),missing()};
    if(probability){ // Wilson score interval, defined even for zero successes.
        const double z=1.959963984540054,n=double(m.n),p=m.mean,d=1+z*z/n;
        const double center=(p+z*z/(2*n))/d;
        const double half=z*std::sqrt(p*(1-p)/n+z*z/(4*n*n))/d;
        return {std::max(0.,center-half),std::min(1.,center+half)};
    }
    if(m.n<2)return {missing(),missing()};
    // Student-t approximation; for n<=30 use an explicit 97.5% t table.
    static const double t[]={0,12.7062,4.30265,3.18245,2.77645,2.57058,2.44691,2.36462,
        2.30600,2.26216,2.22814,2.20099,2.17881,2.16037,2.14479,2.13145,2.11991,
        2.10982,2.10092,2.09302,2.08596,2.07961,2.07387,2.06866,2.06390,2.05954,
        2.05553,2.05183,2.04841,2.04523,2.04227};
    const double df=double(m.n-1),critical=df<=30?t[size_t(df)]:1.959963984540054+2.372271/double(df);
    return {m.mean-critical*m.se(),m.mean+critical*m.se()};
}
inline void number(std::ostream& out,double x){if(std::isfinite(x))out<<std::setprecision(12)<<x;}
inline std::ofstream open(const std::filesystem::path& path){
    std::ofstream out(path,std::ios::binary);
    if(!out)throw std::runtime_error("Cannot write study file: "+pathToUtf8(path));
    out.exceptions(std::ios::badbit|std::ios::failbit);
    out.imbue(std::locale::classic());
    out<<std::setprecision(12);return out;
}
struct Metric {
    std::string name; bool probability=false,conditional=false;
    std::function<double(const RunResult&)> get;
};
inline const std::vector<Metric>& metrics(){
    static const std::vector<Metric> m={
        {"reach_probability",true,false,[](const RunResult&r){return r.reached?1.:0.;}},
        {"height_0_10m_probability",true,false,[](const RunResult&r){return r.inBand?1.:0.;}},
        {"collector_speed_all_runs_m_s",false,false,[](const RunResult&r){return r.reached?r.detectorSpeed:0.;}},
        {"collector_energy_all_runs_J",false,false,[](const RunResult&r){return r.reached?r.detectorEnergy:0.;}},
        {"collector_speed_given_reached_m_s",false,true,[](const RunResult&r){return r.reached?r.detectorSpeed:missing();}},
        {"collector_energy_given_reached_J",false,true,[](const RunResult&r){return r.reached?r.detectorEnergy:missing();}},
        {"max_speed_m_s",false,false,[](const RunResult&r){return r.maxSpeed;}},
        {"max_energy_J",false,false,[](const RunResult&r){return r.maxEnergy;}},
        {"runout_x_m",false,false,[](const RunResult&r){return r.runoutX;}},
        {"max_bounce_height_m",false,false,[](const RunResult&r){return r.maxBounceHeight;}},
        {"duration_s",false,false,[](const RunResult&r){return r.finalState.time;}}
    };return m;
}
inline Moments aggregate(const BatchResult& b,const Metric& m){Moments s;for(const auto&r:b.runs)s.add(m.get(r));return s;}
inline Moments paired(const BatchResult& a,const BatchResult& b,const Metric& m){
    Moments s;if(m.conditional)return s;
    if(a.runs.size()!=b.runs.size())throw std::runtime_error("Paired study requires equal run counts");
    for(size_t i=0;i<a.runs.size();++i)s.add(m.get(a.runs[i])-m.get(b.runs[i]));
    return s;
}
inline std::pair<double,double> effectInterval(const BatchResult& a,const BatchResult& b,const Metric& m){
    const auto effect=paired(a,b,m);
    if(!m.probability)return interval(effect,false);
    if(&a==&b)return {0.,0.};
    if(!effect.n)return {missing(),missing()};
    size_t positive=0,negative=0;
    for(size_t i=0;i<a.runs.size();++i){const double d=m.get(a.runs[i])-m.get(b.runs[i]);positive+=d>0;negative+=d<0;}
    // Two 97.5% Wilson intervals and a Bonferroni union bound. This avoids a
    // misleading zero-width normal interval when no discordant pair was seen.
    auto bounds=[&](size_t k){const double z=2.241402727604947,n=double(effect.n),p=double(k)/n,den=1+z*z/n;
        const double center=(p+z*z/(2*n))/den,half=z*std::sqrt(p*(1-p)/n+z*z/(4*n*n))/den;
        return std::pair<double,double>{std::max(0.,center-half),std::min(1.,center+half)};};
    const auto p=bounds(positive),n=bounds(negative);return {p.first-n.second,p.second-n.first};
}
inline size_t censored(const BatchResult& b){size_t n=0;for(const auto&r:b.runs)n+=r.status=="domain_exit"||r.status=="time_limit"||r.status=="timeout";return n;}
inline void writeConfig(const Config& c,int n,uint32_t seed,const std::filesystem::path& dir){
    std::filesystem::create_directories(dir);
    auto out=open(dir/"config.txt");
    out<<"program_version="<<ROCK_VERSION<<"\nruns="<<n<<"\nseed="<<seed<<"\ncv="<<c.cv<<"\nrandom_surface="<<c.randomSurface
       <<"\ninitial_speed_m_s="<<c.initialSpeed<<"\nradius_m="<<c.radius<<"\ndensity_kg_m3="<<c.density
       <<"\ndt_s="<<c.dt<<"\nmax_time_s="<<c.maxTime<<"\nsample_interval_s="<<c.sampleInterval
       <<"\ncollector_x_m="<<(c.detectorX<0?c.terrain.toeX():c.detectorX)
       <<"\ncollector_reference_z_m="<<c.terrain.toeZ()<<"\nbin_height_m="<<c.binHeight
       <<"\nterrain_width_scale="<<c.terrainWidthScale<<"\nterrain_height_scale="<<c.terrainHeightScale
       <<"\ngrid_width_m="<<c.gridWidth<<"\nshape="<<(c.shape==Shape::Sphere?"sphere":"polygon")
       <<"\npolygon_sides="<<c.sides<<"\nshape_phase_rad="<<c.shapePhase
       <<"\nroll_resistance="<<c.material.rollResistance<<"\nfriction="<<c.material.friction
       <<"\nrestitution="<<c.material.restitution<<"\ndrag="<<c.material.drag
       <<"\nbarrier_kind="<<int(c.barrier.kind)<<"\nbarrier_x_m="<<c.barrier.x
       <<"\nbarrier_height_m="<<c.barrier.height<<"\nbarrier_restitution="<<c.barrier.restitution
       <<"\nbarrier_stiffness_N_m="<<c.barrier.stiffness<<"\nbarrier_damping_N_s_m="<<c.barrier.damping
       <<"\nbarrier_max_deflection_m="<<c.barrier.maxDeflection<<"\nbarrier_mass_kg="<<c.barrier.mass<<"\n";
    auto terrain=open(dir/"terrain.csv");terrain<<"point_index,x_m,z_m,outgoing_roll_resistance,outgoing_friction,outgoing_restitution,outgoing_drag\n";
    for(size_t i=0;i<c.terrain.points.size();++i){
        const auto&p=c.terrain.points[i];terrain<<i<<','<<p.x<<','<<p.z;
        if(i+1<c.terrain.points.size()){
            const auto&m=c.terrain.surfaces.empty()?c.material:c.terrain.surfaces.at(i);
            terrain<<','<<m.rollResistance<<','<<m.friction<<','<<m.restitution<<','<<m.drag;
        }else terrain<<",,,,";
        terrain<<'\n';
    }
}
inline void writeRunMetrics(const BatchResult& b,const std::filesystem::path& dir){
    auto out=open(dir/"run_metrics.csv");
    out<<"run_id,initial_speed_m_s,radius_m,roll_resistance,friction,restitution,drag,reached,in_band,collector_height_m,collector_speed_m_s,collector_energy_J,runout_x_m,max_speed_m_s,max_energy_J,max_bounce_height_m,duration_s,status\n";
    for(size_t i=0;i<b.runs.size();++i){
        const auto&r=b.runs[i];const auto&c=r.config;
        out<<i+1<<','<<c.initialSpeed<<','<<c.radius<<','<<c.material.rollResistance<<','<<c.material.friction<<','
           <<c.material.restitution<<','<<c.material.drag<<','<<r.reached<<','<<r.inBand<<',';
        number(out,r.reached?r.detectorHeight:missing());out<<',';number(out,r.reached?r.detectorSpeed:missing());out<<',';
        number(out,r.reached?r.detectorEnergy:missing());out<<','<<r.runoutX<<','<<r.maxSpeed<<','<<r.maxEnergy<<','
           <<r.maxBounceHeight<<','<<r.finalState.time<<','<<r.status<<'\n';
    }
    // Record sampled segment properties too: otherwise heterogeneous terrain cannot be reproduced.
    auto surf=open(dir/"sampled_surfaces.csv");
    surf<<"run_id,segment,roll_resistance,friction,restitution,drag\n";
    for(size_t i=0;i<b.runs.size();++i){const auto&c=b.runs[i].config;
        for(size_t j=0;j+1<c.terrain.points.size();++j){const auto&m=c.terrain.surfaces.empty()?c.material:c.terrain.surfaces.at(j);
            surf<<i+1<<','<<j<<','<<m.rollResistance<<','<<m.friction<<','<<m.restitution<<','<<m.drag<<'\n';}}
}
inline void stripTrajectories(BatchResult& b){for(auto&r:b.runs){r.trajectory.clear();r.trajectory.shrink_to_fit();}b.routeCells.clear();}
inline void summaryHeader(std::ostream& out){out<<"factor,level,scenario,runs,metric,valid_samples,mean,sd,standard_error,ci95_low,ci95_high,baseline_mean,paired_samples,effect_mean,effect_standard_error,effect_ci95_low,effect_ci95_high\n";}
inline void summaryRows(std::ostream& out,const std::string& factor,const std::string& level,const std::string& scenario,const BatchResult& b,const BatchResult& reference){
    for(const auto&m:metrics()){
        const auto s=aggregate(b,m),ref=aggregate(reference,m),effect=paired(b,reference,m);
        const auto ci=interval(s,m.probability),ei=effectInterval(b,reference,m);
        out<<factor<<','<<level<<','<<scenario<<','<<b.runs.size()<<','<<m.name<<','<<s.n<<',';
        number(out,s.value());out<<',';number(out,s.sd());out<<',';number(out,s.se());out<<',';
        number(out,ci.first);out<<',';number(out,ci.second);out<<',';number(out,ref.value());out<<','<<effect.n<<',';
        number(out,effect.value());out<<',';number(out,effect.se());out<<',';number(out,ei.first);out<<',';number(out,ei.second);out<<'\n';
    }
}
inline void describeEffects(std::ostream& out,const std::string& name,const BatchResult& b,const BatchResult& baseline){
    out<<"\n场景："<<name<<"（到达时间或计算域边界的试验："<<censored(b)<<"/"<<b.runs.size()<<"）\n";
    for(const auto&m:metrics()){
        if(m.conditional)continue;
        if(m.name!="height_0_10m_probability"&&m.name!="collector_speed_all_runs_m_s"&&
           m.name!="collector_energy_all_runs_J"&&m.name!="max_speed_m_s"&&m.name!="runout_x_m")continue;
        const auto s=aggregate(b,m),e=paired(b,baseline,m);const auto ci=effectInterval(b,baseline,m);
        out<<"  "<<m.name<<"：均值=";number(out,s.value());out<<"，相对基准差=";number(out,e.value());
        out<<"，配对差的95%区间=[";number(out,ci.first);out<<", ";number(out,ci.second);out<<"]。";
        if(m.name=="runout_x_m"&&(censored(b)||censored(baseline)))out<<"轨迹范围被截断，不能据此判断最终停积距离的增减。";
        else if(!std::isfinite(ci.first))out<<"有效样本不足，无法估计不确定度。";
        else if(ci.first>0)out<<"本次试验中该指标更高。";
        else if(ci.second<0)out<<"本次试验中该指标更低。";
        else out<<"区间包含零，本次样本尚不能确定增减方向。";
        const auto reference=aggregate(baseline,m);
        if(!m.probability&&std::isfinite(reference.value())&&std::abs(reference.value())>1e-12&&
           std::isfinite(e.value())&&std::abs(e.value())/std::abs(reference.value())<.001&&
           (ci.first>0||ci.second<0)&&!(m.name=="runout_x_m"&&(censored(b)||censored(baseline))))
            out<<"变化幅度小于基准的0.1%，需做步长收敛复核，暂不据此概括物理规律。";
        out<<'\n';
    }
}
inline void commonNotes(std::ostream& out,bool terrainStudy=false){
    out<<"参数对应：initial_speed=初速度，roll_resistance=滚动阻力，friction=摩擦系数，restitution=恢复系数，drag=空气阻力系数。\n"
       <<"概率分母始终为全部试验数，未到达断面的试验没有被剔除。断面指标取石心首次向右通过时的数值。\n"
       <<"all_runs 的断面速度/能量对未到达试验计零，是综合暴露指标，不是条件到达速度。\n"
       <<"given_reached 指标仅统计已到达试验，valid_samples 为对应有效样本数；不提供不同条件样本之间的配对效应。\n"
       <<"概率区间使用 Wilson 95%；概率差使用配对不一致事件的 Wilson/Bonferroni 保守95%区间。\n"
       <<"其他均值及配对差使用近似 Student-t 95%区间；区间衡量蒙特卡洛抽样误差，不包含模型误差。\n"
       <<"各场景使用相同种子和试验数，以共同随机数进行配对比较。\n"
       <<"正态抽样目标标准差=CV×均值；非物理值重采样。截断后实际均值和变异系数可能偏离目标。\n"
       <<"domain_exit/time_limit 表示轨迹被计算域或最长时长截断；runout_x_m 为已观测最远位置，不必然是真实停止位置。\n"
       <<"本分析为探索性多重比较，未进行整体多重检验校正，不能把每个95%区间都当作独立的确认性结论。\n"
       <<"这些结果仅来自教学近似模型，未经现场数据标定，不是工程防护设计依据。\n";
    if(!terrainStudy)out<<"参数研究仅改变指定参数的均值，其余参数均值不变。\n";
}
}

inline int runParameterStudy(const Config& input,int count,uint32_t seed,const std::string& outputDir){
    if(count<2||count>5000)throw std::invalid_argument("Studies require 2 to 5000 runs per scenario");
    Config base=input;base.randomSurface=true;validate(base);
    struct Factor{std::string name;std::function<void(Config&,double)> scale;};
    const std::vector<Factor> factors={
        {"initial_speed",[](Config&c,double f){c.initialSpeed*=f;}},
        {"roll_resistance",[](Config&c,double f){c.material.rollResistance*=f;for(auto&m:c.terrain.surfaces)m.rollResistance*=f;}},
        {"friction",[](Config&c,double f){c.material.friction*=f;for(auto&m:c.terrain.surfaces)m.friction*=f;}},
        {"restitution",[](Config&c,double f){c.material.restitution*=f;for(auto&m:c.terrain.surfaces)m.restitution*=f;}},
        {"drag",[](Config&c,double f){c.material.drag*=f;for(auto&m:c.terrain.surfaces)m.drag*=f;}}
    };
    // Validate every requested mean before creating files: an invalid 1.2x
    // scenario must not leave a report that looks like a completed study.
    for(const auto& factor:factors)for(int percent:{80,100,120}){
        Config candidate=base;factor.scale(candidate,percent/100.);validate(candidate);
    }
    if(count*base.maxTime/base.sampleInterval>4e6)
        throw std::invalid_argument("Too many saved trajectory samples; reduce runs or increase sample interval");
    if(count*std::ceil(base.maxTime/nominalIntegrationStep(base))>2e8)
        throw std::invalid_argument("Too many integration steps; reduce runs or simulation duration");
    const std::filesystem::path root=pathFromUtf8(outputDir);
    if(std::filesystem::exists(root)&&(!std::filesystem::is_directory(root)||!std::filesystem::is_empty(root)))
        throw std::invalid_argument("Study output folder must be new or empty: "+outputDir);
    std::filesystem::create_directories(root);
    auto csv=study::open(root/"parameter_study.csv");study::summaryHeader(csv);
    auto report=study::open(root/"parameter_findings.txt");
    report<<"参数影响研究：单独改变一个参数均值，水平为基准的 0.8 / 1.0 / 1.2 倍\n";study::commonNotes(report);
    study::writeConfig(base,count,seed,root/"baseline");
    auto baseline=runBatch(base,count,seed);study::writeRunMetrics(baseline,root/"baseline");study::stripTrajectories(baseline);
    int scenarios=0;
    for(const auto&factor:factors){for(int percent:{80,100,120}){
        const std::string level=percent==80?"0.8":percent==100?"1.0":"1.2";
        const std::string name=factor.name+"_"+std::to_string(percent);
        if(percent==100){study::summaryRows(csv,factor.name,level,"baseline",baseline,baseline);++scenarios;continue;}
        Config c=base;factor.scale(c,percent/100.);validate(c);
        const auto dir=root/name;study::writeConfig(c,count,seed,dir);
        auto b=runBatch(c,count,seed);study::writeRunMetrics(b,dir);
        study::summaryRows(csv,factor.name,level,name,b,baseline);study::describeEffects(report,name,b,baseline);++scenarios;
    }}
    report<<"\nCSV 中水平1.0的各行均引用同一份已保存的基准试验，不重复运行相同场景。\n";
    csv.flush();report.flush();
    if(!csv||!report)throw std::runtime_error("Failed to flush study results");
    return scenarios;
}

#if ROCK_STAGE >= 4
inline int runTerrainStudy(const Config& input,int count,uint32_t seed,const std::string& outputDir){
    if(count<2||count>5000)throw std::invalid_argument("Studies require 2 to 5000 runs per scenario");
    // Coordinates are already scaled by the caller. The study compares edits
    // relative to those coordinates and keeps the collector at each new toe.
    Config base=input;base.randomSurface=true;base.detectorX=-1.;validate(base);
    if(base.terrain.points.size()<4)
        throw std::invalid_argument("Terrain studies require at least four vertices");
    struct Scenario{std::string factor,level,name;Config config;};
    std::vector<Scenario> scenarios;
    for(int percent:{80,120}){
        const double f=percent/100.;Config width=base,height=base;
        const double x0=base.terrain.points.front().x,z0=base.terrain.toeZ();
        for(auto& p:width.terrain.points)p.x=x0+(p.x-x0)*f;
        width.terrainWidthScale*=f;
        if(width.barrier.kind!=BarrierKind::None)width.barrier.x=x0+(width.barrier.x-x0)*f;
        for(auto& p:height.terrain.points)p.z=z0+(p.z-z0)*f;
        height.terrainHeightScale*=f;
        scenarios.push_back({"horizontal_extent",percent==80?"0.8":"1.2","width_"+std::to_string(percent),width});
        scenarios.push_back({"total_relief",percent==80?"0.8":"1.2","height_"+std::to_string(percent),height});
    }
    Config bench=base;
    const size_t a=std::max(size_t(1),base.terrain.points.size()/2-1),b=a+1;
    const double benchZ=.5*(bench.terrain.points[a].z+bench.terrain.points[b].z);
    bench.terrain.points[a].z=benchZ;bench.terrain.points[b].z=benchZ;
    scenarios.push_back({"local_bench","flat","middle_bench",bench});
    const size_t k=std::max(size_t(1),std::min(base.terrain.points.size()-2,base.terrain.points.size()/3));
    for(int percent:{30,70}){
        Config local=base;
        const double hi=local.terrain.points[k-1].z,lo=local.terrain.points[k+1].z;
        local.terrain.points[k].z=hi-(percent/100.)*(hi-lo);
        scenarios.push_back({"local_drop_share",percent==30?"0.3":"0.7","upper_drop_"+std::to_string(percent),local});
    }
    // Keep these illustrative presets identical to the CLI/UI presets. A
    // uniform terrain uses ONE sampled material per run, like the default
    // baseline; assigning one material per segment changes spatial randomness.
    const Material rockMaterial{.02,.45,.45,1.},soilMaterial{.08,.60,.25,1.},grassMaterial{.12,.75,.20,1.};
    for(const auto& preset:std::vector<std::pair<std::string,Material>>{{"rock",rockMaterial},{"soil",soilMaterial},{"grass",grassMaterial}}){
        Config surface=base;surface.material=preset.second;surface.terrain.surfaces.clear();
        scenarios.push_back({"surface",preset.first,"surface_"+preset.first,surface});
    }
    Config layered=base;layered.terrain.surfaces.clear();
    const size_t segments=layered.terrain.points.size()-1;
    for(size_t i=0;i<segments;++i)
        layered.terrain.surfaces.push_back(i*3<segments?rockMaterial:i*3<2*segments?soilMaterial:grassMaterial);
    scenarios.push_back({"surface","layered","surface_layered",layered});

    // Validate every candidate and the per-batch resource budgets before any
    // directories are created, including geometry overflow after scaling.
    for(const auto& scenario:scenarios)validate(scenario.config);
    if(count*base.maxTime/base.sampleInterval>4e6)
        throw std::invalid_argument("Too many saved trajectory samples; reduce runs or increase sample interval");
    if(count*std::ceil(base.maxTime/nominalIntegrationStep(base))>2e8)
        throw std::invalid_argument("Too many integration steps; reduce runs or simulation duration");
    const std::filesystem::path root=pathFromUtf8(outputDir);
    if(std::filesystem::exists(root)&&(!std::filesystem::is_directory(root)||!std::filesystem::is_empty(root)))
        throw std::invalid_argument("Study output folder must be new or empty: "+outputDir);
    std::filesystem::create_directories(root);
    auto csv=study::open(root/"terrain_study.csv");study::summaryHeader(csv);
    auto report=study::open(root/"terrain_findings.txt");
    report<<"地形影响研究：基准、水平范围、总高差、局部平台、局部坡度与地表材质，共12个场景\n";
    study::commonNotes(report,true);
    report<<"各场景均以当前传入地形为基准进行相对调整；统计断面跟随各自坡脚，覆盖用户指定的固定断面位置。\n"
          <<"水平范围0.8/1.2倍以起点x为中心缩放，同时改变传播距离和坡角；高差0.8/1.2倍以坡脚高程为基准缩放，同时改变初始势能和坡角。\n"
          <<"中部平台将两相邻内部顶点高程设为均值；局部落差比例0.3/0.7调整一个内部顶点，重新分配相邻两段落差。端点不变，局部坡角相互耦合。\n"
          <<"surface_rock/soil/grass采用与命令行相同的示意参数，每次试验全坡共享一组随机材质。\n"
          <<"surface_layered按坡段编号约分三组，依次为岩石、土、草，逐段独立抽样；分组并非按实际长度等分。\n"
          <<"若基准自带逐段材质，统一材质场景同时改变空间抽样结构；分层场景也同时改变材质分布和空间相关性，不能把差异仅归因于某个系数。\n"
          <<"材质预设未经现场标定。若设置拦挡，宽度场景同时缩放其水平位置；拦挡高度与弹性参数保持输入值。\n"
          <<"这些场景为相对比较，不应把几何耦合造成的变化解释成单独坡角的因果效应。\n";
    study::writeConfig(base,count,seed,root/"baseline");
    auto baseline=runBatch(base,count,seed);study::writeRunMetrics(baseline,root/"baseline");study::stripTrajectories(baseline);
    study::summaryRows(csv,"baseline","1.0","baseline",baseline,baseline);
    for(const auto& scenario:scenarios){
        const auto dir=root/scenario.name;study::writeConfig(scenario.config,count,seed,dir);
        auto batch=runBatch(scenario.config,count,seed);study::writeRunMetrics(batch,dir);
        study::summaryRows(csv,scenario.factor,scenario.level,scenario.name,batch,baseline);
        study::describeEffects(report,scenario.name,batch,baseline);
    }
    report<<"\n每个场景的实际坐标与材质均值已保存至terrain.csv；实际随机材质见sampled_surfaces.csv，试验指标见run_metrics.csv。\n"
          <<"CSV包含基准和11个调整场景；水平范围与高差的1.0倍共用同一份基准试验。\n";
    csv.flush();report.flush();
    if(!csv||!report)throw std::runtime_error("Failed to flush study results");
    return int(scenarios.size()+1);
}
#endif
}

#endif

// ==================== cli.hpp ====================

#include <chrono>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>

namespace rock {
inline std::string newResultPath(const std::string& kind){
    auto now=std::chrono::system_clock::now();std::time_t t=std::chrono::system_clock::to_time_t(now);
    std::tm stamp{};
#ifdef _WIN32
    localtime_s(&stamp,&t);
#else
    localtime_r(&t,&stamp);
#endif
    char buffer[64];std::strftime(buffer,sizeof(buffer),"%Y%m%d_%H%M%S",&stamp);
    std::string p="results/"+kind+"_"+buffer;int suffix=1;
    while(std::filesystem::exists(pathFromUtf8(p)))p="results/"+kind+"_"+buffer+"_"+std::to_string(suffix++);
    return p;
}
inline void cliHelp(){
    std::cout<<"Rollingstone "<<ROCK_VERSION<<" | SI units | local educational simulator\n"
      <<"No arguments: interactive canvas.\n"
      <<"--batch --runs 200 --seed 42 [--out results/name]\n"
      <<"--study : paired parameter experiments (5 parameters x 3 mean levels)\n"
      <<"--speed 1.5 --radius 2 --cv .1 --dt .002 --max-time 60 --sample-interval .02\n"
      <<"--roll .02 --friction .45 --restitution .45 --drag 1\n"
      <<"--material rock|soil|grass (illustrative surface presets)\n"
      <<"--detector-x 120 --bin-height .5 --grid-width 1\n"
      <<"CV=0 disables randomness. First forward detector crossing; 0..10m above toe.\n"
      <<"All simulation, statistics and graphics are C++.\n";
#if ROCK_STAGE >= 4
    std::cout<<"--terrain-study : paired comparison of 12 terrain/surface scenarios\n"
      <<"--width-scale 1 --height-scale 1 --terrain-file terrain.csv\n";
#endif
#if ROCK_STAGE >= 5
    std::cout<<"--shape sphere|polygon --sides 3..12 --phase radians\n";
#endif
#if ROCK_STAGE >= 6
    std::cout<<"--barrier none|rigid --barrier-x 112 --barrier-height 6 --barrier-restitution .15\n";
#endif
#if ROCK_STAGE >= 7
    std::cout<<"--barrier elastic --stiffness 2000000 --damping 200000 --max-deflection 3 --barrier-mass 1000\n";
#endif
}
inline int runCommandLine(int argc,char** argv){
    try{
        Config c;
#if ROCK_STAGE >= 3
        c.randomSurface=true;
#endif
        int count=200;uint32_t seed=42;std::string mode="batch",output,terrainFile;
        for(int i=1;i<argc;++i){
            std::string key=argv[i];
            if(key=="--help"||key=="-h"){cliHelp();return 0;}
            if(key=="--batch"){mode="batch";continue;}
            if(key=="--study"){mode="study";continue;}
#if ROCK_STAGE >= 4
            if(key=="--terrain-study"){mode="terrain";continue;}
#endif
            const std::vector<std::string> allowed={"--runs","--seed","--out","--speed","--radius","--cv","--dt","--max-time","--sample-interval","--roll","--friction","--restitution","--drag","--detector-x","--bin-height","--grid-width","--material"
#if ROCK_STAGE >= 4
                ,"--width-scale","--height-scale","--terrain-file"
#endif
#if ROCK_STAGE >= 5
                ,"--shape","--sides","--phase"
#endif
#if ROCK_STAGE >= 6
                ,"--barrier","--barrier-x","--barrier-height","--barrier-restitution"
#endif
#if ROCK_STAGE >= 7
                ,"--stiffness","--damping","--max-deflection","--barrier-mass"
#endif
            };
            if(std::find(allowed.begin(),allowed.end(),key)==allowed.end())throw std::invalid_argument("Unknown option in this version: "+key);
            if(i+1>=argc)throw std::invalid_argument("Missing value for "+key);
            std::string value=argv[++i];
            auto number=[&](){size_t end=0;double v=std::stod(value,&end);if(end!=value.size()||!std::isfinite(v))throw std::invalid_argument("Invalid number for "+key);return v;};
            auto integer=[&](){double n=number();if(n!=std::floor(n)||n<0||n>4294967295.)throw std::invalid_argument("Invalid integer for "+key);return static_cast<uint32_t>(n);};
            if(key=="--runs"){auto n=integer();if(n>5000)throw std::invalid_argument("runs <= 5000 required");count=int(n);}
            else if(key=="--seed")seed=integer();
            else if(key=="--out")output=value;
            else if(key=="--speed")c.initialSpeed=number();
            else if(key=="--radius")c.radius=number();
            else if(key=="--cv")c.cv=number();
            else if(key=="--dt")c.dt=number();
            else if(key=="--max-time")c.maxTime=number();
            else if(key=="--sample-interval")c.sampleInterval=number();
            else if(key=="--roll")c.material.rollResistance=number();
            else if(key=="--friction")c.material.friction=number();
            else if(key=="--restitution")c.material.restitution=number();
            else if(key=="--drag")c.material.drag=number();
            else if(key=="--width-scale")c.terrainWidthScale=number();
            else if(key=="--height-scale")c.terrainHeightScale=number();
            else if(key=="--detector-x")c.detectorX=number();
            else if(key=="--bin-height")c.binHeight=number();
            else if(key=="--grid-width")c.gridWidth=number();
            else if(key=="--terrain-file")terrainFile=value;
            else if(key=="--material"){
                if(value=="rock")c.material={.02,.45,.45,1.};
                else if(value=="soil")c.material={.08,.60,.25,1.};
                else if(value=="grass")c.material={.12,.75,.20,1.};
                else throw std::invalid_argument("Unknown material");
            }
            else if(key=="--shape"){
                if(value=="sphere")c.shape=Shape::Sphere;
                else if(value=="polygon")c.shape=Shape::Polygon;
                else throw std::invalid_argument("shape must be sphere or polygon");
            }
            else if(key=="--sides"){auto n=integer();if(n<3||n>12)throw std::invalid_argument("sides must be 3..12");c.sides=int(n);}
            else if(key=="--phase")c.shapePhase=number();
            else if(key=="--barrier"){
                if(value=="none")c.barrier.kind=BarrierKind::None;
                else if(value=="rigid")c.barrier.kind=BarrierKind::Rigid;
                else if(value=="elastic")c.barrier.kind=BarrierKind::Elastic;
                else throw std::invalid_argument("Unknown barrier");
            }
            else if(key=="--barrier-x")c.barrier.x=number();
            else if(key=="--barrier-height")c.barrier.height=number();
            else if(key=="--barrier-restitution")c.barrier.restitution=number();
            else if(key=="--stiffness")c.barrier.stiffness=number();
            else if(key=="--damping")c.barrier.damping=number();
            else if(key=="--max-deflection")c.barrier.maxDeflection=number();
            else if(key=="--barrier-mass")c.barrier.mass=number();
            else throw std::invalid_argument("Unknown option: "+key);
        }
#if ROCK_STAGE < 3
        if(mode=="study")throw std::invalid_argument("Parameter studies start at v1.3");
#endif
#if ROCK_STAGE < 4
        if(mode=="terrain"||!terrainFile.empty()||c.terrainWidthScale!=1||c.terrainHeightScale!=1)
            throw std::invalid_argument("Terrain studies start at v1.4");
#endif
        if(!terrainFile.empty()){
            std::ifstream f(pathFromUtf8(terrainFile));if(!f)throw std::invalid_argument("Cannot read terrain CSV");
            c.terrain.points.clear();std::string line;bool first=true,hasIndex=false;
            while(std::getline(f,line)){
                if(line.empty())continue;
                if(first && line.compare(0,3,"\xEF\xBB\xBF")==0)line.erase(0,3);
                if(first && line.rfind("point_index,",0)==0){hasIndex=true;first=false;continue;}
                if(first && (line[0]=='x'||line[0]=='X')){first=false;continue;}first=false;
                std::replace(line.begin(),line.end(),',',' ');std::istringstream row(line);Vec2 p;std::string extra;
                if(hasIndex){int index;if(!(row>>index))throw std::invalid_argument("Invalid terrain point index");}
                if(!(row>>p.x>>p.z)||(row>>extra))throw std::invalid_argument("Terrain rows require exactly x,z");
                c.terrain.points.push_back(p);
            }
        }
        if(!(c.terrainWidthScale>0 && c.terrainWidthScale<=5 && c.terrainHeightScale>0 && c.terrainHeightScale<=5))
            throw std::invalid_argument("Terrain scales must be in (0,5]");
        c.terrain.validate();
        const double toeZ=c.terrain.toeZ(),originX=c.terrain.points.front().x;
        for(auto& p:c.terrain.points){p.x=originX+(p.x-originX)*c.terrainWidthScale;p.z=toeZ+(p.z-toeZ)*c.terrainHeightScale;}
        validate(c);
#if ROCK_STAGE >= 3
        if(mode=="study"){if(output.empty())output=newResultPath("parameters");runParameterStudy(c,count,seed,output);std::cout<<"Saved: "<<output<<'\n';return 0;}
#endif
#if ROCK_STAGE >= 4
        if(mode=="terrain"){if(output.empty())output=newResultPath("terrain");runTerrainStudy(c,count,seed,output);std::cout<<"Saved: "<<output<<'\n';return 0;}
#endif
        if(!output.empty() && std::filesystem::exists(pathFromUtf8(output)) && !std::filesystem::is_empty(pathFromUtf8(output)))throw std::invalid_argument("Output folder must be new or empty: "+output);
        auto b=runBatch(c,count,seed);std::cout<<batchSummary(b);
#if ROCK_STAGE >= 2
        if(output.empty())output=newResultPath("batch");exportBatch(b,output);
        std::cout<<"Saved: "<<pathToUtf8(std::filesystem::absolute(pathFromUtf8(output)))<<'\n';
#else
        if(!output.empty())throw std::invalid_argument("CSV export starts at v1.2");
#endif
        return 0;
    }catch(const std::exception& e){std::cerr<<"Error: "<<e.what()<<'\n';return 2;}
}
}


// ==================== viewer.hpp ====================

// Native EasyX front end. The renderer reads the very same samples used by CSV export.
#ifndef NOMINMAX
#define NOMINMAX
#endif
#ifndef UNICODE
#define UNICODE
#endif
#ifndef _UNICODE
#define _UNICODE
#endif
#include <easyx.h>
#include <chrono>
#include <filesystem>
#include <future>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <ctime>

namespace rock { namespace ui {
constexpr int W=1300,H=850,SIDE=342;
const COLORREF BG=RGB(13,23,36),PANEL=RGB(22,36,52),INK=RGB(236,243,250),MUTED=RGB(178,195,211);
const COLORREF ACCENT=RGB(69,211,175),SKYBLUE=RGB(80,169,244),ORANGE=RGB(255,174,83);
struct Rect {int l,t,r,b; bool contains(int x,int y)const{return x>=l&&x<=r&&y>=t&&y<=b;}};
inline std::wstring wide(const std::string& s){
    if(s.empty())return L"";
    int n=MultiByteToWideChar(CP_UTF8,0,s.data(),static_cast<int>(s.size()),nullptr,0);
    if(n<=0)return L"无法显示的文本";
    std::wstring w(static_cast<size_t>(n),L' ');
    MultiByteToWideChar(CP_UTF8,0,s.data(),static_cast<int>(s.size()),&w[0],n);return w;
}
inline std::wstring num(double v,int precision=2){std::wostringstream s;s<<std::fixed<<std::setprecision(precision)<<v;return s.str();}
inline void fontStyle(int size){
    LOGFONT font{};
    // Negative height specifies the glyph height, rather than a cell with internal leading.
    font.lfHeight=-(std::max)(16,size);
    font.lfWeight=size>=19?FW_MEDIUM:FW_NORMAL;
    font.lfCharSet=DEFAULT_CHARSET;
    font.lfOutPrecision=OUT_TT_PRECIS;
    font.lfClipPrecision=CLIP_DEFAULT_PRECIS;
    font.lfQuality=ANTIALIASED_QUALITY;
    font.lfPitchAndFamily=DEFAULT_PITCH|FF_DONTCARE;
    wcscpy_s(font.lfFaceName,L"Microsoft YaHei UI");
    // Retain the font-object reuse optimization, including when switching IMAGE targets.
    LOGFONT current{};gettextstyle(&current);
    if(current.lfHeight==font.lfHeight&&current.lfWeight==font.lfWeight&&
       current.lfQuality==font.lfQuality&&current.lfCharSet==font.lfCharSet&&
       wcscmp(current.lfFaceName,font.lfFaceName)==0)return;
    settextstyle(&font);
}
inline void text(int x,int y,const std::wstring& s,int size=17,COLORREF color=INK){fontStyle(size);settextcolor(color);setbkmode(TRANSPARENT);outtextxy(x,y,s.c_str());}
inline void fittedText(int x,int y,int width,const std::wstring& s,int size=17,COLORREF color=INK){
    fontStyle(size);std::wstring visible=s;
    if(textwidth(visible.c_str())>width){
        while(!visible.empty()&&textwidth((visible+L"…").c_str())>width)visible.pop_back();
        if(!visible.empty()&&visible.back()>=0xD800&&visible.back()<=0xDBFF)visible.pop_back();
        visible+=L"…";
    }
    settextcolor(color);setbkmode(TRANSPARENT);outtextxy(x,y,visible.c_str());
}
inline void box(Rect r,COLORREF color){setfillcolor(color);solidrectangle(r.l,r.t,r.r,r.b);}
inline void lineAt(int x1,int y1,int x2,int y2,COLORREF c,int thickness=1){setlinecolor(c);setlinestyle(PS_SOLID,thickness);line(x1,y1,x2,y2);}
inline void button(Rect r,const std::wstring& label,bool selected=false,bool disabled=false){
    box(r,disabled?RGB(32,44,58):(selected?RGB(30,97,88):RGB(40,58,77)));
    fontStyle(16);settextcolor(disabled?MUTED:(selected?ACCENT:INK));setbkmode(TRANSPARENT);
    outtextxy(r.l+(r.r-r.l-textwidth(label.c_str()))/2,r.t+(r.b-r.t-textheight(label.c_str()))/2,label.c_str());
}
inline std::wstring stopName(const RunResult& r){
    if(r.status=="time_limit")return L"达到时间上限";
    if(r.status=="domain_exit")return L"离开计算域";
    if(r.status=="stopped")return L"自然停稳";
    if(r.status=="intercepted")return L"被拦挡截停";
    return wide(r.status);
}
struct Row {std::wstring name,value;};
struct StudyRow {
    std::string factor,level,metric,scenario;
    int valid=0;
    double mean=0,low=0,high=0;
    bool available=false;
};
struct StudyData {
    Config base;
    int count=0;
    uint32_t seed=0;
    std::string directory;
    std::vector<StudyRow> rows;
    std::map<std::string,std::vector<Vec2>> terrains;
    const StudyRow* find(const std::string& factor,const std::string& level,const std::string& metric)const{
        for(const auto& r:rows)if(r.factor==factor&&r.level==level&&r.metric==metric)return &r;
        return nullptr;
    }
};
inline StudyData readStudy(const Config& c,int count,uint32_t seed,const std::string& directory,bool terrain=false){
    StudyData data;data.base=c;data.count=count;data.seed=seed;data.directory=directory;
    std::ifstream file(pathFromUtf8(directory)/(terrain?"terrain_study.csv":"parameter_study.csv"));
    if(!file)throw std::runtime_error("Cannot read study CSV");
    std::string line;std::getline(file,line);
    while(std::getline(file,line)){
        std::vector<std::string> fields;std::istringstream input(line);std::string field;
        while(std::getline(input,field,','))fields.push_back(field);
        if(fields.size()<11)continue;
        StudyRow row;row.factor=fields[0];row.level=fields[1];row.scenario=fields[2];row.metric=fields[4];
        row.valid=std::stoi(fields[5]);row.available=!fields[6].empty()&&!fields[9].empty()&&!fields[10].empty();
        if(row.available){row.mean=std::stod(fields[6]);row.low=std::stod(fields[9]);row.high=std::stod(fields[10]);}
        data.rows.push_back(row);
    }
    // 行数校验不再硬编码 132/165：汇总表恒为“每个 (factor,level) 场景 × 每个指标”
    // 一行，场景数直接由同批 CSV 的行集合推导；日后增删指标或场景不会误报。
    {
        std::set<std::pair<std::string,std::string>> levels;
        for(const auto& row:data.rows)levels.insert({row.factor,row.level});
        const size_t expected=levels.size()*study::metrics().size();
        if(expected==0||data.rows.size()!=expected)throw std::runtime_error("Incomplete study table");
    }
    if(terrain)for(const auto& row:data.rows)if(data.terrains.find(row.scenario)==data.terrains.end()){
        std::ifstream coordinates(pathFromUtf8(directory)/row.scenario/"terrain.csv");
        if(!coordinates)throw std::runtime_error("Cannot read study terrain coordinates");
        std::getline(coordinates,line);std::vector<Vec2> points;
        while(std::getline(coordinates,line)){std::istringstream input(line);std::string index,x,z;
            if(std::getline(input,index,',')&&std::getline(input,x,',')&&std::getline(input,z,','))points.push_back({std::stod(x),std::stod(z)});}
        data.terrains.emplace(row.scenario,std::move(points));
    }
    return data;
}
struct View {
    Config pending;
    BatchResult batch;
    int count=200,page=0,chosenRow=0,trial=0,speedIndex=0,factorIndex=0;
    uint32_t seed=42;
    bool playing=true,dirty=false,studyDirty=false,analysisMode=false,terrainMode=false,terrainDirty=false,heat=false,overlay=true,scrubbing=false;
    bool energyMode=false,multiReplay=false;
    int terrainFactor=0;
    double time=0;
    int materialPreset=1;
    std::wstring note=L"默认运行 200 次随机试验，CV=0.1；画布与 CSV 使用同一批轨迹。";
    std::filesystem::path lastOutput;
    std::future<BatchResult> future;
    StudyData study;
    std::future<StudyData> studyFuture;
    StudyData terrainStudy;
    std::future<StudyData> terrainFuture;
    IMAGE sceneCache;
    bool sceneCacheValid=false,cachedOverlay=false,cachedHeat=false;
    const RunResult* cachedRuns=nullptr;
    size_t cachedRunCount=0;
    std::chrono::steady_clock::time_point jobStart;
    Rect plot{370,138,1266,518},timeline{376,564,1260,588};
    double worldLeft=0,worldRight=180,worldBottom=-4,worldTop=65;
    const double speeds[5]={.25,.5,1,2,4};
    View(){pending.randomSurface=true; speedIndex=2;}
    bool busy()const{return future.valid()||studyFuture.valid()||terrainFuture.valid();}
    void invalidateScene(){sceneCacheValid=false;}
    void setView(int mode){analysisMode=mode==1;terrainMode=mode==2;energyMode=mode==3;}
    int activeView()const{return energyMode?3:terrainMode?2:analysisMode?1:0;}
    static constexpr int viewCount(){return ROCK_STAGE>=8?4:3;}
    static constexpr int pageCount(){return ROCK_STAGE>=6?5:ROCK_STAGE>=5?4:3;}
    bool hasData()const{return !batch.runs.empty();}
    const RunResult& current()const{return batch.runs[static_cast<size_t>(trial)];}
    const Config& shown()const{return hasData()?batch.base:pending;}
    int replayCount()const{
#if ROCK_STAGE >= 8
        if(multiReplay)return (std::min)(20,static_cast<int>(batch.runs.size()));
#endif
        return hasData()?1:0;
    }
    double replayEnd()const{
        double end=0;for(int i=0;i<replayCount();++i){const auto& samples=batch.runs[(static_cast<size_t>(trial)+i)%batch.runs.size()].trajectory;
            if(!samples.empty())end=(std::max)(end,samples.back().t);}return end;
    }
    void toggleMulti(){multiReplay=!multiReplay;time=0;playing=true;setView(0);}
    Config config()const{
        Config c=pending;
#if ROCK_STAGE >= 4
        const double origin=c.terrain.points.front().x,toe=c.terrain.toeZ();
        for(auto& p:c.terrain.points){p.x=origin+(p.x-origin)*c.terrainWidthScale;p.z=toe+(p.z-toe)*c.terrainHeightScale;}
#endif
        c.terrain.surfaces.clear(); // The selected material applies uniformly to this scenario.
        return c;
    }
    std::vector<Row> rows()const{
        const auto& c=pending;
        if(page==0)return {{L"试验次数 N",num(count,0)},{L"随机种子",num(seed,0)},{L"初始速度 m/s",num(c.initialSpeed)},{L"外接半径 m",num(c.radius)},{L"变异系数 CV",num(c.cv)},{L"最长时长 s",num(c.maxTime,0)},{L"积分步长 ms",num(c.dt*1000,1)},{L"播放倍速",num(speeds[speedIndex])+L"×"}};
        if(page==1)return {{L"地表预设",materialPreset==1?L"岩石":materialPreset==2?L"土层":materialPreset==3?L"草地":L"自定义"},{L"滚动阻力",num(c.material.rollResistance,3)},{L"滑动摩擦",num(c.material.friction)},{L"恢复系数",num(c.material.restitution)},{L"空气阻力 Cd",num(c.material.drag)},{L"材质随机化",c.randomSurface?L"开":L"关"}};
#if ROCK_STAGE >= 5
        if(page==3)return {{L"几何形状",c.shape==Shape::Sphere?L"球体":L"凸多边形"},{L"多边形边数",num(c.sides,0)},{L"初始转角 °",num(c.shapePhase*180/PI,0)}};
#endif
#if ROCK_STAGE >= 6
        if(page==4){std::vector<Row> result={{L"拦挡类型",c.barrier.kind==BarrierKind::None?L"关闭":c.barrier.kind==BarrierKind::Rigid?L"刚性":L"弹性"},{L"位置 x / m",num(c.barrier.x,1)},{L"高度 m",num(c.barrier.height,1)},{L"恢复系数",num(c.barrier.restitution)}};
#if ROCK_STAGE >= 7
            if(c.barrier.kind==BarrierKind::Elastic){result[3]={L"等效质量 kg",num(c.barrier.mass,0)};
                result.push_back({L"刚度 kN/m",num(c.barrier.stiffness/1000,0)});result.push_back({L"阻尼 kN·s/m",num(c.barrier.damping/1000,0)});result.push_back({L"失效位移 m",num(c.barrier.maxDeflection,2)});}
#endif
            return result;}
#endif
#if ROCK_STAGE >= 4
        return {{L"水平宽度倍率",num(c.terrainWidthScale)},{L"垂直高差倍率",num(c.terrainHeightScale)},{L"高度分箱 m",num(c.binHeight)},{L"路径网格宽 m",num(c.gridWidth)}};
#else
        return {{L"高度分箱 m",num(c.binHeight)},{L"路径网格宽 m",num(c.gridWidth)}};
#endif
    }
    void adjust(int dir){
        if(busy())return;
        auto clamp=[](double v,double lo,double hi){return std::clamp(v,lo,hi);};
        auto& c=pending;
        if(page==0){switch(chosenRow){
            case 0:count=std::clamp(count+dir*25,1,1000);break;
            case 1:seed=dir<0?(seed>0?seed-1:0):seed+1;break;
            case 2:c.initialSpeed=clamp(c.initialSpeed+dir*.5,0,50);break;
            case 3:c.radius=clamp(c.radius+dir*.1,.2,3);break;
            case 4:c.cv=clamp(c.cv+dir*.01,0,.5);break;
            case 5:c.maxTime=clamp(c.maxTime+dir*10,10,180);break;
            case 6:c.dt=clamp(c.dt+dir*.0005,.0005,.005);break;
            case 7:speedIndex=std::clamp(speedIndex+dir,0,4);return;}}
        if(page==1){switch(chosenRow){
            case 0:materialPreset=(materialPreset+dir+4)%4;
                if(materialPreset==1)c.material={.02,.45,.45,1.};
                if(materialPreset==2)c.material={.08,.60,.25,1.};
                if(materialPreset==3)c.material={.12,.75,.20,1.};break;
            case 1:c.material.rollResistance=clamp(c.material.rollResistance+dir*.005,0,.3);materialPreset=0;break;
            case 2:c.material.friction=clamp(c.material.friction+dir*.025,0,1.5);materialPreset=0;break;
            case 3:c.material.restitution=clamp(c.material.restitution+dir*.025,0,.95);materialPreset=0;break;
            case 4:c.material.drag=clamp(c.material.drag+dir*.1,0,3);materialPreset=0;break;
            case 5:c.randomSurface=!c.randomSurface;break;}}
        if(page==2){switch(chosenRow){
#if ROCK_STAGE >= 4
            case 0:c.terrainWidthScale=clamp(c.terrainWidthScale+dir*.1,.5,2);break;
            case 1:c.terrainHeightScale=clamp(c.terrainHeightScale+dir*.1,.5,2);break;
            case 2:c.binHeight=clamp(c.binHeight+dir*.25,.25,2);break;
            case 3:c.gridWidth=clamp(c.gridWidth+dir*.5,.5,5);break;
#else
            case 0:c.binHeight=clamp(c.binHeight+dir*.25,.25,2);break;
            case 1:c.gridWidth=clamp(c.gridWidth+dir*.5,.5,5);break;
#endif
        }}
#if ROCK_STAGE >= 5
        if(page==3){switch(chosenRow){
            case 0:c.shape=c.shape==Shape::Sphere?Shape::Polygon:Shape::Sphere;break;
            case 1:c.sides=std::clamp(c.sides+dir,3,12);c.shape=Shape::Polygon;break;
            case 2:c.shapePhase=std::fmod(c.shapePhase+dir*PI/12+2*PI,2*PI);break;}}
#endif
#if ROCK_STAGE >= 6
        if(page==4){switch(chosenRow){
            case 0:c.barrier.kind=static_cast<BarrierKind>((static_cast<int>(c.barrier.kind)+dir+(ROCK_STAGE>=7?3:2))%(ROCK_STAGE>=7?3:2));break;
            case 1:c.barrier.x=clamp(c.barrier.x+dir*2.,5,400);break;
            case 2:c.barrier.height=clamp(c.barrier.height+dir*.5,.5,30);break;
            case 3:if(c.barrier.kind==BarrierKind::Elastic)c.barrier.mass=clamp(c.barrier.mass+dir*250.,250,10000);else c.barrier.restitution=clamp(c.barrier.restitution+dir*.05,0,1);break;
            case 4:c.barrier.stiffness=clamp(c.barrier.stiffness+dir*250000.,250000,10000000);break;
            case 5:c.barrier.damping=clamp(c.barrier.damping+dir*25000.,0,1000000);break;
            case 6:c.barrier.maxDeflection=clamp(c.barrier.maxDeflection+dir*.25,.25,20);break;}}
#endif
        dirty=true;studyDirty=true;terrainDirty=true;if(!hasData())invalidateScene();note=L"设置已修改，当前图为旧结果；B 重算轨迹 / A 参数研究 / G 地形研究。";
    }
    void begin(bool many){
        if(busy())return;
        Config c=config();if(!many)c.cv=0;
        try{validate(c);}catch(const std::exception& e){note=L"参数错误："+wide(e.what());return;}
        const int n=many?count:1;const uint32_t s=seed;
        jobStart=std::chrono::steady_clock::now();
        future=std::async(std::launch::async,[c,n,s]{return runBatch(c,n,s);});
        setView(0);
        note=many?L"正在批量计算；可继续查看现有回放。":L"正在计算固定参数单次轨迹。";
    }
    void beginStudy(){
        if(busy())return;
        if(count<2){note=L"参数分析每组至少需要 2 次试验；请增加 N。";return;}
        Config c=config();c.randomSurface=true;
        try{validate(c);}catch(const std::exception& e){note=L"参数错误："+wide(e.what());return;}
        const int n=count;const uint32_t s=seed;
        const std::string out=newResultPath("parameters_viewer");
        jobStart=std::chrono::steady_clock::now();
        studyFuture=std::async(std::launch::async,[c,n,s,out]{
            runParameterStudy(c,n,s,out);return readStudy(c,n,s,out);
        });
        setView(1);note=L"正在比较 5 个参数 × 3 个水平；完成后显示均值与 95% 置信区间。";
    }
    void beginTerrain(){
#if ROCK_STAGE >= 4
        if(busy())return;
        if(count<2){note=L"地形分析每组至少需要 2 次试验。";return;}
        Config c=config();c.randomSurface=true;c.detectorX=-1.;
        const int n=count;const uint32_t s=seed;const std::string out=newResultPath("terrain_viewer");
        jobStart=std::chrono::steady_clock::now();
        terrainFuture=std::async(std::launch::async,[c,n,s,out]{runTerrainStudy(c,n,s,out);return readStudy(c,n,s,out,true);});
        setView(2);note=L"正在比较地形；蓝色坡脚检测器随各场景坡脚移动，研究完成后显示图表。";
#endif
    }
    void poll(){
        if(future.valid()&&future.wait_for(std::chrono::milliseconds(0))==std::future_status::ready){
            try{batch=future.get();invalidateScene();trial=0;time=0;playing=true;dirty=false;note=L"计算完成。PageUp / PageDown 选择试验；拖动时间轴回放。";}
            catch(const std::exception& e){note=L"计算失败："+wide(e.what());}
        }
        if(studyFuture.valid()&&studyFuture.wait_for(std::chrono::milliseconds(0))==std::future_status::ready){
            try{study=studyFuture.get();studyDirty=false;setView(1);
                lastOutput=std::filesystem::absolute(pathFromUtf8(study.directory));
                note=L"参数对照完成；选择上方参数查看 0.8 / 1.0 / 1.2 倍结果，CSV 已自动保存。";}
            catch(const std::exception& e){note=L"参数分析失败："+wide(e.what());}
        }
        if(terrainFuture.valid()&&terrainFuture.wait_for(std::chrono::milliseconds(0))==std::future_status::ready){
            try{terrainStudy=terrainFuture.get();terrainDirty=false;setView(2);lastOutput=std::filesystem::absolute(pathFromUtf8(terrainStudy.directory));
                note=L"地形研究完成。曲线中的坡形直接读取同批 CSV；选择宽度、高差或局部坡形对照。";}
            catch(const std::exception& e){note=L"地形研究失败："+wide(e.what());}
        }
    }
    void exportData(){
        if(!hasData()||busy())return;
        try{
            std::string out=newResultPath("viewer");
            exportBatch(batch,out);lastOutput=std::filesystem::absolute(pathFromUtf8(out));
            note=L"当前画布对应的批量 CSV 已导出；点击“输出目录 O”查看文件。";
        }catch(const std::exception& e){note=L"导出失败："+wide(e.what());}
    }
    void openOutput(){if(!lastOutput.empty())ShellExecuteW(nullptr,L"open",lastOutput.c_str(),nullptr,nullptr,SW_SHOWNORMAL);}
    void selectTrial(int dir){if(!hasData())return;trial=(trial+dir+static_cast<int>(batch.runs.size()))%static_cast<int>(batch.runs.size());time=0;playing=false;}
    Sample sampleAt(const RunResult& result,double at)const{
        if(result.trajectory.empty())return {};
        const auto& a=result.trajectory;
        auto it=std::lower_bound(a.begin(),a.end(),at,[](const Sample& p,double t){return p.t<t;});
        if(it==a.begin())return *it;if(it==a.end())return a.back();if(it->t==at)return *it;
        const auto& p=*(it-1);const auto& q=*it;double f=(at-p.t)/(q.t-p.t);
        Sample r=p;double* dst[]={&r.x,&r.z,&r.vx,&r.vz,&r.angle,&r.omega,&r.energy,&r.deflection};
        const double pv[]={p.x,p.z,p.vx,p.vz,p.angle,p.omega,p.energy,p.deflection};
        const double qv[]={q.x,q.z,q.vx,q.vz,q.angle,q.omega,q.energy,q.deflection};
        for(int i=0;i<8;++i)*dst[i]=pv[i]+f*(qv[i]-pv[i]);
#if ROCK_STAGE >= 6
        r.barrierVelocity=p.barrierVelocity+f*(q.barrierVelocity-p.barrierVelocity);
        r.barrierForce=p.barrierForce+f*(q.barrierForce-p.barrierForce);
        r.barrierDissipation=p.barrierDissipation+f*(q.barrierDissipation-p.barrierDissipation);
        r.contactImpulse=p.contactImpulse+f*(q.contactImpulse-p.contactImpulse);
#endif
        r.t=at;return r;
    }
    Sample sample()const{return hasData()?sampleAt(current(),time):Sample{};}
    int px(double x)const{return plot.l+static_cast<int>((x-worldLeft)/(worldRight-worldLeft)*(plot.r-plot.l));}
    int py(double z)const{return plot.b-static_cast<int>((z-worldBottom)/(worldTop-worldBottom)*(plot.b-plot.t));}
    void fit(){
        const auto& c=shown();worldLeft=c.terrain.points.front().x-5;worldRight=c.terrain.toeX()+35;worldBottom=c.terrain.toeZ()-4;worldTop=c.terrain.points.front().z+10;
        if(hasData())for(const auto& r:batch.runs){worldRight=(std::max)(worldRight,r.runoutX+5);for(const auto& p:r.trajectory)worldTop=(std::max)(worldTop,p.z+3);}
#if ROCK_STAGE >= 6
        if(c.barrier.kind!=BarrierKind::None){worldRight=(std::max)(worldRight,c.barrier.x+c.barrier.maxDeflection+5);worldTop=(std::max)(worldTop,c.terrain.height(c.barrier.x)+c.barrier.height+5);}
#endif
        // One metre has the same pixel length on both axes, preserving rock geometry.
        const double scale=(std::min)((plot.r-plot.l)/(worldRight-worldLeft),(plot.b-plot.t)/(worldTop-worldBottom));
        worldRight=worldLeft+(plot.r-plot.l)/scale;worldTop=worldBottom+(plot.b-plot.t)/scale;
    }
    void path(const RunResult& r,COLORREF color,double until=-1,int thickness=1){
        setlinecolor(color);setlinestyle(PS_SOLID,thickness);
        const auto& a=r.trajectory;
        for(size_t i=1;i<a.size();++i){if(until>=0&&a[i-1].t>until)break;line(px(a[i-1].x),py(a[i-1].z),px(a[i].x),py(a[i].z));}
    }
    void drawBody(const Config& rc,const Sample& p,COLORREF color){
#if ROCK_STAGE >= 5
        if(rc.shape==Shape::Polygon){
            std::vector<POINT> body;for(const auto& v:bodyVertices(rc,p.angle))body.push_back({px(p.x+v.x),py(p.z+v.z)});
            setfillcolor(color);setlinecolor(ORANGE);setlinestyle(PS_SOLID,1);fillpolygon(body.data(),static_cast<int>(body.size()));
            lineAt(px(p.x),py(p.z),px(p.x+rc.radius*std::cos(p.angle)),py(p.z+rc.radius*std::sin(p.angle)),BG,1);return;
        }
#endif
        const int rx=(std::max)(3,px(p.x+rc.radius)-px(p.x));
        const int ry=(std::max)(3,py(p.z)-py(p.z+rc.radius));
        setfillcolor(color);solidellipse(px(p.x)-rx,py(p.z)-ry,px(p.x)+rx,py(p.z)+ry);
        lineAt(px(p.x),py(p.z),px(p.x+rc.radius*std::cos(p.angle)),py(p.z+rc.radius*std::sin(p.angle)),BG,2);
    }
    void drawBarrier(const Config& c,const Sample& p){
#if ROCK_STAGE >= 6
        if(c.barrier.kind==BarrierKind::None)return;
        const auto& barrier=c.barrier;const bool elastic=barrier.kind==BarrierKind::Elastic;
        const double ground=c.terrain.height(barrier.x),displacement=elastic?p.deflection:0;
        const int base=py(ground),top=py(ground+barrier.height),anchor=px(barrier.x),moved=px(barrier.x+displacement);
        const COLORREF color=p.barrierFailed?RGB(245,91,91):elastic?SKYBLUE:ORANGE;
        if(elastic){
            lineAt(anchor,base,anchor,top,RGB(81,103,123),1);
            const int middle=(base+top)/2;
            for(int i=0;i<12;++i){const int xa=anchor+(moved-anchor)*i/12,xb=anchor+(moved-anchor)*(i+1)/12;
                const int ya=middle+(i==0?0:(i%2?5:-5)),yb=middle+(i==11?0:((i+1)%2?5:-5));lineAt(xa,ya,xb,yb,color,2);}
        }
        if(p.barrierFailed){const int middle=(base+top)/2;lineAt(moved,base,moved,middle+7,color,5);lineAt(moved+5,middle-7,moved+5,top,color,5);}
        else lineAt(moved,base,moved,top,color,5);
        const int labelX=(std::max)(plot.l+12,moved-150),labelY=(std::max)(plot.t+16,top-64);
        lineAt(labelX+83,labelY+18,moved,top,color,1);
        text(labelX,labelY,p.barrierFailed?L"拦挡已失效":elastic?L"弹性拦挡":L"刚性拦挡",13,color);
        if(elastic)text(plot.l+72,plot.t+11,L"当前挠度 "+num(p.deflection,3)+L" m  /  容量 "+num(barrier.maxDeflection,2)+L" m",13,color);
#else
        (void)c;(void)p;
#endif
    }
    void renderStaticScene(){
        const auto& c=shown();fit();box(plot,RGB(17,30,46));
        // Axes use world metres; all panels share the core's upward-positive elevation convention.
        for(double x=0;x<worldRight;x+=20){lineAt(px(x),plot.t,px(x),plot.b,RGB(33,49,65));text(px(x)-7,plot.b+1,num(x,0),13,MUTED);}
        for(double z=std::ceil(worldBottom/10)*10;z<worldTop;z+=10){lineAt(plot.l,py(z),plot.r,py(z),RGB(33,49,65));if(py(z)<plot.b-24)text(plot.l+4,py(z)+1,num(z,0),12,MUTED);}
        // Heat maps show P(a trial visits this cell), not dwell-time density.
        if(heat&&hasData())for(const auto& [cell,n]:batch.routeCells){
            const double p=static_cast<double>(n)/batch.runs.size();
            double x=cell.first*c.gridWidth,z=c.terrain.toeZ()+cell.second*c.binHeight;
            int left=(std::max)(plot.l,px(x)),right=(std::min)(plot.r,px(x+c.gridWidth));
            int top=(std::max)(plot.t,py(z+c.binHeight)),bottom=(std::min)(plot.b,py(z));
            if(left<right&&top<bottom)box({left,top,right,bottom},RGB(static_cast<int>(50+200*p),static_cast<int>(85+80*(1-p)),static_cast<int>(100*(1-p))));
        }
        const auto& pts=c.terrain.points;
        std::vector<POINT> land;land.push_back({px(worldLeft),plot.b});land.push_back({px(worldLeft),py(pts.front().z)});
        for(const auto& p:pts)land.push_back({px(p.x),py(p.z)});
        land.push_back({px(worldRight),py(c.terrain.toeZ())});land.push_back({px(worldRight),plot.b});
        setfillcolor(RGB(56,65,67));solidpolygon(land.data(),static_cast<int>(land.size()));
        for(size_t i=1;i<pts.size();++i)lineAt(px(pts[i-1].x),py(pts[i-1].z),px(pts[i].x),py(pts[i].z),ACCENT,3);
        lineAt(px(c.terrain.toeX()),py(c.terrain.toeZ()),plot.r,py(c.terrain.toeZ()),ACCENT,3);
        const double detector=c.detectorX<0?c.terrain.toeX():c.detectorX;
        lineAt(px(detector),py(c.terrain.toeZ()),px(detector),py(c.terrain.toeZ()+10),SKYBLUE,3);
        text(px(detector)-48,py(c.terrain.toeZ()+10)-24,L"坡脚 0–10 m",14,SKYBLUE);
        if(overlay&&hasData()){
            // 分层叠加：旧版对每条试验绘制全部线段（默认 N=200、0.02 s 采样约
            // 60 万条 line() 调用），重建场景缓存需要数秒。这里跳过两端点落在同一
            // 屏幕像素的线段（路径形状不变），超出上限的试验按等间隔抽样；
            // 当前试验仍由 renderScene 完整绘制，热力图与 CSV 统计不受影响。
            const size_t total=batch.runs.size();
            const size_t cap=80;
            const size_t stride=total<=cap?1:(total+cap-1)/cap;
            setlinecolor(RGB(160,85,78));setlinestyle(PS_SOLID,1);
            for(size_t runIndex=0;runIndex<total;runIndex+=stride){
                const auto& samples=batch.runs[runIndex].trajectory;
                int lastX=0,lastY=0;bool have=false;
                for(const auto& p:samples){
                    const int cx=px(p.x),cy=py(p.z);
                    if(have&&(cx!=lastX||cy!=lastY))line(lastX,lastY,cx,cy);
                    lastX=cx;lastY=cy;have=true;
                }
            }
        }
        text(plot.l+4,plot.t+5,L"高程 m",12,MUTED);
        text(plot.r-151,plot.t+9,L"水平距离 x / m",13,MUTED);
    }
    void renderScene(){
        // Completed trajectories and their fitted world coordinates stay fixed during playback.
        // Only a new batch or a display-layer change rebuilds this GDI image.
        if(!sceneCacheValid||cachedRuns!=batch.runs.data()||cachedRunCount!=batch.runs.size()||cachedOverlay!=overlay||cachedHeat!=heat){
            if(sceneCache.getwidth()!=W||sceneCache.getheight()!=H)sceneCache.Resize(W,H);
            IMAGE* previous=GetWorkingImage();SetWorkingImage(&sceneCache);
            try{setbkcolor(BG);cleardevice();renderStaticScene();}
            catch(...){SetWorkingImage(previous);throw;}
            SetWorkingImage(previous);
            cachedRuns=batch.runs.data();cachedRunCount=batch.runs.size();cachedOverlay=overlay;cachedHeat=heat;sceneCacheValid=true;
        }
        putimage(plot.l,plot.t,plot.r-plot.l+1,plot.b-plot.t+26,&sceneCache,plot.l,plot.t);
        if(hasData()){
            const auto p=sample();path(current(),ORANGE,time,2);
            const auto& rc=current().config;
#if ROCK_STAGE >= 8
            if(multiReplay){for(int i=replayCount()-1;i>0;--i){const auto& run=batch.runs[(static_cast<size_t>(trial)+i)%batch.runs.size()];
                    const COLORREF color=RGB(90+(i*43)%145,120+(i*29)%115,120+(i*71)%115);
                    drawBody(run.config,sampleAt(run,time),color);}
                text(plot.l+72,plot.t+38,L"同步回放 "+num(replayCount(),0)+L" / "+num(static_cast<double>(batch.runs.size()),0)+L" 条独立试验（不含颗粒间碰撞）",13,ACCENT);
                if(rc.barrier.kind!=BarrierKind::None)text(plot.l+72,plot.t+62,L"拦挡状态对应白色高亮试验；各试验相互独立",12,MUTED);
            }
#endif
            drawBody(rc,p,RGB(224,229,233));drawBarrier(rc,p);
        }
    }
    void sidebar(){
        box({14,94,328,806},PANEL);text(28,105,L"试验设置",19);
        const wchar_t* tabs[]={L"实验",L"地表",ROCK_STAGE>=4?L"地形":L"统计",L"形状",L"拦挡"};
        const int tabWidth=300/pageCount();
        for(int i=0;i<pageCount();++i)button({22+i*tabWidth,139,16+(i+1)*tabWidth,171},tabs[i],page==i,busy());
        const auto list=rows();
        for(int i=0;i<static_cast<int>(list.size());++i){int y=186+i*43;
            if(i==chosenRow)box({22,y-1,320,y+38},RGB(31,53,70));
            text(29,y+8,list[i].name,15);text(171,y+8,list[i].value,15,ACCENT);
            button({265,y+2,289,y+32},L"−",false,busy());button({292,y+2,316,y+32},L"+",false,busy());}
        button({24,555,116,593},L"单次 Enter",false,busy());button({121,555,214,593},L"批量 B",true,busy());button({219,555,316,593},L"导出 E",false,busy()||!hasData());
        button({24,605,166,646},L"参数研究 A",false,busy());
#if ROCK_STAGE >= 4
        button({174,605,316,646},L"地形研究 G",false,busy());
#endif
        text(27,660,L"↑↓ 选参数   ←→ 调整   Tab 切页",14,MUTED);
        text(27,685,L"PageUp / PageDown 选择试验",14,MUTED);
        text(27,710,L"空格 暂停   J / L 前后定位",14,MUTED);
        text(27,738,ROCK_STAGE>=8?L"M 多颗同步   F 能量曲线":L"批量按 CV 抽样；单次使用均值",14,MUTED);
        text(27,763,L"V 切换视图   O 输出目录",14,MUTED);
        text(27,784,L"画布、CSV、研究使用同一计算核心",12,ACCENT);
        if(page==2){text(29,398,L"蓝线：坡脚石心高度 0–10 m",14,SKYBLUE);
            text(29,428,L"分箱概率 = 本箱次数 / 全部 N 次",14,MUTED);
            text(29,458,L"每条轨迹仅统计首次向右通过",14,MUTED);
            text(29,488,L"调倍率后按 B；对照研究按 G",14,MUTED);}
#if ROCK_STAGE >= 5
        if(page==3){text(29,365,L"边数 3–12；轮廓即碰撞几何",14,ACCENT);
            text(29,395,L"初始转角也会影响接触与弹跳",14,MUTED);
            text(29,425,L"外接半径在“实验”页调整",14,MUTED);
            text(29,455,L"改变形状后按 B 重新模拟",14,MUTED);}
#endif
#if ROCK_STAGE >= 6
        if(page==4&&pending.barrier.kind!=BarrierKind::Elastic){text(29,405,L"有限高度：可拦截，也可能越顶",14,ACCENT);
            text(29,435,L"拦挡位置为实际世界坐标 x",14,MUTED);
            text(29,465,L"与滚石轮廓接触时改变速度",14,MUTED);}
#endif
    }
    void histogram(){
        Rect area{370,617,826,806};box(area,PANEL);text(385,627,L"首次通过坡脚断面的路线概率",17);
        if(!hasData())return;
        const int n=static_cast<int>(batch.heightCounts.size());if(n==0)return;
        int left=400,top=666,bottom=763,right=810;double maxP=.05;
        for(int v:batch.heightCounts)maxP=(std::max)(maxP,static_cast<double>(v)/batch.runs.size());
        for(int i=0;i<n;++i){double p=static_cast<double>(batch.heightCounts[i])/batch.runs.size();
            int x0=left+static_cast<int>((std::min)(10.,i*batch.base.binHeight)/10*(right-left));
            int x1=left+static_cast<int>((std::min)(10.,(i+1)*batch.base.binHeight)/10*(right-left))-2;
            if(p>0)box({x0,bottom-static_cast<int>((bottom-top)*p/maxP),(std::max)(x0+1,x1),bottom},SKYBLUE);}
        lineAt(left,bottom,right,bottom,MUTED);text(379,652,num(maxP*100,1)+L"%",12,MUTED);
        for(int i=0;i<=5;++i)text(left+i*(right-left)/5-5,bottom+1,num(i*2,0),12,MUTED);
        text(570,782,L"高度 m · 概率分母 N",12,MUTED);
    }
    void statistics(){
        box({842,617,1266,806},PANEL);text(858,627,L"本次试验记录",17);
        if(!hasData())return;
        const auto& r=current();const size_t n=batch.runs.size();
        text(858,658,L"0–10 m 内通过："+num(static_cast<double>(batch.inBand),0)+L" / "+num(static_cast<double>(n),0)+L"  ("+num(100.*batch.inBand/n,1)+L"%)",15,ACCENT);
        text(858,683,L"未到达 "+num(static_cast<double>(batch.notReached),0)+L"    超高 "+num(static_cast<double>(batch.aboveBand),0)+L"    低于范围 "+num(static_cast<double>(batch.belowBand),0),14,MUTED);
        text(858,708,L"当前试验："+stopName(r)+L"，最大速度 "+num(r.maxSpeed,1)+L" m/s",14);
#if ROCK_STAGE >= 6
        if(batch.base.barrier.kind!=BarrierKind::None){
            int hits=0,stops=0,failures=0;for(const auto& run:batch.runs){hits+=run.finalState.barrierHit?1:0;stops+=run.status=="intercepted"?1:0;failures+=run.finalState.barrierFailed?1:0;}
            text(858,733,L"全部试验：命中 "+num(hits,0)+L"   截停 "+num(stops,0)+L"   失效 "+num(failures,0),14,ACCENT);
            if(batch.base.barrier.kind==BarrierKind::Rigid)text(858,758,L"选中试验：累计接触冲量 "+num(r.finalState.contactImpulse/1000.,1)+L" kN·s",13);
            else text(858,758,L"最大挠度 "+num(r.finalState.maxBarrierDeflection,2)+L" m / 峰值支承力 "+num(r.finalState.peakBarrierForce/1000.,0)+L" kN",13);
        }else
#endif
        {text(858,733,r.reached?L"通过高度 "+num(r.detectorHeight)+L" m    速度 "+num(r.detectorSpeed,1)+L" m/s":L"当前试验未到达坡脚检测断面",14);
        text(858,758,r.reached?L"坡脚动能 "+num(r.detectorEnergy/1000000.,3)+L" MJ    首次通过 t="+num(r.detectorTime,1)+L" s":L"条件速度和动能不计入未到达试验",14);}
        text(858,781,L"seed="+num(batch.seed,0)+L"   CV="+num(batch.base.cv)+L"   N="+num(static_cast<double>(n),0),12,MUTED);
    }
    int comparisonCount(bool terrain)const{return terrain&&terrainFactor==2?2:3;}
    std::wstring comparisonLabel(int index,bool terrain)const{
        if(!terrain||terrainFactor<2)return index==0?L"0.8×":index==1?L"原坡形":L"1.2×";
        if(terrainFactor==2)return index==0?L"原坡形":L"平台";
        return index==0?L"落差30%":index==1?L"原坡形":L"落差70%";
    }
    const StudyRow* comparisonRow(int index,int metricIndex,bool terrain)const{
        const char* factors[]={"initial_speed","roll_resistance","friction","restitution","drag"};
        const char* metrics[]={"collector_speed_given_reached_m_s","collector_energy_given_reached_J","height_0_10m_probability"};
        const char* levels[]={"0.8","1.0","1.2"};
        if(!terrain)return study.find(factors[factorIndex],levels[index],metrics[metricIndex]);
        if((terrainFactor!=2&&index==1)||(terrainFactor==2&&index==0))return terrainStudy.find("baseline","1.0",metrics[metricIndex]);
        if(terrainFactor<2)return terrainStudy.find(terrainFactor==0?"horizontal_extent":"total_relief",levels[index],metrics[metricIndex]);
        if(terrainFactor==2)return terrainStudy.find("local_bench","flat",metrics[metricIndex]);
        return terrainStudy.find("local_drop_share",index==0?"0.3":"0.7",metrics[metricIndex]);
    }
    void studyChart(Rect area,int metricIndex,bool terrain=false){
        const wchar_t* titles[]={L"到达样本：坡脚速度",L"到达样本：坡脚动能",L"全部试验：0–10 m 通过率"};
        const wchar_t* units[]={L"m/s",L"MJ",L"%"};
        const double multiplier=metricIndex==1?1e-6:metricIndex==2?100.:1.;
        box(area,PANEL);text(area.l+12,area.t+12,titles[metricIndex],14);
        double maximum=metricIndex==2?100.:.001;
        const int amount=comparisonCount(terrain);
        for(int j=0;j<amount;++j){const auto* r=comparisonRow(j,metricIndex,terrain);
            if(metricIndex!=2&&r&&r->available)maximum=(std::max)(maximum,(std::max)(r->mean,r->high)*multiplier*1.2);}
        const int top=area.t+94,bottom=area.b-54;
        auto yy=[&](double value){return bottom-static_cast<int>((bottom-top)*std::clamp(value/maximum,0.,1.));};
        text(area.l+12,area.t+41,L"0 – "+num(maximum,metricIndex==2?0:2)+L" "+units[metricIndex],12,MUTED);
        lineAt(area.l+18,bottom,area.r-14,bottom,MUTED);
        for(int j=0;j<amount;++j){int center=area.l+(area.r-area.l)*(j+1)/(amount+1);
            const auto* r=comparisonRow(j,metricIndex,terrain);
            if(r&&r->available){
                const double value=r->mean*multiplier;
                box({center-21,yy(value),center+21,bottom},j==1?ACCENT:j==0?SKYBLUE:ORANGE);
                int high=yy(r->high*multiplier),low=yy(r->low*multiplier);
                lineAt(center,high,center,low,INK,2);lineAt(center-7,high,center+7,high,INK,2);lineAt(center-7,low,center+7,low,INK,2);
                text(center-30,(std::max)(top-24,high-25),num(value,metricIndex==2?1:2),13);
            }else text(center-22,top+60,L"无到达",12,MUTED);
            text(center-28,bottom+10,!terrain?(j==0?L"0.8×":j==1?L"1.0×":L"1.2×"):comparisonLabel(j,true),12,j==1?ACCENT:INK);
            if(r)text(center-20,bottom+31,L"n="+num(r->valid,0),10,MUTED);
        }
    }
    void renderTerrain(){
        box({370,138,1266,798},RGB(17,30,46));
        const wchar_t* names[]={L"水平宽度",L"垂直高差",L"中部平台",L"局部坡角"};
        for(int i=0;i<4;++i)button({380+i*220,149,590+i*220,184},names[i],terrainFactor==i);
        if(terrainStudy.rows.empty()){
            text(410,280,terrainFuture.valid()?L"正在计算地形对照…":L"点击“地形研究 G”开始",23,ACCENT);
            text(410,330,L"宽度 / 高差取 0.8、1.0、1.2 倍，并比较平台及局部坡角。",17,MUTED);
            text(410,370,L"每个场景的检测器跟随坡脚；同一组随机种子和试验次数。",17,MUTED);return;
        }
        text(385,197,L"地形研究   每组 N="+num(terrainStudy.count,0)+L"   seed="+num(terrainStudy.seed,0)+L"   CV="+num(terrainStudy.base.cv)+L"   竖线：95% 区间",15,MUTED);
        text(385,223,terrainDirty?L"设置已修改：以下为上次地形研究，请按 G 重算。":L"各场景检测断面随坡脚移动；比较的是地形变化的综合影响。",14,terrainDirty?ORANGE:ACCENT);
        studyChart({380,253,665,549},0,true);studyChart({674,253,959,549},1,true);studyChart({968,253,1256,549},2,true);
        text(385,561,L"实际坡形对照（读取各场景 terrain.csv，横纵坐标使用米）",14);
        const Rect profile{395,596,1240,727};double minX=0,maxX=1,minZ=0,maxZ=1;
        for(int j=0;j<comparisonCount(true);++j){const auto* row=comparisonRow(j,2,true);if(!row)continue;
            const auto it=terrainStudy.terrains.find(row->scenario);if(it==terrainStudy.terrains.end())continue;
            for(const auto& point:it->second){minX=(std::min)(minX,point.x);maxX=(std::max)(maxX,point.x);minZ=(std::min)(minZ,point.z);maxZ=(std::max)(maxZ,point.z);}}
        auto xx=[&](double x){return profile.l+static_cast<int>((x-minX)/(maxX-minX)*(profile.r-profile.l));};
        auto zz=[&](double z){return profile.b-static_cast<int>((z-minZ)/(maxZ-minZ)*(profile.b-profile.t));};
        for(int j=0;j<comparisonCount(true);++j){const auto* row=comparisonRow(j,2,true);if(!row)continue;
            const auto it=terrainStudy.terrains.find(row->scenario);if(it==terrainStudy.terrains.end())continue;
            COLORREF color=j==0?SKYBLUE:j==1?ACCENT:ORANGE;const auto& points=it->second;
            for(size_t k=1;k<points.size();++k)lineAt(xx(points[k-1].x),zz(points[k-1].z),xx(points[k].x),zz(points[k].z),color,2);
            text(685+j*180,582,comparisonLabel(j,true),12,color);
        }
        text(393,734,L"0 m",11,MUTED);text(1160,734,num(maxX,1)+L" m",11,MUTED);
        text(385,758,L"速度和动能只统计到达样本；通过率分母为全部 N 次。",13,MUTED);
        text(385,779,L"宽度、高差也会改变坡角；平台和局部坡角会联动邻段。",12,MUTED);
    }
    void renderStudy(){
        box({370,138,1266,798},RGB(17,30,46));
        const wchar_t* names[]={L"初始速度",L"滚动阻力",L"滑动摩擦",L"恢复系数",L"空气阻力"};
        const char* factors[]={"initial_speed","roll_resistance","friction","restitution","drag"};
        for(int i=0;i<5;++i)button({380+i*175,149,548+i*175,184},names[i],factorIndex==i);
        if(study.rows.empty()){
            text(410,280,studyFuture.valid()?L"正在计算参数对照… 画布仍可操作":L"点击左侧“参数对照 A”开始研究",23,ACCENT);
            text(410,330,L"逐个参数取 0.8 / 1.0 / 1.2 倍均值；其他参数保持基准。",17,MUTED);
            text(410,368,L"每组 N 次，使用相同种子；完成后显示均值和 95% 置信区间。",17,MUTED);
            text(410,406,L"统计 CSV 与图表来自同一参数研究核心。",17,MUTED);return;
        }
        text(385,197,L"每组 N="+num(study.count,0)+L"   seed="+num(study.seed,0)+L"   CV="+num(study.base.cv)+L"   竖线：95% 置信区间",15,MUTED);
        text(385,223,studyDirty?L"设置已修改：以下仍是上次研究，请按 A 重新计算。":L"0.8× / 1.0× / 1.2× 对照；仅改变当前参数的均值。",14,studyDirty?ORANGE:ACCENT);
        studyChart({380,253,665,549},0);studyChart({674,253,959,549},1);studyChart({968,253,1256,549},2);
        const double bases[]={study.base.initialSpeed,study.base.material.rollResistance,study.base.material.friction,study.base.material.restitution,study.base.material.drag};
        text(386,563,std::wstring(names[factorIndex])+L"的基准均值 = "+num(bases[factorIndex],3)+L"；下表为均值 [95% 区间]",15);
        text(387,596,L"水平 / 均值",13,MUTED);text(535,596,L"速度 m/s",13,MUTED);text(789,596,L"动能 MJ",13,MUTED);text(1040,596,L"通过率 %",13,MUTED);
        const char* levels[]={"0.8","1.0","1.2"};
        const char* metrics[]={"collector_speed_given_reached_m_s","collector_energy_given_reached_J","height_0_10m_probability"};
        for(int j=0;j<3;++j){int y=624+j*31;
            text(387,y,num(.8+j*.2,1)+L" / "+num(bases[factorIndex]*(.8+j*.2),3),13,j==1?ACCENT:INK);
            for(int k=0;k<3;++k){const auto* r=study.find(factors[factorIndex],levels[j],metrics[k]);const double scale=k==1?1e-6:k==2?100.:1.;
                std::wstring value=r&&r->available?num(r->mean*scale,k==2?1:2)+L" ["+num(r->low*scale,k==2?1:2)+L", "+num(r->high*scale,k==2?1:2)+L"]":L"无到达样本";
                text(535+k*252,y,value,12,j==1?ACCENT:INK);}
        }
        text(385,727,L"速度 / 动能仅统计到达者；通过率分母为全部 N 次（未到达计为未通过）。",13,MUTED);
        text(385,751,L"正态抽样，越界重采样；CV = 标准差 / 均值。完整数据及配对效应已写入 CSV。",13,MUTED);
        text(385,775,L"点击上方五个参数切换对照图；V 返回批量轨迹。",13,ACCENT);
    }
    void renderEnergy(){
#if ROCK_STAGE >= 8
        box({370,138,1266,798},RGB(17,30,46));text(385,151,L"当前试验的能量随时间变化",20);
        if(!hasData()||current().trajectory.empty()){text(405,290,L"先运行 B 生成轨迹，再查看能量曲线。",20,MUTED);return;}
        const auto& run=current();const auto& samples=run.trajectory;const double mass=run.finalState.mass,toe=run.config.terrain.toeZ();
        const Rect graph{425,220,1238,496};const double end=(std::max)(.001,samples.back().t);
        double lower=0,upper=1;
        for(const auto& point:samples){const double potential=mass*G*(point.z-toe)/1e6,total=point.energy/1e6+potential;
            lower=(std::min)(lower,potential);upper=(std::max)(upper,(std::max)(point.energy/1e6,total));}
        upper*=1.08;lower=lower<0?lower*1.08:0;
        auto xx=[&](double t){return graph.l+static_cast<int>(t/end*(graph.r-graph.l));};
        auto yy=[&](double e){return graph.b-static_cast<int>((e-lower)/(upper-lower)*(graph.b-graph.t));};
        for(int i=0;i<=5;++i){const double t=end*i/5.,e=lower+(upper-lower)*i/5.;
            lineAt(xx(t),graph.t,xx(t),graph.b,RGB(40,56,74));text(xx(t)-12,graph.b+6,num(t,1),11,MUTED);
            lineAt(graph.l,yy(e),graph.r,yy(e),RGB(40,56,74));text(375,yy(e)-7,num(e,1),11,MUTED);}
        const COLORREF colors[]={SKYBLUE,ORANGE,ACCENT};const wchar_t* names[]={L"动能",L"重力势能",L"滚石机械能"};
        for(int type=0;type<3;++type){text(600+type*180,180,names[type],14,colors[type]);
            for(size_t i=1;i<samples.size();++i){auto energy=[&](const Sample& p){double potential=mass*G*(p.z-toe)/1e6;return type==0?p.energy/1e6:type==1?potential:p.energy/1e6+potential;};
                lineAt(xx(samples[i-1].t),yy(energy(samples[i-1])),xx(samples[i].t),yy(energy(samples[i])),colors[type],2);}}
        text(385,191,L"能量 / MJ",13,MUTED);text(1184,516,L"时间 / s",12,MUTED);
        const double at=(std::min)(time,end);lineAt(xx(at),graph.t,xx(at),graph.b,INK,1);
        const auto p=sample();const double kinetic=p.energy/1e6,potential=mass*G*(p.z-toe)/1e6;
        text(385,539,L"试验 "+num(trial+1,0)+L"   t="+num(at)+L" s   速度 "+num(std::hypot(p.vx,p.vz),2)+L" m/s   结局："+stopName(run),15);
        box(timeline,RGB(36,55,74));const double playbackEnd=replayEnd();const int cursor=timeline.l+static_cast<int>((timeline.r-timeline.l)*(playbackEnd>0?time/playbackEnd:0));
        box({timeline.l,timeline.t,cursor,timeline.b},RGB(38,117,113));setfillcolor(ACCENT);solidcircle(cursor,(timeline.t+timeline.b)/2,7);
        text(385,598,L"拖动时间轴定位；PageUp / PageDown 选择试验。所有值取自同一条 CSV 轨迹。",13,MUTED);
        box({382,628,1253,786},PANEL);
        text(400,643,L"当前动能 "+num(kinetic,3)+L" MJ",18,SKYBLUE);text(680,643,L"当前势能 "+num(potential,3)+L" MJ",18,ORANGE);
        text(961,643,L"机械能 "+num(kinetic+potential,3)+L" MJ",18,ACCENT);
        text(400,681,L"势能以坡脚高程为零点；机械能 = 平动与转动动能 + 重力势能。",14,MUTED);
        if(run.config.barrier.kind==BarrierKind::Elastic){const auto& b=run.config.barrier;
            const double stored=p.barrierFailed?0:(.5*b.stiffness*p.deflection*p.deflection+.5*b.mass*p.barrierVelocity*p.barrierVelocity)/1e6;
            text(400,715,L"弹性板储能 "+num(stored,3)+L" MJ    拦挡耗散 "+num(p.barrierDissipation/1e6,3)+L" MJ    当前力 "+num(p.barrierForce/1000.,1)+L" kN",14);
        }else text(400,715,L"接触碰撞、摩擦与空气阻力会耗能；曲线显示计算结果，不假设机械能恒定。",14,MUTED);
        text(400,751,L"F 进入能量视图；M 返回多颗同步回放；V 切换全部视图。",14,ACCENT);
#endif
    }
    void render(){
        setbkcolor(BG);cleardevice();
        text(24,20,L"山区滚石 · 随机试验",29);text(26,60,L"路线概率 / CSV 导出 / 参数研究 / 地形研究",15,MUTED);
#if ROCK_STAGE >= 8
        button({757,58,911,89},L"多颗同步 M",multiReplay);button({921,58,1071,89},L"能量曲线 F",energyMode);
#endif
        if(busy()){
            double elapsed=std::chrono::duration<double>(std::chrono::steady_clock::now()-jobStart).count();
            text(650,28,(terrainFuture.valid()?L"地形分析中… ":studyFuture.valid()?L"参数分析中… ":L"批量计算中… ")+num(elapsed,1)+L" s",19,ORANGE);
        }else {const bool stale=terrainMode?terrainDirty:analysisMode?studyDirty:dirty;
            text(650,28,stale?L"当前显示旧结果，设置尚未重算":L"图表与 CSV 使用同一计算结果",16,stale?ORANGE:ACCENT);}
        button({1070,20,1265,58},L"输出目录 O",false,lastOutput.empty());
        sidebar();
        button({370,94,476,128},playing?L"暂停 Space":L"播放 Space",playing);
        button({484,94,579,128},L"重播 R");button({587,94,696,128},L"轨迹叠加 T",overlay);
        button({704,94,815,128},L"热力图 H",heat);
        button({825,94,906,128},L"上一条");button({914,94,996,128},L"下一条");
        const wchar_t* views[]={L"轨迹",L"参数",L"地形",L"能量"};const int viewWidth=259/viewCount();
        for(int i=0;i<viewCount();++i)button({1007+i*viewWidth,94,1007+(i+1)*viewWidth-6,128},views[i],activeView()==i);
        if(energyMode)renderEnergy();else if(terrainMode)renderTerrain();else if(analysisMode)renderStudy();else{
        renderScene();
        if(hasData()){
            auto p=sample();const double end=replayEnd();
            text(376,543,L"试验 "+num(trial+1,0)+L"/"+num(static_cast<double>(batch.runs.size()),0)+L"   t="+num(time)+L"/"+num(end)+L" s   速度 "+num(std::hypot(p.vx,p.vz))+L" m/s   动能 "+num(p.energy/1000.,1)+L" kJ",15);
            box(timeline,RGB(36,55,74));int xx=timeline.l+static_cast<int>((timeline.r-timeline.l)*(end>0?time/end:0));box({timeline.l,timeline.t,xx,timeline.b},RGB(38,117,113));
            setfillcolor(ACCENT);solidcircle(xx,(timeline.t+timeline.b)/2,7);text(376,593,L"拖动时间轴逐帧查看；热力图为每条轨迹是否访问网格的概率。",12,MUTED);
        }
        histogram();statistics();}
        fittedText(24,818,W-48,note,15,busy()?ORANGE:MUTED);
        FlushBatchDraw();
    }
    void scrub(int x){if(!hasData()||current().trajectory.empty())return;double f=std::clamp(static_cast<double>(x-timeline.l)/(timeline.r-timeline.l),0.,1.);time=f*replayEnd();playing=false;}
    void mouse(const ExMessage& m){
        if(m.message==WM_LBUTTONUP){scrubbing=false;return;}
        if(m.message==WM_MOUSEMOVE&&scrubbing){scrub(m.x);return;}
        if(m.message!=WM_LBUTTONDOWN)return;
        int x=m.x,y=m.y;
        if((activeView()==0||energyMode)&&timeline.contains(x,y)){scrubbing=true;scrub(x);return;}
        const int viewWidth=259/viewCount();
        for(int i=0;i<viewCount();++i)if(Rect{1007+i*viewWidth,94,1007+(i+1)*viewWidth-6,128}.contains(x,y)){setView(i);return;}
#if ROCK_STAGE >= 8
        if(Rect{757,58,911,89}.contains(x,y)){toggleMulti();return;}
        if(Rect{921,58,1071,89}.contains(x,y)){setView(3);return;}
#endif
        if(analysisMode)for(int i=0;i<5;++i)if(Rect{380+i*175,149,548+i*175,184}.contains(x,y)){factorIndex=i;return;}
        if(terrainMode)for(int i=0;i<4;++i)if(Rect{380+i*220,149,590+i*220,184}.contains(x,y)){terrainFactor=i;return;}
        if(Rect{370,94,476,128}.contains(x,y)){playing=!playing;return;}
        if(Rect{484,94,579,128}.contains(x,y)){time=0;playing=true;return;}
        if(Rect{587,94,696,128}.contains(x,y)){overlay=!overlay;return;}
        if(Rect{704,94,815,128}.contains(x,y)){heat=!heat;return;}
        if(Rect{825,94,906,128}.contains(x,y)){selectTrial(-1);return;}
        if(Rect{914,94,996,128}.contains(x,y)){selectTrial(1);return;}
        if(Rect{1070,20,1265,58}.contains(x,y)){openOutput();return;}
        if(busy())return;
        const int tabWidth=300/pageCount();
        for(int i=0;i<pageCount();++i)if(Rect{22+i*tabWidth,139,16+(i+1)*tabWidth,171}.contains(x,y)){page=i;chosenRow=0;return;}
        const auto list=rows();for(int i=0;i<static_cast<int>(list.size());++i){int yy=186+i*43;
            if(Rect{22,yy,320,yy+38}.contains(x,y)){chosenRow=i;if(x>=292)adjust(1);else if(x>=265)adjust(-1);return;}}
        if(Rect{24,555,116,593}.contains(x,y))begin(false);
        if(Rect{121,555,214,593}.contains(x,y))begin(true);
        if(Rect{219,555,316,593}.contains(x,y))exportData();
        if(Rect{24,605,166,646}.contains(x,y))beginStudy();
        if(Rect{174,605,316,646}.contains(x,y))beginTerrain();
    }
    void key(int k){
        if(k=='V'){setView((activeView()+1)%viewCount());return;}
#if ROCK_STAGE >= 8
        if(k=='M'){toggleMulti();return;}if(k=='F'){setView(3);return;}
#endif
        if(k=='O'){openOutput();return;}
        if(analysisMode&&k>='1'&&k<='5'){factorIndex=k-'1';return;}
        if(terrainMode&&k>='1'&&k<='4'){terrainFactor=k-'1';return;}
        if(k==VK_SPACE){playing=!playing;return;}
        if(k=='R'){time=0;playing=true;return;}
        if(k=='T'){overlay=!overlay;return;}if(k=='H'){heat=!heat;return;}
        if(k==VK_PRIOR){selectTrial(-1);return;}if(k==VK_NEXT){selectTrial(1);return;}
        if((k=='J'||k=='L')&&hasData()&&!current().trajectory.empty()){time=std::clamp(time+(k=='J'?-1.:1.),0.,replayEnd());playing=false;return;}
        if(busy())return;
        if(k==VK_RETURN)begin(false);else if(k=='B')begin(true);else if(k=='E')exportData();else if(k=='A')beginStudy();else if(k=='G')beginTerrain();
        else if(k==VK_TAB){page=(page+1)%pageCount();chosenRow=0;}
        else if(k==VK_UP)chosenRow=(chosenRow+static_cast<int>(rows().size())-1)%static_cast<int>(rows().size());
        else if(k==VK_DOWN)chosenRow=(chosenRow+1)%static_cast<int>(rows().size());
        else if(k==VK_LEFT)adjust(-1);else if(k==VK_RIGHT)adjust(1);
    }
};
} // namespace ui

inline int runViewer(){
    HWND window=initgraph(ui::W,ui::H);SetWindowTextW(window,(L"山区滚石模拟 — v1."+std::to_wstring(ROCK_STAGE)).c_str());BeginBatchDraw();
    ui::View view;view.begin(true);auto last=std::chrono::steady_clock::now();bool quit=false;
    while(!quit&&IsWindow(window)){
        ExMessage m;
        while(peekmessage(&m,EX_MOUSE|EX_KEY|EX_WINDOW,true)){
            if(m.message==WM_CLOSE){quit=true;break;}
            if(m.message==WM_KEYDOWN){if(m.vkcode==VK_ESCAPE){quit=true;break;}view.key(m.vkcode);}
            else view.mouse(m);
        }
        view.poll();auto now=std::chrono::steady_clock::now();double elapsed=std::chrono::duration<double>(now-last).count();last=now;
        if(view.playing&&view.hasData()&&!view.current().trajectory.empty()){
            view.time+=(std::min)(elapsed,.1)*view.speeds[view.speedIndex];
            if(view.time>=view.replayEnd()){view.time=view.replayEnd();view.playing=false;}}
        if(!quit&&IsWindow(window))view.render();Sleep(16);
    }
    // A running calculation owns no UI handles. Joining it before destruction avoids a detached thread.
    if(view.future.valid())view.future.wait();
    if(view.studyFuture.valid())view.studyFuture.wait();
    if(view.terrainFuture.valid())view.terrainFuture.wait();
    if(IsWindow(window)){EndBatchDraw();closegraph();}return 0;
}

// Optional development-only screenshot hook. No screenshot is written during normal use.
inline int runViewerSnapshotMode(const wchar_t* filename,int mode){
    HWND window=initgraph(ui::W,ui::H);SetWindowTextW(window,L"山区滚石模拟 — 预览");BeginBatchDraw();
    ui::View view;
#if ROCK_STAGE >= 5
    if(mode==3){view.pending.shape=Shape::Polygon;view.pending.sides=5;view.pending.shapePhase=PI/12;view.page=3;}
#endif
#if ROCK_STAGE >= 6
    if(mode==4){view.pending.barrier.kind=BarrierKind::Rigid;view.page=4;}
#endif
#if ROCK_STAGE >= 7
    if(mode==5){view.pending.barrier.kind=BarrierKind::Elastic;view.pending.barrier.maxDeflection=10;view.page=4;}
#endif
    view.batch=runBatch(view.config(),view.count,42);view.playing=false;
    if(!view.batch.runs.empty()&&!view.current().trajectory.empty())view.time=view.current().trajectory.back().t*.65;
#if ROCK_STAGE >= 8
    if(mode==6){view.multiReplay=true;view.time=view.replayEnd()*.55;}
    if(mode==7){view.setView(3);view.time=view.replayEnd()*.55;}
#endif
#if ROCK_STAGE >= 7
    if(mode==5&&!view.current().trajectory.empty()){
        const auto& samples=view.current().trajectory;
        view.time=std::max_element(samples.begin(),samples.end(),[](const Sample& a,const Sample& b){return std::abs(a.deflection)<std::abs(b.deflection);})->t;
    }
#endif
    if(mode==1){const std::string out=newResultPath("preview_parameters");
        runParameterStudy(view.config(),view.count,42,out);view.study=ui::readStudy(view.config(),view.count,42,out);view.setView(1);view.factorIndex=1;}
#if ROCK_STAGE >= 4
    if(mode==2){const std::string out=newResultPath("preview_terrain");
        runTerrainStudy(view.config(),view.count,42,out);view.terrainStudy=ui::readStudy(view.config(),view.count,42,out,true);view.setView(2);view.page=2;}
#endif
    view.render();saveimage(filename);EndBatchDraw();closegraph();return 0;
}
inline int runViewerSnapshot(const wchar_t* filename,bool analysisMode=false){return runViewerSnapshotMode(filename,analysisMode?1:0);}
} // namespace rock

int wmain(int argc,wchar_t** argv) {
    SetConsoleOutputCP(CP_UTF8);
    try {
        if(argc<=1)return rock::runViewer();
        std::vector<std::string> arguments;
        for(int i=0;i<argc;++i){
            int n=WideCharToMultiByte(CP_UTF8,0,argv[i],-1,nullptr,0,nullptr,nullptr);
            if(n<=0)throw std::runtime_error("Invalid Unicode argument");
            std::string s(static_cast<size_t>(n),'\0');
            WideCharToMultiByte(CP_UTF8,0,argv[i],-1,&s[0],n,nullptr,nullptr);
            s.pop_back();arguments.push_back(std::move(s));
        }
        std::vector<char*> pointers;
        for(auto& s:arguments)pointers.push_back(s.data());
        return rock::runCommandLine(static_cast<int>(pointers.size()),pointers.data());
    }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 2;}
}
