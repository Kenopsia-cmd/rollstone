// 山区滚石模拟：任务1、2、3的单文件 C++ 版本。
// 功能：随机模拟/坡脚0-10m路线概率、CSV输出、CV=0.1参数影响分析。
// 本文件已内联全部项目自定义头文件，不需要复制任何 .hpp。
// 外部依赖：Windows、Visual Studio C++17、EasyX（本机已经安装）。
// Visual Studio：新建“空项目/控制台应用”，将本文件作为唯一有入口的.cpp加入，
// 选择 x64，支持 ISO C++17 / C++20；源文件使用 UTF-8 编码。
// 默认运行显示中文菜单；命令行示例：
//   rollingstone.exe --batch --runs 500 --seed 42
//   rollingstone.exe --study --runs 200 --seed 42 --cv 0.1
// 输出位于运行目录的 results 子目录，不覆盖已有数据。
// 说明：原有动画保留作为演示；作业统计和研究使用新批量核心。
// 教学二维近似模型；不是经现场标定的工程防护设计工具。
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
// ==================== rollingstone.cpp ====================
// ==================== cli.hpp ====================
// ==================== version.hpp ====================
#define ROCK_STAGE 3
#define ROCK_VERSION "v1.3.2"

// ==================== simulation.hpp ====================
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
    bool barrierHit=false,barrierFailed=false;
    std::string stopReason;
};
struct Sample {
    double t=0,x=0,z=0,vx=0,vz=0,angle=0,omega=0,energy=0,deflection=0;
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
    if(c.shape!=Shape::Sphere||c.barrier.kind!=BarrierKind::None)throw std::invalid_argument("This version supports spheres without barriers");
}
inline std::vector<Vec2> bodyVertices(const Config& c,double angle){
    std::vector<Vec2> vertices;
    for(int i=0;i<48;++i){double a=angle+2*PI*i/48;vertices.push_back({c.radius*std::cos(a),c.radius*std::sin(a)});}
    return vertices;
}
inline double supportRadius(const Config& c,double,Vec2){return c.radius;}
inline double kineticEnergy(const State& s){return .5*s.mass*dot(s.v,s.v)+.5*s.inertia*s.omega*s.omega;}
inline State initialState(const Config& c){
    validate(c);State s;
    s.mass=c.density*4./3.*PI*c.radius*c.radius*c.radius;
    s.inertia=.4*s.mass*c.radius*c.radius;
    s.angle=c.shapePhase;
    double x=c.terrain.points.front().x,th=c.terrain.angle(x);
    s.p={x,c.terrain.height(x)+c.radius/std::cos(th)};
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
inline BatchResult runBatch(const Config& base,int count,uint32_t seed) {
    validate(base);
    if(count<1 || count>5000)throw std::invalid_argument("runs must be in [1,5000]");
    if(count*std::ceil(base.maxTime/base.dt)>2e8)
        throw std::invalid_argument("Requested integration work exceeds 200 million steps; reduce runs or time horizon");
    if(count*base.maxTime/base.sampleInterval>4e6)
        throw std::invalid_argument("Too many saved trajectory samples; reduce runs or increase --sample-interval");
    BatchResult b;b.base=base;b.seed=seed;
    b.heightCounts.assign(static_cast<size_t>(std::ceil(10/base.binHeight)),0);
    for(int i=0;i<count;++i){
        std::seed_seq trialSeed{seed,static_cast<uint32_t>(i)};
        std::mt19937 rng(trialSeed);
        auto r=simulate(randomize(base,rng));
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
            if(h<0 || h>10 || x<base.terrain.points.front().x || x>std::max(base.terrain.toeX()+80.,base.detectorX+30.))return;
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
    return s.str();
}
}

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
    runs<<"trial_id,seed,initial_speed_m_s,radius_m,roll_resistance,friction,restitution,drag,shape,sides,reached,in_0_10m,detector_height_m,detector_speed_m_s,detector_energy_J,detector_time_s,runout_x_m,max_speed_m_s,max_energy_J,max_bounce_height_m,duration_s,status,barrier_hit,barrier_failed,max_deflection_m,peak_force_N,dissipated_energy_J\n";
    auto trajectory=csvFile(out/"trajectories.csv");
    trajectory<<"trial_id,t_s,x_m,z_m,height_above_toe_m,vx_m_s,vz_m_s,angle_rad,omega_rad_s,kinetic_energy_J,barrier_deflection_m\n";
    auto surfaces=csvFile(out/"sampled_surfaces.csv");
    surfaces<<"trial_id,segment_index,roll_resistance,friction,restitution,drag\n";
    for(size_t i=0;i<n;++i){const auto& r=b.runs[i];const auto& m=r.config.material;const auto& s=r.finalState;
        runs<<i+1<<','<<b.seed<<','<<r.config.initialSpeed<<','<<r.config.radius<<','<<m.rollResistance<<','<<m.friction<<','<<m.restitution<<','<<m.drag<<','<<shapeName(r.config.shape)<<','<<r.config.sides<<','<<r.reached<<','<<r.inBand<<',';
        if(r.reached)runs<<r.detectorHeight<<','<<r.detectorSpeed<<','<<r.detectorEnergy<<','<<r.detectorTime;
        else runs<<",,,";
        runs<<','<<r.runoutX<<','<<r.maxSpeed<<','<<r.maxEnergy<<','<<r.maxBounceHeight<<','<<s.time<<','<<r.status<<','<<s.barrierHit<<','<<s.barrierFailed<<','<<s.maxBarrierDeflection<<','<<s.peakBarrierForce<<','<<s.barrierDissipation<<'\n';
        for(const auto& t:r.trajectory)trajectory<<i+1<<','<<t.t<<','<<t.x<<','<<t.z<<','<<t.z-c.terrain.toeZ()<<','<<t.vx<<','<<t.vz<<','<<t.angle<<','<<t.omega<<','<<t.energy<<','<<t.deflection<<'\n';
        for(size_t j=0;j+1<r.config.terrain.points.size();++j){
            const auto& sm=r.config.terrain.surfaces.empty()?m:r.config.terrain.surfaces[j];
            surfaces<<i+1<<','<<j<<','<<sm.rollResistance<<','<<sm.friction<<','<<sm.restitution<<','<<sm.drag<<'\n';
        }
    }
    auto terrain=csvFile(out/"terrain.csv");terrain<<"point_index,x_m,z_m\n";
    for(size_t i=0;i<c.terrain.points.size();++i)terrain<<i<<','<<c.terrain.points[i].x<<','<<c.terrain.points[i].z<<'\n';
    auto metadata=csvFile(out/"config.csv");metadata<<"parameter,value\n"
        <<"schema_version,1\nprogram_version,"<<ROCK_VERSION<<"\nseed,"<<b.seed<<"\nruns,"<<n<<"\ndt_s,"<<c.dt<<"\nmax_time_s,"<<c.maxTime
        <<"\nsample_interval_s,"<<c.sampleInterval<<"\ncoefficient_of_variation,"<<c.cv<<"\nrandom_surface,"<<c.randomSurface
        <<"\nradius_m,"<<c.radius<<"\ninitial_speed_m_s,"<<c.initialSpeed<<"\ndensity_kg_m3,"<<c.density
        <<"\nroll_resistance,"<<c.material.rollResistance<<"\nfriction,"<<c.material.friction<<"\nrestitution,"<<c.material.restitution<<"\ndrag,"<<c.material.drag
        <<"\ndetector_x_m,"<<(c.detectorX<0?c.terrain.toeX():c.detectorX)<<"\ndetector_reference_z_m,"<<c.terrain.toeZ()
        <<"\nbin_height_m,"<<c.binHeight<<"\ngrid_width_m,"<<c.gridWidth<<"\nshape,"<<shapeName(c.shape)<<"\nsides,"<<c.sides<<"\nshape_phase_rad,"<<c.shapePhase
        <<"\nbarrier,"<<barrierName(c.barrier.kind)<<"\nbarrier_x_m,"<<c.barrier.x<<"\nbarrier_height_m,"<<c.barrier.height
        <<"\nbarrier_restitution,"<<c.barrier.restitution<<"\nstiffness_N_m,"<<c.barrier.stiffness<<"\ndamping_N_s_m,"<<c.barrier.damping<<"\nmax_deflection_m,"<<c.barrier.maxDeflection
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
       <<"This is a 2-D educational model, not a calibrated engineering design or barrier certification tool.\n";
    for(auto* f:{&bins,&grid,&runs,&trajectory,&surfaces,&terrain,&metadata,&summary,&readme}){
        f->flush();if(!*f)throw std::runtime_error("CSV output failed (disk full or I/O error)");
    }
}
}

#endif
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
       <<"\ngrid_width_m="<<c.gridWidth<<"\nshape="<<(c.shape==Shape::Sphere?"sphere":"polygon")
       <<"\npolygon_sides="<<c.sides<<"\nshape_phase_rad="<<c.shapePhase
       <<"\nroll_resistance="<<c.material.rollResistance<<"\nfriction="<<c.material.friction
       <<"\nrestitution="<<c.material.restitution<<"\ndrag="<<c.material.drag
       <<"\nbarrier_kind="<<int(c.barrier.kind)<<"\nbarrier_x_m="<<c.barrier.x
       <<"\nbarrier_height_m="<<c.barrier.height<<"\nbarrier_restitution="<<c.barrier.restitution
       <<"\nbarrier_stiffness_N_m="<<c.barrier.stiffness<<"\nbarrier_damping_N_s_m="<<c.barrier.damping
       <<"\nbarrier_max_deflection_m="<<c.barrier.maxDeflection<<"\n";
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
inline void commonNotes(std::ostream& out){
    out<<"参数对应：initial_speed=初速度，roll_resistance=滚动阻力，friction=摩擦系数，restitution=恢复系数，drag=空气阻力系数。\n"
       <<"概率分母始终为全部试验数，未到达断面的试验没有被剔除。断面指标取石心首次向右通过时的数值。\n"
       <<"all_runs 的断面速度/能量对未到达试验计零，是综合暴露指标，不是条件到达速度。\n"
       <<"given_reached 指标仅统计已到达试验，valid_samples 为对应有效样本数；不提供不同条件样本之间的配对效应。\n"
       <<"概率区间使用 Wilson 95%；概率差使用配对不一致事件的 Wilson/Bonferroni 保守95%区间。\n"
       <<"其他均值及配对差使用近似 Student-t 95%区间；区间衡量蒙特卡洛抽样误差，不包含模型误差。\n"
       <<"各场景使用相同种子和试验数，以共同随机数进行配对比较。仅改变指定参数的均值，其余参数均值不变。\n"
       <<"正态抽样目标标准差=CV×均值；非物理值重采样。截断后实际均值和变异系数可能偏离目标。\n"
       <<"domain_exit/time_limit 表示轨迹被计算域或最长时长截断；runout_x_m 为已观测最远位置，不必然是真实停止位置。\n"
       <<"本分析为探索性多重比较，未进行整体多重检验校正，不能把每个95%区间都当作独立的确认性结论。\n"
       <<"这些结果仅来自教学近似模型，未经现场数据标定，不是工程防护设计依据。\n";
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
    if(count*std::ceil(base.maxTime/base.dt)>2e8)
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
}

#endif
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
      <<"No arguments: Chinese menu (tasks 1-3).\n"
      <<"--batch --runs 200 --seed 42 [--out results/name]\n"
      <<"--study : paired parameter experiments (5 parameters x 3 mean levels)\n"
      <<"--speed 1.5 --radius 2 --cv .1 --dt .002 --max-time 60 --sample-interval .02\n"
      <<"--roll .02 --friction .45 --restitution .45 --drag 1\n"
      <<"--material rock|soil|grass (illustrative surface presets)\n"
      <<"--detector-x 120 --bin-height .5 --grid-width 1\n"
      <<"CV=0 disables randomness. First forward detector crossing; 0..10m above toe.\n"
      <<"This delivery implements tasks 1-3 only. All functionality is C++.\n";
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
            const std::vector<std::string> allowed={"--runs","--seed","--out","--speed","--radius","--cv","--dt","--max-time","--sample-interval","--roll","--friction","--restitution","--drag","--detector-x","--bin-height","--grid-width","--material"};
            if(std::find(allowed.begin(),allowed.end(),key)==allowed.end())throw std::invalid_argument("Unknown option in tasks 1-3: "+key);
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
            else if(key=="--sides"){auto n=integer();if(n>32)throw std::invalid_argument("sides <= 32 required");c.sides=int(n);}
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
            c.terrain.points.clear();std::string line;bool first=true;
            while(std::getline(f,line)){
                if(line.empty())continue;
                if(first && (line[0]=='x'||line[0]=='X')){first=false;continue;}first=false;
                std::replace(line.begin(),line.end(),',',' ');std::istringstream row(line);Vec2 p;std::string extra;
                if(!(row>>p.x>>p.z)||(row>>extra))throw std::invalid_argument("Terrain rows require exactly x,z");
                c.terrain.points.push_back(p);
            }
        }
        if(!(c.terrainWidthScale>0 && c.terrainWidthScale<=5 && c.terrainHeightScale>0 && c.terrainHeightScale<=5))
            throw std::invalid_argument("Terrain scales must be in (0,5]");
        c.terrain.validate();
        const double toeZ=c.terrain.toeZ();
        for(auto& p:c.terrain.points){p.x*=c.terrainWidthScale;p.z=toeZ+(p.z-toeZ)*c.terrainHeightScale;}
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

using std::min; using std::max;
#include <easyx.h>
// Compatibility for the original drawing code; no private graphics.h is required.
static void legacySetColor(COLORREF color) { setlinecolor(color); settextcolor(color); }
static void legacyFillPolygon(int count,const int* xy) {
    std::vector<POINT> points(static_cast<size_t>(count));
    for(int i=0;i<count;++i) { points[i].x=xy[2*i]; points[i].y=xy[2*i+1]; }
    fillpolygon(points.data(),count);
}
#include <conio.h>
#include <math.h>
#include <cstring>
#include <cmath>
#include <tchar.h>
#include <windows.h>
#include <cstdlib>
#include <cstdio>

// ============================================================================
//  山区滚石模拟（米制版）
//  ---------------------------------------------------------------
//  物理量全部用国际单位制：位置 m、速度 m/s、时间 s、质量 kg、能量 J
//  只有在“画图”这一步才换算成像素，换算比例 PX_PER_M
//  物理模型：
//    · 触地滚动：a = (5/7) g sinθ − c_rr g cosθ − 空气阻力      （实心球纯滚动）
//    · 陡坡打滑：tanθ > μ 时摩擦不足以维持纯滚动，改用 a = g(sinθ − μcosθ)
//    · 空中飞行：竖直重力 + 二次空气阻力 ½ρCdAv²
//    · 碰撞弹跳：法向恢复系数 e，切向按库仑冲量损失
//    · 停止：速度足够小且坡度小于摩擦角（不规则块体被地面抓住）
// ============================================================================

// ------------------------------------------------------------------ 画面
#define WIDTH 1000          // 窗口宽度（像素）
#define HEIGHT 600          // 窗口高度（像素）
#define PX_PER_M 8.0        // 显示比例：1 米 = 8 像素（窗口即 125 m × 75 m）
#define ORIGIN_Y 490.0      // 高程 z=0 对应的屏幕 y 坐标
#define GROUND_Y 490        // 坡脚基准线（屏幕 y）

// ------------------------------------------------------------------ 物理常数
#define GRAVITY 9.81        // 重力加速度 (m/s²)
#define DT 0.002            // 物理步长 (s)
#define STEPS_PER_FRAME 8   // 每帧 8 步 × 2 ms = 16 ms ≈ 实时
#define ROCK_DENSITY 2600.0 // 岩石密度 (kg/m³)
#define AIR_DENSITY 1.225   // 空气密度 (kg/m³)
#define MIN_RADIUS 0.2      // 半径可调下限 (m)，只用来防止误操作，不是默认半径
#define MAX_RADIUS 3.0      // 半径可调上限 (m)
#define BOUNCE_MIN 0.15     // 法向速度小于此值不再弹跳，按接触处理 (m/s)
#define REST_SPEED 0.05     // 判定“停住”的速度阈值 (m/s)

// 可调参数（程序运行时通过菜单修改）
double g_rollResist = 0.02;      // 滚动阻力系数 c_rr（0.02~0.20）
double g_slideFriction = 0.45;   // 地面摩擦系数 μ（tan 24°，坡面摩擦角约 24°）
double g_restitution = 0.45;     // 法向恢复系数 e
double g_dragCoef = 1.0;         // 空气阻力系数 Cd（0 表示真空）
int g_playbackFactor = 1;        // 1x,2x,4x,8x

// ------------------------------------------------------------------ 坐标换算
double toPxX(double x) { return x * PX_PER_M; }
double toPxY(double z) { return ORIGIN_Y - z * PX_PER_M; }

// 地形点（x 水平距离 m，z 高程 m，向上为正）
struct Point {
    double x, z;
    Point(double x = 0, double z = 0) : x(x), z(z) {}
};

// 时间-速度数据结构（替代 std::pair<double,double>）
struct TimeSpeed {
    double t;
    double s;
    TimeSpeed(double tt = 0.0, double ss = 0.0) : t(tt), s(ss) {}
};

// 简单动态数组模板（替代 std::vector）
template<typename T>
class DynamicArray {
private:
    T* data;
    size_t cnt;
    size_t cap;

    void grow() {
        size_t newCap = cap == 0 ? 4 : cap * 2;
        T* newData = new T[newCap];
        if (data) {
            for (size_t i = 0; i < cnt; ++i) newData[i] = data[i];
            delete[] data;
        }
        data = newData;
        cap = newCap;
    }

public:
    DynamicArray() : data(nullptr), cnt(0), cap(0) {}
    ~DynamicArray() { delete[] data; }
    DynamicArray(const DynamicArray& other) : data(nullptr), cnt(0), cap(0) {
        if (other.cnt > 0) {
            data = new T[other.cnt];
            for (size_t i = 0; i < other.cnt; ++i) data[i] = other.data[i];
            cnt = other.cnt; cap = other.cnt;
        }
    }
    DynamicArray& operator=(const DynamicArray& other) {
        if (this == &other) return *this;
        delete[] data;
        data = nullptr; cnt = 0; cap = 0;
        if (other.cnt > 0) {
            data = new T[other.cnt];
            for (size_t i = 0; i < other.cnt; ++i) data[i] = other.data[i];
            cnt = other.cnt; cap = other.cnt;
        }
        return *this;
    }

    void push_back(const T& v) {
        if (cnt >= cap) grow();
        data[cnt++] = v;
    }

    void erase_front() {
        if (cnt == 0) return;
        for (size_t i = 1; i < cnt; ++i) data[i - 1] = data[i];
        --cnt;
    }

    size_t size() const { return cnt; }
    bool empty() const { return cnt == 0; }
    T& operator[](size_t i) { return data[i]; }
    const T& operator[](size_t i) const { return data[i]; }
    void clear() { cnt = 0; }
};

// ------------------------------------------------------------------ 山坡地形
// 13 个控制点，x 从 0 到 120 m，高程从 55.6 m 降到 0 m（坡脚）
class HillSlope {
private:
    DynamicArray<Point> points;  // 地形控制点（米）
    int startX, endX;            // 只用于快速判断，单位仍是米（取整）

public:
    HillSlope() {
        // 一条长 120 m、落差 55.6 m 的山区坡面：上部陡、中部有缓台、坡脚有一段平台
        points.push_back(Point(0.0, 55.6));    // 坡顶（危岩位置）
        points.push_back(Point(10.0, 51.9));
        points.push_back(Point(11.25, 50.6));
        points.push_back(Point(30.0, 35.6));   // 陡坡段
        points.push_back(Point(40.0, 34.4));   // 缓台（约 7°）
        points.push_back(Point(50.0, 23.1));   // 陡坡段
        points.push_back(Point(60.0, 20.6));
        points.push_back(Point(70.0, 15.6));
        points.push_back(Point(80.0, 10.0));
        points.push_back(Point(90.0, 8.75));   // 缓台
        points.push_back(Point(100.0, 2.88));
        points.push_back(Point(120.0, 0.0));   // 坡脚平台（道路位置）

        startX = (int)points[0].x;
        endX = (int)points[points.size() - 1].x;
    }

    // 指定 x 处的地面高程 z (m)
    double getHeight(double x) const {
        if (x <= points[0].x) return points[0].z;
        if (x >= points[points.size() - 1].x) return points[points.size() - 1].z;

        for (size_t i = 0; i < points.size() - 1; i++) {
            if (x >= points[i].x && x <= points[i + 1].x) {
                double t = (x - points[i].x) / (points[i + 1].x - points[i].x);
                return points[i].z + t * (points[i + 1].z - points[i].z);
            }
        }
        return points[0].z;
    }

    // 指定 x 处的坡角 θ（弧度，>0 表示顺坡向下）
    double getSlopeAngle(double x) const {
        if (x <= points[0].x) {
            double dz = points[1].z - points[0].z;
            double dx = points[1].x - points[0].x;
            return -atan2(dz, dx);
        }
        if (x >= points[points.size() - 1].x) {
            double dz = points[points.size() - 1].z - points[points.size() - 2].z;
            double dx = points[points.size() - 1].x - points[points.size() - 2].x;
            return -atan2(dz, dx);
        }

        for (size_t i = 0; i < points.size() - 1; i++) {
            if (x >= points[i].x && x <= points[i + 1].x) {
                double dz = points[i + 1].z - points[i].z;
                double dx = points[i + 1].x - points[i].x;
                return -atan2(dz, dx);
            }
        }
        return 0;
    }

    // 绘制山坡（唯一需要把米换成像素的地方）
    void draw() const {
        for (size_t i = 0; i < points.size() - 1; i++) {
            int x1 = (int)toPxX(points[i].x);
            int y1 = (int)toPxY(points[i].z);
            int x2 = (int)toPxX(points[i + 1].x);
            int y2 = (int)toPxY(points[i + 1].z);

            legacySetColor(RGB(34, 139, 34));
            setlinestyle(PS_SOLID, 3);
            line(x1, y1, x2, y2);

            setfillcolor(RGB(139, 90, 43));
            int pts[8] = { x1, y1, x2, y2, x2, HEIGHT, x1, HEIGHT };
            legacyFillPolygon(4, pts);
        }

        // 坡脚基准线（z = 0）
        legacySetColor(RGB(160, 82, 45));
        setlinestyle(PS_SOLID, 1);
        line(0, (int)toPxY(0), WIDTH, (int)toPxY(0));
    }

    int getStartX() const { return startX; }
    int getEndX() const { return endX; }

    int getPeakX() const {
        if (points.empty()) return startX;
        double bestZ = points[0].z;
        int bestX = (int)points[0].x;
        for (size_t i = 1; i < points.size(); ++i) {
            if (points[i].z > bestZ) {   // 高程最大者即坡顶
                bestZ = points[i].z;
                bestX = (int)points[i].x;
            }
        }
        return bestX;
    }

    int getValleyX() const {
        if (points.empty()) return startX;
        double lowZ = points[0].z;
        int lowX = (int)points[0].x;
        for (size_t i = 1; i < points.size(); ++i) {
            if (points[i].z < lowZ) {    // 高程最小者即坡脚
                lowZ = points[i].z;
                lowX = (int)points[i].x;
            }
        }
        return lowX;
    }
};

// ------------------------------------------------------------------ 滚石
class RollingStone {
private:
    double x, z;        // 球心位置（m）
    double vx, vz;      // 速度（m/s）
    double radius;      // 半径（m）
    bool isRolling;     // 是否还在运动
    double energy;      // 动能（J）
    double omega;       // 角速度（rad/s），纯滚动时 ω = v/r
    double mass;        // 质量（kg），由半径和密度算出
    DynamicArray<Point> trail;
    int trailMax;

    void updateMass() {
        double vol = 4.0 / 3.0 * 3.14159265358979 * radius * radius * radius;
        mass = ROCK_DENSITY * vol;
    }

public:
    // 半径没有默认值：必须由调用者（main）给定，对象只负责保存和使用
    RollingStone(double startX, double startZ, double r, double initVx)
        : x(startX), z(startZ), vx(initVx), vz(0), radius(r), isRolling(true),
          energy(0), omega(0), mass(0), trailMax(40) {
        if (radius < MIN_RADIUS) radius = MIN_RADIUS;
        if (radius > MAX_RADIUS) radius = MAX_RADIUS;
        updateMass();
        // 让球心正好贴在坡面上：沿法向抬高一个半径（真正的坡角在 adjustToGround 里再修正）
        z = startZ + radius;
    }

    // 大气阻力系数：a_drag = kd·v²，kd = ½ρCdA/m
    double dragK() const {
        double A = 3.14159265358979 * radius * radius;
        return 0.5 * AIR_DENSITY * g_dragCoef * A / mass;
    }

    // 物理更新：一步 dt
    void update(const HillSlope& hill, double dt) {
        if (!isRolling) return;

        double zg = hill.getHeight(x);
        double th = hill.getSlopeAngle(x);        // 坡角，>0 顺坡向下
        double ct = cos(th), st = sin(th);
        double nx = st, nz = ct;                  // 坡面外法向
        double dx = ct, dz = -st;                 // 顺坡切向（下坡方向）

        // 球心到坡面的法向距离：竖直高差投影到法向
        double dist = (z - zg) * ct;
        double vt = vx * dx + vz * dz;            // 切向速度（下坡为正）
        double vn = vx * nx + vz * nz;            // 法向速度（离开坡面为正）

        if (dist <= radius + 1e-3) {
            // ---------------- 接触坡面 ----------------
            z = zg + radius / ct;                 // 沿法向推回，恰好在表面

            if (vn < -BOUNCE_MIN) {
                // 弹跳：法向按恢复系数反弹，切向按库仑冲量损失
                double vnOut = -g_restitution * vn;
                double dvt = g_slideFriction * (1.0 + g_restitution) * fabs(vn);
                double sgn = (vt >= 0) ? 1.0 : -1.0;
                double vtOut = fabs(vt) - dvt;
                if (vtOut < 0) vtOut = 0;
                vt = sgn * vtOut;
                vn = vnOut;
                omega = vt / radius;
            }
            else {
                // 贴地运动：纯滚动 or 打滑
                vn = 0;
                double a;
                if (tan(th) > g_slideFriction) {
                    // 摩擦不足以维持纯滚动 → 打滑，转动不再同步
                    a = GRAVITY * (st - g_slideFriction * ct);
                    omega = 0;
                }
                else {
                    // 实心球纯滚动：转动惯量 I = 2/5 mR²，只有 5/7 的重力分量用于平动
                    a = (5.0 / 7.0) * GRAVITY * st - g_rollResist * GRAVITY * ct;
                    omega = vt / radius;
                }
                vt += a * dt;
                // 空气阻力（沿切向，与速度反向）
                double kd = dragK();
                vt -= kd * fabs(vt) * vt * dt;

                // 速度很小且坡度小于摩擦角 → 被地面抓住，停住
                if (fabs(vt) < REST_SPEED && tan(th) <= g_slideFriction) {
                    vt = 0;
                    isRolling = false;
                }
                else if (vt < 0) {
                    vt = 0;
                }
            }
            vx = vt * dx + vn * nx;
            vz = vt * dz + vn * nz;
        }
        else {
            // ---------------- 空中飞行 ----------------
            vz -= GRAVITY * dt;
            double sp = sqrt(vx * vx + vz * vz);
            double kd = dragK();
            if (sp > 1e-9) {
                vx -= kd * sp * vx * dt;
                vz -= kd * sp * vz * dt;
            }
        }

        // 位置更新
        x += vx * dt;
        z += vz * dt;

        // 边界：坡顶以上、坡脚以外
        if (x < hill.getStartX()) {
            x = hill.getStartX();
            vx = 0;
        }
        if (x >= hill.getEndX()) {
            x = hill.getEndX();
            vx = 0;
            isRolling = false;
        }

        // 防止穿过坡面：位置夹回表面，并去掉指向坡内的法向速度
        zg = hill.getHeight(x);
        th = hill.getSlopeAngle(x);
        {
            double ct2 = cos(th), st2 = sin(th);
            if ((z - zg) * ct2 < radius) {
                z = zg + radius / ct2;
                double vn2 = vx * st2 + vz * ct2;
                if (vn2 < 0) {
                    vx -= vn2 * st2;
                    vz -= vn2 * ct2;
                }
            }
        }

        // 动能：滚动时含转动动能 E = ½mv² + ½Iω² = 0.7mv²；其他情况只有平动
        double v2 = vx * vx + vz * vz;
        if (fabs(omega) > 1e-9) energy = 0.5 * mass * v2 + 0.2 * mass * radius * radius * omega * omega;
        else energy = 0.5 * mass * v2;

        // 记录轨迹（世界坐标，米）
        trail.push_back(Point(x, z));
        if ((int)trail.size() > trailMax) trail.erase_front();
    }

    // 绘制石头与轨迹（米 → 像素在这里完成）
    void draw() const {
        // 轨迹
        for (size_t i = 0; i < trail.size(); i++) {
            int alpha = 100 - (int)(i * 100 / trail.size());
            (void)alpha;
            legacySetColor(RGB(255, 100, 100));
            setfillstyle(BS_SOLID);
            int size = 3 - (int)(i * 2 / trail.size());
            if (size < 1) size = 1;
            solidcircle((int)toPxX(trail[i].x), (int)toPxY(trail[i].z), size);
        }

        // 抗锯齿圆掩码（按屏幕半径缓存）
        struct CircleMask {
            int r;
            int sz;
            unsigned char* alpha;
            CircleMask() :r(0), sz(0), alpha(nullptr) {}
            ~CircleMask() { delete[] alpha; }
            void generate(int radiusPx, int ss = 4) {
                if (radiusPx <= 0) return;
                if (r == radiusPx && alpha) return;
                delete[] alpha; alpha = nullptr;
                r = radiusPx;
                sz = 2 * r + 1;
                alpha = new unsigned char[sz * sz];
                int samples = ss * ss;
                double rr = (double)r;
                for (int iy = -r; iy <= r; ++iy) {
                    for (int ix = -r; ix <= r; ++ix) {
                        int inside = 0;
                        for (int sy = 0; sy < ss; ++sy) {
                            for (int sx = 0; sx < ss; ++sx) {
                                double fx = ix + (sx + 0.5) / ss;
                                double fy = iy + (sy + 0.5) / ss;
                                if (fx * fx + fy * fy <= rr * rr) inside++;
                            }
                        }
                        unsigned char a = (unsigned char)((inside * 255) / samples);
                        alpha[(iy + r) * sz + (ix + r)] = a;
                    }
                }
            }
        };

        int ir = (int)(radius * PX_PER_M);
        if (ir < 4) ir = 4;              // 显示时给一个最小半径，否则 1 m 的石头只有 8 px

        static CircleMask mask;
        mask.generate(ir, 3);

        COLORREF stoneColor = RGB(120, 120, 120);
        int mr = mask.r; int msz = mask.sz;

        // 阴影
        for (int yi = 0; yi < msz; ++yi) {
            for (int xi = 0; xi < msz; ++xi) {
                unsigned char a = mask.alpha[yi * msz + xi];
                if (a == 0) continue;
                int px = (int)toPxX(x) - mr + xi + 2;
                int py = (int)toPxY(z) - mr + yi + 3;
                if (px < 0 || px >= WIDTH || py < 0 || py >= HEIGHT) continue;
                COLORREF bg = getpixel(px, py);
                COLORREF sc = RGB(60, 60, 60);
                int alpha = a;
                int inv = 255 - alpha;
                int br = (GetRValue(sc) * alpha + GetRValue(bg) * inv) / 255;
                int bgc = (GetGValue(sc) * alpha + GetGValue(bg) * inv) / 255;
                int bb = (GetBValue(sc) * alpha + GetBValue(bg) * inv) / 255;
                putpixel(px, py, RGB(br, bgc, bb));
            }
        }

        // 主体
        for (int yi = 0; yi < msz; ++yi) {
            for (int xi = 0; xi < msz; ++xi) {
                unsigned char a = mask.alpha[yi * msz + xi];
                if (a == 0) continue;
                int px = (int)toPxX(x) - mr + xi;
                int py = (int)toPxY(z) - mr + yi;
                if (px < 0 || px >= WIDTH || py < 0 || py >= HEIGHT) continue;
                COLORREF bg = getpixel(px, py);
                int alpha = a;
                int inv = 255 - alpha;
                int rr = (GetRValue(stoneColor) * alpha + GetRValue(bg) * inv) / 255;
                int gg = (GetGValue(stoneColor) * alpha + GetGValue(bg) * inv) / 255;
                int bb = (GetBValue(stoneColor) * alpha + GetBValue(bg) * inv) / 255;
                putpixel(px, py, RGB(rr, gg, bb));
            }
        }

        // 高光
        int hr = max(1, (int)(ir / 3));
        static CircleMask hmask;
        hmask.generate(hr, 3);
        COLORREF highlight = RGB(180, 180, 180);
        int hmr = hmask.r; int hmsz = hmask.sz;
        int hx = (int)(toPxX(x) - ir / 3);
        int hy = (int)(toPxY(z) - ir / 3);
        for (int yi = 0; yi < hmsz; ++yi) {
            for (int xi = 0; xi < hmsz; ++xi) {
                unsigned char a = hmask.alpha[yi * hmsz + xi];
                if (a == 0) continue;
                int px = hx - hmr + xi;
                int py = hy - hmr + yi;
                if (px < 0 || px >= WIDTH || py < 0 || py >= HEIGHT) continue;
                COLORREF bg = getpixel(px, py);
                int alpha = a;
                int inv = 255 - alpha;
                int rr = (GetRValue(highlight) * alpha + GetRValue(bg) * inv) / 255;
                int gg = (GetGValue(highlight) * alpha + GetGValue(bg) * inv) / 255;
                int bb = (GetBValue(highlight) * alpha + GetBValue(bg) * inv) / 255;
                putpixel(px, py, RGB(rr, gg, bb));
            }
        }
    }

    bool isMoving() const { return isRolling; }
    double getX() const { return x; }
    double getZ() const { return z; }
    double getSpeed() const { return sqrt(vx * vx + vz * vz); }
    double getEnergy() const { return energy; }
    double getRadius() const { return radius; }
    double getMass() const { return mass; }
    void setRadius(double r) {
        if (r < MIN_RADIUS) r = MIN_RADIUS;
        if (r > MAX_RADIUS) r = MAX_RADIUS;
        radius = r;
        updateMass();
    }
    void adjustToGround(const HillSlope& hill) {
        double zg = hill.getHeight(x);
        z = zg + radius / cos(hill.getSlopeAngle(x));
        vz = 0;
    }
    // 顺坡释放：初速度沿坡面切向，而不是水平（更接近“石块从坡面脱离”）
    void startAlongSlope(const HillSlope& hill, double v0) {
        double th = hill.getSlopeAngle(x);
        vx = v0 * cos(th);
        vz = -v0 * sin(th);
        omega = v0 / radius;
    }
};

// ------------------------------------------------------------------ 主程序
int originalViewer() {
    initgraph(WIDTH, HEIGHT);
    setbkcolor(RGB(135, 206, 235));
    BeginBatchDraw();

    HillSlope hill;
    int initX = hill.getPeakX();
    double initZ = hill.getHeight((double)initX);

    // ---- 滚石参数：在 main 里给定，构造时交给滚石对象（对象内部只保存、不再有默认值）----
    double stoneRadius = 2.0;   // 滚石半径 (m)
    double stoneSpeed = 1.5;    // 顺坡释放的初速度 (m/s)

    RollingStone stone((double)initX, initZ, stoneRadius, 0.0);
    stone.adjustToGround(hill);
    stone.startAlongSlope(hill, stoneSpeed);
    int valleyX = hill.getValleyX();
    bool simEnded = false;

    bool paused = false;
    double stoppedDuration = 0.0;
    bool stoppedTriggered = false;
    bool inMenu = false;
    int menuSel = 0;                       // 0: c_rr 1: μ 2: e 3: Cd 4: 倍速
    const int menuRows = 4;
    const int playbackOptions[4] = { 1, 2, 4, 8 };
    int playbackIndex = 0;
    bool inPauseMenu = false;
    int pauseMenuSel = 0;

    char inputBuf[4][32];
    sprintf_s(inputBuf[0], sizeof(inputBuf[0]), "%.3f", g_rollResist);
    sprintf_s(inputBuf[1], sizeof(inputBuf[1]), "%.3f", g_slideFriction);
    sprintf_s(inputBuf[2], sizeof(inputBuf[2]), "%.3f", g_restitution);
    sprintf_s(inputBuf[3], sizeof(inputBuf[3]), "%.3f", g_dragCoef);

    DynamicArray< TimeSpeed > speedData;
    double simTime = 0.0;
    const double graphWindow = 10.0;

    settextstyle(16, 0, _T("Consolas"));
    legacySetColor(RGB(50, 50, 50));
    outtextxy(10, 20, _T("Rolling Stone Simulator - SI Units (m, s, kg, J)"));
    outtextxy(10, 45, _T("Controls: SPACE=Pause, M=Parameters, P=Pause menu, R=Reset, ESC=Exit"));

    while (true) {
        if (_kbhit()) {
            int key = _getch();
            if (key == 27) break;
            if (key == 32) {
                paused = !paused;
            }
            if (key == 'p' || key == 'P') {
                inPauseMenu = !inPauseMenu;
                pauseMenuSel = paused ? 1 : 0;
            }
            if (inPauseMenu) {
                if (key == 'w' || key == 'W' || key == 's' || key == 'S') {
                    pauseMenuSel = 1 - pauseMenuSel;
                }
                if (key == 13) {
                    if (pauseMenuSel == 0) paused = true; else paused = false;
                    inPauseMenu = false;
                }
            }
            if (key == 'r' || key == 'R') {
                double keepRadius = stone.getRadius();   // 保留用户用 +/- 改过的半径
                stone = RollingStone((double)initX, initZ, keepRadius, 0.0);
                stone.adjustToGround(hill);
                stone.startAlongSlope(hill, stoneSpeed);
                paused = false;
                simEnded = false;
                stoppedDuration = 0.0;
                stoppedTriggered = false;
                speedData.clear();
                simTime = 0.0;
                cleardevice();
            }
            if (key == '+' || key == '=') {
                stone.setRadius(stone.getRadius() + 0.1);
                stone.adjustToGround(hill);
            }
            if (key == '-' || key == '_') {
                stone.setRadius(stone.getRadius() - 0.1);
                stone.adjustToGround(hill);
            }
            if (key == 'm' || key == 'M') {
                inMenu = !inMenu;
                playbackIndex = 0;
                for (int i = 0; i < 4; i++) if (playbackOptions[i] == g_playbackFactor) playbackIndex = i;
            }
            if (inMenu) {
                if (key == 'w' || key == 'W') { menuSel = (menuSel + menuRows) % (menuRows + 1); }
                if (key == 's' || key == 'S') { menuSel = (menuSel + 1) % (menuRows + 1); }
                if (key == 'a' || key == 'A') {
                    if (menuSel == 0) g_rollResist = max(0.0, g_rollResist - 0.01);
                    if (menuSel == 1) g_slideFriction = max(0.0, g_slideFriction - 0.05);
                    if (menuSel == 2) g_restitution = max(0.0, g_restitution - 0.05);
                    if (menuSel == 3) g_dragCoef = max(0.0, g_dragCoef - 0.1);
                    if (menuSel == 4) { playbackIndex = max(0, playbackIndex - 1); g_playbackFactor = playbackOptions[playbackIndex]; }
                }
                if (key == 'd' || key == 'D') {
                    if (menuSel == 0) g_rollResist = min(0.30, g_rollResist + 0.01);
                    if (menuSel == 1) g_slideFriction = min(1.00, g_slideFriction + 0.05);
                    if (menuSel == 2) g_restitution = min(0.90, g_restitution + 0.05);
                    if (menuSel == 3) g_dragCoef = min(2.00, g_dragCoef + 0.1);
                    if (menuSel == 4) { playbackIndex = min(3, playbackIndex + 1); g_playbackFactor = playbackOptions[playbackIndex]; }
                }
                if (key == 13) { inMenu = false; }
            }
        }

        if (!paused) {
            // 每帧 8 个子步（2 ms × 8 = 16 ms），倍速就是“每帧多跑几步”
            int steps = STEPS_PER_FRAME * g_playbackFactor;
            for (int i = 0; i < steps; ++i) {
                stone.update(hill, DT);
                simTime += DT;
            }
            double curSpeed = stone.getSpeed();
            speedData.push_back(TimeSpeed(simTime, curSpeed));
            while (!speedData.empty() && speedData[0].t < simTime - graphWindow) {
                speedData.erase_front();
            }

            if (!simEnded) {
                if (stone.getX() >= hill.getEndX() - 0.5) {
                    MessageBox(NULL, _T("滚石到达坡脚，模拟结束"), _T("提示"), MB_OK | MB_ICONINFORMATION);
                    EndBatchDraw();
                    closegraph();
                    return 0;
                }
            }

            if (!stoppedTriggered && !simEnded) {
                if (curSpeed < REST_SPEED) stoppedDuration += STEPS_PER_FRAME * DT * g_playbackFactor;
                else stoppedDuration = 0.0;
                if (stoppedDuration >= 3.0) {
                    MessageBox(NULL, _T("滚石停止运动，模拟结束"), _T("提示"), MB_OK | MB_ICONINFORMATION);
                    EndBatchDraw();
                    closegraph();
                    return 0;
                }
            }

            cleardevice();
            setbkcolor(RGB(135, 206, 235));

            if (inMenu) {
                setfillcolor(RGB(240, 240, 240));
                solidrectangle(70, 45, 500, 440);
                legacySetColor(RGB(0, 0, 0));
                settextstyle(16, 0, _T("Consolas"));
                outtextxy(90, 55, _T("Parameters (W/S select, A/D change, Enter apply)"));

                TCHAR bufMenu[160];
                for (int i = 0; i < menuRows; i++) {
                    int y = 100 + i * 55;
                    if (i == menuSel) legacySetColor(RGB(255, 0, 0)); else legacySetColor(RGB(0, 0, 0));
                    if (i == 0) _stprintf_s(bufMenu, 160, _T("滚动阻力 c_rr (0.02-0.20):"));
                    else if (i == 1) _stprintf_s(bufMenu, 160, _T("摩擦系数 mu (0.3-0.8):"));
                    else if (i == 2) _stprintf_s(bufMenu, 160, _T("恢复系数 e (0.2-0.6):"));
                    else _stprintf_s(bufMenu, 160, _T("空气阻力 Cd (0=真空):"));
                    outtextxy(90, y, bufMenu);

                    legacySetColor(RGB(0, 0, 0));
                    TCHAR tbuf[64];
                    size_t outLen = 0;
                    mbstowcs_s(&outLen, tbuf, sizeof(tbuf) / sizeof(TCHAR), inputBuf[i], sizeof(inputBuf[i]));
                    outtextxy(400, y, tbuf);
                }
                int y = 100 + menuRows * 55;
                if (menuSel == menuRows) legacySetColor(RGB(255, 0, 0)); else legacySetColor(RGB(0, 0, 0));
                _stprintf_s(bufMenu, 160, _T("播放倍速: x1 x2 x4 x8 (当前 x%d)"), g_playbackFactor);
                outtextxy(90, y, bufMenu);

                setfillcolor(RGB(100, 200, 100));
                solidrectangle(250, 380, 380, 415);
                legacySetColor(RGB(0, 0, 0));
                outtextxy(275, 388, _T("启动模拟"));

                if (GetAsyncKeyState(VK_RETURN) & 0x8000) {
                    double v0 = atof(inputBuf[0]);
                    double v1 = atof(inputBuf[1]);
                    double v2 = atof(inputBuf[2]);
                    double v3 = atof(inputBuf[3]);
                    if (v0 < 0.0) v0 = 0.0; if (v0 > 0.30) v0 = 0.30;
                    if (v1 < 0.0) v1 = 0.0; if (v1 > 1.00) v1 = 1.00;
                    if (v2 < 0.0) v2 = 0.0; if (v2 > 0.90) v2 = 0.90;
                    if (v3 < 0.0) v3 = 0.0; if (v3 > 2.00) v3 = 2.00;
                    g_rollResist = v0; g_slideFriction = v1;
                    g_restitution = v2; g_dragCoef = v3;
                    inMenu = false; paused = false; simEnded = false;
                    stoppedDuration = 0.0; stoppedTriggered = false;
                    cleardevice();
                }
                if (menuSel >= 0 && menuSel <= 3) {
                    for (int c = 0; c < 256; c++) {
                        if (GetAsyncKeyState(c) & 0x8000) {
                            if (c == VK_BACK) {
                                size_t len = strlen(inputBuf[menuSel]);
                                if (len > 0) inputBuf[menuSel][len - 1] = 0;
                            }
                            else if ((c >= '0' && c <= '9') || c == '.' || c == '-') {
                                size_t len = strlen(inputBuf[menuSel]);
                                if (len + 1 < sizeof(inputBuf[menuSel])) {
                                    inputBuf[menuSel][len] = (char)c;
                                    inputBuf[menuSel][len + 1] = 0;
                                }
                            }
                        }
                    }
                }

                legacySetColor(RGB(80, 80, 80));
                outtextxy(90, 360, _T("按 R 重置模拟，按 M 退出菜单"));
            }

            if (inPauseMenu) {
                setfillcolor(RGB(240, 240, 240));
                solidrectangle(WIDTH / 2 - 120, HEIGHT / 2 - 40, WIDTH / 2 + 120, HEIGHT / 2 + 40);
                legacySetColor(RGB(0, 0, 0));
                settextstyle(16, 0, _T("Consolas"));
                outtextxy(WIDTH / 2 - 60, HEIGHT / 2 - 24, _T("Pause Menu (W/S select, Enter)"));
                if (pauseMenuSel == 0) legacySetColor(RGB(255, 0, 0)); else legacySetColor(RGB(0, 0, 0));
                outtextxy(WIDTH / 2 - 40, HEIGHT / 2 - 2, _T("Pause"));
                if (pauseMenuSel == 1) legacySetColor(RGB(255, 0, 0)); else legacySetColor(RGB(0, 0, 0));
                outtextxy(WIDTH / 2 - 40, HEIGHT / 2 + 18, _T("Continue"));
            }

            hill.draw();
            stone.draw();

            // -------- HUD（米制） --------
            TCHAR info[200];
            legacySetColor(RGB(255, 255, 255));
            settextstyle(14, 0, _T("Consolas"));

            _stprintf_s(info, 200, _T("Position: %.1f m"), stone.getX());
            outtextxy(10, HEIGHT - 100, info);

            _stprintf_s(info, 200, _T("Height: %.1f m"), stone.getZ());
            outtextxy(10, HEIGHT - 80, info);

            _stprintf_s(info, 200, _T("Speed: %.1f m/s"), stone.getSpeed());
            outtextxy(10, HEIGHT - 60, info);

            _stprintf_s(info, 200, _T("Energy: %.1f kJ"), stone.getEnergy() / 1000.0);
            outtextxy(10, HEIGHT - 40, info);

            _stprintf_s(info, 200, _T("Diameter: %.2f m   Mass: %.0f kg"), stone.getRadius() * 2.0, stone.getMass());
            outtextxy(10, HEIGHT - 20, info);

            _stprintf_s(info, 200, _T("Status: %s"), stone.isMoving() ? _T("ROLLING") : _T("STOPPED"));
            legacySetColor(stone.isMoving() ? RGB(255, 200, 0) : RGB(100, 255, 100));
            outtextxy(10, HEIGHT - 120, info);

            double slope = hill.getSlopeAngle(stone.getX()) * 180 / 3.14159;
            legacySetColor(RGB(200, 200, 200));
            _stprintf_s(info, 200, _T("Slope: %.1f deg"), slope);
            outtextxy(WIDTH - 130, HEIGHT - 60, info);

            TCHAR paramBuf[240];
            _stprintf_s(paramBuf, 240, _T("c_rr:%.2f mu:%.2f e:%.2f Cd:%.1f x%d"),
                        g_rollResist, g_slideFriction, g_restitution, g_dragCoef, g_playbackFactor);
            legacySetColor(RGB(180, 180, 180));
            settextstyle(12, 0, _T("Consolas"));
            outtextxy(WIDTH - 330, 10, paramBuf);

            // 速度条（0~20 m/s）
            int speedBar = (int)(stone.getSpeed() / 20.0 * 100);
            if (speedBar > 100) speedBar = 100;
            setfillcolor(RGB(255, 100, 100));
            solidrectangle(WIDTH - 120, HEIGHT - 40, WIDTH - 120 + speedBar, HEIGHT - 30);
            legacySetColor(RGB(150, 150, 150));
            rectangle(WIDTH - 120, HEIGHT - 40, WIDTH - 20, HEIGHT - 30);

            // -------- 速度-时间曲线 --------
            const int gx = WIDTH - 260;
            const int gy = 20;
            const int gw = 240;
            const int gh = 120;

            setfillcolor(RGB(135, 206, 235));
            solidrectangle(gx, gy, gx + gw, gy + gh);
            legacySetColor(RGB(200, 200, 200));
            rectangle(gx, gy, gx + gw, gy + gh);

            const int marginL = 36;
            const int marginR = 8;
            const int marginTop = 6;
            const int marginBottom = 20;
            int plotX = gx + marginL;
            int plotY = gy + marginTop;
            int plotW = gw - marginL - marginR;
            int plotH = gh - marginTop - marginBottom;

            double t0 = simTime - graphWindow;
            if (t0 < 0) t0 = 0;

            legacySetColor(RGB(180, 180, 180));
            line(plotX, plotY, plotX, plotY + plotH);
            line(plotX, plotY + plotH, plotX + plotW, plotY + plotH);

            double gridStart = floor(t0);
            double gridEnd = ceil(simTime);
            legacySetColor(RGB(220, 220, 220));
            setlinestyle(PS_DASH, 1);
            for (int g = (int)gridStart; g <= (int)gridEnd; ++g) {
                double gxTime = (double)g;
                double gxp = plotX + (gxTime - (simTime - graphWindow)) / graphWindow * plotW;
                if (gxp >= plotX && gxp <= plotX + plotW) {
                    line((int)gxp, plotY, (int)gxp, plotY + plotH);
                }
            }
            setlinestyle(PS_SOLID, 1);

            double yMax = 8.0;   // 纵轴最小量程 8 m/s
            for (size_t _i = 0; _i < speedData.size(); ++_i) {
                auto& p = speedData[_i];
                if (p.s > yMax) yMax = p.s;
            }
            yMax *= 1.1;

            settextstyle(12, 0, _T("Consolas"));
            legacySetColor(RGB(80, 80, 80));
            TCHAR buf[64];
            _stprintf_s(buf, 64, _T("%.1f"), yMax);
            outtextxy(plotX - 32, plotY - 8, buf);
            _stprintf_s(buf, 64, _T("%.1f"), 0.0);
            outtextxy(plotX - 28, plotY + plotH - 14, buf);

            _stprintf_s(buf, 64, _T("%.1fs"), t0);
            outtextxy(plotX + 2, plotY + plotH + 2, buf);
            _stprintf_s(buf, 64, _T("%.1fs"), simTime);
            outtextxy(plotX + plotW - 40, plotY + plotH + 2, buf);

            _stprintf_s(buf, 64, _T("Speed (m/s)"));
            outtextxy(plotX + 2, plotY - 18, buf);
            _stprintf_s(buf, 64, _T("Time (s)"));
            outtextxy(plotX + plotW / 2 - 22, plotY + plotH + 2, buf);

            if (speedData.size() >= 2) {
                legacySetColor(RGB(255, 0, 0));
                setlinestyle(PS_SOLID, 2);
                for (size_t i = 1; i < speedData.size(); ++i) {
                    double tx1 = speedData[i - 1].t;
                    double ty1 = speedData[i - 1].s;
                    double tx2 = speedData[i].t;
                    double ty2 = speedData[i].s;

                    double nx1 = plotX + (tx1 - (simTime - graphWindow)) / graphWindow * plotW;
                    double nx2 = plotX + (tx2 - (simTime - graphWindow)) / graphWindow * plotW;
                    double ny1 = plotY + plotH - (ty1 / yMax) * plotH;
                    double ny2 = plotY + plotH - (ty2 / yMax) * plotH;

                    if (nx1 < plotX) nx1 = plotX;
                    if (nx1 > plotX + plotW) nx1 = plotX + plotW;
                    if (nx2 < plotX) nx2 = plotX;
                    if (nx2 > plotX + plotW) nx2 = plotX + plotW;
                    if (ny1 < plotY) ny1 = plotY;
                    if (ny1 > plotY + plotH) ny1 = plotY + plotH;
                    if (ny2 < plotY) ny2 = plotY;
                    if (ny2 > plotY + plotH) ny2 = plotY + plotH;

                    line((int)nx1, (int)ny1, (int)nx2, (int)ny2);
                }
                setlinestyle(PS_SOLID, 1);
            }

            FlushBatchDraw();
            Sleep(16);
        }
        else {
            FlushBatchDraw();
            Sleep(50);
        }

        if (!stone.isMoving() && !paused) {
            legacySetColor(RGB(0, 255, 0));
            settextstyle(16, 0, _T("Consolas"));
            outtextxy(WIDTH / 2 - 80, HEIGHT - 150, _T("Stone Stopped! Press R to restart"));
        }
    }

    closegraph();
    return 0;
}

// ==================== entry.hpp ====================

#ifdef _WIN32
#include <windows.h>
namespace rock {
inline int windowsEntry(int argc,wchar_t** argv,int (*animation)()) {
    SetConsoleOutputCP(CP_UTF8);
    std::vector<std::string> arguments;
    for(int i=0;i<argc;++i){
        const int bytes=WideCharToMultiByte(CP_UTF8,0,argv[i],-1,nullptr,0,nullptr,nullptr);
        if(bytes<=0){std::cerr<<"Invalid Unicode argument\n";return 2;}
        std::string a(static_cast<size_t>(bytes),'\0');
        WideCharToMultiByte(CP_UTF8,0,argv[i],-1,&a[0],bytes,nullptr,nullptr);
        a.pop_back();arguments.push_back(std::move(a));
    }
    auto execute=[](std::vector<std::string>& args){std::vector<char*> pointers;for(auto& a:args)pointers.push_back(a.data());return runCommandLine(static_cast<int>(pointers.size()),pointers.data());};
    if(argc>1)return execute(arguments);
    for(;;){
        std::cout<<"\n山区滚石模拟 - 前三项任务（全部 C++）\n"
          <<"1. 随机模拟 + 坡脚路线概率 + CSV导出\n"
          <<"2. 参数影响分析（初速、滚动阻力、摩擦、恢复系数、空气阻力）\n"
          <<"3. 查看原有单次动画演示\n0. 退出\n请选择：";
        std::string choice;if(!std::getline(std::cin,choice)||choice=="0")return 0;
        if(choice=="3"){if(animation)return animation();std::cout<<"此构建不含图形演示。\n";continue;}
        if(choice!="1"&&choice!="2"){std::cout<<"请输入 0、1、2 或 3。\n";continue;}
        std::string n,seed,cv;
        std::cout<<"每组试验次数 [默认200]：";if(!std::getline(std::cin,n))return 0;
        std::cout<<"随机种子 [默认42]：";if(!std::getline(std::cin,seed))return 0;
        std::cout<<"变异系数 [默认0.1]：";if(!std::getline(std::cin,cv))return 0;
        std::vector<std::string> args={arguments[0],choice=="1"?"--batch":"--study","--runs",n.empty()?"200":n,"--seed",seed.empty()?"42":seed,"--cv",cv.empty()?"0.1":cv};
        execute(args);
    }
}
}
#endif

int wmain(int argc,wchar_t** argv) { return rock::windowsEntry(argc,argv,originalViewer); }
