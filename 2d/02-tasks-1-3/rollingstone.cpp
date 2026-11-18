// 山区滚石模拟：任务1、2、3的单文件 C++ 版本。
// 功能：随机模拟/坡脚0-10m路线概率、CSV输出、CV=0.1参数影响分析。
// 本文件已内联全部项目自定义头文件，不需要复制任何 .hpp。
// 外部依赖：Windows、Visual Studio C++17、EasyX（本机已经安装）。
// Visual Studio：新建“空项目/控制台应用”，将本文件作为唯一有入口的.cpp加入，
// 选择 x64，支持 ISO C++17 / C++20；源文件使用 UTF-8 编码。
// 默认运行直接打开统计画布；命令行示例：
//   rollingstone.exe --batch --runs 500 --seed 42
//   rollingstone.exe --study --runs 200 --seed 42 --cv 0.1
// 输出位于运行目录的 results 子目录，不覆盖已有数据。
// 画布、轨迹回放、概率统计、CSV及参数分析使用同一批量计算核心。
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
#define ROCK_VERSION "v1.3.3"

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
      <<"No arguments: interactive canvas (tasks 1-3).\n"
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


// ==================== interactive viewer ====================

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
const COLORREF BG=RGB(13,23,36),PANEL=RGB(22,36,52),INK=RGB(226,236,244),MUTED=RGB(145,165,184);
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
inline void text(int x,int y,const std::wstring& s,int size=17,COLORREF color=INK){settextstyle(static_cast<int>(size*1.3),0,L"Microsoft YaHei");settextcolor(color);setbkmode(TRANSPARENT);outtextxy(x,y,s.c_str());}
inline void box(Rect r,COLORREF color){setfillcolor(color);solidrectangle(r.l,r.t,r.r,r.b);}
inline void lineAt(int x1,int y1,int x2,int y2,COLORREF c,int thickness=1){setlinecolor(c);setlinestyle(PS_SOLID,thickness);line(x1,y1,x2,y2);}
inline void button(Rect r,const std::wstring& label,bool selected=false,bool disabled=false){box(r,disabled?RGB(32,44,58):(selected?RGB(30,97,88):RGB(40,58,77)));text(r.l+9,r.t+7,label,16,disabled?MUTED:(selected?ACCENT:INK));}
inline std::wstring stopName(const RunResult& r){
    if(r.status=="time_limit")return L"达到时间上限";
    if(r.status=="domain_exit")return L"离开计算域";
    if(r.status=="stopped")return L"自然停稳";
    return wide(r.status);
}
struct Row {std::wstring name,value;};
struct StudyRow {
    std::string factor,level,metric;
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
    const StudyRow* find(const std::string& factor,const std::string& level,const std::string& metric)const{
        for(const auto& r:rows)if(r.factor==factor&&r.level==level&&r.metric==metric)return &r;
        return nullptr;
    }
};
inline StudyData readStudy(const Config& c,int count,uint32_t seed,const std::string& directory){
    StudyData data;data.base=c;data.count=count;data.seed=seed;data.directory=directory;
    std::ifstream file(pathFromUtf8(directory)/"parameter_study.csv");
    if(!file)throw std::runtime_error("Cannot read parameter_study.csv");
    std::string line;std::getline(file,line);
    while(std::getline(file,line)){
        std::vector<std::string> fields;std::istringstream input(line);std::string field;
        while(std::getline(input,field,','))fields.push_back(field);
        if(fields.size()<11)continue;
        StudyRow row;row.factor=fields[0];row.level=fields[1];row.metric=fields[4];
        row.valid=std::stoi(fields[5]);row.available=!fields[6].empty()&&!fields[9].empty()&&!fields[10].empty();
        if(row.available){row.mean=std::stod(fields[6]);row.low=std::stod(fields[9]);row.high=std::stod(fields[10]);}
        data.rows.push_back(row);
    }
    if(data.rows.size()!=165)throw std::runtime_error("Incomplete parameter study table");
    return data;
}
struct View {
    Config pending;
    BatchResult batch;
    int count=200,page=0,chosenRow=0,trial=0,speedIndex=0,factorIndex=0;
    uint32_t seed=42;
    bool playing=true,dirty=false,studyDirty=false,analysisMode=false,heat=false,overlay=true,scrubbing=false;
    double time=0;
    int materialPreset=1;
    std::wstring note=L"默认运行 200 次随机试验，CV=0.1；画布与 CSV 使用同一批轨迹。";
    std::filesystem::path lastOutput;
    std::future<BatchResult> future;
    StudyData study;
    std::future<StudyData> studyFuture;
    std::chrono::steady_clock::time_point jobStart;
    Rect plot{370,138,1266,518},timeline{376,564,1260,588};
    double worldLeft=0,worldRight=180,worldBottom=-4,worldTop=65;
    const double speeds[5]={.25,.5,1,2,4};
    View(){pending.randomSurface=true; speedIndex=2;}
    bool busy()const{return future.valid()||studyFuture.valid();}
    bool hasData()const{return !batch.runs.empty();}
    const RunResult& current()const{return batch.runs[static_cast<size_t>(trial)];}
    const Config& shown()const{return hasData()?batch.base:pending;}
    Config config()const{
        Config c=pending;
        c.terrain.surfaces.clear(); // The selected material applies uniformly to this scenario.
        return c;
    }
    std::vector<Row> rows()const{
        const auto& c=pending;
        if(page==0)return {{L"试验次数 N",num(count,0)},{L"随机种子",num(seed,0)},{L"初始速度 m/s",num(c.initialSpeed)},{L"等效半径 m",num(c.radius)},{L"变异系数 CV",num(c.cv)},{L"最长时长 s",num(c.maxTime,0)},{L"积分步长 ms",num(c.dt*1000,1)},{L"播放倍速",num(speeds[speedIndex])+L"×"}};
        if(page==1)return {{L"地表预设",materialPreset==1?L"岩石":materialPreset==2?L"土层":materialPreset==3?L"草地":L"自定义"},{L"滚动阻力",num(c.material.rollResistance,3)},{L"滑动摩擦",num(c.material.friction)},{L"恢复系数",num(c.material.restitution)},{L"空气阻力 Cd",num(c.material.drag)},{L"材质随机化",c.randomSurface?L"开":L"关"}};
        return {{L"高度分箱 m",num(c.binHeight)},{L"路径网格宽 m",num(c.gridWidth)}};
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
            case 0:c.binHeight=clamp(c.binHeight+dir*.25,.25,2);break;
            case 1:c.gridWidth=clamp(c.gridWidth+dir*.5,.5,5);break;}}
        dirty=true;studyDirty=true;note=L"设置已修改：当前图仍为旧结果；按 B 重算轨迹，按 A 重算参数对照。";
    }
    void begin(bool many){
        if(busy())return;
        Config c=config();if(!many)c.cv=0;
        try{validate(c);}catch(const std::exception& e){note=L"参数错误："+wide(e.what());return;}
        const int n=many?count:1;const uint32_t s=seed;
        jobStart=std::chrono::steady_clock::now();
        future=std::async(std::launch::async,[c,n,s]{return runBatch(c,n,s);});
        analysisMode=false;
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
        analysisMode=true;note=L"正在比较 5 个参数 × 3 个水平；完成后显示均值与 95% 置信区间。";
    }
    void poll(){
        if(future.valid()&&future.wait_for(std::chrono::milliseconds(0))==std::future_status::ready){
            try{batch=future.get();trial=0;time=0;playing=true;dirty=false;note=L"计算完成。PageUp / PageDown 选择试验；拖动时间轴回放。";}
            catch(const std::exception& e){note=L"计算失败："+wide(e.what());}
        }
        if(studyFuture.valid()&&studyFuture.wait_for(std::chrono::milliseconds(0))==std::future_status::ready){
            try{study=studyFuture.get();studyDirty=false;analysisMode=true;
                lastOutput=std::filesystem::absolute(pathFromUtf8(study.directory));
                note=L"参数对照完成；选择上方参数查看 0.8 / 1.0 / 1.2 倍结果，CSV 已自动保存。";}
            catch(const std::exception& e){note=L"参数分析失败："+wide(e.what());}
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
    Sample sample()const{
        if(!hasData()||current().trajectory.empty())return {};
        const auto& a=current().trajectory;
        auto it=std::lower_bound(a.begin(),a.end(),time,[](const Sample& p,double t){return p.t<t;});
        if(it==a.begin())return *it;if(it==a.end())return a.back();
        const auto& p=*(it-1);const auto& q=*it;double f=(time-p.t)/(q.t-p.t);
        Sample r;double* dst[]={&r.x,&r.z,&r.vx,&r.vz,&r.angle,&r.omega,&r.energy,&r.deflection};
        const double pv[]={p.x,p.z,p.vx,p.vz,p.angle,p.omega,p.energy,p.deflection};
        const double qv[]={q.x,q.z,q.vx,q.vz,q.angle,q.omega,q.energy,q.deflection};
        for(int i=0;i<8;++i)*dst[i]=pv[i]+f*(qv[i]-pv[i]);r.t=time;return r;
    }
    int px(double x)const{return plot.l+static_cast<int>((x-worldLeft)/(worldRight-worldLeft)*(plot.r-plot.l));}
    int py(double z)const{return plot.b-static_cast<int>((z-worldBottom)/(worldTop-worldBottom)*(plot.b-plot.t));}
    void fit(){
        const auto& c=shown();worldLeft=c.terrain.points.front().x-5;worldRight=c.terrain.toeX()+35;worldBottom=c.terrain.toeZ()-4;worldTop=c.terrain.points.front().z+10;
        if(hasData())for(const auto& r:batch.runs){worldRight=(std::max)(worldRight,r.runoutX+5);for(const auto& p:r.trajectory)worldTop=(std::max)(worldTop,p.z+3);}
        // One metre has the same pixel length on both axes, preserving rock geometry.
        const double scale=(std::min)((plot.r-plot.l)/(worldRight-worldLeft),(plot.b-plot.t)/(worldTop-worldBottom));
        worldRight=worldLeft+(plot.r-plot.l)/scale;worldTop=worldBottom+(plot.b-plot.t)/scale;
    }
    void path(const RunResult& r,COLORREF color,double until=-1,int thickness=1){
        setlinecolor(color);setlinestyle(PS_SOLID,thickness);
        const auto& a=r.trajectory;
        for(size_t i=1;i<a.size();++i){if(until>=0&&a[i-1].t>until)break;line(px(a[i-1].x),py(a[i-1].z),px(a[i].x),py(a[i].z));}
    }
    void renderScene(){
        const auto& c=shown();fit();box(plot,RGB(17,30,46));
        // Axes use world metres; all panels share the core's upward-positive elevation convention.
        for(double x=0;x<worldRight;x+=20){lineAt(px(x),plot.t,px(x),plot.b,RGB(33,49,65));text(px(x)-7,plot.b+6,num(x,0),13,MUTED);}
        for(double z=std::ceil(worldBottom/10)*10;z<worldTop;z+=10){lineAt(plot.l,py(z),plot.r,py(z),RGB(33,49,65));text(plot.l+4,py(z)+1,num(z,0)+L" m",12,MUTED);}
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
            for(const auto& run:batch.runs)path(run,RGB(160,85,78));
        }
        if(hasData()){
            const auto p=sample();path(current(),ORANGE,time,2);
            const auto& rc=current().config;
            const int rx=(std::max)(3,px(p.x+rc.radius)-px(p.x));
            const int ry=(std::max)(3,py(p.z)-py(p.z+rc.radius));
            setfillcolor(RGB(224,229,233));solidellipse(px(p.x)-rx,py(p.z)-ry,px(p.x)+rx,py(p.z)+ry);
            lineAt(px(p.x),py(p.z),px(p.x+rc.radius*std::cos(p.angle)),py(p.z+rc.radius*std::sin(p.angle)),BG,2);
        }
        text(plot.r-151,plot.t+9,L"水平距离 x / m",13,MUTED);
    }
    void sidebar(){
        box({14,94,328,806},PANEL);text(28,105,L"试验设置",19);
        const wchar_t* tabs[]={L"实验",L"地表",L"统计"};
        for(int i=0;i<3;++i)button({22+i*100,139,116+i*100,171},tabs[i],page==i,busy());
        const auto list=rows();
        for(int i=0;i<static_cast<int>(list.size());++i){int y=186+i*43;
            if(i==chosenRow)box({22,y-1,320,y+38},RGB(31,53,70));
            text(29,y+8,list[i].name,15);text(171,y+8,list[i].value,15,ACCENT);
            button({265,y+2,289,y+32},L"−",false,busy());button({292,y+2,316,y+32},L"+",false,busy());}
        button({24,555,116,593},L"单次 Enter",false,busy());button({121,555,214,593},L"批量 B",true,busy());button({219,555,316,593},L"导出 E",false,busy()||!hasData());
        button({24,605,316,646},L"参数对照 A  ·  5 参数 × 3 水平",false,busy());
        text(27,660,L"↑↓ 选参数   ←→ 调整   Tab 切页",14,MUTED);
        text(27,685,L"PageUp / PageDown 选择试验",14,MUTED);
        text(27,710,L"空格 暂停   J / L 前后定位",14,MUTED);
        text(27,738,L"批量按 CV 抽样；单次使用均值",14,MUTED);
        text(27,763,L"V 切换图表   O 输出目录",14,MUTED);
        text(27,784,L"球体 · 固定地形 · 前三项任务",12,ACCENT);
        if(page==2){text(29,300,L"蓝线：坡脚石心高度 0–10 m",14,SKYBLUE);
            text(29,330,L"分箱概率 = 本箱次数 / 全部 N 次",14,MUTED);
            text(29,360,L"每条轨迹仅统计首次向右通过",14,MUTED);
            text(29,390,L"热力图：轨迹访问网格的概率",14,MUTED);}
    }
    void histogram(){
        Rect area{370,617,826,798};box(area,PANEL);text(385,627,L"首次通过坡脚断面的路线概率",17);
        if(!hasData())return;
        const int n=static_cast<int>(batch.heightCounts.size());if(n==0)return;
        int left=400,top=666,bottom=763,right=810;double maxP=.05;
        for(int v:batch.heightCounts)maxP=(std::max)(maxP,static_cast<double>(v)/batch.runs.size());
        for(int i=0;i<n;++i){double p=static_cast<double>(batch.heightCounts[i])/batch.runs.size();
            int x0=left+static_cast<int>((std::min)(10.,i*batch.base.binHeight)/10*(right-left));
            int x1=left+static_cast<int>((std::min)(10.,(i+1)*batch.base.binHeight)/10*(right-left))-2;
            if(p>0)box({x0,bottom-static_cast<int>((bottom-top)*p/maxP),(std::max)(x0+1,x1),bottom},SKYBLUE);}
        lineAt(left,bottom,right,bottom,MUTED);text(379,652,num(maxP*100,1)+L"%",12,MUTED);
        for(int i=0;i<=5;++i)text(left+i*(right-left)/5-5,bottom+6,num(i*2,0),12,MUTED);
        text(605,782,L"石心高度 m；分母为全部 N 次",12,MUTED);
    }
    void statistics(){
        box({842,617,1266,798},PANEL);text(858,627,L"本次试验记录",17);
        if(!hasData())return;
        const auto& r=current();const size_t n=batch.runs.size();
        text(858,658,L"0–10 m 内通过："+num(static_cast<double>(batch.inBand),0)+L" / "+num(static_cast<double>(n),0)+L"  ("+num(100.*batch.inBand/n,1)+L"%)",15,ACCENT);
        text(858,683,L"未到达 "+num(static_cast<double>(batch.notReached),0)+L"    超高 "+num(static_cast<double>(batch.aboveBand),0)+L"    低于范围 "+num(static_cast<double>(batch.belowBand),0),14,MUTED);
        text(858,708,L"当前试验："+stopName(r)+L"，最大速度 "+num(r.maxSpeed,1)+L" m/s",14);
        text(858,733,r.reached?L"通过高度 "+num(r.detectorHeight)+L" m    速度 "+num(r.detectorSpeed,1)+L" m/s":L"当前试验未到达坡脚检测断面",14);
        text(858,758,r.reached?L"坡脚动能 "+num(r.detectorEnergy/1000000.,3)+L" MJ    首次通过 t="+num(r.detectorTime,1)+L" s":L"条件速度和动能不计入未到达试验",14);
        text(858,781,L"seed="+num(batch.seed,0)+L"   CV="+num(batch.base.cv)+L"   N="+num(static_cast<double>(n),0),12,MUTED);
    }
    void studyChart(Rect area,int metricIndex){
        const char* factors[]={"initial_speed","roll_resistance","friction","restitution","drag"};
        const char* metrics[]={"collector_speed_given_reached_m_s","collector_energy_given_reached_J","height_0_10m_probability"};
        const wchar_t* titles[]={L"到达样本：坡脚速度",L"到达样本：坡脚动能",L"全部试验：0–10 m 通过率"};
        const wchar_t* units[]={L"m/s",L"MJ",L"%"};
        const char* levels[]={"0.8","1.0","1.2"};
        const double multiplier=metricIndex==1?1e-6:metricIndex==2?100.:1.;
        box(area,PANEL);text(area.l+12,area.t+12,titles[metricIndex],14);
        double maximum=metricIndex==2?100.:.001;
        for(int j=0;j<3;++j){const auto* r=study.find(factors[factorIndex],levels[j],metrics[metricIndex]);
            if(metricIndex!=2&&r&&r->available)maximum=(std::max)(maximum,(std::max)(r->mean,r->high)*multiplier*1.2);}
        const int top=area.t+75,bottom=area.b-48;
        auto yy=[&](double value){return bottom-static_cast<int>((bottom-top)*std::clamp(value/maximum,0.,1.));};
        text(area.l+12,area.t+41,L"0 – "+num(maximum,metricIndex==2?0:2)+L" "+units[metricIndex],12,MUTED);
        lineAt(area.l+18,bottom,area.r-14,bottom,MUTED);
        for(int j=0;j<3;++j){int center=area.l+54+j*84;
            const auto* r=study.find(factors[factorIndex],levels[j],metrics[metricIndex]);
            if(r&&r->available){
                const double value=r->mean*multiplier;
                box({center-21,yy(value),center+21,bottom},j==1?ACCENT:j==0?SKYBLUE:ORANGE);
                int high=yy(r->high*multiplier),low=yy(r->low*multiplier);
                lineAt(center,high,center,low,INK,2);lineAt(center-7,high,center+7,high,INK,2);lineAt(center-7,low,center+7,low,INK,2);
                text(center-30,(std::max)(top-24,high-25),num(value,metricIndex==2?1:2),13);
            }else text(center-22,top+60,L"无到达",12,MUTED);
            text(center-23,bottom+10,j==0?L"0.8×":j==1?L"1.0×":L"1.2×",14,j==1?ACCENT:INK);
            if(r)text(center-20,bottom+31,L"n="+num(r->valid,0),10,MUTED);
        }
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
    void render(){
        setbkcolor(BG);cleardevice();
        text(24,20,L"山区滚石 · 随机试验",29);text(26,60,L"任务 1 路线概率   /   任务 2 CSV 导出   /   任务 3 参数分析",15,MUTED);
        if(busy()){
            double elapsed=std::chrono::duration<double>(std::chrono::steady_clock::now()-jobStart).count();
            text(650,28,(studyFuture.valid()?L"参数分析中… ":L"批量计算中… ")+num(elapsed,1)+L" s",19,ORANGE);
        }else {const bool stale=analysisMode?studyDirty:dirty;
            text(650,28,stale?L"当前显示旧结果，设置尚未重算":L"图表与 CSV 使用同一计算结果",16,stale?ORANGE:ACCENT);}
        button({1070,20,1265,58},L"输出目录 O",false,lastOutput.empty());
        sidebar();
        button({370,94,476,128},playing?L"暂停 Space":L"播放 Space",playing);
        button({484,94,579,128},L"重播 R");button({587,94,696,128},L"轨迹叠加 T",overlay);
        button({704,94,815,128},L"热力图 H",heat);
        button({825,94,906,128},L"上一条");button({914,94,996,128},L"下一条");
        button({1007,94,1122,128},L"轨迹视图",!analysisMode);button({1130,94,1266,128},L"参数对照 V",analysisMode);
        if(analysisMode)renderStudy();else{
        renderScene();
        if(hasData()){
            auto p=sample();const auto& a=current().trajectory;double end=a.empty()?0:a.back().t;
            text(376,539,L"试验 "+num(trial+1,0)+L"/"+num(static_cast<double>(batch.runs.size()),0)+L"   t="+num(time)+L"/"+num(end)+L" s   速度 "+num(std::hypot(p.vx,p.vz))+L" m/s   动能 "+num(p.energy/1000.,1)+L" kJ",15);
            box(timeline,RGB(36,55,74));int xx=timeline.l+static_cast<int>((timeline.r-timeline.l)*(end>0?time/end:0));box({timeline.l,timeline.t,xx,timeline.b},RGB(38,117,113));
            setfillcolor(ACCENT);solidcircle(xx,(timeline.t+timeline.b)/2,7);text(376,593,L"拖动时间轴逐帧查看；热力图为每条轨迹是否访问网格的概率。",12,MUTED);
        }
        histogram();statistics();}
        text(24,818,note.substr(0,90),15,busy()?ORANGE:MUTED);
        FlushBatchDraw();
    }
    void scrub(int x){if(!hasData()||current().trajectory.empty())return;double f=std::clamp(static_cast<double>(x-timeline.l)/(timeline.r-timeline.l),0.,1.);time=f*current().trajectory.back().t;playing=false;}
    void mouse(const ExMessage& m){
        if(m.message==WM_LBUTTONUP){scrubbing=false;return;}
        if(m.message==WM_MOUSEMOVE&&scrubbing){scrub(m.x);return;}
        if(m.message!=WM_LBUTTONDOWN)return;
        int x=m.x,y=m.y;
        if(!analysisMode&&timeline.contains(x,y)){scrubbing=true;scrub(x);return;}
        if(Rect{1007,94,1122,128}.contains(x,y)){analysisMode=false;return;}
        if(Rect{1130,94,1266,128}.contains(x,y)){analysisMode=true;return;}
        if(analysisMode)for(int i=0;i<5;++i)if(Rect{380+i*175,149,548+i*175,184}.contains(x,y)){factorIndex=i;return;}
        if(Rect{370,94,476,128}.contains(x,y)){playing=!playing;return;}
        if(Rect{484,94,579,128}.contains(x,y)){time=0;playing=true;return;}
        if(Rect{587,94,696,128}.contains(x,y)){overlay=!overlay;return;}
        if(Rect{704,94,815,128}.contains(x,y)){heat=!heat;return;}
        if(Rect{825,94,906,128}.contains(x,y)){selectTrial(-1);return;}
        if(Rect{914,94,996,128}.contains(x,y)){selectTrial(1);return;}
        if(Rect{1070,20,1265,58}.contains(x,y)){openOutput();return;}
        if(busy())return;
        for(int i=0;i<3;++i)if(Rect{22+i*100,139,116+i*100,171}.contains(x,y)){page=i;chosenRow=0;return;}
        const auto list=rows();for(int i=0;i<static_cast<int>(list.size());++i){int yy=186+i*43;
            if(Rect{22,yy,320,yy+38}.contains(x,y)){chosenRow=i;if(x>=292)adjust(1);else if(x>=265)adjust(-1);return;}}
        if(Rect{24,555,116,593}.contains(x,y))begin(false);
        if(Rect{121,555,214,593}.contains(x,y))begin(true);
        if(Rect{219,555,316,593}.contains(x,y))exportData();
        if(Rect{24,605,316,646}.contains(x,y))beginStudy();
    }
    void key(int k){
        if(k=='V'){analysisMode=!analysisMode;return;}
        if(k=='O'){openOutput();return;}
        if(analysisMode&&k>='1'&&k<='5'){factorIndex=k-'1';return;}
        if(k==VK_SPACE){playing=!playing;return;}
        if(k=='R'){time=0;playing=true;return;}
        if(k=='T'){overlay=!overlay;return;}if(k=='H'){heat=!heat;return;}
        if(k==VK_PRIOR){selectTrial(-1);return;}if(k==VK_NEXT){selectTrial(1);return;}
        if((k=='J'||k=='L')&&hasData()&&!current().trajectory.empty()){time=std::clamp(time+(k=='J'?-1.:1.),0.,current().trajectory.back().t);playing=false;return;}
        if(busy())return;
        if(k==VK_RETURN)begin(false);else if(k=='B')begin(true);else if(k=='E')exportData();else if(k=='A')beginStudy();
        else if(k==VK_TAB){page=(page+1)%3;chosenRow=0;}
        else if(k==VK_UP)chosenRow=(chosenRow+static_cast<int>(rows().size())-1)%static_cast<int>(rows().size());
        else if(k==VK_DOWN)chosenRow=(chosenRow+1)%static_cast<int>(rows().size());
        else if(k==VK_LEFT)adjust(-1);else if(k==VK_RIGHT)adjust(1);
    }
};
} // namespace ui

inline int runViewer(){
    HWND window=initgraph(ui::W,ui::H);SetWindowTextW(window,L"山区滚石模拟 — 前三项任务 · v1.3.3");BeginBatchDraw();
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
            if(view.time>=view.current().trajectory.back().t){view.time=view.current().trajectory.back().t;view.playing=false;}}
        if(!quit&&IsWindow(window))view.render();Sleep(16);
    }
    // A running calculation owns no UI handles. Joining it before destruction avoids a detached thread.
    if(view.future.valid())view.future.wait();
    if(view.studyFuture.valid())view.studyFuture.wait();
    if(IsWindow(window)){EndBatchDraw();closegraph();}return 0;
}

// Optional development-only screenshot hook. No screenshot is written during normal use.
inline int runViewerSnapshot(const wchar_t* filename,bool analysisMode=false){
    HWND window=initgraph(ui::W,ui::H);SetWindowTextW(window,L"山区滚石模拟 — 预览");BeginBatchDraw();
    ui::View view;view.batch=runBatch(view.config(),view.count,42);view.playing=false;
    if(!view.batch.runs.empty()&&!view.current().trajectory.empty())view.time=view.current().trajectory.back().t*.65;
    if(analysisMode){const std::string out=newResultPath("preview_parameters");
        runParameterStudy(view.config(),view.count,42,out);view.study=ui::readStudy(view.config(),view.count,42,out);view.analysisMode=true;view.factorIndex=1;}
    view.render();saveimage(filename);EndBatchDraw();closegraph();return 0;
}
} // namespace rock


// Unicode command-line entry; no arguments opens the interactive canvas.
int wmain(int argc, wchar_t** argv) {
    SetConsoleOutputCP(CP_UTF8);
    try {
        if (argc <= 1) return rock::runViewer();
        std::vector<std::string> arguments;
        for (int i=0; i<argc; ++i) {
            const int bytes=WideCharToMultiByte(CP_UTF8,0,argv[i],-1,nullptr,0,nullptr,nullptr);
            if(bytes<=0) throw std::runtime_error("Invalid Unicode argument");
            std::string a(static_cast<size_t>(bytes),'\0');
            WideCharToMultiByte(CP_UTF8,0,argv[i],-1,&a[0],bytes,nullptr,nullptr);
            a.pop_back(); arguments.push_back(std::move(a));
        }
        std::vector<char*> pointers;
        for(auto& a:arguments) pointers.push_back(a.data());
        return rock::runCommandLine(static_cast<int>(pointers.size()),pointers.data());
    } catch(const std::exception& e) {
        std::cerr << "Error: " << e.what() << '\n';
        return 2;
    }
}
