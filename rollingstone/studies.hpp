#pragma once
#include "model.hpp"
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <limits>
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
    if(!out)throw std::runtime_error("Cannot write study file: "+path.string());
    out.exceptions(std::ios::badbit|std::ios::failbit);
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
    out<<"runs="<<n<<"\nseed="<<seed<<"\ncv="<<c.cv<<"\nrandom_surface="<<c.randomSurface
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
    out<<"\nScenario: "<<name<<" (censored runs="<<censored(b)<<"/"<<b.runs.size()<<")\n";
    for(const auto&m:metrics()){
        if(m.conditional)continue;
        if(m.name!="height_0_10m_probability"&&m.name!="collector_energy_all_runs_J"&&m.name!="runout_x_m")continue;
        const auto s=aggregate(b,m),e=paired(b,baseline,m);const auto ci=effectInterval(b,baseline,m);
        out<<"  "<<m.name<<": mean=";number(out,s.value());out<<", difference from baseline=";number(out,e.value());
        out<<", paired 95% CI=[";number(out,ci.first);out<<", ";number(out,ci.second);out<<"]. ";
        if(m.name=="runout_x_m"&&(censored(b)||censored(baseline)))out<<"CENSORED trajectory extent: do not infer final runout ordering.";
        else if(!std::isfinite(ci.first))out<<"Too few samples to estimate uncertainty.";
        else if(ci.first>0)out<<"Higher in this experiment.";
        else if(ci.second<0)out<<"Lower in this experiment.";
        else out<<"Direction is not resolved by this interval.";
        out<<'\n';
    }
}
inline void commonNotes(std::ostream& out){
    out<<"All probabilities use all requested runs as denominator; missing the collector is not discarded.\n"
       <<"Collector measurements refer to the FIRST forward crossing by the stone center.\n"
       <<"all_runs collector speed/energy assign zero to non-arrivals; they are exposure metrics.\n"
       <<"given_reached metrics condition on arrival, have their own valid_samples count, and no paired effect estimate.\n"
       <<"Probability intervals: Wilson 95%; probability-difference intervals: conservative paired-discordance Wilson/Bonferroni 95%.\n"
       <<"Other mean/effect intervals: approximate Student-t 95%.\n"
       <<"Mean/effect intervals assume independent runs; they quantify Monte Carlo sampling error, NOT model error.\n"
       <<"Every scenario uses the same seed and number of runs for common-random-number comparisons.\n"
       <<"Normal inputs have target sigma=CV*mean; nonphysical draws are rejected. Bounded distributions have altered realized moments.\n"
       <<"Only the named scenario factor is changed; other input means remain at their base values.\n"
       <<"Finite max_time may censor runout. Inspect run_metrics.csv status and duration before interpreting results.\n"
       <<"Exploratory multiple comparisons are not corrected; do not treat every 95% interval as a confirmatory test.\n"
       <<"These are educational-model outputs, not validated design values for protective structures.\n";
}

} // namespace study
inline int runTerrainStudy(const Config& input,int count,uint32_t seed,const std::string& outputDir){
    if(count<2)throw std::invalid_argument("Studies require at least 2 runs per scenario");
    Config base=input;base.randomSurface=true;base.detectorX=-1.;validate(base);
    if(base.terrain.points.size()<4)throw std::invalid_argument("Terrain study needs at least 4 control points");
    const std::filesystem::path root=std::filesystem::u8path(outputDir);
    std::filesystem::create_directories(root);
    auto csv=study::open(root/"terrain_study.csv");study::summaryHeader(csv);
    auto report=study::open(root/"terrain_findings.txt");
    report<<"TERRAIN STUDY: horizontal extent, total relief, bench, local slope and surface materials\n";
    study::commonNotes(report);
    report<<"For ALL terrain scenarios, the collector follows the terrain toe (detector_x=-1).\n"
          <<"Width scaling changes travel distance AND slope angles. Height scaling changes potential energy AND slope angles.\n"
          <<"Bench and local-angle edits redistribute neighboring slopes while keeping endpoints fixed. These are coupled geometry effects.\n"
          <<"Material labels are illustrative parameter sets, NOT calibrated site measurements.\n";
    study::writeConfig(base,count,seed,root/"baseline");
    auto baseline=runBatch(base,count,seed);study::writeRunMetrics(baseline,root/"baseline");study::stripTrajectories(baseline);
    study::summaryRows(csv,"baseline","1.0","baseline",baseline,baseline);
    struct Scenario{std::string factor,level,name;Config config;};
    std::vector<Scenario> scenarios;
    for(int percent:{80,120}){
        const double f=percent/100.;Config width=base,height=base;
        const double x0=base.terrain.points.front().x,z0=base.terrain.toeZ();
        for(auto&p:width.terrain.points)p.x=x0+(p.x-x0)*f;
        if(width.barrier.kind!=BarrierKind::None)width.barrier.x=x0+(width.barrier.x-x0)*f;
        for(auto&p:height.terrain.points)p.z=z0+(p.z-z0)*f;
        scenarios.push_back({"horizontal_extent",percent==80?"0.8":"1.2","width_"+std::to_string(percent),width});
        scenarios.push_back({"total_relief",percent==80?"0.8":"1.2","height_"+std::to_string(percent),height});
    }
    Config bench=base;const size_t a=std::max(size_t(1),base.terrain.points.size()/2-1),b=a+1;
    const double benchZ=.5*(bench.terrain.points[a].z+bench.terrain.points[b].z);
    bench.terrain.points[a].z=benchZ;bench.terrain.points[b].z=benchZ;
    scenarios.push_back({"local_bench","flat","middle_bench",bench});
    const size_t k=std::max(size_t(1),std::min(base.terrain.points.size()-2,base.terrain.points.size()/3));
    for(double share:{.3,.7}){
        Config local=base;
        const double hi=local.terrain.points[k-1].z,lo=local.terrain.points[k+1].z;
        local.terrain.points[k].z=hi-share*(hi-lo);
        const std::string level=share<.5?"0.3":"0.7";
        scenarios.push_back({"local_drop_share",level,share<.5?"upper_drop_30":"upper_drop_70",local});
    }
    const Material rockMaterial{.015,.35,.55,1.},soilMaterial{.06,.55,.25,1.},grassMaterial{.10,.65,.20,1.};
    for(const auto&pair:std::vector<std::pair<std::string,Material>>{{"rock",rockMaterial},{"soil",soilMaterial},{"grass",grassMaterial}}){
        Config surface=base;surface.material=pair.second;surface.terrain.surfaces.assign(surface.terrain.points.size()-1,pair.second);
        scenarios.push_back({"surface",pair.first,"surface_"+pair.first,surface});
    }
    Config layered=base;layered.terrain.surfaces.clear();
    const size_t segments=layered.terrain.points.size()-1;
    for(size_t i=0;i<segments;++i)layered.terrain.surfaces.push_back(i*3<segments?rockMaterial:i*3<2*segments?soilMaterial:grassMaterial);
    scenarios.push_back({"surface","layered","surface_layered",layered});
    for(const auto&scenario:scenarios){
        validate(scenario.config);const auto dir=root/scenario.name;
        study::writeConfig(scenario.config,count,seed,dir);
        auto batch=runBatch(scenario.config,count,seed);study::writeRunMetrics(batch,dir);
        study::summaryRows(csv,scenario.factor,scenario.level,scenario.name,batch,baseline);
        study::describeEffects(report,scenario.name,batch,baseline);
    }
    report<<"\nActual edited coordinates and each segment's materials are saved in the scenario's terrain.csv.\n";
    return int(scenarios.size()+1);
}

inline int runParameterStudy(const Config& input,int count,uint32_t seed,const std::string& outputDir){
    if(count<2)throw std::invalid_argument("Studies require at least 2 runs per scenario");
    Config base=input;base.randomSurface=true;validate(base);
    const std::filesystem::path root=std::filesystem::u8path(outputDir);
    std::filesystem::create_directories(root);
    auto csv=study::open(root/"parameter_study.csv");study::summaryHeader(csv);
    auto report=study::open(root/"parameter_findings.txt");
    report<<"PARAMETER STUDY: one mean at a time, levels 0.8 / 1.0 / 1.2\n";study::commonNotes(report);
    study::writeConfig(base,count,seed,root/"baseline");
    auto baseline=runBatch(base,count,seed);study::writeRunMetrics(baseline,root/"baseline");study::stripTrajectories(baseline);
    struct Factor{std::string name;std::function<void(Config&,double)> scale;};
    const std::vector<Factor> factors={
        {"initial_speed",[](Config&c,double f){c.initialSpeed*=f;}},
        {"roll_resistance",[](Config&c,double f){c.material.rollResistance*=f;for(auto&m:c.terrain.surfaces)m.rollResistance*=f;}},
        {"friction",[](Config&c,double f){c.material.friction*=f;for(auto&m:c.terrain.surfaces)m.friction*=f;}},
        {"restitution",[](Config&c,double f){c.material.restitution*=f;for(auto&m:c.terrain.surfaces)m.restitution*=f;}},
        {"drag",[](Config&c,double f){c.material.drag*=f;for(auto&m:c.terrain.surfaces)m.drag*=f;}}
    };
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
    report<<"\nCSV rows at level 1.0 reference the SAME saved baseline run set; no redundant simulations were performed.\n";
    return scenarios;
}
}
