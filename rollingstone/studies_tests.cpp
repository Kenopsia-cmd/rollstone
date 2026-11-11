#include "simulation.hpp"
#include "batch.hpp"
#include "studies.hpp"
#include <chrono>
#include <iostream>
#include <set>

namespace fs=std::filesystem;
static void require(bool value,const std::string& message){
    if(!value)throw std::runtime_error("FAIL: "+message);
}
static std::string read(const fs::path& path){
    std::ifstream in(path,std::ios::binary);
    require(bool(in),"open "+path.u8string());
    return {std::istreambuf_iterator<char>(in),std::istreambuf_iterator<char>()};
}
static std::vector<std::vector<std::string>> csv(const fs::path& path){
    std::istringstream in(read(path));
    std::vector<std::vector<std::string>> table;std::string line;
    while(std::getline(in,line)){
        if(!line.empty()&&line.back()=='\r')line.pop_back();
        std::istringstream row(line);std::vector<std::string> cells;std::string cell;
        while(std::getline(row,cell,','))cells.push_back(cell);
        if(!line.empty()&&line.back()==',')cells.push_back("");
        table.push_back(cells);
    }
    return table;
}
static bool near(double a,double b){return std::abs(a-b)<1e-9;}
int main(int argc,char** argv){
    try{
        const auto tick=std::chrono::high_resolution_clock::now().time_since_epoch().count();
        const fs::path root=argc>1?fs::u8path(argv[1]):fs::path("results")/("study_tests_"+std::to_string(tick));
        require(!fs::exists(root),"test output must be new");
        rock::study::Moments moments;
        for(double v:{1.,2.,3.})moments.add(v);
        require(near(moments.value(),2.)&&near(moments.sd(),1.),"sample mean and unbiased standard deviation");
        rock::study::Moments zeros,ones;
        for(int i=0;i<12;++i){zeros.add(0);ones.add(1);}
        auto zeroCI=rock::study::interval(zeros,true),oneCI=rock::study::interval(ones,true);
        require(zeroCI.first>=0&&zeroCI.second>0&&zeroCI.second<1,"zero successes retain Wilson uncertainty");
        require(oneCI.first>0&&oneCI.first<1&&oneCI.second<=1,"all successes retain Wilson uncertainty");
        rock::BatchResult tinyBase,tinyChanged;
        for(int i=0;i<12;++i){rock::RunResult a,b;a.maxSpeed=10;b.maxSpeed=10.000001;
            a.status=b.status="stopped";tinyBase.runs.push_back(a);tinyChanged.runs.push_back(b);}
        std::ostringstream tinyReport;
        rock::study::describeEffects(tinyReport,"tiny_change",tinyChanged,tinyBase);
        require(tinyReport.str().find("步长收敛复核")!=std::string::npos,"tiny statistical differences flag numerical convergence checks");

        rock::Config config;config.maxTime=30;config.sampleInterval=.05;
        require(rock::runParameterStudy(config,12,42,(root/"parameters").u8string())==15,"five factors times three levels");
        const auto table=csv(root/"parameters"/"parameter_study.csv");
        require(table.size()==166,"15 groups times 11 metrics plus header");
        std::set<std::string> factors,groups;size_t baselineRows=0;
        for(size_t i=1;i<table.size();++i){
            const auto& r=table[i];require(r.size()==17,"consistent summary columns");
            factors.insert(r[0]);groups.insert(r[0]+":"+r[1]);
            require(r[3]=="12","requested count recorded");
            if(r[4].find("probability")!=std::string::npos){
                const double mean=std::stod(r[6]);
                require(r[5]=="12"&&mean>=0&&mean<=1,"unconditional probability denominator");
                require(std::stod(r[9])>=0&&std::stod(r[10])<=1,"probability interval in [0,1]");
            }
            if(r[2]=="baseline"){
                ++baselineRows;
                if(r[4].find("given_reached")==std::string::npos){
                    require(r[12]=="12"&&near(std::stod(r[13]),0),"baseline paired effect is zero");
                    require(near(std::stod(r[15]),0)&&near(std::stod(r[16]),0),"identical baseline comparison interval is zero");
                }
            }
        }
        require(factors.size()==5&&groups.size()==15&&baselineRows==55,"complete factor/level coverage");
        const auto runs=csv(root/"parameters"/"baseline"/"run_metrics.csv");
        require(runs.size()==13,"one saved row per independent run");
        bool variedSpeed=false,variedFriction=false;
        for(size_t i=1;i<runs.size();++i){
            require(runs[i].size()==18,"run schema");
            variedSpeed|=!near(std::stod(runs[i][1]),config.initialSpeed);
            variedFriction|=!near(std::stod(runs[i][4]),config.material.friction);
            require(std::stod(runs[i][1])>=0&&std::stod(runs[i][5])>=0&&std::stod(runs[i][5])<=1,"sampled inputs respect physical range");
        }
        require(variedSpeed&&variedFriction,"initial speed and surface parameters both sampled");
        const auto surfaces=csv(root/"parameters"/"baseline"/"sampled_surfaces.csv");
        require(surfaces.size()==1+12*(config.terrain.points.size()-1),"all sampled segment inputs saved");
        for(size_t i=1;i<surfaces.size();++i){
            const auto& row=surfaces[i];require(row.size()==6,"sampled surface schema");
            const auto& run=runs.at(std::stoul(row[0]));
            for(size_t k=2;k<6;++k)require(near(std::stod(row[k]),std::stod(run[k+1])),"saved segment inputs match sampled global material");
        }
        const auto report=read(root/"parameters"/"parameter_findings.txt");
        require(report.find("参数影响研究")!=std::string::npos&&report.find("方向")!=std::string::npos,"UTF-8 Chinese findings generated from metrics");

        const auto before=read(root/"parameters"/"parameter_study.csv");
        bool rejected=false;
        try{rock::runParameterStudy(config,12,42,(root/"parameters").u8string());}catch(const std::exception&){rejected=true;}
        require(rejected&&before==read(root/"parameters"/"parameter_study.csv"),"repeated output refused without overwriting");
        fs::create_directories(root/"occupied");
        {std::ofstream sentinel(root/"occupied"/"keep.txt");sentinel<<"original";}
        rejected=false;
        try{rock::runParameterStudy(config,12,42,(root/"occupied").u8string());}catch(const std::exception&){rejected=true;}
        require(rejected&&read(root/"occupied"/"keep.txt")=="original"&&!fs::exists(root/"occupied"/"parameter_study.csv"),"unrelated nonempty folder preserved");
        auto invalid=config;invalid.initialSpeed=200;
        rejected=false;
        try{rock::runParameterStudy(invalid,12,42,(root/"invalid_mean").u8string());}catch(const std::invalid_argument&){rejected=true;}
        require(rejected&&!fs::exists(root/"invalid_mean"),"1.2x mean validated before writing any files");
        rejected=false;
        try{rock::runParameterStudy(config,1,42,(root/"invalid_count").u8string());}catch(const std::invalid_argument&){rejected=true;}
        require(rejected&&!fs::exists(root/"invalid_count"),"study count checked before writing");

        auto expensive=config;expensive.dt=.00005;expensive.maxTime=600;expensive.sampleInterval=600;
        rejected=false;
        try{rock::runParameterStudy(expensive,20,42,(root/"invalid_budget").u8string());}catch(const std::invalid_argument&){rejected=true;}
        require(rejected&&!fs::exists(root/"invalid_budget"),"integration budget checked before writing");
        expensive.sampleInterval=.00005;
        rejected=false;
        try{rock::runParameterStudy(expensive,2,42,(root/"invalid_samples").u8string());}catch(const std::invalid_argument&){rejected=true;}
        require(rejected&&!fs::exists(root/"invalid_samples"),"sample budget checked before writing");

        config.maxTime=.002;
        require(rock::runParameterStudy(config,2,42,(root/"no_arrivals").u8string())==15,"all-not-reached study remains valid");
        const auto missing=csv(root/"no_arrivals"/"parameter_study.csv");
        for(size_t i=1;i<missing.size();++i){const auto& row=missing[i];
            if(row[4].find("probability")!=std::string::npos)require(row[5]=="2"&&near(std::stod(row[6]),0),"non-arrivals remain in denominator");
            if(row[4].find("given_reached")!=std::string::npos)require(row[5]=="0"&&row[6].empty()&&row[9].empty()&&row[10].empty(),"no-arrival conditional statistics stay missing");
        }
        std::cout<<"Task 3 tests passed: 15 groups, paired baseline, sampling, intervals, no-arrivals, Chinese report, output protection.\n";
        return 0;
    }catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}
}
