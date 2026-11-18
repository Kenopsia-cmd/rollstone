#include "simulation.hpp"
#include "batch.hpp"
#include "studies.hpp"
#include <cassert>
#include <iostream>
#include <iterator>
#include <set>

static std::vector<std::vector<std::string>> csv(const std::filesystem::path& path){
    std::ifstream f(path);assert(f);std::vector<std::vector<std::string>> rows;std::string line;
    while(std::getline(f,line)){std::vector<std::string> row;std::istringstream s(line);std::string cell;
        while(std::getline(s,cell,','))row.push_back(cell);
        if(!line.empty()&&line.back()==',')row.push_back("");
        rows.push_back(row);}
    return rows;
}
int main(int argc,char**argv){
    const std::filesystem::path out=argc>1?std::filesystem::u8path(argv[1]):std::filesystem::path("study_test_output");
    rock::Config c;c.maxTime=30.;c.sampleInterval=.05;
    rock::study::Moments known;for(double v:{1.,2.,3.})known.add(v);
    assert(std::abs(known.value()-2)<1e-12&&std::abs(known.sd()-1)<1e-12);
    rock::study::Moments zero;for(int i=0;i<10;++i)zero.add(0);
    auto w=rock::study::interval(zero,true);assert(w.first>=0&&w.second>0&&w.second<1);
    assert(rock::runParameterStudy(c,12,42,(out/"parameters").string())==15);
    assert(rock::runTerrainStudy(c,12,42,(out/"terrain").string())==12);
    auto p=csv(out/"parameters"/"parameter_study.csv");assert(p.size()==166);
    auto t=csv(out/"terrain"/"terrain_study.csv");assert(t.size()==133);
    for(const auto*table:{&p,&t})for(size_t i=1;i<table->size();++i){const auto&r=table->at(i);assert(r.size()==17);
        if(r[4].find("probability")!=std::string::npos){const double m=std::stod(r[6]);assert(m>=0&&m<=1);assert(std::stoul(r[5])==12);}
        if(r[2]=="baseline"&&r[4].find("given_reached")==std::string::npos)assert(std::abs(std::stod(r[13]))<1e-12);}
    auto base=csv(out/"terrain"/"baseline"/"terrain.csv"),wide=csv(out/"terrain"/"width_120"/"terrain.csv");
    assert(std::abs(std::stod(wide.back()[1])-std::stod(base.back()[1])*1.2)<1e-9);
    assert(std::abs(std::stod(wide.back()[2])-std::stod(base.back()[2]))<1e-9);
    auto grass=csv(out/"terrain"/"surface_grass"/"terrain.csv");assert(std::abs(std::stod(grass[1][3])-.1)<1e-12);
    auto bench=csv(out/"terrain"/"middle_bench"/"terrain.csv");const size_t a=std::max(size_t(1),c.terrain.points.size()/2-1);
    assert(std::abs(std::stod(bench[a+1][2])-std::stod(bench[a+2][2]))<1e-12);
    auto runs=csv(out/"parameters"/"baseline"/"run_metrics.csv");assert(runs.size()==13);
    auto surfaces=csv(out/"terrain"/"surface_layered"/"sampled_surfaces.csv");assert(surfaces.size()==1+12*(c.terrain.points.size()-1));
    bool threw=false;try{rock::runParameterStudy(c,1,42,(out/"invalid").string());}catch(const std::invalid_argument&){threw=true;}assert(threw);
    std::cout<<"Study tests passed: statistics, paired baselines, 15 parameter groups, 12 terrain cases, saved sampled inputs.\n";
}
