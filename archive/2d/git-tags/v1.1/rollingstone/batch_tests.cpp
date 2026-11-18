#include "cli.hpp"
#include <cassert>
#include <iostream>
#if ROCK_STAGE >= 2
#include <fstream>
#endif
int main(){
    rock::Config c;c.maxTime=30;c.sampleInterval=.04;
#if ROCK_STAGE >= 3
    c.randomSurface=true;
#endif
    auto a=rock::runBatch(c,30,42),b=rock::runBatch(c,30,42);
    assert(a.runs.size()==30 && a.heightCounts==b.heightCounts);
    assert(a.inBand+a.aboveBand+a.belowBand+a.notReached==30);
    int total=0;for(auto n:a.heightCounts)total+=n;assert(total==int(a.inBand));
    for(const auto& cell:a.routeCells)assert(cell.second>0 && cell.second<=30);
    for(size_t i=0;i<a.runs.size();++i){
        assert(a.runs[i].detectorHeight==b.runs[i].detectorHeight);
        assert(a.runs[i].config.initialSpeed==b.runs[i].config.initialSpeed);
        for(const auto& s:a.runs[i].trajectory)assert(std::isfinite(s.energy)&&std::isfinite(s.x)&&std::isfinite(s.z));
    }
    auto other=rock::runBatch(c,3,7);assert(other.runs[0].config.initialSpeed!=a.runs[0].config.initialSpeed);
    c.cv=0;auto fixed=rock::runBatch(c,3,1);assert(fixed.runs[0].detectorHeight==fixed.runs[2].detectorHeight);
    bool caught=false;try{rock::runBatch(c,0,0);}catch(const std::invalid_argument&){caught=true;}assert(caught);
#if ROCK_STAGE >= 2
    std::string dir=rock::newResultPath("csv_test");rock::exportBatch(a,dir);
    std::ifstream f(std::filesystem::u8path(dir)/"runs.csv");std::string row;size_t lines=0,cols=0;
    while(std::getline(f,row)){const size_t n=std::count(row.begin(),row.end(),',');if(lines==0)cols=n;else assert(n==cols);++lines;}
    assert(lines==31 && cols==26);
    caught=false;try{rock::exportBatch(a,dir);}catch(const std::runtime_error&){caught=true;}assert(caught);
    c.maxTime=.002;auto noArrival=rock::runBatch(c,2,2);assert(noArrival.notReached==2);
    const auto emptyDir=rock::newResultPath("csv_missing");rock::exportBatch(noArrival,emptyDir);
    std::ifstream mf(std::filesystem::u8path(emptyDir)/"runs.csv");lines=0;
    while(std::getline(mf,row)){assert(std::count(row.begin(),row.end(),',')==26);++lines;}assert(lines==3);
#endif
    std::cout<<"Batch checks passed: reproducibility, probability accounting, finite trajectories, CSV schema.\n";
}
