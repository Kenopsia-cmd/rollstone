// Acceptance tests for the three-task canvas. Include after core.inc + viewer.hpp.
// Deliberately calls the same key handlers used by the visible application.
#include <set>
#include <thread>

namespace rock::ui::acceptance {
inline void require(bool value,const std::string& message){
    if(!value)throw std::runtime_error("UI acceptance: "+message);
}
inline bool nearlyEqual(double a,double b,double tolerance=1e-8){
    return std::abs(a-b)<=tolerance*(1+std::abs(a)+std::abs(b));
}
using Table=std::vector<std::vector<std::string>>;
inline Table readCsv(const std::filesystem::path& file){
    std::ifstream input(file);require(bool(input),"missing CSV: "+pathToUtf8(file));
    Table table;std::string line;
    while(std::getline(input,line)){
        if(!line.empty()&&line.back()=='\r')line.pop_back();
        std::vector<std::string> row;std::istringstream stream(line);std::string cell;
        while(std::getline(stream,cell,','))row.push_back(cell);
        if(!line.empty()&&line.back()==',')row.emplace_back();
        table.push_back(std::move(row));
    }
    require(!table.empty(),"empty CSV: "+pathToUtf8(file));return table;
}
inline size_t column(const Table& table,const std::string& name){
    const auto it=std::find(table.front().begin(),table.front().end(),name);
    require(it!=table.front().end(),"missing CSV column: "+name);
    return static_cast<size_t>(it-table.front().begin());
}
inline std::string metadata(const Table& table,const std::string& name){
    for(size_t i=1;i<table.size();++i)if(table[i].size()==2&&table[i][0]==name)return table[i][1];
    throw std::runtime_error("UI acceptance: missing metadata: "+name);
}
inline void awaitWork(View& view,const char* operation){
    const auto deadline=std::chrono::steady_clock::now()+std::chrono::seconds(60);
    while(view.busy()){
        view.poll();
        require(std::chrono::steady_clock::now()<deadline,std::string(operation)+" timed out");
        std::this_thread::sleep_for(std::chrono::milliseconds(1));
    }
}
inline void checkBatch(const View& view,size_t expected){
    const auto& batch=view.batch;
    require(batch.runs.size()==expected,"displayed batch must contain all requested trials");
    require(&view.shown()==&batch.base,"canvas must use the displayed batch configuration");
    size_t reached=0,inBand=0,above=0,below=0;
    std::vector<int> bins(batch.heightCounts.size(),0);
    for(const auto& run:batch.runs){
        require(run.config.shape==Shape::Sphere&&run.config.barrier.kind==BarrierKind::None,
                "three-task batch must retain the original sphere without a barrier");
        require(run.trajectory.size()>2,"visible trajectories must come from actual sampled runs");
        if(run.reached){
            ++reached;
            if(run.inBand){++inBand;++bins[(std::min)(bins.size()-1,
                static_cast<size_t>((std::max)(0.,run.detectorHeight)/batch.base.binHeight))];}
            else if(run.detectorHeight>10)++above;
            else ++below;
        }
    }
    require(batch.reached==reached&&batch.inBand==inBand&&batch.aboveBand==above&&batch.belowBand==below,
            "canvas category counts must match the individual displayed trajectories");
    require(batch.notReached+reached==expected,"all trials must be accounted for in the denominator");
    require(batch.heightCounts==bins,"height chart must match this same batch");
    for(const auto& entry:batch.routeCells)require(entry.second>0&&entry.second<=static_cast<int>(expected),
            "route-cell probabilities count each trial at most once");
}
inline void checkExport(const View& view,const std::filesystem::path& folder){
    const auto& batch=view.batch;const auto n=batch.runs.size();
    const auto runs=readCsv(folder/"runs.csv"),bins=readCsv(folder/"height_probability.csv");
    const auto grid=readCsv(folder/"route_probability.csv"),config=readCsv(folder/"config.csv");
    require(runs.size()==n+1,"E must export every displayed trial");
    require(bins.size()==batch.heightCounts.size()+1,"E must export the displayed height bins");
    require(std::stoull(metadata(config,"runs"))==n,"export denominator must use displayed batch N");
    require(nearlyEqual(std::stod(metadata(config,"coefficient_of_variation")),batch.base.cv),"export CV mismatch");
    require(nearlyEqual(std::stod(metadata(config,"initial_speed_m_s")),batch.base.initialSpeed),"export must not mix in pending input");
    require(metadata(config,"shape")=="sphere"&&metadata(config,"barrier")=="none","export expanded task scope");
    const auto trial=column(runs,"trial_id"),speed=column(runs,"initial_speed_m_s"),reached=column(runs,"reached");
    for(size_t i=0;i<n;++i){
        require(std::stoull(runs[i+1][trial])==i+1,"exported trial order mismatch");
        require(nearlyEqual(std::stod(runs[i+1][speed]),batch.runs[i].config.initialSpeed),"CSV and displayed samples differ");
        require((std::stoi(runs[i+1][reached])!=0)==batch.runs[i].reached,"CSV and displayed crossing status differ");
    }
    const auto binCount=column(bins,"trial_count"),binN=column(bins,"total_trials"),binP=column(bins,"probability");
    for(size_t i=0;i<batch.heightCounts.size();++i){
        require(std::stoi(bins[i+1][binCount])==batch.heightCounts[i],"CSV and displayed histogram differ");
        require(std::stoull(bins[i+1][binN])==n,"histogram denominator changed");
        require(nearlyEqual(std::stod(bins[i+1][binP]),double(batch.heightCounts[i])/n),"histogram probability mismatch");
    }
    require(grid.size()==batch.routeCells.size()+1,"CSV and displayed route cells differ");
    const auto cellCount=column(grid,"visited_trials"),cellN=column(grid,"total_trials"),cellP=column(grid,"probability");
    size_t row=1;
    for(const auto& entry:batch.routeCells){
        require(std::stoi(grid[row][cellCount])==entry.second&&std::stoull(grid[row][cellN])==n,"route cell count mismatch");
        require(nearlyEqual(std::stod(grid[row][cellP]),double(entry.second)/n),"route cell probability mismatch");++row;
    }
    const auto tracks=readCsv(folder/"trajectories.csv");
    size_t samples=0;for(const auto& run:batch.runs)samples+=run.trajectory.size();
    require(tracks.size()==samples+1,"CSV must export precisely the displayed trajectories");
    const auto x=column(tracks,"x_m"),z=column(tracks,"z_m"),t=column(tracks,"t_s");
    row=1;
    for(const auto& run:batch.runs)for(const auto& point:run.trajectory){
        require(nearlyEqual(std::stod(tracks[row][x]),point.x)&&nearlyEqual(std::stod(tracks[row][z]),point.z)
            &&nearlyEqual(std::stod(tracks[row][t]),point.t),"playback sample differs from exported trajectory");++row;
    }
}
inline void checkStudy(const View& view,int expectedCount,double expectedSpeed){
    const auto& study=view.study;
    require(study.count==expectedCount&&study.seed==view.seed,"analysis must retain its sampled N and seed");
    require(nearlyEqual(study.base.initialSpeed,expectedSpeed)&&nearlyEqual(study.base.cv,.1),"analysis configuration mismatch");
    require(study.base.randomSurface,"analysis must randomize surface parameters");
    require(study.rows.size()==165,"analysis must load all 5 factors x 3 levels x 11 metrics");
    const auto data=readCsv(pathFromUtf8(study.directory)/"parameter_study.csv");
    require(data.size()==166,"analysis output has missing scenarios");
    const auto factor=column(data,"factor"),level=column(data,"level"),metric=column(data,"metric");
    const auto count=column(data,"runs"),valid=column(data,"valid_samples"),mean=column(data,"mean");
    const auto low=column(data,"ci95_low"),high=column(data,"ci95_high");
    std::map<std::string,std::set<std::string>> levels;
    for(size_t i=1;i<data.size();++i){
        const auto& row=data[i];
        require(std::stoi(row[count])==expectedCount,"scenario sample count mismatch");
        const auto* chart=study.find(row[factor],row[level],row[metric]);
        require(chart!=nullptr,"CSV row not loaded into chart data");
        require(chart->valid==std::stoi(row[valid]),"chart valid-sample count differs from CSV");
        const bool available=!row[mean].empty()&&!row[low].empty()&&!row[high].empty();
        require(chart->available==available,"chart must preserve unavailable conditional metrics");
        if(available)require(nearlyEqual(chart->mean,std::stod(row[mean]))&&nearlyEqual(chart->low,std::stod(row[low]))
                &&nearlyEqual(chart->high,std::stod(row[high])),"chart values/uncertainty differ from CSV");
        levels[row[factor]].insert(row[level]);
    }
    const std::set<std::string> expectedLevels={"0.8","1.0","1.2"};
    const std::set<std::string> expectedFactors={"initial_speed","roll_resistance","friction","restitution","drag"};
    std::set<std::string> actualFactors;
    for(const auto& entry:levels){actualFactors.insert(entry.first);require(entry.second==expectedLevels,"factor lacks one of its three levels");}
    require(actualFactors==expectedFactors,"analysis changed task scope or omitted a factor");
    require(study.find("absent","1.0","reach_probability")==nullptr,"missing chart data must not silently select another factor");
}

inline int run(){
    View view;
    require(view.count==200&&view.seed==42&&nearlyEqual(view.pending.cv,.1)&&view.pending.randomSurface,
            "startup must request 200 randomized trials with CV 0.1 and seed 42");
    require(view.overlay,"startup must show the trajectory ensemble");
    view.key('B');require(view.busy(),"B must start batch computation");awaitWork(view,"B batch");
    require(!view.dirty,"fresh batch cannot be marked stale");checkBatch(view,200);
    bool differentSpeed=false;
    for(const auto& run:view.batch.runs)differentSpeed|=!nearlyEqual(run.config.initialSpeed,view.batch.base.initialSpeed);
    require(differentSpeed,"B must sample a random ensemble, not repeat a deterministic run");
    view.scrub(view.timeline.r);
    require(nearlyEqual(view.sample().x,view.current().trajectory.back().x)&&!view.playing,"timeline must replay this batch's endpoint");
    view.selectTrial(-1);require(view.trial==199,"trial selection must cover all requested runs");
    view.selectTrial(1);require(view.trial==0,"trial selection must wrap back to the first run");

    const double oldSpeed=view.batch.base.initialSpeed;
    const auto oldHistogram=view.batch.heightCounts;
    view.page=0;view.chosenRow=2;view.adjust(1);
    require(view.dirty&&view.studyDirty,"changing a sampled parameter must invalidate both results");
    require(nearlyEqual(view.shown().initialSpeed,oldSpeed)&&!nearlyEqual(view.pending.initialSpeed,oldSpeed),
            "editing pending inputs must preserve the displayed batch snapshot");
    require(view.batch.heightCounts==oldHistogram,"editing inputs must not modify the displayed chart");
    view.key('E');require(!view.lastOutput.empty(),"E must create an output directory");
    checkExport(view,view.lastOutput);

    view.count=4;
    const double firstStudySpeed=view.pending.initialSpeed;
    view.key('A');require(view.busy(),"A must compute actual parameter study data");awaitWork(view,"A analysis");
    require(view.analysisMode&&!view.studyDirty,"A must select the completed analysis chart");
    require(view.dirty&&view.batch.runs.size()==200&&nearlyEqual(view.batch.base.initialSpeed,oldSpeed),
            "parameter analysis must not replace the separate displayed trajectory batch");
    checkStudy(view,4,firstStudySpeed);
    const std::string firstStudyFolder=view.study.directory;
    const double firstStudyMean=view.study.rows.front().mean;
    view.page=0;view.chosenRow=2;view.adjust(1);
    require(view.studyDirty&&nearlyEqual(view.study.base.initialSpeed,firstStudySpeed)
            &&nearlyEqual(view.study.rows.front().mean,firstStudyMean),"pending edits must retain and mark old analysis data");

    view.key('B');require(view.busy(),"B must recalculate pending input");awaitWork(view,"updated B batch");
    checkBatch(view,4);
    require(!view.dirty&&nearlyEqual(view.batch.base.initialSpeed,view.pending.initialSpeed),"B must publish the new snapshot");
    require(view.studyDirty&&view.study.directory==firstStudyFolder,"B must not mark the older study current");
    view.key('A');require(view.busy(),"A must recalculate stale study");awaitWork(view,"updated A analysis");
    require(!view.studyDirty&&view.study.directory!=firstStudyFolder,"A must publish a new analysis output");
    checkStudy(view,4,view.pending.initialSpeed);

    View controls;const auto terrain=controls.pending.terrain.points;
    for(int page=0;page<3;++page){
        require(controls.page==page,"settings tabs must have exactly three pages");
        const auto rows=controls.rows();require(!rows.empty(),"settings page is empty");
        for(size_t i=0;i<rows.size();++i){
            controls.chosenRow=static_cast<int>(i);controls.adjust(1);validate(controls.config());
            controls.adjust(-1);validate(controls.config());
            require(controls.pending.shape==Shape::Sphere&&controls.pending.barrier.kind==BarrierKind::None,
                    "task 4/5 or barrier settings leaked into the three-task interface");
        }
        controls.key(VK_TAB);
    }
    require(controls.page==0,"Tab must wrap after the third page");
    require(controls.pending.terrain.points.size()==terrain.size(),"controls must retain the original slope geometry");
    for(size_t i=0;i<terrain.size();++i)require(nearlyEqual(controls.pending.terrain.points[i].x,terrain[i].x)
            &&nearlyEqual(controls.pending.terrain.points[i].z,terrain[i].z),"task 4 terrain geometry edits are out of scope");
    std::cout<<"UI acceptance passed: default N=200/CV=0.1; B ensemble and counters; E matching CSV; "
        "A five factors x three levels with CSV-backed means/CI; independent stale snapshots; three-page scope.\n";
    return 0;
}
} // namespace rock::ui::acceptance

inline int runVisibleThreeTests(){return rock::ui::acceptance::run();}
