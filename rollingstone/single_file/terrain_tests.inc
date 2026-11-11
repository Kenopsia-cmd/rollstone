// Include after the assembled core. No main; call rock::runTerrainStudyTests
// with a new/empty temporary directory. These checks use short trajectories.
#if ROCK_STAGE >= 4
namespace rock {
namespace terrain_tests {
inline void require(bool ok,const std::string& message){
    if(!ok)throw std::runtime_error("Terrain study regression: "+message);
}
inline std::string read(const std::filesystem::path& path){
    std::ifstream file(path,std::ios::binary);require(bool(file),"missing file "+pathToUtf8(path));
    std::ostringstream out;out<<file.rdbuf();return out.str();
}
inline std::vector<std::vector<std::string>> table(const std::filesystem::path& path){
    std::istringstream input(read(path));std::string line;
    std::vector<std::vector<std::string>> rows;
    while(std::getline(input,line)){
        if(!line.empty()&&line.back()=='\r')line.pop_back();
        std::vector<std::string> cells;size_t start=0;
        for(;;){const size_t comma=line.find(',',start);
            cells.push_back(line.substr(start,comma==std::string::npos?comma:comma-start));
            if(comma==std::string::npos)break;
            start=comma+1;
        }
        rows.push_back(std::move(cells));
    }
    return rows;
}
inline double setting(const std::filesystem::path& path,const std::string& name){
    std::istringstream input(read(path));std::string line;
    while(std::getline(input,line))if(line.compare(0,name.size()+1,name+"=")==0)
        return std::stod(line.substr(name.size()+1));
    throw std::runtime_error("Missing config field: "+name);
}
inline bool near(double a,double b){return std::abs(a-b)<1e-9*std::max(1.,std::max(std::abs(a),std::abs(b)));}
}
inline void runTerrainStudyTests(const std::filesystem::path& testRoot){
    using namespace terrain_tests;
    require(!std::filesystem::exists(testRoot)||
        (std::filesystem::is_directory(testRoot)&&std::filesystem::is_empty(testRoot)),"test root must be new or empty");
    std::filesystem::create_directories(testRoot);
    Config input;input.dt=.01;input.maxTime=.03;input.sampleInterval=.02;
    input.detectorX=150.;input.terrainWidthScale=1.5;input.terrainHeightScale=1.25;
    for(auto& p:input.terrain.points){p.x+=10.;p.z+=7.;}
    const auto before=input.terrain.points;
    const auto root=testRoot/pathFromUtf8("地形 回归结果");
    require(runTerrainStudy(input,2,42,pathToUtf8(root))==12,"must generate 12 scenarios");
    require(input.detectorX==150.&&!input.randomSurface&&input.terrain.surfaces.empty(),"input config changed");
    for(size_t i=0;i<before.size();++i)
        require(input.terrain.points[i].x==before[i].x&&input.terrain.points[i].z==before[i].z,"input geometry changed");
    const auto rows=table(root/"terrain_study.csv");
    require(rows.size()==1+12*study::metrics().size(),"unexpected summary row count");
    for(const auto& row:rows)require(row.size()==17,"summary must have 17 columns");
    for(size_t i=1;i<rows.size();++i){const auto& row=rows[i];
        require(row[3]=="2","incorrect run count");
        if(row[4]=="reach_probability"||row[4]=="height_0_10m_probability"){
            const double p=std::stod(row[6]);require(p>=0&&p<=1&&row[5]=="2","probability denominator changed");
        }
        if(row[2]=="baseline"&&row[4].find("given_reached")==std::string::npos){
            require(near(std::stod(row[13]),0.),"baseline effect must be zero");
            if(row[4].find("probability")!=std::string::npos)
                require(near(std::stod(row[15]),0.)&&near(std::stod(row[16]),0.),"baseline probability interval must be zero");
        }
    }
    const auto base=table(root/"baseline"/"terrain.csv"),width=table(root/"width_120"/"terrain.csv"),height=table(root/"height_80"/"terrain.csv");
    const double startX=before.front().x,toeZ=before.back().z;
    require(base.size()==before.size()+1,"wrong terrain point count");
    for(size_t i=0;i<before.size();++i){
        require(near(std::stod(width[i+1][1]),startX+(before[i].x-startX)*1.2),"width rescaled wrong origin or twice");
        require(near(std::stod(width[i+1][2]),before[i].z),"width changed height");
        require(near(std::stod(height[i+1][2]),toeZ+(before[i].z-toeZ)*.8),"height changed toe reference");
    }
    require(near(setting(root/"width_120"/"config.txt","terrain_width_scale"),1.8),"width metadata not cumulative");
    require(near(setting(root/"height_80"/"config.txt","terrain_height_scale"),1.),"height metadata not cumulative");
    require(near(setting(root/"width_120"/"config.txt","collector_x_m"),std::stod(width.back()[1])),"collector must follow edited toe");
    const auto bench=table(root/"middle_bench"/"terrain.csv");
    const size_t a=std::max(size_t(1),before.size()/2-1);
    require(near(std::stod(bench[a+1][2]),std::stod(bench[a+2][2])),"middle bench is not flat");
    for(int percent:{30,70}){
        const auto drop=table(root/("upper_drop_"+std::to_string(percent))/"terrain.csv");
        const size_t k=std::max(size_t(1),std::min(before.size()-2,before.size()/3));
        require(near(std::stod(drop[k+1][2]),before[k-1].z-percent/100.*(before[k-1].z-before[k+1].z)),"local drop share incorrect");
    }
    const auto grass=table(root/"surface_grass"/"terrain.csv"),soil=table(root/"surface_soil"/"terrain.csv");
    require(near(std::stod(grass[1][3]),.12)&&near(std::stod(grass[1][4]),.75)&&near(std::stod(grass[1][5]),.20),"grass preset differs from CLI");
    require(near(std::stod(soil[1][3]),.08)&&near(std::stod(soil[1][4]),.60)&&near(std::stod(soil[1][5]),.25),"soil preset differs from CLI");
    require(read(root/"baseline"/"run_metrics.csv")==read(root/"surface_rock"/"run_metrics.csv"),"uniform material changed sampling structure");
    const auto sampled=table(root/"surface_grass"/"sampled_surfaces.csv");
    const size_t segments=before.size()-1;
    require(sampled.size()==1+2*segments,"missing sampled segment materials");
    for(size_t trial=0;trial<2;++trial)for(size_t j=1;j<segments;++j)for(size_t col=2;col<6;++col)
        require(sampled[1+trial*segments+j][col]==sampled[1+trial*segments][col],"uniform material must be shared across all segments per run");
    const auto layered=table(root/"surface_layered"/"terrain.csv");
    for(size_t i=0;i<segments;++i)
        require(near(std::stod(layered[i+1][3]),i*3<segments?.02:i*3<2*segments?.08:.12),"layered material thirds incorrect");
    const auto original=read(root/"terrain_study.csv");bool rejected=false;
    try{runTerrainStudy(input,2,42,pathToUtf8(root));}catch(const std::invalid_argument&){rejected=true;}
    require(rejected&&read(root/"terrain_study.csv")==original,"existing report overwritten");
    auto rejectBeforeFiles=[&](Config bad,int runs,const std::string& name){
        const auto path=testRoot/name;bool failed=false;
        try{runTerrainStudy(bad,runs,42,pathToUtf8(path));}catch(const std::invalid_argument&){failed=true;}
        require(failed&&!std::filesystem::exists(path),name+" should fail before creating output");
    };
    rejectBeforeFiles(input,1,"bad_count");
    Config few=input;few.terrain.points.resize(3);rejectBeforeFiles(few,2,"few_vertices");
    Config invalid=input;invalid.terrain.points[1].x=invalid.terrain.points[0].x;rejectBeforeFiles(invalid,2,"bad_geometry");
    Config overflow=input;overflow.terrain.points={{0,3},{5e307,2},{1e308,1},{1.6e308,0}};rejectBeforeFiles(overflow,2,"invalid_candidate");
    Config tooManySteps=input;tooManySteps.maxTime=600.;tooManySteps.dt=.00005;tooManySteps.sampleInterval=1.;
    rejectBeforeFiles(tooManySteps,5000,"integration_budget");
    Config tooManySamples=input;tooManySamples.maxTime=600.;tooManySamples.dt=.01;tooManySamples.sampleInterval=.01;
    rejectBeforeFiles(tooManySamples,5000,"sample_budget");
    require(read(root/"terrain_findings.txt").find("地形影响研究")!=std::string::npos,"Chinese report missing");
}
}
#endif
