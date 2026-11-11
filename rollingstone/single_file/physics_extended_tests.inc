#include <iostream>
#include <iomanip>
namespace rock::extended_tests {
inline void require(bool v,const char* message){if(!v)throw std::runtime_error(std::string("Extended physics: ")+message);}
inline bool close(double a,double b,double tolerance=1e-8){return std::abs(a-b)<=tolerance*(1+std::abs(a)+std::abs(b));}
inline Config flat(){Config c;c.terrain.points={{0,0},{120,0}};return c;}
inline double total(const State& s,const Config& c){return kineticEnergy(s)+s.mass*G*s.p.z+detail::panelStored(s,c)+s.barrierDissipation;}
inline void ticks(State& s,const Config& c,int n){for(int i=0;i<n&&s.moving;++i)step(s,c);}
inline int run(){
    Config c;auto original=initialState(c),wrapped=original;
    for(int i=0;i<500;++i){stepSphereLegacy(original,c);step(wrapped,c);}
    require(original.p.x==wrapped.p.x&&original.p.z==wrapped.p.z&&original.v.x==wrapped.v.x&&original.omega==wrapped.omega,
            "sphere/no-barrier results must preserve the first-three core exactly");
#if ROCK_STAGE >= 5
    c=flat();c.shape=Shape::Polygon;c.radius=1;c.sides=4;c.material={0,0,1,0};
    auto s=initialState(c);require(close(s.inertia,s.mass/3),"square lamina inertia");
    s.v={10,-5};s.omega=0;
    std::vector<detail::Contact> contacts={{{-.5,-.5},{0,1},0,0,0,0,c.material},{{.5,-.5},{0,1},0,0,0,0,c.material}};
    const double before=kineticEnergy(s);detail::groundImpulse(s,contacts);
    require(kineticEnergy(s)<=before*(1+1e-10),"multiple contacts must not add kinetic energy");
    require(s.v.z>4.99&&std::abs(s.omega)<1e-5,"symmetric elastic face impact must retain symmetry");
    c.terrain.points={{0,0},{9,0},{10,1},{11,0},{120,0}};c.shapePhase=PI/4;s=initialState(c);s.p={10,1.2};s.v={0,0};s.omega=0;
    double correction=0;auto ridge=detail::groundContacts(s,c,correction);
    require(!ridge.empty()&&correction>.50,"terrain corner against polygon edge must be detected");
    step(s,c);require(s.p.z-std::sqrt(.5)>=1-1e-7,"ridge must not penetrate the polygon lower edge");
    // Ordinary outer dt changes share internal substeps. This checks consistent
    // stepping, not convergence of a long, repeatedly impacting trajectory.
    c=Config();c.shape=Shape::Polygon;c.sides=7;c.initialSpeed=8;c.maxTime=8;
    c.terrain.points={{0,12},{30,0}};c.detectorX=15;c.material={.01,.4,.2,0};
    c.dt=.002;const auto coarse=simulate(c);c.dt=.001;const auto fine=simulate(c);c.dt=.0005;const auto finer=simulate(c);
    require(coarse.reached&&fine.reached&&finer.reached,"controlled polygon reaches collector");
    require(std::abs(coarse.detectorSpeed-finer.detectorSpeed)<.05&&std::abs(fine.detectorSpeed-finer.detectorSpeed)<.05,
            "polygon collector speed must agree when internal substeps are identical");
    require(std::abs(coarse.detectorHeight-finer.detectorHeight)<.02,"outer-step collector interpolation consistency");
    const double e0=coarse.trajectory.front().energy+coarse.finalState.mass*G*coarse.trajectory.front().z;
    for(const auto& p:coarse.trajectory)require(p.energy+coarse.finalState.mass*G*p.z<=e0*(1+1e-8),"passive polygon travel energy bound");
    auto drop=[](double h){Config dropConfig=flat();dropConfig.shape=Shape::Polygon;dropConfig.sides=4;dropConfig.radius=1;
        dropConfig.shapePhase=PI/4;dropConfig.material={0,0,.6,0};dropConfig.dt=h;
        State falling=initialState(dropConfig);falling.p={10,5};falling.v={0,-1};falling.omega=0;
        for(int i=0;i<static_cast<int>(1.2/h+.5);++i){step(falling,dropConfig);
            for(const auto& vertex:bodyVertices(dropConfig,falling.angle))require(falling.p.z+vertex.z>=-1e-7,"controlled drop must not penetrate ground");}
        return falling;};
    const auto drop1=drop(.0005),drop2=drop(.00025),drop3=drop(.000125);
    require(std::abs(drop1.v.z-drop3.v.z)<.02&&std::abs(drop2.v.z-drop3.v.z)<.01,"single polygon impact velocity refinement");
    require(std::abs(drop1.p.z-drop3.p.z)<.01&&std::abs(drop3.omega)<1e-5,"single polygon impact pose refinement");
    std::cout<<"stage 5: polygon corner/edge contact, inertia, passive impacts, outer-step consistency and single-impact refinement passed\n";
#else
    c.shape=Shape::Polygon;bool rejected=false;try{validate(c);}catch(const std::invalid_argument&){rejected=true;}require(rejected,"polygon stage gate");
#endif
#if ROCK_STAGE >= 6
    c=flat();c.radius=.05;c.dt=.01;c.material={0,0,0,0};c.barrier.kind=BarrierKind::Rigid;c.barrier.x=10;c.barrier.height=1;
    auto fast=initialState(c);fast.p={9.9,.95};fast.v={100,40};fast.omega=0;step(fast,c);
    require(fast.barrierHit&&fast.contactImpulse>0,"high-speed corner crossing must not miss finite barrier");
    require(fast.peakBarrierForce==0&&fast.barrierForce==0,"rigid impact must report impulse rather than fake j/dt peak force");
    c=flat();c.shape=Shape::Polygon;c.radius=1;c.sides=4;c.barrier.kind=BarrierKind::Rigid;c.barrier.x=10;c.barrier.height=1;
    auto tip=initialState(c);tip.p={9.5,1.5};tip.angle=0;
    const auto cap=detail::barrierContacts(tip,c,-1);require(!cap.empty(),"polygon edge at panel cap must contact");
    bool upward=false;for(const auto& p:cap)upward|=p.normal.z>.5;
    require(upward,"polygon tip impact must use edge normal with vertical impulse");
    c=flat();c.radius=1;c.material.drag=0;c.barrier.kind=BarrierKind::Rigid;c.barrier.x=10;c.barrier.height=4;
    auto reflected=initialState(c);reflected.p={8.9,2};reflected.v={20,0};reflected.omega=0;const double energy=total(reflected,c);
    ticks(reflected,c,20);require(reflected.barrierHit&&reflected.v.x<0,"rigid face impact reflects incoming rock");
    require(reflected.p.x+c.radius<=c.barrier.x+1e-6&&total(reflected,c)<=energy*(1+1e-8),"rigid collision geometry and energy");
    auto over=initialState(c);over.p={8.9,8};over.v={20,0};over.omega=0;ticks(over,c,100);
    require(!over.barrierHit&&over.barrierPassed&&over.barrierOvertopped&&!over.barrierIntercepted,"clean overtopping must differ from interception");
    c=Config();c.barrier.kind=BarrierKind::Rigid;auto stopped=simulate(c);
    require(stopped.finalState.barrierIntercepted&&stopped.status=="intercepted","stable rigid interception has its own status");
    std::cout<<"stage 6: high-speed barrier hit, cap normal, impulse/energy, overtopping and rigid interception passed\n";
#else
    c=flat();c.barrier.kind=BarrierKind::Rigid;bool rejectedRigid=false;try{validate(c);}catch(const std::invalid_argument&){rejectedRigid=true;}require(rejectedRigid,"rigid stage gate");
#endif
#if ROCK_STAGE >= 7
    c=Config();c.barrier.kind=BarrierKind::Elastic;c.barrier.maxDeflection=5;
    c.dt=.002;const auto elastic=simulate(c);c.dt=.000125;const auto elasticFine=simulate(c);
    require(elastic.finalState.barrierHit&&!elastic.finalState.barrierFailed&&elastic.finalState.barrierIntercepted,"adequate elastic barrier intercepts");
    require(elastic.finalState.maxBarrierDeflection>.01&&elastic.finalState.barrierDissipation>0,"elastic deformation and dissipation");
    require(std::abs(elastic.finalState.maxBarrierDeflection-elasticFine.finalState.maxBarrierDeflection)<.03,"elastic displacement dt convergence");
    require(std::abs(elastic.finalState.peakBarrierForce-elasticFine.finalState.peakBarrierForce)/elasticFine.finalState.peakBarrierForce<.02,"spring/damper peak force must not diverge as j/dt");
    const double initial=elastic.trajectory.front().energy+elastic.finalState.mass*G*elastic.trajectory.front().z;
    for(const auto& p:elastic.trajectory){const double stored=p.barrierFailed?0:.5*c.barrier.stiffness*p.deflection*p.deflection+.5*c.barrier.mass*p.barrierVelocity*p.barrierVelocity;
        require(p.energy+elastic.finalState.mass*G*p.z+stored+p.barrierDissipation<=initial*(1+1e-8),"coupled panel energy bound");}
    require(!elastic.trajectory.front().barrierHit&&!elastic.trajectory.front().barrierFailed,"frame flags describe frame time, not the final state");
    require(elastic.trajectory.back().barrierHit&&elastic.trajectory.back().contactImpulse>0,"sample contains barrier telemetry");
    c.barrier.maxDeflection=.005;const auto broken=simulate(c);
    require(broken.finalState.barrierFailed&&broken.finalState.barrierPassed&&!broken.finalState.barrierIntercepted,"capacity failure must allow passage and not count as interception");
    std::cout<<"stage 7: elastic displacement/force convergence, coupled energy, frame telemetry and failure passage passed\n";
#else
    c=flat();c.barrier.kind=BarrierKind::Elastic;bool rejectedElastic=false;try{validate(c);}catch(const std::invalid_argument&){rejectedElastic=true;}require(rejectedElastic,"elastic stage gate");
#endif
    return 0;
}
}
inline int runExtendedPhysicsTests(){return rock::extended_tests::run();}
