#pragma once
#include "cli.hpp"
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
