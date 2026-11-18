#include "cli.hpp"
#ifdef _WIN32
#include "entry.hpp"
int wmain(int argc,wchar_t** argv){return rock::windowsEntry(argc,argv,nullptr);}
#else
int main(int argc,char** argv){return rock::runCommandLine(argc,argv);}
#endif
