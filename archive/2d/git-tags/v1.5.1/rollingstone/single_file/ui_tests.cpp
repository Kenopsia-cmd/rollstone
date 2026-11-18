#define wmain rockApplicationMain
#include "rollingstone_tasks_1_3.cpp"
#undef wmain
#include "ui_tests.inc"
int wmain() { try { return runVisibleThreeTests(); } catch(const std::exception& e) { std::cerr << e.what() << std::endl; return 1; } }
