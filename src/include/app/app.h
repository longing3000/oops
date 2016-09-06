#ifndef APP_H
#define APP_H

#include "include/oops.h"

string PROJECT_PATH;
string LOG_PATH;
string INPUT_PATH;
string OUTPUT_PATH;
string CONFIG_PATH;
string DEBUG_PATH;

cSPINDATA SPIN_DATABASE=cSPINDATA();

clock_t global_start,global_end;
int global_bug_node;

#endif
