#ifndef FILEIO_H
#define FILEIO_H

#include "classdef.h"

#include <string>

// Truncate a name the way the fixed MAX_NAME_LEN fields in classdef.h store it.
//
// Every external identifier that is later compared against one of those fields
// has to pass through here first, so that both sides of the comparison are cut
// at the same point. The --control list and the --skip-gene values are compared
// against ITEM_STRUCT::name and GROUP_STRUCT::name respectively; matching them
// at full length against a stored name that was truncated loses the match
// silently, which makes a long control sgRNA count as absent and lets a long
// skipped gene be permuted anyway.
std::string nameAsStored(const std::string &name);

//Read input file. File Format: <item id> <group id> <list id> <value>. Return 1 if success, -1 if failure
int ReadFile(char *fileName, GROUP_STRUCT*& groups, int maxGroupNum, int *groupNum, LIST_STRUCT *lists, int maxListNum, int *listNum);

//Save group information to output file. Format <group id> <number of items in the group> <lo-value> <false discovery rate>
int SaveGroupInfo(char *fileName, GROUP_STRUCT *groups, int groupNum);





#endif
