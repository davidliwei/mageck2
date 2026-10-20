#ifndef CLASSDEF_H
#define CLASSDEF_H



#include <stdio.h>
#include <stdlib.h>
#include <memory.h>

#define NDEBUG
#include <assert.h>

// Size of the fixed name buffer held inline by all three structs below, so it
// multiplies by the number of sgRNAs and genes rather than being paid once. At
// 10000 it made GROUP_STRUCT 10,056 bytes and the groups array 1.01 GB. The
// longest names RRA sees are MSigDB gene set ids, which run to roughly 200
// characters; longer ones are truncated with a warning by copyName() in
// fileio.cpp rather than overflowing, which is what the unbounded strcpy calls
// it replaced would have done.
#define MAX_NAME_LEN 512             //maximum length of item name, group name or list name
#define CDF_MAX_ERROR 1E-10        //maximum error in Cumulative Distribution Function estimation in beta statistics

// Initial capacity of the groups array, not a limit: getGroupListNum() doubles
// it as needed while counting the input. It used to be allocated at
// MAX_GROUP_NUM up front, which cost 1.01 GB whatever the input contained.
#define INIT_GROUP_NUM 4096        //initial number of groups; grown on demand
#define MAX_LIST_NUM 1000          //maximum number of list 
#define RAND_PASS_NUM 100          //number of passes in random simulation for computing FDR

#define MAX_WORD_NUM 1000        //maximum number of word

//WL: define boolean variables
//typedef int bool;
//#define true 1
//#define false 0


typedef struct // item definition; i.e., sgRNA
{
	char name[MAX_NAME_LEN];       //name of the item
	int listIndex;                 //index of list storing the item
	double value;                  //value of measurement
	double percentile;             //percentile in the list
	double prob;			//The probability of each sgRNA; added by Wei
  int isChosen;  //whether this sgRNA should be considered in calculation
} ITEM_STRUCT;

typedef struct // group definition; i.e., gene
{
	char name[MAX_NAME_LEN];       //name of the group
	ITEM_STRUCT *items;            //items in the group
	int itemNum;                   //number of items in the group
  int maxItemNum;                // max number of items
	double loValue;                //lo-value in RRA
  double pvalue;                //p value for permutation
	double fdr;                    //false discovery rate
  int isbad;                    //if the lovalue is too low (i.e., higher than the given percentile)
  int goodsgrnas;               //sgRNAs with significant changes
  int controlsgs;               //the numbers of control sgs
} GROUP_STRUCT;

typedef struct //list definition; i.e., gene groups
{
	char name[MAX_NAME_LEN];       //name of the list
	double *values;                //values of items in the list, used for sorting
	int itemNum;                   //number of items in the list
  int maxItemNum;               //max item number
} LIST_STRUCT;


#endif
