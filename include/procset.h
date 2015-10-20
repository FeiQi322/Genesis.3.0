// $Id: procset.h,v 1.2 2010/09/23 14:42:51 mraymond Exp $

#ifndef __PROCSET_H
#define __PROCSET_H

//
// Copyright 2006, Silicon Graphics, Inc.
// ALL RIGHTS RESERVED
//

#include <stdint.h>
#include <unistd.h>
#include <sys/types.h>


typedef enum {
	PSET_CLASS_ASH 	= 0,		// used by Array Services
	PSET_CLASS_USER1 = 1,		// for users
	PSET_CLASS_USER2 = 2,		// for users
	PSET_NUM_CLASSES 
} procset_class_t;

typedef uint64_t procset_t;


//
// Function Prototypes
//

int procset_proc_join_newset(procset_class_t pclass, procset_t pset, 
	uid_t uid, procset_t *r_pset);
	//
	// Create a new Process Set "pset" in the DB for Process Set
	// Class pclass.  Associates user id "uid" with this new Process Set.
	// Associates calling process with the new Process Set.
	//
	// The following error checks are made:
	//	o Process Set "pset" must be an empty set (not existent)
	//	  upon entry.
	//	o The calling process must belong to user "uid".
	//

int procset_proc_join_set(procset_class_t pclass, procset_t pset,
        uid_t uid, procset_t *r_pset);
	//
	// Assigns calling process to Process Set "pset" in the DB for 
	// Process Set Class pclass.  If the Process Set does not 
	// yet exist, create it and associate it with user "uid".
	//

int procset_proc_moveinto_set(procset_class_t pclass, procset_t pset,
        procset_t *);
	//
	// Same as procset_proc_join_set ?
	//

int procset_proc_get_set(procset_class_t pclass, pid_t pid, procset_t *r_pset);
	//
	// Returns in *r_pset the Process Set to which "pid" belongs.
	//

int procset_get_setcount(procset_class_t pclass, int *r_cnt);
	//
	// Returns the number of active Process Sets for Process Set Class
	// pclass.
	//

int procset_get_setlist(procset_class_t pclass, procset_t *set_array,
        int array_size, int *r_cnt);
	//
	// Returns in "set_array" the list of active Process Sets for Process 
	// Set Class pclass.  Argument "array_size" is the size of array 
	// set_array.  On return, *r_cnt contains the actual number of
	// Process Sets returned.
	//

int procset_set_get_pidcount(procset_class_t pclass, procset_t pset,
	int *r_cnt);
	//
	// Returns in *r_cnt the number of processes that are currently
	// members of Process Set "pset".
	//

int procset_set_get_pidlist(procset_class_t, procset_t,
        pid_t *pid_array, int array_size, int *r_cnt);
	//
	// Returns in "pid_array" the list of processes that are
	// members of Process Set "pset". Upon return, *r_cnt is
	// assigned the number of processes written to "pid_array".
	// 

#endif // ! __PROCSET_H
