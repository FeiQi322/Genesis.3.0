#ifndef DSPEXTAPI_H
#define DSPEXTAPI_H

#include <espMemHeap.h>
#include <sscErrors.h>


#define DSPAPI 


typedef void* TabEntryHandle;   /* Table Entry Handle      */
typedef void* ObjCompHandle;    /* Object Component Handle */
typedef void* scopeHandle;      /* Scope Handle            */


/*========== Virtual Machine Run-time Execution Block ==========*/

typedef struct {

    long          *sp; 		/* stack pointer  */
    memHeapHandle  pageHeap;	/* page execution level heap */
    memHeapHandle  scopeHeap;	/* page scope level heap */
    sscErrorHandle err;		/* error handle */
    streamHandle   out;         /* standard DSP output stream */
    
} dspVMExec;


/*========== DSP run-time function pointer type ==========*/

typedef void FuncPointer(dspVMExec *vm);

/*========== Context extension functions (definition time) ==========*/

/*  Variable definition   */

TabEntryHandle DSPAPI dspDefineVariable (scopeHandle scope, 
					 const char *spec, const char *varname, 
                                         long value);

#define DISPL(f,s) ((int)(&((s*)0)->f)) /* bias of field f in struct s */

/*  Object type definition   */

TabEntryHandle DSPAPI dspDefineObjectType (scopeHandle scope, 
					   const char *objname);

/* Add components to construct the object */

ObjCompHandle DSPAPI dspAddField        (scopeHandle, TabEntryHandle, 
    					 const char *spec, const char *fieldname, 
                                         int displ);
ObjCompHandle DSPAPI dspAddMethod       (scopeHandle, TabEntryHandle, 
					 const char *spec, const char *methodname, 
                                         FuncPointer funcpointer);
ObjCompHandle DSPAPI dspAddProperty     (scopeHandle, TabEntryHandle, 
					 const char *spec, const char *propname, 
                                         FuncPointer getpointer, FuncPointer setpointer);
ObjCompHandle DSPAPI dspAddIndexedAccess(scopeHandle, TabEntryHandle, 
					 const char *spec,
                                         FuncPointer getpointer, FuncPointer setpointer);

/*  Function definition   */

TabEntryHandle DSPAPI dspDefineFunction (scopeHandle scope, 
					 const char *spec, const char *funcname,
                                         FuncPointer funcpointer);

/*========== Run-time callable functions ==========*/

void DSPAPI dspException(dspVMExec *vm, const char* errmess, ...);

typedef void shutdownFP(void *data);
void DSPAPI dspOnPageShutdown   (dspVMExec *vm, shutdownFP *shutdown, void *data);
void DSPAPI dspNoPageShutdown   (dspVMExec *vm, shutdownFP *shutdown, void *data);
void DSPAPI dspOnSessionShutdown(dspVMExec *vm, shutdownFP *shutdown, void *data);
void DSPAPI dspNoSessionShutdown(dspVMExec *vm, shutdownFP *shutdown, void *data);

void DSPAPI dspLock(void);
void DSPAPI dspUnlock(void);

void Dprintf(const char *fmt, ...);

/*========== Functions defined by DSP plugin header ==========*/

char *cloneString(memHeapHandle h, const char *sstr);
 
/*========== DSO function defined by plugin developer 
             and called by plugin header code           ==========*/

void defineDSPExtension(scopeHandle scope);


#endif
